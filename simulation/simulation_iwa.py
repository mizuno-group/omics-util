# -*- coding: utf-8 -*-
"""
Created on 2024-11-27 (Wed) 21:52:13

@author: I.Azuma
"""

import os
import random
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from tqdm import tqdm

# adataのlog判定で使用、スパースかどうかの確認
from scipy.sparse import issparse

class BaseSimulator:
    def __init__(self, adata, sample_size, method='dirichlet'):
        self.sample_size = sample_size
        self.method = method # 現在はディリクレ分布と一様分布のみ
        self.summary_df = None
        self.cell_idx_dict = None
        self.adata = None

        # 各組織で以下を定義
        self.basic_cells = [] # 免疫細胞をはじめとした共通細胞
        self.tissue_specific_cells = [] # 組織特異的な実質細胞および常在性の免疫細胞

    def _ensure_counts(self):
        """adata.Xがlog1pならカウントに戻し、負の値を0にする処理"""
        if self.adata is None: return

        # 1. ログ変換の解除
        if 'log1p' in self.adata.uns:
            print("Log-transformation detected. Reverting to linear scale...")
            # 直接adata.Xを上書き
            if issparse(self.adata.X):
                self.adata.X.data = np.expm1(self.adata.X.data)
            else:
                self.adata.X = np.expm1(self.adata.X)
            del self.adata.uns['log1p'] # 対数変換した記録を削除

        # 2. 負の値を0に丸める (バッチ補正後の微小な負値を排除)
        if issparse(self.adata.X):
            self.adata.X.data = np.maximum(self.adata.X.data, 0)
        else:
            self.adata.X = np.maximum(self.adata.X, 0)
        
        print("Data check: Negative values clipped to 0. This adata contains count data.")

    # initで指定したmethodでassign
    # 後で下のやつと統合する
    def assign(self, tissue_w=None, **kwargs):
        method_name = f"_assign_{self.method}"
        if not hasattr(self, method_name):
            raise ValueError(f"Method {method_name} not found.")
        
        func = getattr(self, method_name)
        
        # 実行
        basic_summary = func(self.basic_cells, **kwargs)
        tissue_summary = func(self.tissue_specific_cells, **kwargs)
    
    # 一様分布で細胞比率を決定
    def _assign_uniform(self, cell_types: list, sparse=True):
        """Generate uniform distribution for cell type proportions."""
        final_res = []
        for idx in range(self.sample_size):
            np.random.seed(seed=idx)
            if sparse:
                # Randomly select subset of cell types
                use_cell_types = np.random.choice(
                    cell_types, 
                    size=np.random.randint(1, len(cell_types) + 1), 
                    replace=False
                )
                p_list = np.random.rand(len(use_cell_types))
                
                # Create full proportion list
                final_p_list = np.zeros(len(cell_types))
                for j, c in enumerate(use_cell_types):
                    final_p_list[cell_types.index(c)] = p_list[j]
            else:
                final_p_list = np.random.rand(len(cell_types))
            
            # Normalize to sum to 1
            norm_p_list = final_p_list / final_p_list.sum() if final_p_list.sum() > 0 else final_p_list
            final_res.append(norm_p_list)
        
        return pd.DataFrame(final_res, columns=cell_types)
    
    # ディリクレ分布で細胞比率を決定
    def _assign_dirichlet(self, cell_types: list, alpha=1.0, do_viz=False):
        """Generate Dirichlet distribution for cell type proportions."""
        alpha_vec = [alpha] * len(cell_types)
        np.random.seed(seed=42)
        data = np.random.dirichlet(alpha_vec, size=self.sample_size)

        if do_viz and len(cell_types) > 1:
            plt.hist(data[:, 1], bins=50, alpha=0.7, color='blue', label=f'alpha={alpha}')
            plt.xlabel('Value')
            plt.ylabel('Frequency')
            plt.title('Distribution of proportion')
            plt.legend()
            plt.show()
        
        return pd.DataFrame(data, columns=cell_types)

    def assign(self, nonim_w=None):
        """Assign cell type proportions using specified method."""
        methods = {
            'uniform_sparse': lambda types: self.assign_uniform(types, sparse=True),
            'uniform': lambda types: self.assign_uniform(types, sparse=False), 
            'dirichlet': lambda types: self.assign_dirichlet(types, alpha=1.0)
        }
        
        if self.method not in methods:
            raise ValueError(f"Method not supported. Choose from: {list(methods.keys())}")
        
        # Generate proportions for immune and non-immune cells
        im_summary = methods[self.method](self.immune_cells)
        non_im_summary = methods[self.method](self.non_immune_cells)
        
        # Normalize to sum to 1
        self.im_summary = im_summary.div(im_summary.sum(axis=1), axis=0)
        self.non_im_summary = non_im_summary.div(non_im_summary.sum(axis=1), axis=0)

        # Apply random weights if specified
        if nonim_w is not None:
            random_w = np.random.uniform(low=nonim_w, high=1.0, size=self.sample_size)
            for i in range(self.sample_size):
                w = random_w[i]
                self.im_summary.iloc[i] *= (1 - w)
                self.non_im_summary.iloc[i] *= w

        # Combine and normalize final result
        summary_df = pd.concat([self.im_summary, self.non_im_summary], axis=1)
        self.summary_df = summary_df.div(summary_df.sum(axis=1), axis=0)

    # 多分ここはさわらなくていい
    def create_ref(self):
        """Create reference matrix from training data."""
        raw_exp = self._get_expression_matrix()
        pooled_exp = []
        
        for cell_type, idx_dict in self.cell_idx_dict.items():
            train_idx = idx_dict['train']
            tmp_mean = raw_exp[train_idx].mean(axis=0)
            pooled_exp.append(tmp_mean)

        ref_df = pd.DataFrame(pooled_exp).T
        ref_df.index = self.adata.var_names
        ref_df.columns = self.cell_idx_dict.keys()
        return ref_df
    
    def _get_expression_matrix(self):
        """Helper method to get expression matrix from adata."""
        if hasattr(self.adata.X, 'todense'):
            return np.array(self.adata.X.todense())
        else:
            return np.array(self.adata.X)
            
    
    def set_data(self, summary_df=None, cell_idx_dict=None, adata=None):
        if summary_df is not None:
            self.summary_df = summary_df
        if cell_idx_dict is not None:
            self.cell_idx_dict = cell_idx_dict
        if adata is not None:
            self.adata = adata

    # ほぼ共通なのでBaseSimulater内部に移動、target_colを外部指定にして汎用性獲得
    def split_cell_idx(self, target_col='cell_type', save_dir='./data/cell_idx', train_ratio=0.7):
        """Split cell indices into train/test sets."""
        if self.adata is None:
            raise ValueError("adata must be set before splitting cell indices")

        # 共通のシード (ループの外に出したほうがよい？)
        random.seed(42)
        np.random.seed(42)

        info_df = self.adata.obs
            
        # immune/non_immuneなどの分割の仕方は後で統一
        cell_types = self.immune_cells + self.non_immune_cells
        cell_idx_dict = {}
        
        for cell in cell_types:
            cellname = self.filter_dict[cell]
            cell_mask = info_df[target_col] == cell # target_colを外部指定にして汎用性を獲得
            target_idx = np.where(cell_mask)[0].tolist() # whereにして位置をしっかりと獲得
            cell_size = len(target_idx)

            # 細胞が検出されなかったときは警告
            if cell_size == 0:
                print(f"Warning: {cell} not found in {target_col}")
                continue
            else:
                print(f"{cell}: {cell_size} cells detected")
            
            # Train/Test split
            shuffle_idx = random.sample(target_idx, cell_size)
            split_point = int(cell_size * train_ratio)
            train_idx = shuffle_idx[:split_point]
            test_idx = shuffle_idx[split_point:]
            cell_idx_dict[cell] = {'train': train_idx, 'test': test_idx}

            # Save indices if directory specified
            if save_dir:
                os.makedirs(save_dir, exist_ok=True)
                pd.to_pickle(train_idx, os.path.join(save_dir, f'{cellname}_train_idx.pkl'))
                pd.to_pickle(test_idx, os.path.join(save_dir, f'{cellname}_test_idx.pkl'))
        
        self.cell_idx_dict = cell_idx_dict