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
        self.adata = adata
        self.filter_dict = None

        # 各組織で3つのグループを定義
        self.parenchymal_cells = []  # 組織の主機能を担う細胞 (P)
        self.supporting_cells = []   # 構造・維持・支持細胞 (S)
        self.immune_cells = [] # 免疫細胞 (I)

    # initで指定したmethodでassign
    # 後で下のやつと統合する
    def assign(self, group_weights=None, **kwargs):
        method_name = f"_assign_{self.method}"
        if not hasattr(self, method_name):
            raise ValueError(f"Method {method_name} not found.")
        func = getattr(self, method_name)
        
        # 各グループ内部での比率を計算
        # 各関数の戻り値は (sample_size, 各グループの細胞種数) の DataFrame
        p_df = func(self.parenchymal_cells, **kwargs) if self.parenchymal_cells else None
        s_df = func(self.supporting_cells, **kwargs) if self.supporting_cells else None
        i_df = func(self.immune_cells, **kwargs) if self.immune_cells else None

        # 特定のグループがない可能性を考慮
        group_dfs = [
            ('P', p_df),
            ('S', s_df),
            ('I', i_df)
            ]
        active_groups = [(name, df) for name, df in group_dfs if df is not None]

        # グループ間の重み付け
        if group_weights is not None:
            # 指定された固定重みを使用
            weights = np.array(group_weights)
            mask = [df is not None for _, df in group_dfs]
            active_w = weights[mask]
            active_w = active_w / active_w.sum()
            active_ratios = np.tile(active_w, (self.sample_size, 1))
        else:
            # 重み指定がない場合、グループ間の比率もランダムに決定 (Dirichlet alpha=1.0)
            n_groups = len(active_groups)
            active_ratios = np.random.dirichlet([1.0] * n_groups, size=self.sample_size)


        # 統合とスケーリング
        summaries = []
        for g_idx, (name, df) in enumerate(active_groups):
            if df is None:
                continue

            df_scaled = df.mul(active_ratios[:, g_idx], axis=0)
            summaries.append(df_scaled)

        # 全細胞種を結合
        self.summary_df = pd.concat(summaries, axis=1)
        # 最後に全体で合計1になるよう微調整 (なくてもいいかも？)
        self.summary_df = self.summary_df.div(self.summary_df.sum(axis=1), axis=0)
        return self.summary_df
                
    
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
        data = np.random.dirichlet(alpha_vec, size=self.sample_size)

        if do_viz and len(cell_types) > 1:
            plt.hist(data[:, 1], bins=50, alpha=0.7, color='blue', label=f'alpha={alpha}')
            plt.xlabel('Value')
            plt.ylabel('Frequency')
            plt.title('Distribution of proportion')
            plt.legend()
            plt.show()
        
        return pd.DataFrame(data, columns=cell_types)

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
        # 細胞×遺伝子のnp.arrayを獲得
        """Helper method to get expression matrix from adata."""
        if hasattr(self.adata.X, 'todense'):
            return np.array(self.adata.X.todense())
        else:
            return np.array(self.adata.X)

    # ほぼ共通なのでBaseSimulater内部に移動、target_colを外部指定にして汎用性獲得
    # trainのみが欲しいときはtrain_ratio=1.0でcreate_sim_bulkでmode='train'にすればOK
    def split_cell_idx(self, target_col='cell_type', save_dir='./data/cell_idx', train_ratio=0.7):
        """Split cell indices into train/test sets."""
        if self.adata is None:
            raise ValueError("adata must be set before splitting cell indices")

        # 共通のシード (ループの外に出したほうがよい？)
        random.seed(42)
        np.random.seed(42)

        info_df = self.adata.obs
            
        # 分割の仕方は一旦これでいいか？
        cell_types = self.parenchymal_cells + self.supporting_cells + self.immune_cells
        cell_idx_dict = {}
        
        for cell in cell_types:
            cellname = self.filter_dict[cell] if self.filter_dict is not None else cell
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

    def set_data(self, summary_df=None, cell_idx_dict=None, adata=None):
        if summary_df is not None:
            self.summary_df = summary_df
        if cell_idx_dict is not None:
            self.cell_idx_dict = cell_idx_dict
        if adata is not None:
            self.adata = adata

    def create_sim_bulk(self, pool_size=500, mode='train', adata_path=None):
        """Create simulated bulk expression data."""
        if self.summary_df is None or self.cell_idx_dict is None:
            raise ValueError("summary_df and cell_idx_dict must be set")
        
        # Load adata if path provided, otherwise use self.adata
        if adata_path:
            adata = sc.read_h5ad(adata_path)
        elif self.adata is not None:
            adata = self.adata
        else:
            raise ValueError("Please specify the correct data path, or register the adata object!")
        
        # adataのQC
        self._ensure_counts()

        # expression matrixをnp.arrayのrawデータとして保持、重複ありの回収に対応
        raw_exp = self._get_expression_matrix()

        total_cells = self.summary_df.columns.tolist()
        pooled_exp = []
        np.random.seed(42)
        
        for idx in tqdm(range(len(self.summary_df))):
            p_list = self.summary_df.iloc[idx].values
            bulk_single = np.zeros(raw_exp.shape[1]) # 初期化を遺伝子数サイズの配列に変更
            
            for j, p in enumerate(p_list):
                cell = total_cells[j]
                cellname = cell
                tmp_size = int(pool_size * p)
                candi_idx = self.cell_idx_dict[cellname][mode] # trainとtestはそれぞれ別で実行

                # Generate reproducible random indices
                rng = np.random.RandomState(42 + idx * len(total_cells) + j)
                select_idx = rng.choice(candi_idx, size=tmp_size, replace=len(candi_idx) < tmp_size) # 細胞数が足りないときは重複を許可
                
                # Get expression data for selected cells
                # unique_colsを廃止して普通の足し合わせにしました。
                # 外部でraw_expを作成しておくことで重複ありでも回収可能
                bulk_single += np.array(raw_exp[select_idx].sum(axis=0)).flatten()
            
            if bulk_single is not None:
                pooled_exp.append(bulk_single.flatten())
        
        # Create DataFrame
        pooled_exp = np.array(pooled_exp)
        bulk_df = pd.DataFrame(pooled_exp.T)
        bulk_df.index = adata.var_names  # gene names
        return bulk_df

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


class Liver_Simulator(BaseSimulator):
    def __init__(self, adata, sample_size=8000, method='dirichlet'):
        super().__init__(adata, sample_size, method)
        self.parenchymal_cells = ['Hepatocytes']
        self.supporting_cells = ['Cholangiocytes', 'Fibroblasts']
        self.immune_cells = ['Neutrophils', 'Macrophages', 'Monocytes', 'NK', 'CD4Tcells', 'CD8Tcells', 'B cells', 'Dendritic']