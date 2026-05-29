import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import torch
from scvi.external import SysVI
from tqdm import tqdm  # 進捗表示用

sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=150, color_map='viridis')
sc.logging.print_versions()

def prep_train(combined_adata, n_top_genes=3000, batch_key="origin"):
    # 全遺伝子で学習すると重すぎる＆ノイズが入るため、通常は上位2000-4000遺伝子に絞ります
    sc.pp.highly_variable_genes(
        combined_adata,
        n_top_genes=n_top_genes,
        batch_key=batch_key, # バッチごとに分散を見て選ぶ
        subset=True          # 学習用にサブセット化
    )

    # モデルセットアップ
    SysVI.setup_anndata(
        combined_adata,
        batch_key=batch_key,           # ここに 'Liver' / 'PBMC' の列を指定
        categorical_covariate_keys=[], # もしドナーIDなど他のバッチがあればここに追加
        layer=None                     # .X (log1p正規化済み) を使用
    )

    # モデル初期化（VampPrior使用）
    model = SysVI(
        combined_adata,
        prior="vamp",
        n_prior_components=10,
    )

    return model


def plot_loss(model, save_path_root="/save/dir"):
    history = model.history

    # グラフの作成（1行4列）横幅を24に拡張
    fig, axes = plt.subplots(1, 4, figsize=(24, 5))

    # 1. トータル Loss (Total Loss / ELBO)
    # --- Train ---
    if 'train_loss' in history:
        history['train_loss'].plot(ax=axes[0], label='Train')
    elif 'elbo_train' in history:
        history['elbo_train'].plot(ax=axes[0], label='Train')
    else:
        print("Warning: トータルLoss(Train)のキーが見つかりません。")
    # --- Validation ---
    if 'validation_loss' in history:
        history['validation_loss'].plot(ax=axes[0], label='Validation')
    elif 'elbo_validation' in history:
        history['elbo_validation'].plot(ax=axes[0], label='Validation')

    axes[0].set_title('Total Loss')
    axes[0].set_xlabel('Epochs')
    axes[0].legend() # 凡例を表示

    # 2. 再構成誤差 (Reconstruction Loss)
    if 'reconstruction_loss_train' in history:
        history['reconstruction_loss_train'].plot(ax=axes[1], label='Train')
    if 'reconstruction_loss_validation' in history:
        history['reconstruction_loss_validation'].plot(ax=axes[1], label='Validation')
    axes[1].set_title('Reconstruction Loss')
    axes[1].set_xlabel('Epochs')
    axes[1].legend()

    # 3. KLダイバージェンス (KL Divergence)
    if 'kl_local_train' in history:
        history['kl_local_train'].plot(ax=axes[2], label='Train')
    if 'kl_local_validation' in history:
        history['kl_local_validation'].plot(ax=axes[2], label='Validation')
    axes[2].set_title('KL Divergence')
    axes[2].set_xlabel('Epochs')
    axes[2].legend()

    # 4. サイクル一貫性損失 (Cycle Loss)
    if 'cycle_loss_train' in history:
        history['cycle_loss_train'].plot(ax=axes[3], label='Train')
    elif 'cycle_loss' in history:
        history['cycle_loss'].plot(ax=axes[3], label='Train')

    if 'cycle_loss_validation' in history:
        history['cycle_loss_validation'].plot(ax=axes[3], label='Validation')
    axes[3].set_title('Cycle Consistency Loss')
    axes[3].set_xlabel('Epochs')
    axes[3].legend()

    plt.tight_layout()
    save_path = os.path.join(save_path_root, "loss.png")
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def get_latent_feature(model, combined_adata):
    # latent取得
    z = model.get_latent_representation()

    # AnnDataに入れる
    adata_latent = sc.AnnData(z)

    # combined_adata.obsをlatent_adataにcopyする
    adata_latent.obs = combined_adata.obs.copy()

    return adata_latent


def umap_and_viz(adata, n_neighbors=15, n_pcs=None, save_path_root="/save/dir", data_type="latent"):
    """
    data_typeはfigの名前分けのみに使用
    基本はlatentとpseude
    """
    # 近傍グラフ & UMAP
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(adata)

    # 可視化（バッチ）
    sc.pl.umap(adata, color="origin", show=False)
    save_path = os.path.join(save_path_root, str(data_type)+"_umap_origin.png")
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()

    # 可視化（細胞型）
    sc.pl.umap(adata, color="cell_type", show=False)
    save_path = os.path.join(save_path_root, str(data_type)+"_umap_cell_type.png")
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()


def get_sysvi_corrected_mini_batch(model, adata, target_batch_name, batch_size=1024, layer_name="X"):
    module = model.module
    module.eval()
    device = list(module.parameters())[0].device
    
    # 1. データの準備
    x_data = adata.layers[layer_name] if layer_name in adata.layers else adata.X
    n_obs = x_data.shape[0]
    
    # バッチIDの準備
    registry = model.adata_manager.get_state_registry("batch")
    mapping = registry.categorical_mapping
    batch_series = adata.obs[registry.original_key].astype("category").cat.set_categories(mapping)
    real_batch_codes = batch_series.cat.codes.values
    
    try:
        target_id = np.where(mapping == target_batch_name)[0][0]
    except IndexError:
        raise ValueError(f"Target batch '{target_batch_name}' not found.")

    corrected_values_list = []

    # 2. ミニバッチ処理
    with torch.no_grad():
        for i in tqdm(range(0, n_obs, batch_size), desc="Decoding in mini-batches"):
            # バッチのスライス
            start, end = i, min(i + batch_size, n_obs)
            
            # 必要な分だけGPUへ送る
            x_chunk = x_data[start:end]
            if hasattr(x_chunk, "toarray"): x_chunk = x_chunk.toarray()
            x_tensor = torch.from_numpy(x_chunk).float().to(device)
            
            rb_indices = torch.from_numpy(real_batch_codes[start:end]).long().to(device).unsqueeze(1)
            tb_indices = torch.full((end - start, 1), target_id, device=device).long()

            # --- Forward Pass ---
            # Encoder
            enc_out = module.encoder(x=x_tensor, batch_index=rb_indices, cat_list=[])
            z = enc_out["q"] # 疑似バルク用なら .mean を使うのもアリ
            
            # Decoder
            dec_out = module.decoder(x=z, batch_index=tb_indices, cat_list=[])
            
            # 結果を抽出して即座にCPUへ（GPUメモリを解放するため）
            if "q_dist" in dec_out:
                val = dec_out["q_dist"].mean.cpu().numpy()
            else:
                val = dec_out["q"].cpu().numpy()
            
            corrected_values_list.append(val)
            
            # 明示的なメモリ解放（念のため）
            del x_tensor, rb_indices, tb_indices, enc_out, dec_out
            if i % (batch_size * 10) == 0:
                torch.cuda.empty_cache()

    # 3. CPU上で結合してDataFrame化
    corrected_values = np.concatenate(corrected_values_list, axis=0)
    return pd.DataFrame(corrected_values, index=adata.obs_names, columns=adata.var_names)


def plot_gene_dist(corrected_df, save_path_root="/save/dir", data_type="decoder_raw"):
    """
    data_typeはfigの名前分けのみに使用
    基本はdecoder_rawとdecoder_fix
    """
    # corrected_df が DataFrame の場合
    for i in range(100):
        # .iloc を使って列（遺伝子）を指定
        show_data = corrected_df.iloc[:, i] 
        plt.hist(show_data, bins=100, alpha=0.3) # alphaを薄くすると重なりが見やすい

    plt.xlabel("Reconstructed Expression (Linear)")
    plt.ylabel("Frequency")
    plt.title("Distribution of 100 Genes")
    save_path = os.path.join(save_path_root, str(data_type)+"_pseud_dist.png")
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()


# =============================================================================
def sysvi_execution(combined_adata, n_top_genes=3000, batch_key="origin",
                    plan_kwargs={
                        "reconstruction_weight": 1.0,      # 再構成誤差を重視
                        "kl_weight": 1.0,                  # KLを弱める
                        "z_distance_cycle_weight": 2.0     # サイクル整合性を強力に
                        },
                    target_batch_name="liver", 
                    save_path_root="/save/dir",
                    max_epochs=300,
                    save_adata=None
                    ):
    # モデルの準備
    model = prep_train(combined_adata, n_top_genes=n_top_genes, batch_key=batch_key)

    # 学習
    model.train(
        max_epochs=max_epochs,
        early_stopping=True,
        # ▼ GPU設定
        accelerator="gpu",
        devices=1,
        batch_size=4096,  
        plan_kwargs=plan_kwargs
    )

    # lossの描画
    plot_loss(model, save_path_root=save_path_root)

    # 潜在空間情報の獲得
    adata_latent = get_latent_feature(model, combined_adata)
    umap_and_viz(adata_latent, n_neighbors=15, save_path_root=save_path_root, data_type="latent") # 潜在空間はKLがかかっているのでスケーリング必要なし

    # 一度きれいにする
    torch.cuda.empty_cache()

    # 全員を "liver" の条件に合わせて補正
    corrected_df = get_sysvi_corrected_mini_batch(
        model, 
        combined_adata, 
        target_batch_name=target_batch_name, 
        layer_name="X"
    )

    # 出力の発現分布を確認する
    plot_gene_dist(corrected_df, save_path_root=save_path_root, data_type="decoder_raw")

    # 0以下の値を0に丸め込む
    corrected_df = corrected_df.clip(lower=0)

    # 0に丸め込んだのちもう一度発現分布を確認する
    plot_gene_dist(corrected_df, save_path_root=save_path_root, data_type="decoder_fix")

    # combined_adata.Xをデコーダの出力の発現量 (対数変換)で上書きする
    combined_adata.X = corrected_df.loc[combined_adata.obs_names, combined_adata.var_names].values

    # UMAP描画用疑似データの準備
    adata_pseudo = combined_adata.copy() # データをコピー
    sc.pp.scale(adata_pseudo, max_value=10) # データをスケール
    sc.pp.pca(adata_pseudo, n_comps=50) # PCAを実行、3000次元から下げる
    umap_and_viz(adata_pseudo, n_pcs=30, save_path_root=save_path_root, data_type="pseudo") # PCAデータを用いるためn_pcsをintで指定

    if save_adata:
        save_adata_path = os.path.join(save_path_root, save_adata)
        # 一旦ここでデータを保存する
        combined_adata.write(save_adata_path)