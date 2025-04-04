#!/usr/bin/env python3
"""
Created on 2025-03-11 (Tue) 14:07:11

@author: I.Azuma
"""
# %%
import pandas as pd


def load_counts(filein='', targets=[]):
    """
    ### Input
    	Depth	Name	                            #Cells
    0	NaN	    Specimen_001_APAP_Ly6G_1_006.fcs	200000
    1	>	    CD45+Dump-	                        122377
    2	> >	    Single Cells	                    113197
    3	> > >	lymphocytes	                        97321
    4	> > >	main_population	                    112431


    ### Output
    Name	main_population	lymphocytes	Single Cells	Others	NK
    Sample					
    Specimen_001_APAP_1_003.fcs	65145	53733	65685	48083	4955
    Specimen_001_APAP_2_004.fcs	55090	47193	55546	42900	3709
    Specimen_001_APAP_3_005.fcs	130694	106052	131663	101058	3830
    Specimen_001_APAP_Ly6G_1_006.fcs	112431	97321	113197	86154	10616
    Specimen_001_APAP_Ly6G_2_007.fcs	119637	91245	120841	79270	15527

    """
    df = pd.read_csv(filein, usecols=["Name", "#Cells", "Depth"])

    # Depth列を処理して階層を抽出
    df["Depth"] = df["Depth"].fillna("").apply(lambda x: len(x)//2)

    # Depthが0の行をサンプルとして抽出
    samples = df[df["Depth"] == 0]

    # 結果を保存するリスト
    aggregated_data = []

    # サンプルごとに処理
    for i in range(len(samples)):
        sample_start_index = samples.iloc[i].name
        if i + 1 < len(samples):
            # 次のサンプルまでのデータを抽出
            sample_end_index = samples.iloc[i + 1].name
            sample_data = df.iloc[sample_start_index + 1: sample_end_index]
        else:
            # 最後のサンプルの場合は最後まで抽出
            sample_data = df.iloc[sample_start_index + 1:]
        
        # サンプル名を取り出し、そのサンプルに関連する細胞データを集約
        # 同名の細胞が複数ある場合は、Depthが大きいものを残す
        sample_name = df.iloc[sample_start_index]["Name"]
        sample_cells = sample_data.sort_values(by=["Name", "Depth"], ascending=False).drop_duplicates(subset=["Name"])

        if len(targets) > 0:
            sample_cells = sample_cells[sample_cells["Name"].isin(targets)]
        
        # サンプル名を列として追加
        sample_cells["Sample"] = sample_name
        
        # 結果をリストに追加
        aggregated_data.append(sample_cells)

    # 集約データをデータフレームに変換
    aggregated_df = pd.concat(aggregated_data, axis=0)

    # サンプルごとに細胞数を列として並べる
    #aggregated_df_pivot = aggregated_df.pivot(index="Sample", columns=["Name", "Depth"], values="#Cells")
    aggregated_df_pivot = aggregated_df.pivot(index="Sample", columns=["Name"], values="#Cells")

    display(aggregated_df_pivot.style.background_gradient(cmap='viridis', axis=0))

    return aggregated_df_pivot

def norm_counts(df, base_col):
    """
    ### Input: load_counts() output

    """
    df_norm = df.div(df[base_col], axis=0)
    display(df_norm.style.background_gradient(cmap='YlOrRd', axis=0))

    return df_norm
