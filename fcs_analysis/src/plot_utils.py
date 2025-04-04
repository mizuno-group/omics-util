# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 10:23:01 2022

@author: I.Azuma
"""
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors

def plot_box(merge_df,immune="Monocyte"):
    """
    Box Plot for detecting FACS immune cell population

    Parameters
    ----------
    merge_df : dataframe
                      Ctrl_1    Ctrl_2    Ctrl_3  ...      CA_3      CA_4       CA5
        abT         0.144218  0.354034  0.128803  ...  0.269444  0.269753  0.216667
        gdT         0.015701  0.044586  0.017606  ...  0.038849  0.035183  0.052812
        NKT         0.207992  0.356555  0.173388  ...  0.034363  0.017509  0.032319
    immune : TYPE, str
        Immune cell number. The default is "Monocyte".

    """
    immune_df = merge_df.T[[immune]]
    variable = [t.split("_")[0] for t in merge_df.columns.tolist()]
    immune_df["variable"]=variable
        
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    sns.boxplot(x='variable', y=immune, data=immune_df, showfliers=False, ax=ax)
    sns.stripplot(x='variable', y=immune, data=immune_df, jitter=True, color='black', ax=ax)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().yaxis.set_ticks_position('left')
    plt.gca().xaxis.set_ticks_position('bottom')
    plt.show()

def plot_bar(merge_df,immune="Monocyte",sample_sorting=True,sort_sample=["Ctrl","APAP","MDA","ANIT","CIV","CA"]):
    """
    Bar Plot for detecting FACS immune cell population

    Parameters
    ----------
    merge_df : dataframe
                      Ctrl_1    Ctrl_2    Ctrl_3  ...      CA_3      CA_4       CA5
        abT         0.144218  0.354034  0.128803  ...  0.269444  0.269753  0.216667
        gdT         0.015701  0.044586  0.017606  ...  0.038849  0.035183  0.052812
        NKT         0.207992  0.356555  0.173388  ...  0.034363  0.017509  0.032319
    immune : TYPE, str
        Immune cell number. The default is "Monocyte".

    """
    immune_df = merge_df.T[[immune]]
    immune_df.index = [t.split("_")[0] for t in immune_df.index.tolist()]
    if sample_sorting:
        immune_df = immune_df.T[sort_sample].T
    else:
        pass
    samples = immune_df.index.unique()
    summary = [immune_df.loc[sample,immune] for sample in samples]
    means = [immune_df.loc[sample,immune].mean() for sample in samples]
    bar = [immune_df.loc[sample,immune].sem() for sample in samples]
    
    fig, ax = plt.subplots()
    # barplot
    ax.bar(range(len(samples)),means, yerr=bar, capsize=12,width=0.6,color="firebrick",alpha=0.8)
    # dot plot
    for (lis,numb) in zip(summary,range(len(summary))):
        x = [numb]*len(lis)
        ax.plot(x, [float(i) for i in lis],marker='o',color='black',markersize=8,linestyle='None')
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=90)
    ax.set_ylabel('Cell Population')
    ax.set_xlabel('Conditions')
    
    #ax.set_xlim(-1,len(samples))
    ax.set_title(immune)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().yaxis.set_ticks_position('left')
    plt.gca().xaxis.set_ticks_position('bottom')
    plt.show()

def plot_group_scatter(merged_df, cell1='neutrophil', cell2='B'):
    """
    Bar Plot for detecting FACS immune cell population

    Parameters
    ----------
    merge_df : dataframe
                      Ctrl_1    Ctrl_2    Ctrl_3  ...      CA_3      CA_4       CA5
        abT         0.144218  0.354034  0.128803  ...  0.269444  0.269753  0.216667
        gdT         0.015701  0.044586  0.017606  ...  0.038849  0.035183  0.052812
        NKT         0.207992  0.356555  0.173388  ...  0.034363  0.017509  0.032319
    immune : TYPE, str
        Immune cell number. The default is "Monocyte".

    """
    merged_df = merged_df.T
    v1 = merged_df[cell1]
    v2 = merged_df[cell2]

    # カラーパレットの準備
    tab_colors = mcolors.TABLEAU_COLORS
    color_list = list(tab_colors.keys())

    merged_df['group'] = [t.split("_")[0] for t in merged_df.index]
    colors = {}
    for i, group in enumerate(merged_df['group'].unique()):
        colors[group] = color_list[i]
    #colors = {'PBS-PBS': 'red', 'APAP-PBS': 'blue', 'APAP-Ly6G': 'green'}

    fig, ax = plt.subplots(figsize=(5, 5))
    for group in colors.keys():
        plt.scatter(v1[merged_df['group'] == group], v2[merged_df['group'] == group], label=group, color=colors[group], alpha=0.8, zorder=2)

    plt.legend()
    plt.xlabel(cell1)
    plt.ylabel(cell2)
    plt.grid(zorder=1)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().yaxis.set_ticks_position('left')
    plt.gca().xaxis.set_ticks_position('bottom')
    plt.show()

# %% Legacy codes FIXME: remove or modify
def plot_traffick(df,immune="Monocyte",sample_sorting=True,sort_sample=["APAP#12hr","APAP#24hr","APAP#48hr"]):
    immune_df = df.T[[immune]]
    immune_df.index = [t.split("_")[0] for t in immune_df.index.tolist()]
    if sample_sorting:
        immune_df = immune_df.T[sort_sample].T
    else:
        pass
    samples = immune_df.index.unique().tolist()
    summary = [immune_df.loc[sample,immune] for sample in samples]
    means = [immune_df.loc[sample,immune].mean() for sample in samples]
    bar = [immune_df.loc[sample,immune].sem() for sample in samples]
    
    fig, ax = plt.subplots()
    # barplot
    ax.bar([0],[1], yerr=[0], capsize=12,width=0.4,color="yellow",alpha=0.8)
    ax.bar(range(1,len(samples)+1),means, yerr=bar, capsize=12,width=0.4,color="firebrick",alpha=0.8)
    
    #ax.bar(range(len(samples)+1),[1]+means, yerr=[0] + bar, capsize=12,width=0.4,color="firebrick",alpha=0.8)
    # dot plot
    for (lis,numb) in zip(summary,range(len(summary))):
        x = [numb+1]*len(lis)
        ax.plot(x, [float(i) for i in lis],marker='o',color='black',markersize=8,linestyle='None')
    ax.set_xticks(range(len(samples)+1))
    ax.set_xticklabels(["Ctrl"]+samples, rotation=90)
    ax.set_ylabel('Cell Population')
    ax.set_xlabel('Conditions')
    
    ax.set_xlim(-0.5,len(samples)+1)
    ax.set_title(immune)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().yaxis.set_ticks_position('left')
    plt.gca().xaxis.set_ticks_position('bottom')
    plt.show()
