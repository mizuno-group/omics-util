# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 12:59:35 2019

Binomial test

@author: tadahaya
"""
import pandas as pd
import numpy as np
import scipy.stats as stats
import statsmodels.stats.multitest as multitest
from scipy.stats import rankdata


class Calculator():
    def __init__(self):
        self.res = pd.DataFrame()


    def calc(self,obj,ref,whole,focus=None,**kwargs):
        self.res = do_binom(obj,ref,whole,focus=None,**kwargs)
        return self.res


    def get_details(self):
        return self.res


def do_binom(obj,ref,whole,focus=None,correction="fdr_bh",mode="greater"):
    """
    conduct Binomial test and obtain p value corrected for multiple tests
    all elements should be given as ID, except for term
    
    Parameters
    ----------
    obj: set
        a set of variables in signature of interest
        
    ref: dict
        a dict of term and member list of datasets

    whole: set
        a set of whole variables adjusted between datasets and interest 
        
    method: str
        indicate method for correcting multiple tests
        depend on "statsmodels.stats.multitest.multipletests"
        
    focus: int
        export results by XX th lowest p value

    mode: str
        indicate the type of significant judging
        "greater", "two-sided", or "less"

    """
    pval = []
    overlap = []
    total = []
    hit = []
    ap = pval.append
    ap2 = overlap.append
    ap3 = total.append
    ap4 = hit.append
    n_whole = len(whole)
    keys = list(ref.keys())
    values = list(ref.values())
    for m in values:
        n_set = len(m)
        pc = n_set/n_whole
        inter = obj & m
        n_inter = len(inter)
        ap(stats.binom_test(n_inter,n_set,pc,alternative=mode))
        ap2(inter)
        ap3(len(m))
        ap4(n_inter)
    if len(pval)==0:
        res = pd.DataFrame(columns=["p value","adjusted p value","overlap","hit No.",
                                    "total No."])
    else:
        res = pd.DataFrame({"p value":pval,"overlap":overlap,
                            "hit No.":hit,"total No.":total},index=keys).sort_values(by="p value")
        fxn = lambda x: len(x) > 0
        res = res[res["overlap"].map(fxn)]
        if res.shape[0]!=0:
            res["adjusted p value"] = multitest.multipletests(res["p value"],alpha=0.05,method=correction)[1]
            res = res.sort_values(by="p value")
            res = res.loc[:,["p value","adjusted p value","overlap","hit No.","total No."]] # sort
    if (focus is None) or (res.shape[0]==0):
        pass
    else:
        res = res.iloc[:focus,:]
    return res