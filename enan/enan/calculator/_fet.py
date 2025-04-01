# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 12:59:35 2019

Fisher's exact test

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
        self.res = do_fet(obj,ref,whole,focus=None,**kwargs)
        return self.res


    def get_details(self):
        return self.res


def do_fet(obj,ref,whole,correction="fdr_bh",focus=None,mode="greater"):
    """
    conduct Fisher's exact test p value corrected for multiple tests
    all elements should be given as ID, except for term
    
    Parameters
    ----------
    obj: set
        a set of variables in signature of interest
        
    ref: dict
        a dict of term and member list of datasets

    whole: set
        a set of whole variables adjusted between datasets and interest 
        
    correction: str
        indicate method for correcting multiple tests
        depend on "statsmodels.stats.multitest.multipletests"
        
    focus: int
        export results by XX th lowest p value
    
    mode: int
        indicate the type of test
        "two-sided", "less", "greater"
        in many case of ontology analysis, "greater" is appropriate

    ### contingency table ###
                      signature(+)  signature(-)  sum
    
    dataset var (+)       n11           n12       n1p
    
    dataset var (-)       n21           n22       n2p
    
        sum               np1           np2       npp

    """
    pval = []
    overlap = []
    total = []
    hit = []
    ap = pval.append
    ap2 = overlap.append
    ap3 = total.append
    ap4 = hit.append
    np1 = len(obj)
    np2 = len(whole) - np1
    keys = list(ref.keys())
    values = list(ref.values())
    for m in values:
        inter = obj & m
        n11 = len(inter)
        n12 = len(m) - n11
        x = np.array([[n11,n12],[np1 - n11,np2 - n12]])
        ap(stats.fisher_exact(x,alternative=mode)[1])
        ap2(inter)
        ap3(len(m))
        ap4(n11)
    if len(pval)==0:
        res = pd.DataFrame(columns=["p value","adjusted p value","overlap","hit No.",
                                    "total No."])
    else:
        res = pd.DataFrame({"p value":pval,"overlap":overlap,
                            "hit No.":hit,"total No.":total},index=keys)
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