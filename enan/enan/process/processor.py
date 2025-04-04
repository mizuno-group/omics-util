# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:05:34 2019

PreProcessor class

Analyzer o-- Processor

@author: tadahaya
"""
import pandas as pd
import numpy as np

class Processor():
    def __init__(self):
        pass

    def __vec2tpl(self,mtx,fold=2.0,method="iqr",nmin=None,nmax=None):
        """ convert dataframe into tuple of tags """
        sample_name = list(mtx.columns)
        n_sample = len(sample_name)
        n_feature = len(mtx.index)
        if nmin is None:
            nmin = 15
        if nmax is None:
            nmax = int(0.01*n_feature)
        if method=="std":
            upper,lower = outlier_std(mtx=mtx,fold=fold,axis=0)
        else:
            upper,lower = outlier_iqr(mtx=mtx,fold=fold,axis=0)
        res = []
        ap = res.append
        for i in range(n_sample):
            temp = mtx.iloc[:,i].sort_values(ascending=False)
            up_val = upper[i]
            low_val = lower[i]
            temp_l = list(temp.index)
            upper_tag = set(temp[temp > up_val].index)
            lower_tag = set(temp[temp < low_val].index)
            n_up = len(upper_tag)
            n_low = len(lower_tag)
            if n_up > nmax:
                upper_tag = set(temp_l[:nmax])
            elif n_up < nmin:
                upper_tag = set(temp_l[:nmin])
            if n_low > nmax:
                lower_tag = set(temp_l[-nmax:])
            elif n_low < nmin:
                lower_tag = set(temp_l[-nmin:])
            ap((upper_tag,lower_tag))
        return res


    def vec2set(self,mtx,fold=2.0,two_sided=True,method="iqr",nmin=None,nmax=None):
        """
        convert dataframe into dict of tags

        Parameters
        ----------
        mtx: dataframe
            feature x sample matrix

        fold: float
            determine threshold of outliers

        two_sided: boolean
            determine whether up/down is kept

        method: str
            "std" or "iqr"

        """
        temp = self.__vec2tpl(mtx,fold,method,nmin,nmax)
        if two_sided:
            dic = dict(zip(list(mtx.columns),temp))
        else:
            temp = [v[0] | v[1] for v in temp]
            dic = dict(zip(list(mtx.columns),temp))
        return dic


def outlier_std(mtx,fold,axis=0):
    """ calculate upper and lower values for outlier detection """
    loc = np.mean(mtx.values,axis=axis)
    scale = np.std(mtx.values,axis=axis)        
    upper = loc + fold*scale
    lower = loc - fold*scale
    return upper,lower    
    
    
def outlier_iqr(mtx,fold,axis=0):
    """ calculate upper and lower values for outlier detection """
    q3 = np.percentile(mtx.values,75,axis=axis)
    q1 = np.percentile(mtx.values,25,axis=axis) 
    scale = q3 - q1        
    upper = q3 + fold*scale
    lower = q1 - fold*scale
    return upper,lower