# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:05:34 2019

Adjuster class

Analyzer o-- DataControl o-- Data o-- Adjuster

@author: tadahaya
"""
import pandas as pd

__all__ = ["SeqAdjuster","SetAdjuster","SetTSAdjuster","VectorAdjuster"]

# delegation
class SeqAdjuster():
    def __init__(self):
        pass

    def adjust(self,data,whole,nmin=3):
        """
        adjust seq (set or list) according to whole

        Parameters
        ----------
        data: set or list
            a set or list indicating feature sets

        nmin: int
            indicate minimum number of each set

        """
        return set(data) & whole


class SetAdjuster():
    def __init__(self):
        pass

    def adjust(self,data,whole,nmin=3):
        """
        adjust set (dict) according to whole

        Parameters
        ----------
        data: dict
            dict indicating feature sets

        nmin: int
            indicate minimum number of each set

        """
        new_keys = []
        new_values = []
        ap = new_keys.append
        ap2 = new_values.append
        for k,v in data.items():
            new = v & whole
            if len(new) >= nmin:
                ap(k)
                ap2(new)
        return dict(zip(new_keys,new_values))



class SetTSAdjuster():
    def __init__(self):
        pass

    def adjust(self,data,whole,nmin=3):
        """
        adjust two-sided set (dict) according to whole

        Parameters
        ----------
        data: dict
            dict indicating feature sets

        nmin: int
            indicate minimum number of each set

        """
        new_keys = []
        new_values = []
        ap = new_keys.append
        ap2 = new_values.append
        for k,v in data.items():
            new1 = v[0] & whole
            new2 = v[1] & whole
            if (len(new1) >= nmin) and (len(new2) >= nmin):
                ap(k)
                ap2((new1,new2))
        return dict(zip(new_keys,new_values))


class VectorAdjuster():
    def __init__(self):
        pass

    def adjust(self,data,whole):
        """
        adjust vector (dataframe) according to whole

        Parameters
        ----------
        data: dict
            dict indicating feature sets

        """
        new = set(data.index) & whole
        return data.loc[new,:].sort_index()