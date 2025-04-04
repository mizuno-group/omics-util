# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:05:34 2019

BT class

@author: tadahaya
"""
import pandas as pd
import numpy as np
from itertools import chain
import random
import string

from .process.processor import Processor
from .analyzer import Analyzer
from .data.data_control import BTDataControl
from .calculator._binom import Calculator
from .plot._plot import PlotFET

# concrete class
class BT(Analyzer):
    def __init__(self):
        self.data = BTDataControl()
        self.__process = Processor()
        self.__calc = Calculator()
        self.__plot = PlotFET()
        self.__whole = set()
        self.__ref = dict()
        self.res = pd.DataFrame()


    ### data processing ###
    def check_ref(self,keyword:str):
        """ check contents of reference data """
        if len(self.__ref)==0:
            raise ValueError("!! fit() before this process !!")
        try:
            temp = self.__ref[keyword]
            print("{0}: {1}".format(keyword,temp))
        except KeyError:
            print("!! Wrong keyword !!")
            hit = {v for v in self.__ref.keys() if keyword in v}
            print("perhaps: {}".format(hit))

    def vector2set(self,data,fold:float=3.0,
                   nmin:int=None,nmax:int=None,**kwargs):
        """
        convert dataframe to the outlier set for reference or object

        Parameters
        ----------
        data: dataframe
            feature x sample dataframe

        fold: float
            indicates fold change determining the outliers

        nmin,nmax: int
            indicate the minimum/maximum number of each set 

        """
        return self.__process.vec2set(mtx=data,fold=fold,nmin=nmin,nmax=nmax,
                                      two_sided=False,**kwargs)


    ### data control ###
    def fit(self,data:dict,keep_whole:bool=False,nmin=None):
        """
        set a reference data instance
        
        Parameters
        ----------
        data: dict
            a dictionary of sets like {"XXXX":{"aa","bb"},"YYYY":{"cc","dd","ee"},...}

        keep_whole: boolean
            whether whole features is conserved when already registered

        nmin: int
            indicates the number of features necessary for each set

        """
        self.data.set_ref(data=data)
        self.__ref = self.data.get_ref()
        if keep_whole:
            if len(self.__whole)==0:
                raise ValueError("!! set_whole() or turn off keep_whole !!")
            else:
                self.data.adjust_ref()
        else:
            self.__whole = set(chain.from_iterable(self.__ref.values()))
            self.data.set_whole(self.__whole)
        if nmin is not None:
            temp = self.__ref.copy()
            for k,v in self.__ref.items():
                if len(v) < nmin:
                    del temp[k]
            self.__ref = temp
            self.data.set_ref(data=self.__ref)


    def set_whole(self,whole:set):
        """
        set whole features
        
        Parameters
        ----------
        whole: set
            indicates whole features
        
        """
        self.data.set_whole(whole)
        self.__whole = self.data.get_whole()
        if len(self.data.get_ref())!=0:
            self.data.adjust_ref()

    def get_ref(self):
        """ get reference data instance """
        return self.__ref

    def get_whole(self):
        return self.__whole


    ### calculation ###
    def calc(self,data:set,**kwargs): # realization
        """
        conduct Binomial test and obtain p value corrected for multiple tests
        all elements should be given as ID, except for term
        
        Parameters
        ----------
        data: set
            features of interest to be anlyzed

        correction: str
            indicate method for correcting multiple tests
            depend on "statsmodels.stats.multitest.multipletests"
            
        focus: int
            export results by XX th lowest p value

        mode: str
            indicate the type of significant judging
            "greater", "two-sided", or "less"
    
        """
        self.data.set_obj(data)
        self.res = self.__calc.calc(obj=self.data.get_obj(),ref=self.__ref,whole=self.__whole,**kwargs)
        return self.res


    ### visualization ###
    def set_res(self,df):
        """ set a result """
        self.res = df


    def plot(self,display_num:int=5,**kwargs): # realization
        """
        visualize a result of enrichment analysis

        Parameters
        ----------
        display_num: int
            determine how many groups are visualized

        fileout: str
            indicate the path for the output image

        dpi: int
            indicate dpi of the output image
            
        thresh: float
            indicate the threshold value of adjusted p value to be colored

        xlabel,ylabel: str
            indicate the name of x and y axes

        title: str
            indicate the title of the plot

        color: str
            indicate the color of the bars

        alpha: float
            indicate transparency of the bars: (0,1)

        height: float
            indicate the height of the bars

        fontsize: float
            indicate the fontsize in the plot

        textsize: float
            indicate the fontsize of the texts in the bars

        figsize: tuple
            indicate the size of the plot

        """
        self.__plot.plot(res=self.res,focus=display_num,**kwargs)


    ### other ###
    def generate_test_data(self):
        """ generate data for test """    
        # prepare reference
        key = [''.join(random.choices(string.ascii_letters,k=5)) for i in range(10)]

        def fxn():
            return {random.randint(0,50) for j in range(random.randint(5,30))}
        
        val = [fxn() for i in range(9)] + [{1,2,3,4,5,6,7}]
        ref = dict(zip(key,val))

        # prepare object
        obj = {1,2,3,4,5,6,7,8,9,100,1000}
        print("generate test data as follows:")
        print('    ref: dict, {"XXXX":{"aa","bb"},"YYYY":{"cc","dd","ee"},...}')
        print('    obj: set, {"aa","bb","cc"}')
        return ref,obj
