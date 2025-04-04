# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:05:34 2019

Connect class

@author: tadahaya
"""
import pandas as pd
import numpy as np
from itertools import chain
import copy
import random
import string

from .process.processor import Processor
from .analyzer import Analyzer
from .data.data_control import ConnectivityDataControl
from .calculator._connectivity import Calculator
from .plot._plot import PlotGSEA

# concrete class
class Connect(Analyzer):
    def __init__(self):
        self.data = ConnectivityDataControl()
        self.__process = Processor()
        self.__calc = Calculator()
        self.__plot = PlotGSEA()
        self.__whole = set()
        self.__ref = dict()
        self.__obj = set()
        self.res = pd.DataFrame()


    ### data processing ###
    def check_ref(self,keyword:str):
        """ check contents of reference data """
        if len(self.__ref)==0:
            raise ValueError("!! fit() before this process !!")
        try:
            temp = self.__ref[keyword]
            print("up: {}".format(temp[0]))
            print("")
            print("down: {}".format(temp[0]))
        except KeyError:
            print("!! Wrong keyword !!")
            hit = {v for v in self.__ref.keys() if keyword in v}
            print("perhaps: {}".format(hit))

    def vector2set(self,data,fold:float=3.0,
                   nmin:int=None,nmax:int=None,**kwargs):
        """
        convert dataframe to the outlier set for reference data

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
                                      two_sided=True,**kwargs)


    ### data control ###
    def fit(self,data:dict,keep_whole:bool=False,nmin=None):
        """
        set a reference data instance
        
        Parameters
        ----------
        data: dict of up-/down-tags
            keys: tag name
            values: tuple of up-/down-gene set
            {tag_name:(up-tag set,down-tag set)}

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
            temp = [v[0]|v[1] for v in self.__ref.values()]
            self.__whole = set(chain.from_iterable(temp))
            self.data.set_whole(self.__whole)
        if nmin is not None:
            temp = self.__ref.copy()
            for k,v in self.__ref.items():
                if (len(v[0]) < nmin) or (len(v[1]) < nmin):
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
        return self.data.get_ref()

    def get_obj(self):
        """ get an object data instance """
        return self.data.get_obj()

    def get_whole(self):
        return self.__whole


    ### calculation ###
    def calc(self,data): # realization
        """
        conduct connectivity analysis

        Parameters
        ----------
        data: dataframe
            feature x sample dataframe

        Returns res
        -------
        res: df
            gene set enrichment score

        """
        self.data.set_obj(data)
        self.__obj = self.data.get_obj()
        temp = self.__obj.copy()
        ref_ins = copy.deepcopy(self.data.ref)
        ref_ins.set_whole(set(temp.index))
        ref_ins.adjust()
        self.res = self.__calc.calc(obj=temp,ref=ref_ins.get_data())
        del ref_ins
        return self.res


    ### visualization ###
    def set_res(self,data):
        """ set a result """
        self.res = data


    def plot(self,sample_name:str=None,highlight:list=[],ylabel:str="connectivity score",**kwargs): # realization
        """
        visualize a result of connectivity score

        Parameters
        ----------
        sample_name: str
            indicate the sample name to be visualized

        highlight: list
            indicate the plots to be highlightened

        fileout: str
            indicate the path for the output image

        dpi: int
            indicate dpi of the output image
            
        ylabel: str
            indicate the name of y axis

        title: str
            indicate the title of the plot

        color: str
            indicate the color of the bars

        fontsize: float
            indicate the fontsize in the plot

        size: float
            indicate the size of the plot

        figsize: tuple
            indicate the size of the plot

        """
        if sample_name is None:
            col = list(self.res.columns)
            for v in col:
                self.__plot.plot(data=self.res[v],highlight=highlight,ylabel=ylabel,**kwargs)
        else:
            focused = self.res[sample_name]
            self.__plot.plot(data=focused,highlight=highlight,ylabel=ylabel,**kwargs)


    ### other ###
    def generate_test_data(self):
        """ generate data for test """    
        # prepare reference
        key = [''.join(random.choices(string.ascii_letters, k=5)) for i in range(10)]
        ref = pd.DataFrame([np.random.randn(10000) for i in range(10)],index=key).T
        temp = ref.iloc[:,0]
        ref = self.vector2set(ref)

        # prepare object
        col = ["xxxx","yyyy"]
        obj = pd.DataFrame([np.random.randn(10000) for i in range(2)],index=col).T
        obj["zzzz"] = temp
        print("generate test data as follows:")
        print('    ref: dict, {"XXXX":({"aa","bb"},{"cc","dd","ee"}),...}')
        print('    obj: dataframe, feature x sample')
        return ref,obj