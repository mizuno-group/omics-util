# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:05:34 2019

ssGSEA class

@author: tadahaya
"""
import pandas as pd
import numpy as np
from itertools import chain
import random
import string
from tqdm import tqdm

from .process.processor import Processor
from .analyzer import Analyzer
from .data.data_control import ssGSEADataControl
from .calculator._gsea import Calculator
from .plot._plot import PlotSsGSEA

# concrete class
class ssGSEA(Analyzer):
    def __init__(self):
        self.data = ssGSEADataControl()
        self.__process = Processor()
        self.__calc = Calculator()
        self.__plot = PlotSsGSEA()
        self.__whole = set()
        self.__ref = dict()
        self.__obj = pd.DataFrame()
        self.__method = ""
        self.alpha = 0.0
        self.__fterm = None
        self.__mode = None
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
        convert dataframe to the outlier set for reference

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
        return self.data.get_ref()

    def get_obj(self):
        """ get an object data instance """
        return self.data.get_obj()

    def get_whole(self):
        return self.__whole


    def get_calculator(self):
        """ get calculater object """
        return self.__calc


    ### calculation ###
    def calc(self,data,fterm=None,method:str="standard",alpha:float=0.25): # realization
        """
        conduct ssGSEA

        Parameters
        -------
        data: dataframe
            feature x sample dataframe

        fterm: str or int
            indicate the term of interest or the corresponding No.

        method: str
            indicate a method for calculating the enrichment score
            "starndard": employed in the original paper Barbie, et al, 2009
            "kuiper": Kuiper test statistics, good when up/down genes are mixed, tail sensitive
            "gsva": GSVA like statistics, good when unidirection (ex. up only)

        alpha: float, (0,1]
            indicate weight of center
            0.25 is employed in the original paper (Barbie,et al.,2009)

        Returns res
        -------
        res: df
            gene set enrichment score

        """
        self.__method = method
        self.__alpha = alpha
        self.__fterm = fterm
        if method=="standard":
            self.__calc.to_standard()
            print("Standard method",flush=True)
        elif method=="kuiper":
            self.__calc.to_kuiper()
            print("Kuiper method",flush=True)
        elif method=="gsva":
            self.__calc.to_gsva()
            print("GSVA method",flush=True)
        else:
            raise ValueError("!! Wrong method: choose 'standard', 'kuiper', or 'gsva' !!")
        self.data.set_obj(data)
        self.__obj = self.data.get_obj()
        temp = self.__obj.copy()
        col = list(temp.columns)
        if fterm is None:
            self.__mode = "exploratory"
            self.__calc.to_expssgsea()
            res = []
            ap = res.append
            for v in tqdm(col):
                ap(self.__calc.calc(obj=temp[v],ref=self.__ref,alpha=alpha))
            res = pd.concat(res,axis=1,join="inner")
            res.columns = col
        else:
            self.__mode = "focused"
            self.__calc.to_ssgsea()
            self.__fterm = fterm
            res = self.__calc.calc(obj=temp,ref=self.__ref,alpha=alpha,fterm=fterm)
        self.res = res
        return res


    def normalize_score(self, data=None):
        """
        normalize enrichment score with maximum value
        Note that this method can be applied
        only when the standard method was employed for calculation
        
        """
        if self.__mode=="focused":
            raise ValueError("!! Only applicable to the exploratory mode !!")
        if data is None:
            n = self.data.obj.data.shape[0] # 230317, n is obtained from the previous data
        else:
            n = data.shape[0] # 230317, n is obtained from the indicated data            
        if self.res.empty:
            raise ValueError("!! No result stored: calc() before this method !!")
        key = list(self.__ref.keys())
        val = list(self.__ref.values())
        dic = dict(zip(key,[len(v) for v in val]))
        idx = list(self.res.index)
        t = np.array([dic[v] for v in idx])
        if type(self.res)==type(pd.DataFrame()):
            kmax = np.c_[1 - t/n]
            val = self.res.values/kmax
        else:
            kmax = 1 - t/n
            val = self.res.values.flatten()/kmax
        return pd.DataFrame(val,columns=self.res.columns,index=idx)


    ### visualization ###
    def set_res(self,data):
        """ set a result """
        self.res = data


    def plot(self,keyword:list=[],fterm:str=None,mode:str=None,**kwargs): # realization
        """
        visualize a result of enrichment analysis

        Parameters
        ----------
        fterm: str
            indicate the group of interest

        mode: str
            indicate ssGSEA mode: 'exploratory' or 'focused'
            if None (default), plot the current result

        keyword: list
            indicate samples to be visualized

        fileout: str
            indicate the path for the output image

        dpi: int
            indicate dpi of the output image
            
        xlabel,ylabel: str
            indicate the name of x and y axes

        title: str
            indicate the title of the plot

        color: str
            indicate the color of swarmplot (if sample size is less than 30)

        palette: list
            indicate color of boxplot

        alpha: float
            indicate transparency of the bars: (0,1)

        size: float
            size of the markers

        fontsize: float
            indicate the fontsize in the plot

        textsize: float
            indicate the fontsize of the texts in the bars

        figsize: tuple
            indicate the size of the plot

        """
        if mode is None: # plot the current result
            if self.__mode=="exploratory":
                # exploratory mode
                if fterm is None:
                    raise ValueError("!! Indicate fterm (focused term) !!")
                else:
                    data = self.res.T[fterm]
                    self.__plot.plot(data=data,keyword=keyword,focus=fterm,**kwargs)
            elif self.__mode=="focused":
                # focused mode
                self.__plot.plot(data=self.res,keyword=keyword,focus=self.__fterm,**kwargs)
            else:
                raise ValueError("!! Indicate mode: 'exploratory' or 'focused' !!")
        elif mode=="exploratory": # plot a loaded result
            # exploratory mode
            if fterm is None:
                raise ValueError("!! Indicate fterm (focused term) !!")
            else:
                data = self.res.T[fterm]
                self.__plot.plot(data=data,keyword=keyword,focus=fterm,**kwargs)
        elif mode=="focused": # plot a loaded result
            # focused mode
            if fterm is None:
                raise ValueError("!! Indicate fterm (focused term) !!")
            else:
                self.__plot.plot(data=self.res,keyword=keyword,focus=fterm,**kwargs)
        else:
            raise ValueError("!! Wrong mode: 'exploratory' or 'focused' !!")


    ### other ###
    def generate_test_data(self):
        """ generate data for test """    
        # prepare reference
        key = [''.join(random.choices(string.ascii_letters, k=5)) for i in range(10)]

        def fxn():
            return {random.randint(0,50) for j in range(random.randint(10,100))}

        val = [fxn() for i in range(9)] + [{1,2,3,4,5,6,7}]
        ref = dict(zip(key,val))

        # prepare object
        col = ["control_{}".format(i + 1) for i in range(5)]
        col += ["treated_{}".format(i + 1) for i in range(5)]
        obj = pd.DataFrame([np.random.randn(10000) for i in range(10)],index=col).T
        print("generate test data as follows:")
        print('    ref: dict, {"XXXX":{"aa","bb"},"YYYY":{"cc","dd","ee"},...}')
        print('    obj: dataframe, feature x sample')
        return ref,obj