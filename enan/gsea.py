# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:05:34 2019

GSEA class

@author: tadahaya
"""
import pandas as pd
import numpy as np
from itertools import chain
import random
import string

from .process.processor import Processor
from .analyzer import Analyzer
from .data.data_control import GSEADataControl
from .calculator._gsea import Calculator
from .plot._plot import PlotGSEA

# concrete class
class GSEA(Analyzer):
    def __init__(self):
        self.data = GSEADataControl()
        self.__process = Processor()
        self.__calc = Calculator()
        self.__plot = PlotGSEA()
        self.__whole = set()
        self.__obj = pd.DataFrame()
        self.__ref = dict()
        self.__method = ""
        self.alpha = 0.0
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
    def calc(self,data,method:str="standard",alpha:float=0.0): # realization
        """
        conduct GSEA

        Parameters
        -------
        data: dataframe
            feature x sample dataframe

        method: str
            indicate a method for calculating the enrichment score
            "starndard": employed in the original paper Barbie, et al, 2009
            "kuiper": Kuiper test statistics, good when up/down genes are mixed, tail sensitive
            "gsva": GSVA like statistics, good when unidirection (ex. up only)

        alpha: float, (0,1]
            indicate weight of center
            0 means no weight and is employed well

        Returns res
        -------
        res: df
            gene set enrichment score

        """
        self.__method = method
        self.__alpha = alpha
        if method=="standard":
            self.__calc.to_standard()
            print("Standard method")
        elif method=="kuiper":
            self.__calc.to_kuiper()
            print("Kuiper method")
        elif method=="gsva":
            self.__calc.to_gsva()
            print("GSVA method")
        else:
            raise ValueError("!! Wrong method: choose 'standard', 'kuiper', or 'gsva' !!")
        self.data.set_obj(data)
        self.__obj = self.data.get_obj()
        temp = self.__obj.copy()
        col = list(temp.columns)
        res = []
        ap = res.append
        for v in col:
            ap(self.__calc.calc(obj=temp[v],ref=self.__ref,alpha=alpha))
        res = pd.concat(res,axis=1,join="inner")
        res.columns = col
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


    def plot(self,sample_name:str=None,highlight:list=[],ylabel:str="enrichment score",**kwargs): # realization
        """
        visualize a result of GSEA
        Parameters
        ----------
        sample_name: str
            indicate the sample name to be visualized

        highlight: list
            indicate the reference keys to be highlightened in the plot

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
        

    def plot_running(
            self,sample_name:str=None,fterm:str="",title:str="",
            barcode_params:dict=None,
            heatmap:bool=True,heatmap_params:dict=None,**kwargs
            ): # realization
        """
        visualize a result of GSEA, running sum plot

        Parameters
        ----------
        sample_name: str
            indicate the sample name to be visualized

        fterm: str or int
            indicate the term of interest or the corresponding No.

        fileout: str
            indicate the path for the output image

        dpi: int
            indicate dpi of the output image
            
        xlabel,ylabel: str
            indicate the name of x and y axes

        title: str
            indicate the title of the plot

        color: str
            indicate the color of the bars

        fontsize: float
            indicate the fontsize in the plot

        figsize: tuple
            indicate the size of the plot

        barcode_params: dict
            indicate the parameters for barcode
            color: str
            size: float

        heatmap_params: dict
            indicate the parameters for heatmap
            cmap: str
            size: float
            alpha: float

        """
        if sample_name is None:
            raise ValueError("!! Indicate sample_name !!")
        data = self.__obj.copy()
        focused = data[sample_name]
        res = self.__calc.calc(obj=focused,ref=self.__ref,alpha=self.__alpha)
        es = self.__calc.es
        keys = self.__calc.keys
        if len(es)==0:
            raise ValueError("!! No Enrichment score !!")
        else:
            if type(fterm)==str:
                focus_key = fterm
                try:
                    focus_num = keys.index(fterm)
                except KeyError:
                    print("!! Wrong key: change fterm !!")
            elif type(fterm)==int:
                if len(es) < fterm:
                    raise ValueError("!! focused number is larger than column No. !!")
                focus_num = fterm
                focus_key = keys[fterm]
            else:
                raise TypeError("!! Wrong type: focus should be str or int !!")
        loc = self.__calc.loc
        if len(title) > 0:
            if heatmap:
                self.__plot.plot_running_heatmap(
                    data=focused.values,es=es,loc=loc,focus=focus_num,title=title,
                    barcode_params=barcode_params,heatmap_params=heatmap_params,**kwargs
                    )
            else:
                self.__plot.plot_running(
                    es=es,loc=loc,focus=focus_num,barcode_params=barcode_params**kwargs
                    )
        else:
            if heatmap:
                self.__plot.plot_running_heatmap(
                    data=focused.values,es=es,loc=loc,focus=focus_num,title=title,
                    barcode_params=barcode_params,heatmap_params=heatmap_params,**kwargs
                    )
            else:
                self.__plot.plot_running(
                    es=es,loc=loc,focus=focus_num,barcode_params=barcode_params,**kwargs
                    )


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
        obj = [np.random.randn(10000) for i in range(3)]
        obj = pd.DataFrame({"xxxx":obj[0],"yyyy":obj[1],"zzzz":obj[2]})
        print("generate test data as follows:")
        print('    ref: dict, {"XXXX":{"aa","bb"},"YYYY":{"cc","dd","ee"},...}')
        print('    obj: dataframe, feature x sample')
        return ref,obj
