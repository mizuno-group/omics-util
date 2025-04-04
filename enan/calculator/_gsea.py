# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 15:48:00 2020

Gene set enrichment analysis

@author: tadahaya
"""
import pandas as pd
import numpy as np

### preprocessing ###
class Process():
    """ prepare data to be analyzed """
    def __init__(self):
        self.__process = GSEAProcess()

    def do(self,data,ref,**kwargs):
        return self.__process.do(data,ref,**kwargs)

    def to_gsea(self):
        self.__process = GSEAProcess()

    def to_ssgsea(self):
        self.__process = ssGSEAProcess()

    def to_expssgsea(self):
        self.__process = ExpssGSEAProcess()


class GSEAProcess():
    def do(self,data,ref):
        """
        prepare target data for GSEA
        
        data: series
            a series of interest data

        """
        return data.sort_values(ascending=False),ref


class ssGSEAProcess():
    def do(self,data,ref,fterm):
        """
        prepare target data for ssGSEA focused mode
        
        data: dataframe
            a dataframe of interest data (feature x sample)
        
        fterm: str
            a term of interest in refrence data set
        
        """
        return _sort_data(data),ref[fterm]


class ExpssGSEAProcess():
    def do(self,data,ref):
        """
        prepare target data for ssGSEA exploratory mode
        
        data: series
            a series of interest data

        """
        return data.sort_values(ascending=False),ref


### scoreing method ###
class Method():
    """ indicate scoreing method of enrichment scores for output """
    def __init__(self):
        self.__method = StandardMethod()

    def calc(self,es,names):
        return self.__method.calc(es,names)

    def to_standard(self):
        self.__method = StandardMethod()

    def to_kuiper(self):
        self.__method = KuiperMethod()

    def to_gsva(self):
        self.__method = GSVAMethod()


class StandardMethod():
    def calc(self,es,names):
        """ calculate Enrichment score from es with the standard method """
        tmp_max = np.max(es,axis=0)
        tmp_min = np.min(es,axis=0)
        if np.abs(tmp_max) >np.abs(tmp_min):
            res = pd.DataFrame(tmp_max,index=names,
                                columns=["ES"]).sort_values(by="ES",ascending=False)
        else:
            res = pd.DataFrame(tmp_min,index=names,
                                columns=["ES"]).sort_values(by="ES",ascending=False)
        return res


class KuiperMethod():
    def calc(self,es,names):
        """ calculate Enrichment score from es with the method in Kuiper method """
        res = pd.DataFrame(np.abs(np.max(es,axis=0)) + np.abs(np.min(es,axis=0)),
                            index=names,columns=["ES"]).sort_values(by="ES",ascending=False)
        return res


class GSVAMethod():
    def calc(self,es,names):
        """ calculate Enrichment score from es with the method in GSVA """
        maxi = np.max(es,axis=0)
        mini = np.min(es,axis=0)
        res = pd.DataFrame(np.where(maxi < 0,0,maxi) - np.where(mini > 0,0,mini),
                            index=names,columns=["ES"]).sort_values(by="ES",ascending=False)
        return res


### Algorithm ###
class Algorithm():
    """ indicate analysis algorithm """
    def __init__(self,):
        self.__algorithm = GSEAAlgorithm()

    def calc(self,data,tag,alpha):
        """
        Returns: (es,keys,loc)
        ----------
        es: 2d array
            enrichment scores of each location (feature x tag/sample)

        keys: list
            tag names or sample names for output

        loc: 2d array
            location of tag-positive elements in data (tag/sample x feature)

        """
        return self.__algorithm.calc(data,tag,alpha)

    def to_gsea(self):
        self.__algorithm = GSEAAlgorithm()

    def to_ssgsea(self):
        self.__algorithm = ssGSEAAlgorithm()

    def to_expssgsea(self):
        self.__algorithm = ExpssGSEAAlgorithm()


class GSEAAlgorithm():
    def calc(self,data,tag,alpha=0):
        """
        standard GSEA algorithm
        
        Parameters
        ----------
        data: series
            a sorted series of interest data in descending order (high values get high ranks)
            
        tag: dict
            keys: group name
            values: members of each group

        """
        keys = list(tag.keys())
        values = list(tag.values())
        sorted_list = list(data.index)
        sorted_val = np.abs(list(data))
        n = len(sorted_list)
        ### grasp where the tagged genes
        loc = _location(sorted_list,values)
        ### calculate posi
        posi0 = np.power(sorted_val,alpha)*loc
        den_p = np.sum(posi0,axis=1)
        posi = _accumulative(posi0,axis=1)/den_p
        ### calculate nega
        nega0 = 1 - loc
        # n_nh = n - np.array([len(v) for v in values])
        n_nh = np.sum(nega0,axis=1) # 211014
        nega = _accumulative(nega0,axis=1)/n_nh
        ### calculate ES
        es = posi - nega
        return es,keys,loc


class ssGSEAAlgorithm():
    def calc(self,data,tag,alpha=0.25):
        """
        ssGSEA focusing on a tag
        
        Parameters
        ----------
        data: dataframe
            a dataframe composed of sorted features (descending; high values get high ranks)
            each column corresponds to each sample
            
        tag: set
            a set of elements constituting a signature        
        
        """
        keys = list(data.columns)
        n = data.shape[0]
        ### grasp where the tagged genes
        tag = list(tag)
        loc = data.replace(tag,1).T
        loc = (loc.values==1)
        ### calculate posi
        full = np.arange(1,n + 1,1)[::-1]
        posi0 = np.power(full,alpha)*loc
        den_p = np.sum(posi0,axis=1)
        posi = _accumulative(posi0,axis=1)/den_p
        ### calculate nega
        nega0 = 1 - loc
        n_nh = np.sum(nega0,axis=1)
        nega = _accumulative(nega0,axis=1)/n_nh
        ### calculate ES
        es = posi - nega
        return es,keys,loc


class ExpssGSEAAlgorithm():
    def calc(self,data,tag,alpha=0.25):
        """
        ssGSEA for exploratory analysis by searching whole tags
        
        Parameters
        ----------
        data: series
            a sorted series of interest data in descending order (high values get high ranks)
            
        tag: dict
            keys: group name
            values: members of each group
        
        """        
        keys = list(tag.keys())
        values = list(tag.values())
        data = list(data.index)
        n = len(data)
        ### grasp where the tagged genes
        loc = _location(data,values)
        ### calculate posi
        full = np.arange(1,n + 1,1)[::-1]
        posi0 = np.power(full,alpha)*loc
        den_p = np.sum(posi0,axis=1)
        posi = _accumulative(posi0,axis=1)/den_p
        ### calculate nega
        nega0 = 1 - loc
        n_nh = np.sum(nega0,axis=1)
        nega = _accumulative(nega0,axis=1)/n_nh
        ### calculate ES
        es = posi - nega
        return es,keys,loc


### Context ###
class Calculator():
    """ GSEA client """
    def __init__(self):
        self.__process = Process()
        self.__algorithm = Algorithm()
        self.__method = Method()
        self.es = np.array([[],[]]) # for plot
        self.loc = np.array([[],[]]) # for plot
        self.keys = []
        self.res = pd.DataFrame()


    def calc(self,obj,ref,alpha,**kwargs):
        """
        conduct GSEA

        Parameters
        ----------
        alpha: float, (0,1]
            indicate weight of center
            0.25 is derived from the original paper of ssGSEA (Barbie,et al.,2009)
        
        fterm: str
            a term of interest in refrence data set

        """
        obj,tag = self.__process.do(obj,ref,**kwargs)
        self.es,self.keys,self.loc = self.__algorithm.calc(obj,tag,alpha)
        self.res = self.__method.calc(self.es,self.keys)
        return self.res


    def get_details(self):
        """
        obtain detailed data for visualization

        Returns (es,loc,keys,res)
        ----------
        es: 2d array
            enrichment scores of each location (feature x tag/sample)

        loc: 2d array
            location of tag-positive elements in data (tag/sample x feature)

        keys: list
            tag names or sample names for output

        """
        return self.es,self.loc,self.keys,self.res


    def to_gsea(self):
        self.__process.to_gsea()
        self.__algorithm.to_gsea()

    def to_ssgsea(self):
        self.__process.to_ssgsea()
        self.__algorithm.to_ssgsea()

    def to_expssgsea(self):
        self.__process.to_expssgsea()
        self.__algorithm.to_expssgsea()

    def to_standard(self):
        self.__method.to_standard()

    def to_kuiper(self):
        self.__method.to_kuiper()

    def to_gsva(self):
        self.__method.to_gsva()


### Common Functions ###
def _location(data,tags):
    """
    determine the location of the interest elements in a sample
    returns tag No. x feature No. matrix
    
    Parameters
    ----------
    data: list
        sorted list by values (descending)

    tags: list
        a list of tag sets

    """
    res = []
    ap = res.append
    for tag in tags:
         ap([(ele in tag) for ele in data])
    return np.array(res)


def _accumulative(X,axis=0):
    """
    calculate accumulative sum
    
    Parameters
    ----------
    axis: 0 or 1
        determine sum direction
    
    """
    if axis==1:
        X = X.T
    x = np.zeros(X.shape[1])
    res = []
    ap = res.append
    for v in X:
        x += v
        ap(x.copy())
    res = np.array(res)
    return res


def _sort_data(df,ascending=False):
    """
    sort dataframe by values in columns and return corresponding indices
    
    Parameters
    ----------
    df: dataframe

    ascending: boolean
        indicates sort direction

    """
    col = list(df.columns)
    res = list(map(lambda x: list(df[x].sort_values(ascending=ascending).index),col))
    return pd.DataFrame(res,index=col).T
