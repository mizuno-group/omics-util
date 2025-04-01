# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 15:59:27 2020

algorithm of connectivity score calculation

@author: tadahaya
"""
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import rankdata
from numpy import random as rnd

class Calculator():
    """ calculate connectivity score """
    def __init__(self):
        self.vj = []
        self.res = pd.DataFrame()


    def calc(self,obj,ref):
        """
        calculate connectivity scores

        Parameters
        ----------
        object : dataframe
            feature x sample dataframe
        
        ref: dict of up-/down-tags
            keys: tag name
            values: tuple of up-/down-gene set
            {tag_name:(up-tag set,down-tag set)}
    
        Returns
        -------
        connectivity score: float

        """
        self.vj = generate_v(obj)
        val = list(ref.values())
        score = []
        ap = score.append
        for t in val:
            temp = []
            ap2 = temp.append
            for v in self.vj:
                ap2(calc_kss(t,v))
            ap(temp)
        self.res = pd.DataFrame(score,columns=list(obj.columns),index=list(ref.keys()))
        return self.res


    def get_details(self):
        """
        obtain detailed data for visualization

        Returns (vj,res)
        ----------
        vj: dict
            dictionary indicating the locations of tag-positive features

        """
        return self.vj,self.res


############### functions #####################################################
def generate_v(target):
    """
    generate a lisf of dictionaries corresponding to V(j) in CMap; Lamb J, Science, 2006

    Parameters
    ----------
    target: dataframe
        target response data
        features x samples

    """
    target = pd.DataFrame(target.sort_index())
    ind = list(target.index)
    if len(target.columns) < 2:
        rank = rankdata(-target.T.values).tolist()
        dic = dict(zip(ind,rank))
        return [dic] # corresponding to V(j) in CMap
    else:
        ranks = [rankdata(-v).tolist() for v in target.T.values]
        dics = [dict(zip(ind,v)) for v in ranks] # list of Vj
        return dics    


def _kss_max(n,tu,td):
    """ to normalize the difference of KS max dependent on t """
    return 2 + 1/n - (tu + td)/n


def _ab(tag,vj):
    """
    calculate a and b in CMap; Lamb J, Science, 2006

    Parameters
    ----------
    tag: list
        a feature list for tag

    vj: dict
        dictionary corresponding to V(j) in CMap; Lamb J, Science, 2006
        key, gene; val, rank

    """
    t = len(tag)
    n = len(vj)
    rank = [vj[v] for v in tag]
    rank.sort()
    lst_a = []
    ap_a = lst_a.append
    lst_b = []
    ap_b = lst_b.append
    for j in range(t):
        a = (j + 1)/t - rank[j]/n
        b = rank[j]/n - j/t 
        ap_a(a)
        ap_b(b)
    a_max = np.max(lst_a)
    b_max = np.max(lst_b)
    return (a_max,b_max)


def calc_kss(tag,vj):
    """
    calculate Kolmogorov-Smirnov statistics as in CMap; Lamb J, Science, 2006

    Parameters
    ----------
    tag: tuple
        tuple of up-/down-gene lists; (up,down)
        sorted with the values in the descending order

    vj: dict
        dictionary corresponding to V(j) in CMap; Lamb J, Science, 2006
        key, gene; val, rank

    """
    a_up,b_up = _ab(tag[0],vj)
    a_dn,b_dn = _ab(tag[1],vj)
    if a_up > b_up:
        ks_up = a_up
    else:
        ks_up = -1*b_up
    if a_dn > b_dn:
        ks_dn = a_dn
    else:
        ks_dn = -1*b_dn
    if ks_up*ks_dn > 0:
        ks = 0
    else:
        ks = ks_up - ks_dn
    n = len(vj)
    tu = len(tag[0])
    td = len(tag[1])
    kssmax = _kss_max(n,tu,td)    
    return ks/kssmax            


# def calc_p(target,tag,kss,permutation=10000):
#     """
#     calculate p value of Kolmogorov-Smirnov statistics based on permutation
#     randomize target data to generate null distribution

#     Parameters
#     ----------
#     target: dataframe or series
#         target response data

#     tag: tuple
#         tuple of up-/down-gene lists; (up,down)
#         sorted with the values in the descending order
#         << Caution >> take a tag, not tags (because it is too time-consuming)

#     kss: default, None
#         KS score of real data

#     Return
#     ----------
#     float: permutation p value
#     dataframe: null distribution

#     """
#     target = pd.DataFrame(target.sort_index())
#     ind = list(target.index)
#     perm_set = int(permutation/1000)
#     if perm_set==0:
#         perm_set = 1
#     res = []
#     ap = res.append
#     print("permutation ({} times)".format(perm_set))
#     for i in range(perm_set):
#         print(i)
#         y = 2*stats.bernoulli.rvs(0.5,loc=0,size=len(ind)) - 1
#         rand = np.array([rnd.permutation(y) for i in range(1000)])
#         rand = pd.DataFrame(target.values.flatten()*rand,columns=target.index).T
#         vjs = generate_v(rand)
#         ap([calc_kss(tag,w) for w in vjs])
#     res = list(chain.from_iterable(res))
#     ks_p = np.count_nonzero(np.abs(np.array(res)) >= np.abs(kss))
#     if ks_p/permutation==1.0:
#         pval = 1.0
#     else:    
#         pval = 2*ks_p/permutation
#     return (pval,res)


# def calc_p_fast(tag,vj,kss,permutation=10000):
#     """
#     calculate p value of Kolmogorov-Smirnov statistics based on permutation
#     random sampling of up-/down-tags

#     Parameters
#     ----------
#     tag: tuple
#         tuple of up-/down-gene lists; (up,down)
#         sorted with the values in the descending order

#     vj: dict
#         dictionary corresponding to V(j) in CMap 2006,Science
#         key, gene; val, rank

#     Return
#     ----------
#     float: permutation p value
#     dataframe: null distribution

#     """
#     n_up = len(tag[0])
#     n_dn = len(tag[1])
#     keys = list(vj.keys())
#     res = []
#     ap = res.append
#     for i in range(permutation):
#         rnd.shuffle(keys)
#         tag_random = (keys[:n_up],keys[n_up:n_up + n_dn])
#         ap(calc_kss(tag_random,vj))
#     ks_p = np.count_nonzero(np.abs(np.array(res)) >= np.abs(kss))
#     pval = 2*ks_p/permutation
#     if pval > 1:
#         pval = 1.0
#     return (pval,res)