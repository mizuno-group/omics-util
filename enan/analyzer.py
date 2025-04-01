# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:05:34 2019

Analyzer class

@author: tadahaya
"""
import pandas as pd
import numpy as np

from .data.data_control import DataControl
from .process.processor import Processor

# abstract class
class Analyzer():
    def __init__(self):
        self.data = DataControl()
        self.__process = Processor()
        self.__calc = None
        self.__plot = None
        self.__whole = set()
        self.__ref = None
        self.__obj = None
        self.res = None


    ### data processing ###
    def check_ref(self):
        """ check contents of reference data """
        raise NotImplementedError

    def vector2set(self):
        """ convert dataframe to the outlier set """
        raise NotImplementedError


    ### data control ###
    def fit(self,data=None):
        """ set a reference data instance """
        raise NotImplementedError

    def set_whole(self,whole=None):
        """ set whole features """
        raise NotImplementedError

    def get_ref(self):
        """ get reference data instance """
        raise NotImplementedError

    def get_whole(self):
        """ get whole set """
        raise NotImplementedError

    def get_res(self):
        """ get result """
        raise NotImplementedError


    ### calculator ###
    # abstract method
    def calc(self):
        """ conduct calculation """
        raise NotImplementedError


    ### visualization ###
    # abstract method
    def set_res(self):
        """ set a result """
        raise NotImplementedError

    def plot(self):
        """ visualize the result """
        raise NotImplementedError