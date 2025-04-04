# -*- coding: utf-8 -*-
"""
Created on Sat May 30 13:05:34 2020

Data Control class

Analyzer o-- DataControl

@author: tadahaya
"""
import pandas as pd
import numpy as np

from .data import *

# abstract factory
class DataControl():
    def __init__(self):
        self.obj = Data()
        self.ref = Data()
        self.whole = set()

    def set_whole(self,whole):
        """ set whole features """
        self.whole = whole
        self.obj.set_whole(self.whole)
        self.ref.set_whole(self.whole)

    def set_obj(self,data):
        """ create object Data instance for analysis """
        self.obj.set_data(data)

    def set_ref(self,data):
        """ create reference Data instance for analysis """
        self.ref.set_data(data)

    def get_obj(self):
        """ get object data instance """
        return self.obj.get_data()

    def get_ref(self):
        """ get reference data instance """
        return self.ref.get_data()

    def get_whole(self):
        """ get whole """
        return self.whole

    def adjust_obj(self):
        """ adjust object """
        self.obj.adjust()

    def adjust_ref(self):
        """ adjust reference """
        self.ref.adjust()


# concrete class
class FETDataControl(DataControl):
    def __init__(self):
        super().__init__()
        self.obj = SeqData()
        self.ref = SetData()


# concrete class
class BTDataControl(DataControl):
    def __init__(self):
        super().__init__()
        self.obj = SeqData()
        self.ref = SetData()


# concrete class
class GSEADataControl(DataControl):
    def __init__(self):
        super().__init__()
        self.obj = VectorData()
        self.ref = SetData()


# concrete class
class ssGSEADataControl(DataControl):
    def __init__(self):
        super().__init__()
        self.obj = VectorData()
        self.ref = SetData()


# concrete class
class ConnectivityDataControl(DataControl):
    def __init__(self):
        super().__init__()
        self.obj = VectorData()
        self.ref = SetTSData()