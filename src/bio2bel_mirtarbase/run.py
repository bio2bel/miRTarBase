# -*- coding: utf-8 -*-
"""depricated. first script worte by Charlie to explain some stuff"""
import pandas as pd

DATA_URL = 'http://mirtarbase.mbc.nctu.edu.tw/cache/download/6.1/miRTarBase_MTI.xlsx'


def get_data():
    """Gets miRTarBase Interactions table

    :rtype: pandas.DataFrame
    """
    return pd.read_excel(DATA_URL)
