# -*- coding: utf-8 -*-

"""Constants for Bio2BEL miRTarBase."""

import os

from bio2bel.utils import get_data_dir

VERSION = '0.2.0'

MODULE_NAME = 'mirtarbase'
DATA_DIR = get_data_dir(MODULE_NAME)

#: Data source
DATA_URL = 'http://mirtarbase.mbc.nctu.edu.tw/cache/download/6.1/miRTarBase_MTI.xlsx'
DATA_FILE_PATH = os.path.join(DATA_DIR, 'miRTarBase_MTI.xlsx')
