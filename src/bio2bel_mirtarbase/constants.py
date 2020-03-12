# -*- coding: utf-8 -*-

"""Constants for Bio2BEL miRTarBase."""

import os

from bio2bel.utils import get_data_dir

VERSION = '0.3.0-dev'

MODULE_NAME = 'mirtarbase'
DATA_DIR = get_data_dir(MODULE_NAME)

MIRTARBASE_VERSION = '7.0'

#: Data source
MIRTARBASE_DATA_URL = f'http://mirtarbase.mbc.nctu.edu.tw/cache/download/{MIRTARBASE_VERSION}/miRTarBase_MTI.xlsx'
VERSIONED_DATA_FOLDER = os.path.join(DATA_DIR, MIRTARBASE_VERSION)
os.makedirs(VERSIONED_DATA_FOLDER, exist_ok=True)
MIRTARBASE_RAW_DATA_PATH = os.path.join(VERSIONED_DATA_FOLDER, 'miRTarBase_MTI.xlsx')
MIRTARBASE_PROCESSED_DATA_PATH = os.path.join(VERSIONED_DATA_FOLDER, 'mirtarbase.tsv')
