# -*- coding: utf-8 -*-

"""Constants for Bio2BEL miRTarBase."""

import os

from bio2bel.utils import get_data_dir

VERSION = '0.3.0-dev'

MODULE_NAME = 'mirtarbase'
DATA_DIR = get_data_dir(MODULE_NAME)

MIRTARBASE_VERSION = '6.1'

#: Data source
DATA_URL = f'http://mirtarbase.mbc.nctu.edu.tw/cache/download/{MIRTARBASE_VERSION}/miRTarBase_MTI.xlsx'
DATA_FILE_PATH = os.path.join(DATA_DIR, 'miRTarBase_MTI.xlsx')
