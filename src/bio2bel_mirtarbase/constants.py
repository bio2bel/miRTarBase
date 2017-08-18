# -*- coding: utf-8 -*-

import os

#: Data source
DATA_URL = 'http://mirtarbase.mbc.nctu.edu.tw/cache/download/6.1/miRTarBase_MTI.xlsx'

MIRTARBASE_DATA_DIR = os.path.join(os.path.expanduser('~'), '.pymirtarbase')

if not os.path.exists(MIRTARBASE_DATA_DIR):
    os.makedirs(MIRTARBASE_DATA_DIR)

MIRTARBASE_DATABASE_NAME = 'miRTarBase.db'
MIRTARBASE_SQLITE_PATH = 'sqlite:///' + os.path.join(MIRTARBASE_DATA_DIR, MIRTARBASE_DATABASE_NAME)

MIRTARBASE_CONFIG_FILE_PATH = os.path.join(MIRTARBASE_DATA_DIR, 'config.ini')