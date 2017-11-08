# -*- coding: utf-8 -*-

import os

#: Data source
DATA_URL = 'http://mirtarbase.mbc.nctu.edu.tw/cache/download/6.1/miRTarBase_MTI.xlsx'

MIRTARBASE_DATA_DIR = os.path.join(os.path.expanduser('~'), '.pybel', 'bio2bel', 'mirtarbase')
os.makedirs(MIRTARBASE_DATA_DIR, exist_ok=True)

MIRTARBASE_DATABASE_NAME = 'mirtarbase.db'
MIRTARBASE_DATABASE_PATH = os.path.join(MIRTARBASE_DATA_DIR, MIRTARBASE_DATABASE_NAME)
MIRTARBASE_SQLITE_PATH = 'sqlite:///' + MIRTARBASE_DATABASE_PATH

MIRTARBASE_CONFIG_FILE_PATH = os.path.join(MIRTARBASE_DATA_DIR, 'config.ini')
