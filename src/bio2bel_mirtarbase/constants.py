# -*- coding: utf-8 -*-

import os

#: Data source
DATA_URL = 'http://mirtarbase.mbc.nctu.edu.tw/cache/download/6.1/miRTarBase_MTI.xlsx'

BIO2BEL_DIR = os.environ.get('BIO2BEL_DIRECTORY', os.path.join(os.path.expanduser('~'), '.pybel', 'bio2bel'))
DATA_DIR = os.path.join(BIO2BEL_DIR, 'mirtarbase')
os.makedirs(DATA_DIR, exist_ok=True)

DEFAULT_CACHE_NAME = 'mirtarbase.db'
DEFAULT_CACHE_PATH = os.path.join(DATA_DIR, DEFAULT_CACHE_NAME)
DEFAULT_CACHE_CONNECTION = os.environ.get('BIO2BEL_MIRTARBASE_DB', 'sqlite:///' + DEFAULT_CACHE_PATH)

MIRTARBASE_CONFIG_FILE_PATH = os.path.join(DATA_DIR, 'config.ini')
