# -*- coding: utf-8 -*-

"""Parsers for miRTarBase resources."""

import logging
import os
from typing import Optional
from urllib.request import urlretrieve

import pandas as pd

from .constants import DATA_FILE_PATH, DATA_URL

log = logging.getLogger(__name__)


def download_data(force_download: bool = False):
    """Download the miRTarBase Excel sheet to a local path.

    :param force_download: If true, don't download the file again if it already exists
    """
    if not os.path.exists(DATA_FILE_PATH) or force_download:
        log.info('downloading %s to %s', DATA_URL, DATA_FILE_PATH)
        urlretrieve(DATA_URL, DATA_FILE_PATH)
    else:
        log.info('using cached data at %s', DATA_FILE_PATH)

    return DATA_FILE_PATH


def get_data(url: Optional[str] = None, cache: bool = True, force_download: bool = False) -> pd.DataFrame:
    """Get miRTarBase Interactions table and exclude rows with NULL values.

    :param url: location that goes into :func:`pandas.read_excel`. Defaults to :data:`DATA_URL`.
    :param cache: If true, the data is downloaded to the file system, else it is loaded from the internet
    :param force_download: If true, overwrites a previously cached file
    """
    if url is None and cache:
        url = download_data(force_download=force_download)

    df = pd.read_excel(url or DATA_URL)

    # find null rows
    null_rows = pd.isnull(df).any(1).nonzero()[0]
    return df.drop(null_rows)
