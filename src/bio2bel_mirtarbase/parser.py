# -*- coding: utf-8 -*-

"""Parsers for miRTarBase resources."""

import pandas as pd
from bio2bel.downloading import make_downloader

from pyobo import get_id_name_mapping, get_name_id_mapping, get_xrefs
from pyobo.cache_utils import cached_df
from pyobo.sources.mirbase import get_mature_id_to_name, get_mature_to_premature
from .constants import MIRTARBASE_DATA_URL, MIRTARBASE_PROCESSED_DATA_PATH, MIRTARBASE_RAW_DATA_PATH

__all__ = [
    'mirtarbase_data_downloader',
    'get_mirtarbase_df',
    'get_preprocessed_df',
]

mirtarbase_data_downloader = make_downloader(MIRTARBASE_DATA_URL, MIRTARBASE_RAW_DATA_PATH)

COLUMNS = [
    'mirtarbase_id',
    'mirna_name',
    'mirna_species_name',
    'target_name',
    'target_entrez_id',
    'target_species_name',
    'experiments',
    'support_type',
    'pubmed_id',
]


def _reverse_dict(d):
    return {v: k for k, v in d.items()}


@cached_df(path=MIRTARBASE_PROCESSED_DATA_PATH, dtype={
    'mirna_species_taxonomy_id': str,
    'target_species_taxonomy_id': str,
    'target_entrez_id': str,
    'target_hgnc_id': str,
    'pubmed_id': str,
})
def get_preprocessed_df() -> pd.DataFrame:
    """Get a preprocessed dataframe."""
    df = get_mirtarbase_df()

    hgnc_id_to_entrez_id = get_xrefs('hgnc', 'entrez')
    entrez_id_to_hgnc_id = _reverse_dict(hgnc_id_to_entrez_id)
    mirbase_id_to_name = get_id_name_mapping('mirbase')
    mirbase_mature_id_to_name = get_mature_id_to_name()
    mirbase_mature_name_to_id = _reverse_dict(mirbase_mature_id_to_name)
    mirbase_mature_id_to_premature_id = get_mature_to_premature()
    taxonomy_name_to_id = get_name_id_mapping('ncbitaxon')

    df['target_hgnc_id'] = df['target_entrez_id'].map(entrez_id_to_hgnc_id.get)

    df['mirbase_mature_id'] = df['mirna_name'].map(mirbase_mature_name_to_id.get)
    df['mirbase_id'] = df['mirbase_mature_id'].map(mirbase_mature_id_to_premature_id.get)
    df['mirbase_name'] = df['mirbase_id'].map(mirbase_id_to_name.get)

    df['mirna_species_taxonomy_id'] = df['mirna_species_name'].map(taxonomy_name_to_id.get)
    df['target_species_taxonomy_id'] = df['target_species_name'].map(taxonomy_name_to_id.get)

    return df


def get_mirtarbase_df() -> pd.DataFrame:
    """Get miRTarBase Interactions table and exclude rows with NULL values."""
    path = mirtarbase_data_downloader()

    df = pd.read_excel(path, names=COLUMNS, dtype={
        'Target Gene (Entrez Gene ID)': str,
        'References (PMID)': str,
    })
    # find null rows
    # null_rows = pd.isnull(df).any(1).nonzero()[0]
    # df= df.drop(null_rows)
    return df
