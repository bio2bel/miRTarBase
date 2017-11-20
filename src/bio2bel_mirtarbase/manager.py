# -*- coding: utf-8 -*-

import configparser
import logging
import os
import time
from urllib.request import urlretrieve

import pandas as pd
import pyhgnc
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from tqdm import tqdm

from bio2bel_mirtarbase.constants import CONFIG_FILE_PATH, DATA_DIR, DATA_URL, DEFAULT_CACHE_CONNECTION
from bio2bel_mirtarbase.models import Base, Evidence, Interaction, Mirna, Species, Target

log = logging.getLogger(__name__)

DATA_FILE_PATH = os.path.join(DATA_DIR, 'miRTarBase_MTI.xlsx')


def download_data(force_download=False):
    if not os.path.exists(DATA_FILE_PATH) or force_download:
        urlretrieve(DATA_URL, DATA_FILE_PATH)


def get_data(source=None):
    """Gets miRTarBase Interactions table and exclude rows with NULL values

    :param str source: location that goes into :func:`pandas.read_excel`. Defaults to :data:`DATA_URL`.
    :rtype: pandas.DataFrame
    """
    df = pd.read_excel(source or DATA_URL)

    # find null rows
    null_rows = pd.isnull(df).any(1).nonzero()[0]
    return df.drop(null_rows)


def build_entrez_map(pyhgnc_connection=None):
    """

    :param Optional[str] pyhgnc_connection:
    :rtype: dict[str,pyhgnc.manager.models.HGNC]
    """
    log.info('getting entrez mapping')

    if isinstance(pyhgnc_connection, pyhgnc.QueryManager):
        pyhgnc_manager = pyhgnc_connection
    else:
        pyhgnc_manager = pyhgnc.QueryManager(connection=pyhgnc_connection)
        log.info('using PyHGNC connection: %s', pyhgnc_manager.connection)

    t = time.time()
    emap = {
        model.entrez: model
        for model in pyhgnc_manager.hgnc()
        if model.entrez
    }
    log.info('got entrez mapping in %.2f seconds', time.time() - t)
    return emap


class Manager(object):
    """Manages the mirTarBase database"""

    def __init__(self, connection=None):
        """
        :param str connection: The connection string
        """
        self.connection = self.get_connection(connection)
        self.engine = create_engine(self.connection)
        self.session_maker = sessionmaker(bind=self.engine, autoflush=False, expire_on_commit=False)
        self.session = self.session_maker()
        self.create_all()

    @staticmethod
    def get_connection(connection=None):
        """Return the SQLAlchemy connection string if it is set

        :param connection: get the SQLAlchemy connection string
        :rtype: str
        """
        if connection:
            return connection

        config = configparser.ConfigParser()

        cfp = CONFIG_FILE_PATH

        if os.path.exists(cfp):
            log.info('fetch database configuration from {}'.format(cfp))
            config.read(cfp)
            connection = config['database']['sqlalchemy_connection_string']
            log.info('load connection string from {}: {}'.format(cfp, connection))
            return connection

        with open(cfp, 'w') as config_file:
            config['database'] = {'sqlalchemy_connection_string': DEFAULT_CACHE_CONNECTION}
            config.write(config_file)
            log.info('create configuration file {}'.format(cfp))

        return DEFAULT_CACHE_CONNECTION

    def create_all(self, check_first=True):
        """Creates all tables"""
        Base.metadata.create_all(self.engine, checkfirst=check_first)

    def drop_all(self, check_first=True):
        """Drops all tables"""
        log.info('dropping tables')
        Base.metadata.drop_all(self.engine, checkfirst=check_first)

    @staticmethod
    def ensure(connection=None):
        """Checks and allows for a Manager to be passed to the function

        :type connection: None or str or Manager
        :rtype: Manager
        """
        if connection is None or isinstance(connection, str):
            return Manager(connection=connection)

        if isinstance(connection, Manager):
            return connection

        raise TypeError

    def populate(self, source=None, update_pyhgnc=False, pyhgnc_connection=None):
        """Populate database with the data from miRTarBase

        :param str source: path or link to data source needed for :func:`get_data`
        :param bool update_pyhgnc: Should HGNC be updated?
        :param Optional[str] pyhgnc_connection: Optional connection string for pyhgnc
        """
        if update_pyhgnc:
            pyhgnc.update(connection=pyhgnc_connection)

        t = time.time()
        log.info('getting data')
        df = get_data(source)
        log.info('got data in %.2f seconds', time.time() - t)

        mirna_set = {}
        target_set = {}
        species_set = {}

        emap = build_entrez_map(pyhgnc_connection=pyhgnc_connection)

        log.info('building models')
        t = time.time()
        for (index, mirtarbase_id, mirtarbase_name, mirna_species, target, entrez, taret_species, exp, sup_type,
             pubmed) in tqdm(df.itertuples(), total=len(df.index)):
            # create new miRNA instance
            if mirtarbase_id not in mirna_set:

                species = species_set.get(mirna_species)

                if species is None:
                    species = species_set[mirna_species] = Species(name=mirna_species)
                    self.session.add(species)

                new_mirna = Mirna(
                    mirtarbase_id=mirtarbase_id,
                    mirtarbase_name=mirtarbase_name,
                    species=species
                )
                self.session.add(new_mirna)
                mirna_set[mirtarbase_id] = new_mirna

            entrez = str(int(entrez))

            # create new target instance
            if entrez not in target_set:
                species = species_set.get(taret_species)

                if species is None:
                    species = species_set[taret_species] = Species(name=taret_species)
                    self.session.add(species)

                new_target = Target(
                    entrez_id=entrez,
                    species=species,
                    target_gene=target,
                )

                if entrez in emap:
                    g_first = emap[entrez]
                    new_target.hgnc_symbol = g_first.symbol
                    new_target.hgnc_id = str(g_first.identifier)

                self.session.add(new_target)
                target_set[entrez] = new_target

            # create new evidence instance
            new_evidence = Evidence(experiment=exp, support=sup_type, reference=pubmed)
            self.session.add(new_evidence)

            # create new interaction instance
            new_interaction = Interaction(mirna=mirna_set[mirtarbase_id], target=target_set[entrez],
                                          evidence=new_evidence)
            self.session.add(new_interaction)

        log.info('built models in %.2f seconds', time.time() - t)

        log.info('committing models')
        t = time.time()

        try:
            self.session.commit()
            log.info('committed after %.2f seconds', time.time() - t)
        except:
            self.session.rollback()
            log.exception('commit failed after %.2f seconds', time.time() - t)

    def map_entrez_to_hgnc(self):
        """Function to map entrez identifiers to HGNC identifiers"""
        raise NotImplementedError

    def query_MTIs(self, query_mir):
        """Find all MTI's for a given miRTarBase identifier

        :param query_mir: miRTarBase identifier of interest
        :return targets: list of all targets of query_mir
        """

        raise NotImplementedError

    def query_mirna_by_mirtarbase_identifier(self, mirtarbase_id):
        """Gets an miRNA by its miRTarBase identifier

        :param str mirtarbase_id: An miRTarBase identifier
        :rtype: Optional[Mirna]
        """
        return self.session.query(Mirna).filter(Mirna.mirtarbase_id == mirtarbase_id).one_or_none()

    def query_targets(self, targets):
        """Find all targets

        :param list[str] targets: list of HGNC names
        """
        raise NotImplementedError

    def query_target_by_entrez_id(self, entrez_id):
        """Query for one target

        :param str entrez_id: Entrez gene identifier
        :rtype: Optional[Target]
        """
        return self.session.query(Target).filter(Target.entrez_id == entrez_id).one_or_none()

    def query_target_by_hgnc_symbol(self, hgnc_symbol):
        """Query for one target

        :param str hgnc_symbol: HGNC gene symbol
        :rtype: Optional[Target]
        """
        return self.session.query(Target).filter(Target.hgnc_symbol == hgnc_symbol).one_or_none()

    def query_target_by_hgnc_identifier(self, hgnc_id):
        """Query for one target

        :param str hgnc_id: HGNC gene identifier
        :rtype: Optional[Target]
        """
        return self.session.query(Target).filter(Target.hgnc_id == hgnc_id).one_or_none()


if __name__ == '__main__':
    logging.basicConfig(level=20)
    log.setLevel(20)
    m = Manager()
    m.drop_all()
    m.populate()
