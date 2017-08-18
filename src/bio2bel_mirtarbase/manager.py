# -*- coding: utf-8 -*-

import configparser
import logging
import os

import pandas as pd
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from .models import Base, Mirna, Target, Evidence, Interaction
from .constants import (
    DATA_URL,
    MIRTARBASE_SQLITE_PATH,
    MIRTARBASE_CONFIG_FILE_PATH,
)

log = logging.getLogger(__name__)

def get_data():
    """Gets miRTarBase Interactions table and exclude rows with NULL values

    :rtype: pandas.DataFrame
    """
    df = pd.read_excel("/home/colin/SCAI/test.xlsx")
    # find null rows
    null_rows = pd.isnull(df).any(1).nonzero()[0]
    return df.drop(null_rows)

class Manager(object):
    def __init__(self, connection=None):
        self.connection = self.get_connection(connection)
        self.engine = create_engine(self.connection)
        self.sessionmake = sessionmaker(bind=self.engine, autoflush=False, expire_on_commit=False)
        self.session = self.sessionmake()
        self.make_tables()

    @staticmethod
    def get_connection(connection=None):
        """Return the SQLAlchemy connection string if it is set
        :param connection: get the SQLAlchemy connection string
        :rtype: str
        """
        if connection:
            return connection

        config = configparser.ConfigParser()

        cfp = MIRTARBASE_CONFIG_FILE_PATH

        if os.path.exists(cfp):
            log.info('fetch database configuration from {}'.format(cfp))
            config.read(cfp)
            connection = config['database']['sqlalchemy_connection_string']
            log.info('load connection string from {}: {}'.format(cfp, connection))
            return connection

        with open(cfp, 'w') as config_file:
            config['database'] = {'sqlalchemy_connection_string': MIRTARBASE_SQLITE_PATH}
            config.write(config_file)
            log.info('create configuration file {}'.format(cfp))

        return MIRTARBASE_SQLITE_PATH

    def make_tables(self, check_first=True):
        """ create tables """
        Base.metadata.create_all(self.engine, checkfirst=check_first)

    def populate(self):
        """Populate database with the data from miRTarBase

        :param session: session object from sqlalchemy
        :return:
        """
        df = get_data()
        mirna_set = {}
        target_set = {}
        # iterate through rows and construct tables from it
        for index, mir_id, mirna, species_mirna, target, entrez, species_target, exp, sup_type, pubmed in df.itertuples():
            print(mir_id)
            # create new miRNA instance
            if mir_id not in mirna_set:
                new_mirna = Mirna(mirtarbase_id=mir_id, mir_name=mirna, species=species_mirna)
                mirna_set[mir_id] = new_mirna

            # create new target instance
            if entrez not in target_set:
                new_target = Target(target_gene=target, entrez_id=int(entrez), species=species_target)
                target_set[entrez] = new_target

            # create new evidence instance
            new_evidence = Evidence(experiment=exp, support=sup_type, reference=int(pubmed))

            # create new interaction instance
            new_interaction = Interaction(mirna=mirna_set[mir_id], target=target_set[entrez], evidence=new_evidence)

            # add instances to session
            self.session.add_all([mirna_set[mir_id], target_set[entrez], new_evidence, new_interaction])
        self.session.commit()