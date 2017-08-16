# -*- coding: utf-8 -*-

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, scoped_session

import pandas as pd

DATA_URL = 'http://mirtarbase.mbc.nctu.edu.tw/cache/download/6.1/miRTarBase_MTI.xlsx'

def get_data():
    """Gets miRTarBase Interactions table and exclude rows with NULL values

    :rtype: pandas.DataFrame
    """
    df = pd.read_excel("test.xlsx")
    # find null rows
    null_rows = pd.isnull(df).any(1).nonzero()[0]
    return df.drop(null_rows)

class Manager(object):
    def __init__(self, connection=None):
        self.connection = 'sqlite:///miRTarBase.db'
        self.engine = create_engine(self.connection)
        self.sessionmake = sessionmaker(bind=self.engine, autoflush=False, expire_on_commit=False)
        self.session = self.sessionmake()
        self.make_tables()


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