# -*- coding: utf-8 -*-

from sqlalchemy import Column, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine
from sqlalchemy.orm import relationship
from sqlalchemy.orm import sessionmaker
import pandas as pd

DATA_URL = 'http://mirtarbase.mbc.nctu.edu.tw/cache/download/6.1/miRTarBase_MTI.xlsx'

# create base class
Base = declarative_base()


class Mirna(Base):
    """create mirna table that stores information about the miRNA"""
    __tablename__ = "mirna"

    id = Column(Integer, primary_key=True)
    mirtarbase_id = Column(String, nullable=False, unique=True, doc="miRTarBase ID")
    mir_name = Column(String, nullable=False, doc="miRNA name")
    species = Column(String, nullable=False, doc="Species associated with miRNA")


class Target(Base):
    """build target table, which stores information about the target gene"""
    __tablename__ = "target"

    id = Column(Integer, primary_key=True)
    target_gene = Column(String, nullable=False, doc="Target gene name")
    entrez_id = Column(Integer, nullable=False, unique=True, doc="Target gene Entrez ID")
    species = Column(String, nullable=False, doc="Species associated with target gene")


class Evidence(Base):
    """build Evidence table used to store MTI's and their evidence"""
    __tablename__ = "evidence"

    id = Column(Integer, primary_key=True)
    experiment = Column(String, nullable=False, doc="Experiments made for evidence")
    support = Column(String, nullable=False, doc="Type and strength of the MTI")
    reference = Column(Integer, nullable=False, doc="Reference PubMed ID")


class Interaction(Base):
    """build Interaction table used to store miRNA and target relations"""
    __tablename__ = "interaction"

    id = Column(Integer, primary_key=True)
    mirna_id = Column(Integer, ForeignKey("mirna.id"))
    mirna = relationship("Mirna", backref="interaction")

    target_id = Column(Integer, ForeignKey("target.id"))
    target = relationship("Target", backref="interaction")

    evidence_id = Column(Integer, ForeignKey("evidence.id"))
    evidence = relationship("Evidence", backref="interaction")


def get_data():
    """Gets miRTarBase Interactions table and exclude rows with NULL values

    :rtype: pandas.DataFrame
    """
    df = pd.read_excel(DATA_URL)
    # find null rows
    null_rows = pd.isnull(df).any(1).nonzero()[0]
    return df.drop(null_rows)


def populate(session):
    """Populate database with the data from miRTarBase

    :param session: session object from sqlalchemy
    :return:
    """
    df = get_data()
    mirna_set = {}
    target_set = {}

    # iterate through rows and construct tables from it
    for index, mir_id, mirna, species_mirna, target, entrez, species_target, exp, sup_type, pubmed in df.itertuples():
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
        session.add_all([mirna_set[mir_id], target_set[entrez], new_evidence, new_interaction])
    session.commit()


if __name__ == '__main__':
    # create engine used to create tables
    engine = create_engine('sqlite:///miRTarBase.db')
    # create tables
    Base.metadata.create_all(engine)
    # create session
    Session = sessionmaker(bind=engine)
    session = Session()
    # populate database
    populate(session)