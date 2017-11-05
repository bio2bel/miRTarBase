# -*- coding: utf-8 -*-

from pybel.constants import FUNCTION, IDENTIFIER, MIRNA, NAME, NAMESPACE, RNA
from sqlalchemy import Column, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

ENTREZ_GENE_ID = 'EGID'
MIRTARBASE_ID = 'MIRTARBASE'

MIRTARBASE_PREFIX = 'mirtarbase'
MIRNA_TABLE_NAME = '{}_mirna'.format(MIRTARBASE_PREFIX)
TARGET_TABLE_NAME = '{}_target'.format(MIRTARBASE_PREFIX)
EVIDENCE_TABLE_NAME = '{}_evidence'.format(MIRTARBASE_PREFIX)
INTERACTION_TABLE_NAME = '{}_interaction'.format(MIRTARBASE_PREFIX)

# create base class
Base = declarative_base()


class Mirna(Base):
    """Create mirna table that stores information about the miRNA"""
    __tablename__ = MIRNA_TABLE_NAME

    id = Column(Integer, primary_key=True)

    mirtarbase_id = Column(String, nullable=False, unique=True, doc="miRTarBase ID")
    mir_name = Column(String, nullable=False, doc="miRNA name")
    species = Column(String, nullable=False, doc="Species associated with miRNA")

    def serialize_to_bel(self):
        """Function to serialize to PyBEL node data dictionary.

        :rtype: dict
        """
        return {
            FUNCTION: MIRNA,
            NAMESPACE: MIRTARBASE_ID,
            IDENTIFIER: self.mirtarbase_id,
            NAME: self.mir_name
        }


class Target(Base):
    """Build target table, which stores information about the target gene"""
    __tablename__ = TARGET_TABLE_NAME

    id = Column(Integer, primary_key=True)

    target_gene = Column(String, nullable=False, doc="Target gene name")
    entrez_id = Column(String, nullable=False, unique=True, doc="Target gene Entrez ID")
    species = Column(String, nullable=False, doc="Species associated with target gene")
    hgnc_symbol = Column(String, nullable=True, doc="HGNC Gene Symbol")
    hgnc_id = Column(String, nullable=True, doc="HGNC Gene Identifier")

    def serialize_to_entrez_node(self):
        """Function to serialize to PyBEL node data dictionary.

        :rtype: dict
        """
        return {
            FUNCTION: RNA,
            NAMESPACE: ENTREZ_GENE_ID,
            IDENTIFIER: self.entrez_id,
            NAME: self.target_gene,
        }

    def serialize_to_hgnc_node(self):
        """Function to serialize to PyBEL node data dictionary.

        :rtype: dict
        """
        return {
            FUNCTION: RNA,
            NAMESPACE: 'HGNC',
            IDENTIFIER: self.hgnc_id,
            NAME: self.hgnc_symbol
        }

    def to_json(self):
        """Returns this object as JSON

        :rtype: dict
        """
        return {
            'id': self.id,
            'species': self.species,
            'HGNC': {
                'symbol': self.hgnc_symbol,
                'id': self.hgnc_id
            },
            'Entrez': {
                'id': self.entrez_id,
                'name': self.target_gene
            }
        }


class Evidence(Base):
    """Build Evidence table used to store MTI's and their evidence"""
    __tablename__ = EVIDENCE_TABLE_NAME

    id = Column(Integer, primary_key=True)

    experiment = Column(String, nullable=False,
                        doc="Experiments made to find miRNA - target interaction. E.g. 'Luciferase reporter assay//qRT-PCR//Western blot'")
    support = Column(String, nullable=False,
                     doc="Type and strength of the miRNA - target interaction. E.g. 'Functional MTI (Weak)'")
    reference = Column(String, nullable=False, doc="Reference PubMed Identifier")


class Interaction(Base):
    """Build Interaction table used to store miRNA and target relations"""
    __tablename__ = INTERACTION_TABLE_NAME

    id = Column(Integer, primary_key=True)

    mirna_id = Column(Integer, ForeignKey("{}.id".format(MIRNA_TABLE_NAME)))
    mirna = relationship("Mirna", backref="interactions")

    target_id = Column(Integer, ForeignKey("{}.id".format(TARGET_TABLE_NAME)))
    target = relationship("Target", backref="interactions")

    evidence_id = Column(Integer, ForeignKey("{}.id".format(EVIDENCE_TABLE_NAME)))
    evidence = relationship("Evidence", backref="interactions")
