# -*- coding: utf-8 -*-

from sqlalchemy import Column, ForeignKey, Index, Integer, String, UniqueConstraint
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

from pybel.constants import DIRECTLY_DECREASES
from pybel.dsl import mirna, rna

ENTREZ_GENE_ID = 'EGID'
MIRTARBASE = 'MIRTARBASE'

MIRTARBASE_PREFIX = 'mirtarbase'
MIRNA_TABLE_NAME = '{}_mirna'.format(MIRTARBASE_PREFIX)
TARGET_TABLE_NAME = '{}_target'.format(MIRTARBASE_PREFIX)
SPECIES_TABLE_NAME = '{}_species'.format(MIRTARBASE_PREFIX)
EVIDENCE_TABLE_NAME = '{}_evidence'.format(MIRTARBASE_PREFIX)
INTERACTION_TABLE_NAME = '{}_interaction'.format(MIRTARBASE_PREFIX)

# create base class
Base = declarative_base()


class Mirna(Base):
    """Create mirna table that stores information about the miRNA"""
    __tablename__ = MIRNA_TABLE_NAME

    id = Column(Integer, primary_key=True)

    name = Column(String, nullable=False, unique=True, index=True, doc="miRTarBase name")

    species_id = Column(Integer, ForeignKey('{}.id'.format(SPECIES_TABLE_NAME)), nullable=False, doc='The host species')
    species = relationship('Species')

    def serialize_to_bel(self):
        """Function to serialize to PyBEL node data dictionary.

        :rtype: dict
        """
        return mirna(
            namespace=MIRTARBASE,
            name=str(self.name)
        )

    def __str__(self):
        return self.mirtarbase_name


class Target(Base):
    """Build target table, which stores information about the target gene"""
    __tablename__ = TARGET_TABLE_NAME

    id = Column(Integer, primary_key=True)

    name = Column(String, nullable=False, index=True, doc="Target gene name")
    entrez_id = Column(String, nullable=False, unique=True, index=True, doc="Entrez gene identifier")

    hgnc_symbol = Column(String, nullable=True, unique=True, index=True, doc="HGNC gene symbol")
    hgnc_id = Column(String, nullable=True, unique=True, index=True, doc="HGNC gene identifier")

    species_id = Column(Integer, ForeignKey('{}.id'.format(SPECIES_TABLE_NAME)), nullable=False, doc='The host species')
    species = relationship('Species')

    def __str__(self):
        return self.name

    def serialize_to_entrez_node(self):
        """Function to serialize to PyBEL node data dictionary.

        :rtype: dict
        """
        return rna(
            namespace=ENTREZ_GENE_ID,
            identifier=str(self.entrez_id),
            name=str(self.name)
        )

    def serialize_to_hgnc_node(self):
        """Function to serialize to PyBEL node data dictionary.

        :rtype: dict
        """
        if self.hgnc_id is None:
            raise ValueError('missing HGNC information for Entrez Gene {}'.format(self.entrez_id))

        return rna(
            namespace='HGNC',
            identifier=str(self.hgnc_id),
            name=str(self.hgnc_symbol)
        )

    def to_json(self, include_id=True):
        """Returns this object as JSON

        :rtype: dict
        """
        rv = {
            'species': self.species.to_json(),
            'identifiers': [
                self.serialize_to_entrez_node(),
                self.serialize_to_hgnc_node()
            ]
        }

        if include_id:
            rv['id'] = self.id

        return rv


class Species(Base):
    """Represents a species"""
    __tablename__ = SPECIES_TABLE_NAME

    id = Column(Integer, primary_key=True)

    name = Column(String, nullable=False, unique=True, index=True, doc='The scientific name for the species')

    def to_json(self, include_id=True):
        rv = {
            'name': str(self.name)
        }

        if include_id:
            rv['id'] = self.id

        return rv

    def __str__(self):
        return self.name


class Interaction(Base):
    """Build Interaction table used to store miRNA and target relations"""
    __tablename__ = INTERACTION_TABLE_NAME

    id = Column(Integer, primary_key=True)

    mirtarbase_id = Column(String, nullable=False, unique=True, index=True,
                           doc="miRTarBase interaction identifier which is unique for a pair of miRNA and RNA targets")

    mirna_id = Column(Integer, ForeignKey("{}.id".format(MIRNA_TABLE_NAME)), nullable=False, index=True,
                      doc='The miRTarBase identifier of the interacting miRNA')
    mirna = relationship("Mirna", backref="interactions")

    target_id = Column(Integer, ForeignKey("{}.id".format(TARGET_TABLE_NAME)), nullable=False, index=True,
                       doc='The Entrez gene identifier of the interacting RNA')
    target = relationship("Target", backref="interactions")

    __table_args__ = (
        UniqueConstraint('mirna_id', 'target_id'),
        Index('interaction_idx', 'mirna_id', 'target_id', unique=True),
    )

    def __str__(self):
        return '{} =| {}'.format(self.mirna.name, self.target.name)


class Evidence(Base):
    """Build Evidence table used to store MTI's and their evidence"""
    __tablename__ = EVIDENCE_TABLE_NAME

    id = Column(Integer, primary_key=True)

    experiment = Column(String, nullable=False,
                        doc="Experiments made to find miRNA - target interaction. E.g. 'Luciferase reporter assay//qRT-PCR//Western blot'")
    support = Column(String, nullable=False,
                     doc="Type and strength of the miRNA - target interaction. E.g. 'Functional MTI (Weak)'")
    reference = Column(String, nullable=False, doc="Reference PubMed Identifier")

    interaction_id = Column(Integer, ForeignKey("{}.id".format(INTERACTION_TABLE_NAME)),
                            doc='The interaction for which this evidence was captured')
    interaction = relationship("Interaction", backref="evidences")

    def __str__(self):
        return '{}: {}'.format(self.reference, self.support)

    def add_to_graph(self, graph):
        try:
            target_node = self.interaction.target.serialize_to_hgnc_node()
        except ValueError:
            target_node = self.interaction.target.serialize_to_entrez_node()

        graph.add_qualified_edge(
            self.interaction.mirna.serialize_to_bel(),
            target_node,
            relation=DIRECTLY_DECREASES,
            evidence=str(self.support),
            citation=str(self.reference),
            annotations={
                'Experiment': str(self.experiment),
                'SupportType': str(self.support),
            }
        )
