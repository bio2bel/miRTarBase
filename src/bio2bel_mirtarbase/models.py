# -*- coding: utf-8 -*-

"""SQLAlchemy models for Bio2BEL miRTarBase."""

from typing import Iterable, Mapping

import pybel.dsl
import pybel.dsl
from pybel import BELGraph
from sqlalchemy import Column, ForeignKey, Index, Integer, String, UniqueConstraint
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

from .constants import MODULE_NAME

NCBIGENE = 'ncbigene'
MIRBASE = 'mirbase'
HGNC = 'hgnc'

MIRNA_TABLE_NAME = f'{MODULE_NAME}_mirna'
TARGET_TABLE_NAME = f'{MODULE_NAME}_target'
SPECIES_TABLE_NAME = f'{MODULE_NAME}_species'
EVIDENCE_TABLE_NAME = f'{MODULE_NAME}_evidence'
INTERACTION_TABLE_NAME = f'{MODULE_NAME}_interaction'

# create base class
Base = declarative_base()


class Species(Base):
    """Represents a species."""

    __tablename__ = SPECIES_TABLE_NAME

    id = Column(Integer, primary_key=True)

    taxonomy_id = Column(String(255), nullable=True, unique=True, index=True, doc='The NCBI taxonomy identifier')
    name = Column(String(255), nullable=False, unique=True, index=True, doc='The scientific name for the species')

    def to_json(self, include_id: bool = True) -> Mapping:
        """Serialize to JSON.

        :param include_id: Include the database identifier?
        """
        rv = {
            'name': str(self.name),
        }

        if include_id:
            rv['id'] = self.id

        return rv

    def __str__(self):  # noqa: D105
        return self.name


class Mirna(Base):
    """Create mirna table that stores information about the miRNA."""

    __tablename__ = MIRNA_TABLE_NAME

    id = Column(Integer, primary_key=True)

    name = Column(String(31), nullable=False, unique=True, index=True, doc="miRTarBase name")
    mirbase_id = Column(String(255), nullable=True, unique=True, index=True, doc="miRBase identifier")
    entrez_id = Column(String(255), nullable=True, unique=True, index=True, doc="Entrez Gene identifier")

    species_id = Column(Integer, ForeignKey(f'{Species.__tablename__}.id'), nullable=False, doc='The host species')
    species = relationship(Species)

    def as_bel(self) -> pybel.dsl.MicroRna:
        """Serialize to a PyBEL node data dictionary."""
        return pybel.dsl.MicroRna(
            namespace=MIRBASE,
            name=self.name,
            identifier=self.mirbase_id,
        )

    @staticmethod
    def filter_name_in(names: Iterable[str]):
        """Build a name filter."""
        return Mirna.name.in_(names)

    def __str__(self):  # noqa: D105
        return self.name


class Target(Base):
    """Represents a target RNA."""

    __tablename__ = TARGET_TABLE_NAME

    id = Column(Integer, primary_key=True)

    name = Column(String(63), nullable=False, index=True, doc="Target gene name")
    entrez_id = Column(String(32), nullable=False, unique=True, index=True, doc="Entrez gene identifier")

    hgnc_symbol = Column(String(32), nullable=True, unique=True, index=True, doc="HGNC gene symbol")
    hgnc_id = Column(String(32), nullable=True, unique=True, index=True, doc="HGNC gene identifier")

    species_id = Column(Integer, ForeignKey(f'{SPECIES_TABLE_NAME}.id'), nullable=False, doc='The host species')
    species = relationship('Species')

    def __str__(self):  # noqa: D105
        return self.name

    def serialize_to_entrez_node(self) -> pybel.dsl.Rna:
        """Serialize to PyBEL node data dictionary."""
        return pybel.dsl.Rna(
            namespace=NCBIGENE,
            identifier=str(self.entrez_id),
            name=str(self.name),
        )

    def serialize_to_hgnc_node(self) -> pybel.dsl.Rna:
        """Serialize to PyBEL node data dictionary."""
        if self.hgnc_id is None:
            raise ValueError(f'missing HGNC information for Entrez Gene {self.entrez_id}')

        return pybel.dsl.Rna(
            namespace=HGNC,
            identifier=str(self.hgnc_id),
            name=str(self.hgnc_symbol)
        )

    def to_json(self, include_id=True) -> Mapping:
        """Return this object as JSON."""
        rv = {
            'species': self.species.to_json(),
            'identifiers': [
                self.serialize_to_entrez_node(),
                self.serialize_to_hgnc_node(),
            ]
        }

        if include_id:
            rv['id'] = self.id

        return rv


class Interaction(Base):
    """Build Interaction table used to store miRNA and target relations."""

    __tablename__ = INTERACTION_TABLE_NAME

    id = Column(Integer, primary_key=True)

    mirtarbase_id = Column(String(64), nullable=False, unique=True, index=True,
                           doc="miRTarBase interaction identifier which is unique for a pair of miRNA and RNA targets")

    mirna_id = Column(Integer, ForeignKey(f'{MIRNA_TABLE_NAME}.id'), nullable=False, index=True,
                      doc='The miRTarBase identifier of the interacting miRNA')
    mirna = relationship(Mirna, backref="interactions")

    target_id = Column(Integer, ForeignKey(f'{TARGET_TABLE_NAME}.id'), nullable=False, index=True,
                       doc='The Entrez gene identifier of the interacting RNA')
    target = relationship(Target, backref="interactions")

    __table_args__ = (
        UniqueConstraint('mirna_id', 'target_id'),
        Index('interaction_idx', 'mirna_id', 'target_id', unique=True),
    )

    def __str__(self):  # noqa: D105
        return f'{self.mirna.name} =| {self.target.name}'


class Evidence(Base):
    """Build Evidence table used to store MTI's and their evidence."""

    __tablename__ = EVIDENCE_TABLE_NAME

    id = Column(Integer, primary_key=True)

    experiment = Column(
        String(255),
        nullable=False,
        doc="Experiments made to find miRNA - target interaction."
            " E.g. 'Luciferase reporter assay//qRT-PCR//Western blot'"
    )
    support = Column(String(255), nullable=False,
                     doc="Type and strength of the miRNA - target interaction. E.g. 'Functional MTI (Weak)'")
    reference = Column(String(255), nullable=False, doc="Reference PubMed Identifier")

    interaction_id = Column(Integer, ForeignKey(f'{INTERACTION_TABLE_NAME}.id'),
                            doc='The interaction for which this evidence was captured')
    interaction = relationship(Interaction, backref="evidences")

    def __str__(self):  # noqa: D105
        return '{}: {}'.format(self.reference, self.support)

    def add_to_graph(self, graph: BELGraph) -> str:
        """Add this edge to the BEL graph and return the ket for that edge."""
        return self._add_to_graph(
            graph,
            self.interaction.mirna.as_bel(),
            self.interaction.target.serialize_to_entrez_node(),
        )

    def _add_to_graph(self, graph: BELGraph, source: pybel.dsl.MicroRna, target: pybel.dsl.Rna) -> str:
        """Add this edge to the BEL graph and return the ket for that edge."""
        return graph.add_directly_decreases(
            source,
            target,
            evidence=str(self.support),
            citation=str(self.reference),
            annotations={
                'Experiment': str(self.experiment),
                'SupportType': str(self.support),
            }
        )
