from sqlalchemy import Column, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

from pybel.constants import FUNCTION, NAMESPACE, NAME, MIRNA, RNA

ENTREZ_GENE_ID = 'EGID'
MIRTARBASE_ID = 'MIRTARBASE'

# create base class
Base = declarative_base()


class Mirna(Base):
    """Create mirna table that stores information about the miRNA"""
    __tablename__ = "mirna"

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
            NAME: self.mirtarbase_id
        }


class Target(Base):
    """Build target table, which stores information about the target gene"""
    __tablename__ = "target"

    id = Column(Integer, primary_key=True)
    target_gene = Column(String, nullable=False, doc="Target gene name")
    entrez_id = Column(Integer, nullable=False, unique=True, doc="Target gene Entrez ID")
    species = Column(String, nullable=False, doc="Species associated with target gene")
    hgnc_symbol = Column(String, nullable=True, doc="")
    hgnc_id = Column(String, nullable=True, doc="")

    def serialize_to_entrez_node(self):
        """Function to serialize to PyBEL node data dictionary.

        :rtype: dict
        """
        return {
            FUNCTION: RNA,
            NAMESPACE: ENTREZ_GENE_ID,
            NAME: self.entrez_id
        }

    def serialize_to_hgnc_node(self):
        """Function to serialize to PyBEL node data dictionary.

        :rtype: dict
        """
        return {
            FUNCTION: RNA,
            NAMESPACE: 'HGNC',
            NAME: self.hgnc_symbol
        }



class Evidence(Base):
    """Build Evidence table used to store MTI's and their evidence"""
    __tablename__ = "evidence"

    id = Column(Integer, primary_key=True)
    experiment = Column(String, nullable=False, doc="Experiments made to find miRNA - target interaction. E.g. 'Luciferase reporter assay//qRT-PCR//Western blot'")
    support = Column(String, nullable=False, doc="Type and strength of the miRNA - target interaction. E.g. 'Functional MTI (Weak)'")
    reference = Column(Integer, nullable=False, doc="Reference PubMed ID")


class Interaction(Base):
    """Build Interaction table used to store miRNA and target relations"""
    __tablename__ = "interaction"

    id = Column(Integer, primary_key=True)
    mirna_id = Column(Integer, ForeignKey("mirna.id"))
    mirna = relationship("Mirna", backref="interactions")

    target_id = Column(Integer, ForeignKey("target.id"))
    target = relationship("Target", backref="interactions")

    evidence_id = Column(Integer, ForeignKey("evidence.id"))
    evidence = relationship("Evidence", backref="interactions")
