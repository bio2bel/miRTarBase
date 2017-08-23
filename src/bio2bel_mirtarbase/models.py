from sqlalchemy import Column, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

# create base class
Base = declarative_base()


class Mirna(Base):
    """Create mirna table that stores information about the miRNA"""
    __tablename__ = "mirna"

    id = Column(Integer, primary_key=True)
    mirtarbase_id = Column(String, nullable=False, unique=True, doc="miRTarBase ID")
    mir_name = Column(String, nullable=False, doc="miRNA name")
    species = Column(String, nullable=False, doc="Species associated with miRNA")


class Target(Base):
    """Build target table, which stores information about the target gene"""
    __tablename__ = "target"

    id = Column(Integer, primary_key=True)
    target_gene = Column(String, nullable=False, doc="Target gene name")
    entrez_id = Column(Integer, nullable=False, unique=True, doc="Target gene Entrez ID")
    species = Column(String, nullable=False, doc="Species associated with target gene")


class Evidence(Base):
    """Build Evidence table used to store MTI's and their evidence"""
    __tablename__ = "evidence"

    id = Column(Integer, primary_key=True)
    experiment = Column(String, nullable=False, doc="Experiments made for evidence")
    support = Column(String, nullable=False, doc="Type and strength of the MTI")
    reference = Column(Integer, nullable=False, doc="Reference PubMed ID")


class Interaction(Base):
    """Build Interaction table used to store miRNA and target relations"""
    __tablename__ = "interaction"

    id = Column(Integer, primary_key=True)
    mirna_id = Column(Integer, ForeignKey("mirna.id"))
    mirna = relationship("Mirna", backref="interaction")

    target_id = Column(Integer, ForeignKey("target.id"))
    target = relationship("Target", backref="interaction")

    evidence_id = Column(Integer, ForeignKey("evidence.id"))
    evidence = relationship("Evidence", backref="interaction")
