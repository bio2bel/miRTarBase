from sqlalchemy import Column, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine
from sqlalchemy.orm import relationship

#create base class
Base = declarative_base()

#create mirna table that stores information about the miRNA
class Mirna(Base):
    __tablename__ = "mirna"
    
    id = Column(Integer, primary_key = True)
    mirtarbase_id = Column(String, nullable = False, unique = True)
    mir_name = Column(String, nullable = False)
    species = Column(String, nullable = False)
    
    def __repr__(self):
        return "<Mirna(id='%s', fullname='%s', species='%s')>" % (
                                self.mirtarbase_id, self.mir_name, self.species)
    
#build target table, which stores information about the target gene
class Target(Base):
    __tablename__ = "target"
    
    id = Column(Integer, primary_key = True)
    target_gene = Column(String, nullable = False)
    entrez_id = Column(String, nullable = False, unique = True)
    species = Column(String, nullable = False)

#build Evidence table used to store MTI's and their evidence
class Evidence(Base):
    __tablename__ = "evidence"
    
    id = Column(Integer, primary_key = True)
    experiment = Column(String, nullable = False)
    support = Column(String, nullable = False)
    reference = Column(String, nullable = False)

#build Interaction table used to store miRNA and target relations
class Interaction(Base):
    __tablename__ = "interaction"
    
    id = Column(Integer, primary_key = True)
    mirna_id = Column(Integer, ForeignKey("mirna.id"))
    mirna = relationship("Mirna", backref = "interaction")
    
    target_id = Column(Integer, ForeignKey("target.id"))
    target = relationship("Target", backref = "interaction")
    
    evidence_id = Column(Integer, ForeignKey("evidence.id"))
    evidence = relationship("Evidence", backref = "interaction")

#create engine used to create tables
engine = create_engine('sqlite:///miRTarBase.db')
#create tables
Base.metadata.create_all(engine)
