from sqlalchemy import Column, String, DateTime, func, Index
from sqlalchemy.dialects.postgresql import UUID
import uuid
from .db import Base

class Molecule(Base):
    __tablename__ = "molecules"
    id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    smiles = Column(String, nullable=False, unique=True)
    name = Column(String, nullable=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    
    __table_args__ = (
        Index('ix_molecules_smiles', 'smiles'),
    )
