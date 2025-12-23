from pydantic import BaseModel, field_validator
from typing import Optional
from uuid import UUID
from rdkit import Chem

class MoleculeBase(BaseModel):
    smiles: str
    name: Optional[str] = None

    @field_validator("smiles")
    @classmethod
    def validate_smiles(cls, v: str) -> str:
        if not Chem.MolFromSmiles(v):
            raise ValueError("Invalid SMILES")
        return v

class MoleculeCreate(MoleculeBase):
    id: Optional[UUID] = None

class MoleculeUpdate(BaseModel):
    smiles: Optional[str] = None
    name: Optional[str] = None

    @field_validator("smiles")
    @classmethod
    def validate_smiles(cls, v: str | None) -> str | None:
        if v is not None and not Chem.MolFromSmiles(v):
            raise ValueError("Invalid SMILES")
        return v

class MoleculeOut(MoleculeBase):
    id: UUID

class MoleculeList(BaseModel):
    total: int
    items: list[MoleculeOut]

class SearchRequest(BaseModel):
    substructure: str

    @field_validator("substructure")
    @classmethod
    def validate_sub(cls, v: str) -> str:
        if not Chem.MolFromSmiles(v):
            raise ValueError("Invalid substructure SMILES")
        return v

class SearchResponse(BaseModel):
    matches: list[MoleculeOut]
