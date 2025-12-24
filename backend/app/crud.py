from sqlalchemy.orm import Session
from sqlalchemy import select, func
from uuid import UUID
from . import models, schemas

def get_molecule_by_smiles(db: Session, smiles: str) -> models.Molecule | None:
    return db.execute(
        select(models.Molecule).where(models.Molecule.smiles == smiles)
    ).scalar_one_or_none()

def create_molecule(db: Session, data: schemas.MoleculeCreate) -> models.Molecule:
    m = models.Molecule(id=data.id, smiles=data.smiles, name=data.name)
    db.add(m)
    db.commit()
    db.refresh(m)
    return m

def get_molecule(db: Session, id_: UUID) -> models.Molecule | None:
    return db.get(models.Molecule, id_)

def update_molecule(db: Session, id_: UUID, data: schemas.MoleculeUpdate) -> models.Molecule | None:
    m = db.get(models.Molecule, id_)
    if not m:
        return None
    if data.smiles is not None:
        m.smiles = data.smiles
    if data.name is not None:
        m.name = data.name
    db.commit()
    db.refresh(m)
    return m

def delete_molecule(db: Session, id_: UUID) -> bool:
    m = db.get(models.Molecule, id_)
    if not m:
        return False
    db.delete(m)
    db.commit()
    return True

def list_molecules(db: Session, offset: int, limit: int) -> tuple[int, list[models.Molecule]]:
    total = db.execute(select(func.count(models.Molecule.id))).scalar_one()
    items = db.execute(select(models.Molecule).offset(offset).limit(limit)).scalars().all()
    return total, items
