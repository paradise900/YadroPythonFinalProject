from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session
from uuid import UUID
from .. import schemas, crud, models
from ..db import get_db
from ..config import get_settings
from ..cache import cache_with_ttl, invalidate_search_cache
from ..search import find_matches
from fastapi import status

router = APIRouter()
settings = get_settings()

@router.post("/molecules", response_model=schemas.MoleculeOut, status_code=201)
def create(data: schemas.MoleculeCreate, db: Session = Depends(get_db)):
    existing = crud.get_molecule_by_smiles(db, data.smiles)
    if existing:
        return existing
    m = crud.create_molecule(db, data)
    invalidate_search_cache()
    return m

@router.get("/molecules/{id}", response_model=schemas.MoleculeOut)
def read(id: UUID, db: Session = Depends(get_db)):
    m = crud.get_molecule(db, id)
    if not m:
        raise HTTPException(status_code=404, detail="Not found")
    return m

@router.put("/molecules/{id}", response_model=schemas.MoleculeOut)
def update(id: UUID, data: schemas.MoleculeUpdate, db: Session = Depends(get_db)):
    m = crud.update_molecule(db, id, data)
    if not m:
        raise HTTPException(status_code=404, detail="Not found")
    invalidate_search_cache()
    return m

@router.delete("/molecules/{id}", status_code=204)
def delete(id: UUID, db: Session = Depends(get_db)):
    ok = crud.delete_molecule(db, id)
    if not ok:
        raise HTTPException(status_code=404, detail="Not found")
    invalidate_search_cache()
    return

@router.get("/molecules", response_model=schemas.MoleculeList)
def list_molecules(
    db: Session = Depends(get_db),
    limit: int = Query(default=settings.PAGINATION_DEFAULT_LIMIT, ge=1, le=settings.PAGINATION_MAX_LIMIT),
    offset: int = Query(default=0, ge=0),
):
    total, items = crud.list_molecules(db, offset=offset, limit=limit)
    return {"total": total, "items": items}

def _search_cache_key(substructure: str, limit: int, offset: int) -> str:
    return f"search:{substructure}:{limit}:{offset}"

@cache_with_ttl(lambda substructure, limit, offset, items: _search_cache_key(substructure, limit, offset))
def _search_impl(substructure: str, limit: int, offset: int, items: list[dict]) -> list[dict]:
    matches = find_matches(items, substructure)
    return matches[offset: offset + limit]

@router.post("/search", response_model=schemas.SearchResponse, status_code=200)
def search(
    req: schemas.SearchRequest,
    db: Session = Depends(get_db),
    limit: int = Query(default=settings.PAGINATION_DEFAULT_LIMIT, ge=1, le=settings.PAGINATION_MAX_LIMIT),
    offset: int = Query(default=0, ge=0),
):
    _, items = crud.list_molecules(db, offset=0, limit=10_000_000)
    objects = [
        {"id": str(m.id), "smiles": m.smiles, "name": m.name}
        for m in items
    ]
    matches = _search_impl(req.substructure, limit, offset, objects)
    return {"matches": matches}
