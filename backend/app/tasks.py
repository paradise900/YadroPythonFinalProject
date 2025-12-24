from .celery_app import celery_app
from .search import find_matches
from rdkit import Chem

@celery_app.task(name="app.tasks.substructure_search")
def substructure_search_task(objects: list[dict], substructure: str) -> list[dict]:
    if not Chem.MolFromSmiles(substructure):
        raise ValueError("Invalid substructure SMILES")
    return find_matches(objects, substructure)
