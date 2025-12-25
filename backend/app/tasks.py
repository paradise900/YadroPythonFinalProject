from .celery_app import celery_app
from .search import find_matches
from rdkit import Chem

@celery_app.task(bind=True, name="app.tasks.substructure_search")
def substructure_search_task(self, objects: list[dict], substructure: str) -> list[dict]:
    if not Chem.MolFromSmiles(substructure):
        raise ValueError("Invalid substructure SMILES")

    total = len(objects)
    chunk = max(1, total // 10)

    results: list[dict] = []
    for idx, obj in enumerate(objects, start=1):
        if idx % chunk == 0 or idx == total:
            self.update_state(
                state="PROGRESS",
                meta={"current": idx, "total": total}
            )
        # фильтрация
        matches = find_matches([obj], substructure)
        if matches:
            results.extend(matches)

    return results
