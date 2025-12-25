from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session
from celery.result import AsyncResult
from ..db import get_db
from .. import crud
from ..schemas import SearchRequest
from ..celery_app import celery_app
from ..config import get_settings

router = APIRouter()
settings = get_settings()

@router.post("/search/async")
def start_async_search(
    req: SearchRequest,
    db: Session = Depends(get_db),
    limit: int = Query(default=settings.PAGINATION_DEFAULT_LIMIT, ge=1, le=settings.PAGINATION_MAX_LIMIT),
    offset: int = Query(default=0, ge=0),
):
    _, items = crud.list_molecules(db, offset=0, limit=10_000_000)
    objects = [{"id": str(m.id), "smiles": m.smiles, "name": m.name} for m in items]

    task = celery_app.send_task(
        "app.tasks.substructure_search",
        args=[objects, req.substructure],
    )
    return {"task_id": task.id}

@router.get("/search/tasks/{task_id}")
def get_search_task_status(task_id: str):
    task_result = AsyncResult(task_id, app=celery_app)

    if task_result.state == "PENDING":
        return {"status": "pending", "progress": None}
    elif task_result.state == "PROGRESS":
        meta = task_result.info or {}
        return {
            "status": "progress",
            "progress": {
                "current": meta.get("current", 0),
                "total": meta.get("total", 0),
            },
        }
    elif task_result.state == "SUCCESS":
        return {
            "status": "success",
            "result": task_result.result,
        }
    elif task_result.state in {"FAILURE", "REVOKED"}:
        return {"status": "failed", "error": str(task_result.info)}
    else:
        return {"status": task_result.state}
