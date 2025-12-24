from celery import Celery
from .config import get_settings

s = get_settings()
celery_app = Celery(
    "yadro_tasks",
    broker=s.CELERY_BROKER_URL,
    backend=s.CELERY_RESULT_BACKEND,
)
celery_app.conf.task_routes = {"app.tasks.*": {"queue": "search"}}
