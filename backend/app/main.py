import logging
from fastapi import FastAPI
from .routers import molecules
from .config import get_settings

logging.basicConfig(level=logging.INFO)
settings = get_settings()
app = FastAPI(title=settings.APP_NAME)
app.include_router(molecules.router, prefix=settings.API_PREFIX)

@app.get("/health")
def health():
    return {"status": "ok"}
