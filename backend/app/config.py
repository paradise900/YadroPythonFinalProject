from pydantic_settings import BaseSettings
from functools import lru_cache

class Settings(BaseSettings):
    APP_NAME: str = "yadro-molecules"
    API_PREFIX: str = "/api"
    DATABASE_URL: str = "postgresql+psycopg2://postgres:postgres@db:5432/molecules"
    REDIS_URL: str = "redis://redis:6379/0"
    CACHE_TTL_SECONDS: int = 600
    CELERY_BROKER_URL: str = "redis://redis:6379/1"
    CELERY_RESULT_BACKEND: str = "redis://redis:6379/2"
    PAGINATION_DEFAULT_LIMIT: int = 20
    PAGINATION_MAX_LIMIT: int = 100

    class Config:
        env_file = ".env"

@lru_cache
def get_settings() -> Settings:
    return Settings()
