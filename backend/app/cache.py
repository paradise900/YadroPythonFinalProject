import json
import logging
import redis
from typing import Callable, Any
from .config import get_settings

log = logging.getLogger(__name__)
settings = get_settings()

def get_redis() -> redis.Redis| None:
    try:
        return redis.from_url(settings.REDIS_URL, decode_responses=True)
    except Exception as e:
        log.warning(f'Redis unavailable: {e}')
        return None

def cache_with_ttl(key_builder: Callable[..., str]):
    def decorator(fn: Callable[..., Any]):
        def wrapper(*args, **kwargs):
            r = get_redis()
            key = key_builder(*args, **kwargs)
            if r:
                try:
                    cached = r.get(key)
                    if cached:
                        return json.loads(cached)
                except Exception as e:
                    log.warning(f'Redis get failed: {e}')
            result = fn(*args, **kwargs)
            if r:
                try:
                    r.setex(key, settings.CACHE_TTL_SECONDS, json.dumps(result))
                except Exception as e:
                    log.warning(f'Redis get failed: {e}')
            return result
        return wrapper
    return decorator

def invalidate_search_cache():
    r = get_redis()
    if not r:
        return
    try:
        cursor = 0
        while True:
            cursor, keys = r.scan(cursor=cursor, match='search:*', count=100)
            if keys:
                r.delete(*keys)
            if cursor == 0:
                break
    except Exception as e:
        log.warning(f"Redis invalidation failed: {e}")
