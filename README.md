# YADRO Molecules Service

Pet‑проект по заданию YADRO: сервис хранения молекул в SMILES‑нотации с субструктурным и similarity‑поиском. Реализован на FastAPI с PostgreSQL, SQLAlchemy, Redis, Celery, Nginx, Docker/Compose, CI и лёгким фронтендом.

## Функциональность

- Хранение молекул:
  - SMILES, название, UUID, дата добавления.
  - Уникальность по SMILES (одна молекула — один SMILES).
- Субструктурный поиск:
  - Синхронный поиск по RDKit через `POST /api/search`.
  - Валидация SMILES (молекул и субструктур) через RDKit.
- Асинхронный поиск:
  - Celery‑таска для субструктурного поиска.
  - API для запуска async‑поиска и проверки статуса задачи (через Redis).
- Similarity‑поиск (базовый):
  - Расчёт сходства Tanimoto по Morgan‑fingerprints (утилита в коде, при необходимости легко оборачивается в эндпоинт).
- UI:
  - HTML/CSS/JS интерфейс для добавления молекул, просмотра списка и субструктурного поиска.
  - Уведомления об ошибках и успешных операциях.

## Архитектура

- **Backend** (`backend/app`):
  - `main.py` — создание FastAPI‑приложения, подключение роутеров под `/api`.
  - `config.py` — конфигурация (DATABASE_URL, REDIS_URL, CELERY_*, лимиты пагинации) через `pydantic-settings`.
  - `db.py` — SQLAlchemy `engine`, `SessionLocal`, `Base`.
  - `models.py` — модель `Molecule`: UUID PK, уникальный `smiles`, индекс по `smiles`.
  - `schemas.py` — Pydantic‑схемы, валидация SMILES/субструктур через RDKit.
  - `crud.py` — CRUD‑функции и пагинация.
  - `search.py` — RDKit‑логика: субструктурный поиск и similarity (Tanimoto).
  - `cache.py` — Redis‑кеш (TTL, инвалидация, graceful‑fallback при недоступности Redis).
  - `celery_app.py`, `tasks.py` — Celery‑приложение и задача асинхронного субструктурного поиска.
  - `routers/molecules.py` — CRUD для молекул + синхронный субструктурный поиск.
  - `routers/search_async.py` — асинхронный поиск и эндпоинт статуса задачи.

- **База данных**:
  - PostgreSQL (сервис `db` в docker-compose).
  - SQLAlchemy ORM.
  - Alembic‑миграции (`backend/alembic.ini`, `backend/migrations/*`), запуск через `alembic upgrade head` при старте backend.

- **Кеш и брокер**:
  - Redis (сервис `redis`):
    - кеширует результаты поиска по ключам `search:{substructure}:{limit}:{offset}` с TTL;
    - используется как брокер и backend для Celery.

- **Celery**:
  - Worker (сервис `worker`) выполняет задачу субструктурного поиска.
  - Flower (сервис `flower`) для мониторинга задач.

- **Nginx**:
  - Отдаёт статику фронтенда из `frontend/`.
  - Проксирует `/api/*` на FastAPI (`backend:8000`).

- **Frontend** (`frontend/index.html` + `frontend/app.js` или inline JS):
  - Формы для добавления молекул, просмотра списка и субструктурного поиска.
  - Обработка 404/422 и других ошибок через уведомления.

- **CI**:
  - GitHub Actions (`.github/workflows/ci.yml`):
    - black, flake8, mypy;
    - pytest с покрытием;
    - сборка Docker‑образа backend.

## Запуск

### Требования

- Docker + Docker Compose.
- Свободен порт `80` (Nginx) и `5555` (Flower, опционально).

### Старт проекта

```bash
docker compose build
docker compose up -d
```

Поднимаются сервисы:

- `db` — PostgreSQL
- `redis` — Redis
- `backend` — FastAPI + Alembic
- `worker` — Celery worker
- `flower` — мониторинг Celery (http://localhost:5555)
- `nginx` — фронт + reverse proxy (http://localhost)

### Миграции (если нужно вручную)

```bash
docker compose exec backend alembic upgrade head
```

## API

### CRUD для молекул

- **Создать молекулу**

  `POST /api/molecules`

  Пример запроса:

  ```json
  {
    "smiles": "O",
    "name": "Water",
    "id": "опционально, UUID"
  }
  ```

  Пример ответа 201:

  ```json
  {
    "id": "313cf931-e67c-47b1-9089-70e3e6d54f41",
    "smiles": "O",
    "name": "Water"
  }
  ```

  Если молекула с таким SMILES уже существует, возвращается существующая (идемпотентность по SMILES).

- **Получить по id**

  `GET /api/molecules/{id}` → 200 или 404.

- **Обновить**

  `PUT /api/molecules/{id}`:

  ```json
  {
    "name": "New name",
    "smiles": "CCO"
  }
  ```

  → 200 или 404.

- **Удалить**

  `DELETE /api/molecules/{id}` → 204 или 404.

- **Список с пагинацией**

  `GET /api/molecules?limit=20&offset=0`

  ```json
  {
    "total": 2,
    "items": [
      {
        "id": "313cf931-e67c-47b1-9089-70e3e6d54f41",
        "smiles": "O",
        "name": "Water"
      }
    ]
  }
  ```

### Синхронный субструктурный поиск

`POST /api/search`

```json
{ "substructure": "c1ccccc1" }
```

Ответ 200:

```json
{
  "matches": [
    {
      "id": "08321ffc-470a-4a1a-b32a-458786a9c253",
      "smiles": "c1ccccc1",
      "name": "Benzene"
    }
  ]
}
```

Особенности:

- Если БД пуста → 404 `{"detail": "No molecules in database"}`.
- Если молекулы есть, но совпадений нет → 404 `{"detail": "No molecules match given substructure"}`.
- Если субструктура — невалидный SMILES → 422 (ошибка валидации Pydantic/RDKit).

Результаты кешируются в Redis с TTL, кеш инвалидации при изменении молекул.

### Асинхронный поиск (Celery)

- **Запуск async‑поиска**

  `POST /api/search/async`

  ```json
  { "substructure": "c1ccccc1" }
  ```

  Ответ:

  ```json
  { "task_id": "celery-task-id" }
  ```

- **Статус задачи**

  `GET /api/search/tasks/{task_id}`

  Возможные ответы:

  ```json
  { "status": "pending" }
  ```

  ```json
  {
    "status": "progress",
    "progress": { "current": 42, "total": 100 }
  }
  ```

  ```json
  {
    "status": "success",
    "result": [
      { "id": "...", "smiles": "...", "name": "..." }
    ]
  }
  ```

  ```json
  {
    "status": "failed",
    "error": "..." 
  }
  ```

## Frontend

Фронтенд находится в папке `frontend/` и разворачивается Nginx'ом на `http://localhost`.

Возможности:

- **Добавление молекулы**:
  - SMILES (обязателен), название, UUID (опционально).
  - Уведомления об успехе/ошибках (включая ошибки SMILES).

- **Список молекул**:
  - Таблица с названием, SMILES, ID (обрезанный), кнопкой удаления.
  - Отображение общего количества и пустое состояние.

- **Субструктурный поиск**:
  - Поле для SMILES субструктуры.
  - Карточки результатов (название, SMILES, ID).
  - Обработка:
    - 404 ("нет молекул в БД" или "нет совпадений") — пустое состояние + уведомление.
    - 422 ("некорректный SMILES") — ошибка в уведомлении.

## Тесты и качество кода

- Тесты (`backend/tests`):
  - unit‑тесты RDKit‑поиска,
  - API‑тесты CRUD и поиска,
  - тесты кеша и инвалидации.

- CI (GitHub Actions):
  - Запуск при `push`/`pull_request`.
  - Линтеры: `black`, `flake8`, `mypy` (опционально).
  - Тесты: `pytest` с покрытием.
  - Сборка Docker‑образа backend.

## Стек технологий

- **Язык**: Python 3.11
- **Backend**: FastAPI, Pydantic v2, SQLAlchemy 2.x, Alembic
- **БД**: PostgreSQL
- **Кеш/брокер**: Redis
- **Async задачи**: Celery
- **Хемоинформатика**: RDKit (SMILES, подструктуры, fingerprints)
- **Инфраструктура**: Docker, Docker Compose, Nginx
- **CI**: GitHub Actions (pytest, coverage, black, flake8, mypy, docker build)

## Структура проекта

```
.
├── backend/
│   ├── app/
│   │   ├── main.py
│   │   ├── config.py
│   │   ├── db.py
│   │   ├── models.py
│   │   ├── schemas.py
│   │   ├── crud.py
│   │   ├── search.py
│   │   ├── cache.py
│   │   ├── celery_app.py
│   │   ├── tasks.py
│   │   └── routers/
│   │       ├── molecules.py
│   │       └── search_async.py
│   ├── tests/
│   ├── alembic.ini
│   ├── migrations/
│   ├── Dockerfile
│   ├── requirements.txt
│   └── entrypoint.sh
├── frontend/
│   ├── index.html
│   └── app.js
├── nginx/
│   └── nginx.conf
├── docker-compose.yml
├── .github/
│   └── workflows/
│       └── ci.yml
└── README.md
```
