#!/bin/bash
set -e

echo "Waiting for postgres..."
while ! nc -z db 5432; do
  sleep 0.1
done
echo "PostgreSQL started"

cd /app

echo "Creating database tables..."
python -m app.init_db

echo "Starting application..."
exec uvicorn app.main:app --host 0.0.0.0 --port 8000
