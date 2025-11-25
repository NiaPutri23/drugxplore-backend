FROM python:3.11-slim

RUN apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt /app/
RUN pip install --no-cache-dir -r requirements.txt

COPY . /app/

ENV PYTHONUNBUFFERED=1

EXPOSE 8000

CMD ["gunicorn", "antimalaria_ml.wsgi:application", "--bind", "0.0.0.0:8000", "--workers", "1", "--timeout", "120"]
