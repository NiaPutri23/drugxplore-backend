FROM python:3.11-slim

# Install required system dependencies
RUN apt-get update && apt-get install -y wget curl build-essential && rm -rf /var/lib/apt/lists/*

# Set the working directory
WORKDIR /app

# Copy requirements and install first (better caching)
COPY requirements.txt /app/
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of the code
COPY . /app/

ENV PYTHONUNBUFFERED 1
ENV DEBUG=False

# Port Gunicorn akan berjalan
EXPOSE 8000

# Optimal Gunicorn command for a multi-service API
# Workers = 4 (default for fast, I/O-bound tasks)
# Timeout = 30 (standard for quick responses)
# Bind = 0.0.0.0:8000 (essential for Docker networking)
CMD ["gunicorn", \
     "--workers", "4", \
     "--timeout", "30", \
     "--bind", "0.0.0.0:8000", \
     "antimalaria_backend.wsgi:application"]