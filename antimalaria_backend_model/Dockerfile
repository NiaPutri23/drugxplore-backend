# Base image
FROM python:3.9-slim

# Install system dependencies
RUN apt-get update && apt-get install -y wget curl && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy project files
COPY . /app/

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Environment variables
ENV PYTHONUNBUFFERED=1
ENV DEBUG=False

# Expose port for the model service
EXPOSE 8001

# Gunicorn command
CMD ["gunicorn", "--workers", "1", "--timeout", "120", "--keep-alive", "10", "antimalaria_ml.wsgi:application"]
