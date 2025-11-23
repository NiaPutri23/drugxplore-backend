FROM python:3.11-slim

# 2. Instal dependensi sistem yang mungkin diperlukan oleh library ML
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# 3. Tetapkan direktori kerja
WORKDIR /app

# 4. Salin requirements.txt dan install (agar lebih efisien dalam caching Docker)
COPY requirements.txt /app/
# Pastikan semua library ML (TensorFlow, Scikit-learn, dll.) ada di sini
RUN pip install --no-cache-dir -r requirements.txt

# 5. Salin sisa kode proyek (termasuk model ML yang dilatih)
COPY . /app/

# Lingkungan
ENV PYTHONUNBUFFERED 1
ENV DEBUG=False

# 6. Expose port yang akan digunakan
EXPOSE 8000

# 7. Command Gunicorn yang stabil untuk ML
CMD ["gunicorn", \
     "--workers", "1", \
     "--timeout", "120", \
     "--keep-alive", "10", \
     "--bind", "0.0.0.0:8000", \
     "antimalaria_ml.wsgi:application"]