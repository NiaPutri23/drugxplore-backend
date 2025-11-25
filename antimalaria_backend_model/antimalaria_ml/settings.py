import os
import environ
import dj_database_url
from pathlib import Path

# Build paths inside the project like this: BASE_DIR / 'subdir'.
BASE_DIR = Path(__file__).resolve().parent.parent

# --- ML Model Path Configuration ---
# Pastikan model ML berada di sub-folder media
MEDIA_ROOT = BASE_DIR / 'media'
ML_MODEL_DIR = MEDIA_ROOT / "ml_models" 

# Initialize environment variables
env = environ.Env()
# Di production (Railway), JANGAN baca .env karena variabel diinject langsung.
# environ.Env.read_env(BASE_DIR / '.env') 


# --- CORE SECURITY SETTINGS (PRODUCTION READY) ---
# SECURITY WARNING: keep the secret key used in production secret!
# Ambil SECRET_KEY dari Environment Variable di Railway
SECRET_KEY = os.getenv('SECRET_KEY')

# SECURITY WARNING: don't run with debug turned on in production!
# Ambil DEBUG dari Environment Variable. Default ke False (Production Mode)
DEBUG = os.getenv("DEBUG", "False").lower() == "true"

# Layanan ML adalah internal, jadi ALLOWED_HOSTS bisa tetap longgar
ALLOWED_HOSTS = ['*']

# Application definition
BUILTIN_APPS = [
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
]

THIRD_PARTY_APPS = [
    'rest_framework',
    'corsheaders',
]

LOCAL_APPS = [
    'api', # Your ML API app
]

INSTALLED_APPS = BUILTIN_APPS + THIRD_PARTY_APPS + LOCAL_APPS

# --- MIDDLEWARE ---
MIDDLEWARE = [
    # Tambahkan CORSHeaders Middleware
    'corsheaders.middleware.CorsMiddleware',
    'django.middleware.security.SecurityMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
]

SITE_ID = 1

# Layanan ML adalah internal, jadi boleh CORS_ALLOW_ALL_ORIGINS=True 
# (karena hanya API Backend yang memanggilnya)
CORS_ALLOW_ALL_ORIGINS = True

ROOT_URLCONF = 'antimalaria_ml.urls'

# REST Framework: Boleh AllowAny karena mengandalkan firewall internal Railway
REST_FRAMEWORK = {
    'DEFAULT_AUTHENTICATION_CLASSES': [
        'rest_framework.authentication.SessionAuthentication',
        # Hapus BasicAuthentication jika tidak diperlukan
    ],
    'DEFAULT_PERMISSION_CLASSES': [
        'rest_framework.permissions.AllowAny',
    ],
}

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
            ],
        },
    },
]

WSGI_APPLICATION = 'antimalaria_ml.wsgi.application'


# --- DATABASE CONFIGURATION (OPTIMIZED for NeonDB/Railway) ---
# Menggunakan dj_database_url untuk parsing URL koneksi tunggal
if os.getenv('DATABASE_URL'):
    DATABASES = {
        'default': dj_database_url.config(
            default=os.getenv('DATABASE_URL'),
            conn_max_age=600,
            # Wajib untuk koneksi ke NeonDB di Production
            ssl_require=True 
        )
    }
else:
    # Fallback untuk pengembangan lokal tanpa DATABASE_URL
    DATABASES = {
        'default': {
            'ENGINE': 'django.db.backends.sqlite3',
            'NAME': BASE_DIR / 'db.sqlite3',
        }
    }


# Password validation
# ... (Biarkan pengaturan validasi password tetap seperti ini)

# Internationalization
LANGUAGE_CODE = 'en-us'
TIME_ZONE = 'Asia/Jakarta'
USE_I18N = True
USE_TZ = True


# Static files (CSS, JavaScript, Images)
# Layanan ML mungkin tidak perlu mengumpulkan statis, tapi biarkan ada untuk Django
STATIC_URL = 'static/'
STATIC_ROOT = BASE_DIR / 'staticfiles'

# Media files: Penting untuk lokasi model ML Anda
MEDIA_URL = '/media/'
# MEDIA_ROOT sudah didefinisikan di awal sebagai BASE_DIR / 'media'

DEFAULT_AUTO_FIELD = 'django.db.models.BigAutoField'