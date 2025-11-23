import os
import environ
import dj_database_url # Import untuk koneksi DB
from pathlib import Path
from datetime import timedelta

BASE_DIR = Path(__file__).resolve().parent.parent

# ML_MODEL_DIR tidak diperlukan di API backend, tetapi jika ada, gunakan:
# ML_MODEL_DIR = BASE_DIR / "media" / "ml_models"

# Inisialisasi Environment
env = environ.Env()
# Di production (Railway), JANGAN baca .env karena variabel di-inject langsung.
# environ.Env.read_env(os.path.join(BASE_DIR, '.env'))

# --- CORE SECURITY & ENVIRONMENT ---

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = os.getenv("SECRET_KEY", 'default-insecure-key-for-local-dev')

# Ambil DEBUG dari ENV. Default ke False jika tidak ada
DEBUG = os.getenv("DEBUG", "False").lower() == "true" 

# Host yang diizinkan
if DEBUG:
    ALLOWED_HOSTS = ['*']
else:
    # Ambil domain publik Railway dan tambahkan host lain jika perlu
    RAILWAY_PUBLIC_DOMAIN = os.getenv('RAILWAY_PUBLIC_DOMAIN')
    ALLOWED_HOSTS = [RAILWAY_PUBLIC_DOMAIN] if RAILWAY_PUBLIC_DOMAIN else []
    ALLOWED_HOSTS.append('.railway.app') # Mengizinkan semua subdomain Railway

# --- ML SERVICE CONFIGURATION ---
# URL internal untuk berkomunikasi dengan layanan ML
ML_SERVICE_URL = os.getenv("ML_SERVICE_URL") 
# Contoh nilai di Railway: http://antimalaria-ml:8000

# --- APP DEFINITION ---

BUILTIN_APPS = [
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django.contrib.sites',
]

THIRD_PARTY_APPS = [
    'rest_framework',
    'rest_framework.authtoken',
    'corsheaders',
    'drf_spectacular'
]

LOCAL_APPS = [
    'api',
    'api.v1.prediction_compounds',
]

INSTALLED_APPS = BUILTIN_APPS + THIRD_PARTY_APPS + LOCAL_APPS

# ... (SOCIALACCOUNT_PROVIDERS, SITE_ID, REST_FRAMEWORK, SPECTACULAR_SETTINGS biarkan)

# --- MIDDLEWARE ---
MIDDLEWARE = [
    'corsheaders.middleware.CorsMiddleware',
    'django.middleware.security.SecurityMiddleware',
    'whitenoise.middleware.WhiteNoiseMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
]

# --- CORS & CSRF (Production Ready) ---

FRONTEND_URL = os.getenv("FRONTEND_URL", "http://localhost:5173")

if DEBUG:
    CORS_ALLOW_ALL_ORIGINS = True
    CSRF_TRUSTED_ORIGINS = [FRONTEND_URL]
    CSRF_COOKIE_SECURE = False
    SESSION_COOKIE_SECURE = False
else:
    # Di Production, batasi CORS ke domain frontend
    CORS_ALLOW_ALL_ORIGINS = False
    CORS_ALLOWED_ORIGINS = [FRONTEND_URL] 
    
    # TRUSTED ORIGINS: Tambahkan domain Railway dan Frontend
    trusted_origins = [FRONTEND_URL]
    if RAILWAY_PUBLIC_DOMAIN:
        trusted_origins.append(f"https://{RAILWAY_PUBLIC_DOMAIN}")
        
    CSRF_TRUSTED_ORIGINS = trusted_origins
    
    # Wajib di Production karena menggunakan HTTPS
    CSRF_COOKIE_SECURE = True
    SESSION_COOKIE_SECURE = True

CORS_ALLOW_CREDENTIALS = True
CSRF_COOKIE_SAMESITE = 'Lax'
SESSION_COOKIE_SAMESITE = 'Lax'
SESSION_COOKIE_HTTP_ONLY = True

# --- STATIC FILES ---
STATIC_URL = '/static/'
STATIC_ROOT = BASE_DIR / 'staticfiles'
MEDIA_URL = '/media/'
MEDIA_ROOT = BASE_DIR / 'media'
STATICFILES_STORAGE = 'whitenoise.storage.CompressedManifestStaticFilesStorage'

# ... (ROOT_URLCONF, TEMPLATES, WSGI_APPLICATION biarkan)


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
    
# Hapus konfigurasi TEST: MIRROR, karena sering bermasalah dengan dj-database-url

AUTH_USER_MODEL = 'api.CustomUser'

# ... (AUTH_PASSWORD_VALIDATORS, I18N, TZ, AUTO_FIELD biarkan)

AUTHENTICATION_BACKENDS = (
    'django.contrib.auth.backends.ModelBackend',
    # 'allauth.account.auth_backends.AuthenticationBackend',
)

LOGIN_REDIRECT_URL = '/'
LOGOUT_REDIRECT_URL = '/'

# ... (LOGGING, EMAIL_BACKEND, MEDIA_URL/ROOT biarkan)