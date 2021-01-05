"""
Django settings for ngsdb project.

Generated by 'django-admin startproject' using Django 3.1.4.

For more information on this file, see
https://docs.djangoproject.com/en/3.1/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/3.1/ref/settings/
"""

from pathlib import Path

# Build paths inside the project like this: BASE_DIR / 'subdir'.
BASE_DIR = Path(__file__).resolve().parent.parent


# data file dir ---------------
import os
from os import path
workdir = os.getcwd()
ngsdb_DIR = path.join(path.dirname(BASE_DIR), "results", "sqlite3")

EXP_DB_PATH = os.path.join(ngsdb_DIR, 'exp', 'exp.sqlite3')
EXP_DB_PATH = os.path.join(workdir, 'results', 'sqlite3', 'exp', 'exp.sqlite3')
GFF_DB_PATH = os.path.join(workdir, 'results', 'sqlite3', 'gff_sqlite3', 'gtf.sqlite3')


#SITE_TITLE = 'Plukenetia volubilis'

#SERVER_NAME = 'http://pvdb.xtbg.ac.cn'

#FTP_NAME = 'ftp://pvdb.xtbg.ac.cn'

#GENOME_FASTA = os.path.join(BASE_DIR, 'PVDB.genome.fa')
# -----------------------------


# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/3.1/howto/deployment/checklist/

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = 's#0s)u3x!r4c0d=xw4m+$z(307vf-&$1%ou+c2#kgd#7=505qx'

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = True

ALLOWED_HOSTS = ['*']


# Application definition

INSTALLED_APPS = [
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'home', # first page of website
    'geneAnno', # gene annotation from nr/nt/pfam/go/kegg
    'geneExpAtlas', # gene expression matrix
    'blastplus', # blast tool
    'igv',
    'wooey',
    'tools',
    'search',
]

MIDDLEWARE = [
    'django.middleware.security.SecurityMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
]

ROOT_URLCONF = 'ngsdb.urls'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
            ],
        },
    },
]

WSGI_APPLICATION = 'ngsdb.wsgi.application'


# Database
# https://docs.djangoproject.com/en/3.1/ref/settings/#databases

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': BASE_DIR / 'db.sqlite3',
    },
    # add new database -------------------------------------------------
    'expDb': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': EXP_DB_PATH
    },
    'gffDb': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': GFF_DB_PATH
    },
    # ------------------------------------------------------------------
}

# use multi-database in django -----------------------------
DATABASE_ROUTERS = ['ngsdb.db_router.DatabaseAppsRouter']
DATABASE_APPS_MAPPING = {
    # example:
    #'app_name':'database_name',
    'geneExpAtlas': 'expDb',
    'geneAnno': 'gffDb',
}
#-----------------------------------------------------------


# Password validation
# https://docs.djangoproject.com/en/3.1/ref/settings/#auth-password-validators

AUTH_PASSWORD_VALIDATORS = [
    {
        'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator',
    },
]


# Internationalization
# https://docs.djangoproject.com/en/3.1/topics/i18n/

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'UTC'

USE_I18N = True

USE_L10N = True

USE_TZ = True


# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/3.1/howto/static-files/

STATIC_URL = '/static/'

# add to collect stastic file-------------------------
import os
# 如下为了关闭开发模式的时候还能找到静态文件
STATICFILES_DIRS = (
    os.path.join(BASE_DIR, "commonstatic"),
    os.path.join(workdir, "results"),
)
# 如下为了关闭开发模式的时候还能找到静态文件
STATIC_ROOT = os.path.join(workdir,'collectstatic')
# ----------------------------------------------------

# allow iframe
X_FRAME_OPTIONS = 'SAMEORIGIN'
