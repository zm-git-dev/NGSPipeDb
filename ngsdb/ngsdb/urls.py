"""ngsdb URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/3.1/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path

from django.urls import include # Add a URL to urlpatterns

urlpatterns = [
    path('admin/', admin.site.urls),
    path(r'', include('home.urls')), # home page
    path(r'blastplus/', include('blastplus.urls')), # blast page
    path(r'geneAnno/', include('geneAnno.urls')), # gborwse page
    path(r'geneExpAtlas/', include('geneExpAtlas.urls')), # exp
    path(r'igv/', include('igv.urls')), # exp
    path(r'tools/', include('wooey.urls')),
]
