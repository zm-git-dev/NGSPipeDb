from django.urls import path

from . import views

urlpatterns = [
    path('', views.index),
    path('gff2json', views.gff2json),
]