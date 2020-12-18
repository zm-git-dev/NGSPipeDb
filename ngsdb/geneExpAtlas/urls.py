from django.urls import path

from . import views

urlpatterns = [
    path('', views.index),
    path('exp_json', views.exp_json),
]