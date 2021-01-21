from django.urls import path

from . import views

app_name = 'home'

urlpatterns = [
    path('about/', views.about),
    path('', views.index, name='home'),
]