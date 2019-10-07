from django.urls import path
from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path('md_analysis/', views.md_analysis, name='md_analysis'),
    path('about/', views.about, name='about'),
    path('applications/', views.applications, name='applications'),
    path('contact/', views.contact, name='contact'),
    path('members/', views.members, name='members'),
    path('papers/', views.papers, name='papers'),
    path('uploaded/', views.uploaded, name='uploaded')
    #path('news/', include('news.urls'))
]