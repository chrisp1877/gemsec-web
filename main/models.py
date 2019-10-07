from django.db import models
from django import forms


class Lab_Member(models.Model):
    user_id=models.IntegerField(default=0, primary_key=True)
    user_name=models.CharField(max_length=16, blank=True, null=True)
    # email

class Project(models.Model):
    project_id=models.IntegerField(default=0, primary_key=True)

class Event(models.Model):
    event_id=models.IntegerField(default=0, primary_key=True)



# Create your models here.
