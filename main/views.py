from django.shortcuts import render

# Create your views here.

from django.http import HttpResponseRedirect
from django.http import HttpResponse
from django.shortcuts import render
from .forms import MD_Visualization_Form
from django.core.files.storage import FileSystemStorage
from django.contrib.auth.models import User

from . import functions

# Imaginary function to handle an uploaded file.
# from .src.Main_MD_analysis import *

# Home page of application portal (where GEMSEC future homepage will live)
def index(request):
    context = {}
    return render(request, 'main/index.html', context)


def md_analysis(request):
    submitted = False
    if request.method == 'POST':
        form = MD_Visualization_Form(request.POST, request.FILES)
        if form.is_valid():
            submitted = True
            analysis = form.cleaned_data["analysis"]
            visualization = form.cleaned_data["viz"]
            upload_file = form.cleaned_data["file"]
            handle_uploaded_file(request.FILES['file'])
            Main_MD_analysis.md_vis(request.FILES['file'].name, analysis, visualization)
            return render(request, 'main/md_analysis.html', locals())
    else:
        form = MD_Visualization_Form()
        return render(request, 'main/md_analysis.html', locals())

def handle_uploaded_file(f):
    with open('media/'+f.name, 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)

def uploaded(request):
    return render(request, 'main/uploaded.html', {})

def about(request):
    context = {}
    return render(request, 'main/about.html', context)

def applications(request):
    context = {}
    return render(request, 'main/applications.html', context)

def contact(request):
    context = {}
    return render(request, 'main/contact.html', context)

def members(request):
    users = User.objects.all()
    context = {}
    return render(request, 'main/members.html', context)

def papers(request):
    context = {}
    return render(request, 'main/papers.html', context)



