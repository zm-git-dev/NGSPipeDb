from django.shortcuts import render

import math
from django.http import JsonResponse, HttpResponse
from .models import Exp
from django.core import serializers
import json

# Create your views here.

def index(request):
    return render(request, 'geneExpAtlas/index.html')

def exp_json(request):
    # https://github.com/rg3915/django-datatables-experiment/issues/1
    # ?start=0&length=10
    exps = Exp.objects.all()
    total = exps.count()
    
    _start = request.GET.get('start')
    _length = request.GET.get('length')
    if _start and _length:
        start = int(_start)
        length = int(_length)
        page = math.ceil(start / length) + 1 # math.ceil函数返回大于或等于一个给定数字的最小整数。
        per_page = length

        exps = exps[start: start+length]

    data = serializers.serialize("json", exps)
    data = json.loads(data)
    response = {
        'data': data,
        'page': page,  # [opcional]
        'per_page': per_page,  # [opcional]
        'recordsTotal': total,
        'recordsFiltered': total,
    }
    
    return JsonResponse(response, content_type='application/json')