import sys
import os

pwd = os.path.dirname(os.path.realpath(__file__)) # script
workdir = os.path.join(os.getcwd(), 'ngsdb')

sys.path.append(workdir)
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'ngsdb.settings')

import django
django.setup()

from wooey.models.core import Script

scripts = Script.objects.all()
print('before:')
print(scripts)
print(len(scripts))

scripts = Script.objects.all().delete()

scripts = Script.objects.all()
print('after delete')
print(scripts)
print(len(scripts))