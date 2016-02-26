import os
import glob
import shutil
from os.path import join, dirname, realpath, isfile, isdir
D = dirname(dirname(realpath(__file__)))
def teardown_module(module):
    def remove(a):
        if isfile(a): os.remove(a)
        elif isdir(a): shutil.rmtree(a)
    for filename in glob.glob('*.exo'):
        remove(filename)
    for filename in glob.glob('*.pvd'):
        remove(filename)
    for filename in glob.glob('*.vtu.d'):
        remove(filename)
    for filename in glob.glob('*.pyc'):
        remove(filename)
