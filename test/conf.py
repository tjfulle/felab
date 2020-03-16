import os
import glob
import shutil
from os.path import join, dirname, realpath, isfile, isdir, basename

D = dirname(dirname(realpath(__file__)))


def teardown_module(module):
    def remove(a):
        if isfile(a):
            os.remove(a)
        elif isdir(a):
            shutil.rmtree(a)

    for (dirname, dirs, files) in os.walk(D):
        d = basename(dirname)
        if d == ".git":
            del dirs[:]
            continue
        if d.endswith(".vtu.d"):
            remove(d)
            del dirs[:]
        for filename in files:
            if filename.endswith((".exo", ".pyc", ".pvd")):
                remove(join(dirname, filename))
            if filename == ".DS_Store":
                remove(join(dirname, filename))
