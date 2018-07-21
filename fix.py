import os
import sys
import glob
import shutil
from subprocess import Popen

x = {
    'pyfem2': 'felab',
     }

for (pat, rep) in x.items():
    command = "gsed -i 's:\<{0}\>:{1}:g' *.py".format(pat, rep)
    for (dirname, dirs, files) in os.walk('./'):
        if os.path.basename(dirname) == '.git':
            del dirs[:]
            continue
        if os.path.realpath(dirname) == os.path.realpath('./'):
            continue
        if not glob.glob(os.path.join(dirname, './*.py')):
            continue
        p = Popen(command, cwd=dirname, shell=True)
        p.wait()

