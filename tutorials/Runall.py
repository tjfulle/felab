import os
import sys
import glob
from subprocess import Popen
from os.path import join, dirname, realpath
D = dirname(dirname(realpath(__file__)))
env = dict(os.environ)
env['PYTHONPATH'] = D
env['NOGRAPHICS'] = '1'
files = [f for f in glob.glob('*.py') if 'Runall' not in f]
for filename in files:
    proc = Popen(['python', filename], env=env)
    proc.wait()
    assert proc.returncode == 0
for filename in glob.glob('*.exo'):
    os.remove(filename)
