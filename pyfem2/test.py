import glob, imp, os, re, shutil, sys
F = os.path.realpath(__file__)
D = os.path.dirname(F)
sys.path.insert(0, D)

def teardown_module(module):
    def remove(a):
        if os.path.isfile(a): os.remove(a)
        elif os.path.isdir(a): shutil.rmtree(a)
    for filename in glob.glob('*.exo'):
        remove(filename)
    for filename in glob.glob('*.pvd'):
        remove(filename)
    for filename in glob.glob('*.vtu.d'):
        remove(filename)
    for filename in glob.glob('*.pyc'):
        remove(filename)

def load():
    """dynamically load a test from a module"""
    # Finally, we retrieve the Class
    found = {}
    regex = re.compile('\s*def test_')
    filenames = [f for f in glob.glob('*.py') if regex.search(open(f).read())
                 and f != 'test.py']
    for filename in filenames:
        name = os.path.splitext(os.path.basename(filename))[0]
        module = __import__(name)
        fname = lambda fun: fun + '_' + module.__name__
        found.update(dict([(fname(fun), getattr(module, fun))
                           for fun in dir(module) if fun.startswith('test_')]))
    globals().update(found)
    return
load()
