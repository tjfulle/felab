import os
import pytest
import shutil


this_path = os.path.dirname(os.path.realpath(__file__))


def remove(path):
    if os.path.isfile(path):
        os.remove(path)
    elif os.path.isdir(path):
        shutil.rmtree(path)


@pytest.fixture(scope="session", autouse=True)
def cleanall():
    yield
    for (dirname, dirs, files) in os.walk("."):
        d = os.path.basename(dirname)
        if dirname.endswith(".git"):
            del dirs[:]
            continue
        elif dirname.endswith( ".vtu.d"):
            remove(d)
            del dirs[:]
            continue
        for name in files:
            if name.endswith((".exo", ".pyc", ".pvd")) or name in (".DS_Store",):
                remove(os.path.join(dirname, name))


@pytest.fixture(scope="function", autouse=False)
def data_path():
    yield os.path.join(this_path, "data")
