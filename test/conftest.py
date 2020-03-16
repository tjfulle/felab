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
def remove_generated_files():
    yield
    for (dirname, dirs, files) in os.walk(this_path):
        d = os.path.basename(dirname)
        if d == ".git":
            del dirs[:]
            continue
        if d.endswith(".vtu.d"):
            remove(d)
            del dirs[:]
        for filename in files:
            if filename.endswith((".exo", ".pyc", ".pvd")):
                remove(os.path.join(dirname, filename))
            if filename == ".DS_Store":
                remove(os.path.join(dirname, filename))


@pytest.fixture(scope="function", autouse=False)
def data_path():
    yield os.path.join(this_path, "data")
