try:
    string_types = (basestring,)
except NameError:
    string_types = (str, bytes)


def is_stringlike(a):
    return hasattr(a, "strip")


def is_listlike(a):
    """Is `a` like a list?"""
    if isinstance(a, string_types):
        return False
    try:
        len(a)
        return True
    except TypeError:
        return False
