import sys
from numpy import *
from math import log as logm
from constants import *

# --------------------------------------------------------------------------- #
# -------------------------- UTILITY FUNCTIONS ------------------------------ #
# --------------------------------------------------------------------------- #
class UserInputError(Exception):
    def __init__(self, message):
        sys.tracebacklimit = -1
        super(UserInputError, self).__init__(message)

def IX(*args):
    out = []
    nd = len(args)
    for k, new in enumerate(args):
        new = asarray(new)
        if new.ndim != 1:
            raise ValueError("Cross index must be 1 dimensional")
        new = new.reshape((1,)*k + (new.size,) + (1,)*(nd-k-1))
        out.append(new)
    return tuple(out)

def Flatten(seq):
    return [a for b in seq for a in b]

def is_stringlike(a):
    return hasattr(a, 'strip')

def is_listlike(a):
    return (not hasattr(a, 'strip') and
            hasattr(a, '__getitem__') or hasattr(a, '__iter__'))

def DigitCount(p, d=1):
    return len([i for i in IntegerDigits(p) if i==d])

def FromDigits(seq):
    return sum([x * 10 ** i for (i,x) in enumerate(seq[::-1])])

def IntegerDigits(n):
    if not n:
        return [0]*MDOF
    m = 10**MDOF
    n += m
    return [n / 10 ** i % 10 for i in range(int(logm(n,10)),-1,-1)][1:]

def sign(a, b=None):
    sgn = {True: -1., False: 1.}[a<0.]
    if b is not None:
        sgn = b * sgn
    return sgn

def count_digits(seq, d=1):
    return len([i for i in seq if i == d])

def normal2d(xp):
    xp = asarray(xp)
    dx, dy = xp[1,:] - xp[0,:]
    n = array([dy, -dx], dtype=float)
    return n / sqrt(dot(n, n))
