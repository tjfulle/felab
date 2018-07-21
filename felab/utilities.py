import sys
import logging
from numpy import *
from numpy.linalg import solve, lstsq, LinAlgError
from math import log as logm
try:
    import scipy.linalg.flapack as flapack
except ImportError:
    flapack = None

from .constants import *

# --------------------------------------------------------------------------- #
# -------------------------- UTILITY FUNCTIONS ------------------------------ #
# --------------------------------------------------------------------------- #
class UserInputError(Exception):
    def __init__(self, message):
        sys.tracebacklimit = 1
        super(UserInputError, self).__init__(message)

def IX(*args):
    out = []
    nd = len(args)
    for k, new in enumerate(args):
        new = asarray(new)
        if new.ndim != 1:
            raise ValueError("CROSS INDEX MUST BE 1 DIMENSIONAL")
        new = new.reshape((1,)*k + (new.size,) + (1,)*(nd-k-1))
        out.append(new)
    return tuple(out)

def Flatten(seq):
    return [a for b in seq for a in b]

try:
  string_types = (basestring,)
except NameError:
  string_types = (str, bytes)

def is_stringlike(a):
    return hasattr(a, 'strip')

def is_listlike(a):
    """Is `a` like a list?"""
    if isinstance(a, string_types):
        return False
    try:
        len(a)
        return True
    except TypeError:
        return False

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

def normal2d(xc):
    xc = asarray(xc)
    if xc.shape[0] == 2:
        # LINEAR ELEMENT
        dx, dy  = xc[1,:] - xc[0,:]
        n = array([dy, -dx], dtype=float)
    elif xc.shape[0] == 3:
        p = poly1d(polyfit(xc[:,0], xc[:,1], 2))
        dp = polyder(p)
        dy = average([dp(xc[1,0]), dp(xc[0,0]), dp(xc[2,0])])
        n = array([dy, -1.], dtype=float)
    return n / sqrt(dot(n, n))

def linsolve(A, b, symmetric=True):
    """Interface to the lapack dposv solve function in scipy.linalg

    Parameters
    ----------
    A : ndarray
        Real, symmetric, positive-definite matrix (the stiffness matrix)
    b : ndarray
        RHS of system of equations

    Returns
    -------
    c : ndarray
    x : ndarray
        Solution to A x = b
    info : int
        info < 0 -> bad input
        info = 0 -> solution is in x
        ifno > 0 -> singular matrix

    Notes
    -----
    dposv solves the system of equations A x = b using lapack's dposv
    procedure. This interface function is used to avoid the overhead of
    calling down in to scipy, converting arrays to fortran order, etc.

    """
    try:
        F = b.asarray()
    except AttributeError:
        F = asarray(b)

    use_np_solve = not symmetric or flapack == None
    x, info = None, 1
    if not use_np_solve:
        c, x, info = flapack.dposv(A, F, lower=0, overwrite_a=0, overwrite_b=0)
        if info < 0:
            raise ValueError("ILLEGAL VALUE IN {0}-TH ARGUMENT OF "
                             "INTERNAL DPOSV".format(-info))
        if info != 0:
            use_np_solve = True

    if use_np_solve:
        try:
            x = solve(A, F)
            info = 0
        except LinAlgError:
            raise RuntimeError('ATTEMPTING TO SOLVE UNDER CONSTRAINED SYSTEM')

    if info > 0:
        logging.warn("LINSOLVE FAILED, USING LEAST SQUARES "
                     "TO SOLVE SYSTEM")
        x = lstsq(A, F)[0]

    return x

def iso_dev_split1(ndir, nshr, numdim, D):
    ntens = ndir + nshr
    D1, D2 = zeros((ntens, ntens)), eye(ntens)
    D1[:ndir,:ndir] = D[0,1]
    return D1, D-D1

def iso_dev_split2(ndir, nshr, numdim, D):
    I = array([1. for i in range(2)] + [0. for i in range(2)])
    a = dot(I, D) / numdim
    Diso = outer(a, I)
    return Diso, D-Diso

iso_dev_split = iso_dev_split1

def get_procname(proc):
    return {STATIC: 'STATIC',
            DYNAMIC: 'DYNAMIC',
            HEAT_TRANSFER: 'HEAT TRANSFER'}[proc]

def get_steptypname(step_type):
    return step_type.upper()

def emptywithlists(n):
    a = zeros(n, dtype=object)
    a[:] = [[] for _ in range(n)]
    return a

def axialv(a):
    """Construct the axial vector associated with a

    w_i = -1/2 e_ijk a_jk
    Parameters
    ----------
    a : ndarray (3, 3)
        Second-order tensors

    Returns
    -------
    w : ndarray (3,)
        Axial vector associated with a
    """
    return .5 * array([a[2, 1] - a[1, 2],
                       a[0, 2] - a[2, 0],
                       a[1, 0] - a[0, 1]])

def axialt(a):
    """Construct the axial tensor associated with a
    W_ij = e_ijk a_k

    Parameters
    ----------
    a : ndarray (3,)
        Vector

    Returns
    -------
    W : ndarray (3,)
        Axial tensor associated with a
    """
    return array([[   0, -a[2], a[1]],
                  [a[2],   0, -a[0]],
                  [-a[1], a[0],   0]])

def asvec(t, ndir=3, nshr=3):
    """Cast the symmetric tensor as an array

    Parameters
    ----------
    t : ndarray (3,3)
        Tensor

    Returns
    -------
    W : ndarray (ndir+nshr,)
        Array

    """
    a = zeros(ndir+nshr)
    a[:ndir] = diag(t)[:ndir]
    a[ndir:ndir+nshr] = [t[0,1], t[1,2], t[0,2]][:nshr]
    return a
