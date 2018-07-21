import os
import logging
from numpy import *
import numpy.linalg as la

from .utilities import *
from .constants import *


def assemble_system(model, u, du, Q, svtab, svars, dltyp, dload, predef,
                    procedure, step_type, time=array([0.,0.]), dtime=1.,
                    period=1., istep=1, iframe=1, nlgeom=False, ninc=None,
                    cflag=STIFF_AND_RHS, disp=0):
    """
    Assembles the global system of equations

    Parameters
    ----------
    u, du : ndarray
        Value of DOFs at beginning of step and increment, respectively.
    Q : ndarray
        Force array due to Neumann boundary condtions
    svtab : ndarray of int
        svtab[iel] are the indices for state variables for element iel
    svars : ndarray
        svtab[:,svtab[iel]] are the state variables for element iel
    dltyp : ndarray
        Distributed load type specifier
    dload : ndarray
        Distributed loads
    predef : ndarray
        Predefined fields
    procedure : symbolic constant
        The procedure specifier
    step_type : symbolic constant
        The step type
    time : ndarray, optional {array([0.,0.])}
        time[0] is the step time, time[1] the total time
    dtime : float {1.}
        Time increment
    period : float {1.}
        Step period
    istep, iframe : int, optional {1, 1}
        Step and frame numbers
    nlgeom : bool, optional {False}
        Nonlinear geometry
    ninc : int, opitional {None}
        Current increment
    cflag : symbolic constant, optional {STIFF_AND_RHS}

    Returns
    -------
    K : ndarray
        The (N,N) global stiffness array, where N is the total number of degrees
        of freedom in the probem.
    F : ndarray
        The (N,0) global RHS array, where N is the total number of degrees
        of freedom in the probem.

    Notes
    -----
    ``assemble`` implements a simplified assembler that adopts the
    following assumptions:

    - nodes are ordered continuously from 0 to :math:`n-1`;
    - there are no multifreedom constraints; and
    - the global stiffness matrix is stored as a full symmetric matrix.

    """
    procname = get_procname(procedure)
    steptypname = get_steptypname(step_type)
    msg  = 'ASSEMBLING GLOBAL SYSTEM OF EQUATIONS\n      '
    msg += 'PROCEDURE: {0}, STEP TYPE: {1}, NLGEOM: {2}\n      '.format(
        procname, steptypname, nlgeom)
    tf = time[-1] + dtime
    msg += 'STEP: {0}, FRAME: {1}, TIME: {2}'.format(istep, iframe, tf)
    if ninc is not None:
        msg += ', INCREMENT: {0}'.format(ninc)
    if not ninc:
        logging.debug(msg)
    elif ninc == 1 and iframe == 1:
        logging.debug(msg)

    if cflag not in CFLAGS:
        raise ValueError('UNKNOWN COMPUTE QUANTITY')

    compute_stiff = cflag in (STIFF_AND_RHS, STIFF_ONLY)
    compute_rhs = cflag in (STIFF_AND_RHS, RHS_ONLY, MASS_AND_RHS)
    compute_mass = cflag in (MASS_AND_RHS, MASS_ONLY)

    if compute_stiff:
        K = zeros((model.numdof, model.numdof))

    if compute_mass:
        M = zeros((model.numdof, model.numdof))

    if compute_rhs:
        fext = zeros(model.numdof)
        fint = zeros(model.numdof)

    # INTERPOLATE FIELD VARIABLES
    fac1 = time[1] / (time[0] + period)
    fac2 = (time[1]+dtime) / (time[0] + period)
    x0 = (1. - fac1) * predef[0] + fac1 * predef[1]
    xf = (1. - fac2) * predef[0] + fac2 * predef[1]
    predef_i = array([x0, xf-x0])

    # COMPUTE THE ELEMENT STIFFNESS AND SCATTER TO GLOBAL ARRAY
    for (ieb, eb) in enumerate(model.mesh.element_blocks):
        for (e, xel) in enumerate(eb.labels):

            # ELEMENT STIFFNESS
            iel = model.mesh.elemap[xel]
            el = model.elements[iel]
            eft = model.eftab[iel]
            response = el.response(u[eft], du[eft], time, dtime, istep,
                                   iframe, svars[:,svtab[iel]],
                                   dltyp[iel], dload[iel],
                                   predef_i[:,:,el.inodes], procedure, nlgeom,
                                   cflag, step_type)

            if cflag == STIFF_AND_RHS:
                K[IX(eft, eft)] += response[0]
                fext[eft] += response[1]
                fint[eft] += response[2]

            elif cflag == MASS_AND_RHS:
                M[IX(eft, eft)] += response[0]
                fext[eft] += response[1]
                fint[eft] += response[2]

            elif cflag == MASS_ONLY:
                M[IX(eft, eft)] += response

            elif cflag == STIFF_ONLY:
                K[IX(eft, eft)] += response

            elif cflag == RHS_ONLY:
                fext[eft] += response[0]
                fint[eft] += response[1]

    if compute_rhs:
        fext += Q

    if compute_mass:
        # LUMPED MASS MATRIX
        M = array([sum(row) for row in M])

    if cflag == STIFF_AND_RHS:
        if disp:
            return K, fext, fint
        return K, fext - fint

    elif cflag == STIFF_ONLY:
        return K

    elif cflag == RHS_ONLY:
        if disp:
            return fext, fint
        return fext - fint

    elif cflag == MASS_ONLY:
        return M

    elif cflag == MASS_AND_RHS:
        if disp:
            return M, fext, fint
        return M, fext - fint
