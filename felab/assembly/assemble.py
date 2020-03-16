import numpy as np

import felab.util.tty as tty
from felab.util.numeric import IX
from felab.constants import (
    STATIC_ITERATIVE,
    STATIC_DIRECT,
    DYNAMIC,
    HEAT_TRANSFER_STEADY_STATE,
    RHS_ONLY,
    MASS_ONLY,
    MASS_AND_RHS,
    STIFF_AND_RHS,
    STIFF_ONLY,
)


def get_procname(proc):
    return {
        STATIC_DIRECT: "DIRECT STATIC",
        STATIC_ITERATIVE: "ITERATIVE STATIC",
        DYNAMIC: "DYNAMIC",
        HEAT_TRANSFER_STEADY_STATE: "STEADY STATE HEAT TRANSFER",
    }[proc]


def assemble_system(
    rhs,
    A,
    svtab,
    svars,
    energy,
    Q,
    u,
    du,
    v,
    a,
    time,
    dtime,
    kstep,
    kinc,
    kiter,
    dltyp,
    dlmag,
    predef,
    lflags,
    ddlmag,
    mdload,
    pnewdt,
    period,
    element_blocks,
    element_map,
    elements,
    element_freedom_table,
):
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
    dlmag : ndarray
        Distributed loads
    predef : ndarray
        Predefined fields
    time : ndarray, optional {np.array([0.,0.])}
        time[0] is the step time, time[1] the total time
    dtime : float {1.}
        Time increment
    period : float {1.}
        step period
    kstep, kinc : int, optional {1, 1}
        step and increment numbers
    nlgeom : bool, optional {False}
        Nonlinear geometry
    kiter : int, opitional {None}
        Current iteration
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
    procname = get_procname(lflags[0])
    nlgeom = lflags[1] == 1
    msg = "ASSEMBLING GLOBAL SYSTEM OF EQUATIONS\n      "
    msg += "PROCEDURE: {0}, NLGEOM: {1}\n      ".format(procname, nlgeom)

    tf = time[-1] + dtime
    msg += "STEP: {0}, INCREMENT: {1}, TIME: {2}".format(kstep, kinc, tf)
    if kiter is not None:
        msg += ", INCREMENT: {0}".format(kiter)

    if not kiter:
        tty.debug(msg)

    elif kiter == 1 and kinc == 1:
        tty.debug(msg)

    if lflags[2] not in (STIFF_AND_RHS, STIFF_ONLY, MASS_ONLY, RHS_ONLY, MASS_AND_RHS):
        raise NotImplementedError

    # INTERPOLATE FIELD VARIABLES
    fac1 = time[1] / (time[0] + period)
    fac2 = (time[1] + dtime) / (time[0] + period)
    x0 = (1.0 - fac1) * predef[0] + fac1 * predef[1]
    xf = (1.0 - fac2) * predef[0] + fac2 * predef[1]
    predef_i = np.array([x0, xf - x0])

    # COMPUTE THE ELEMENT STIFFNESS AND SCATTER TO GLOBAL ARRAY
    for (ieb, eb) in enumerate(element_blocks):
        for (e, xel) in enumerate(eb.labels):

            # ELEMENT STIFFNESS
            iel = element_map[xel]
            el = elements[iel]
            eft = element_freedom_table[iel]
            n = len(eft)
            A_e = np.zeros((n, n))
            rhs_e = np.zeros(n)
            el.response(
                rhs_e,
                A_e,
                svars[:, svtab[iel]],
                energy,
                u[eft],
                du[eft],
                v,
                a,
                time,
                dtime,
                kstep,
                kinc,
                dltyp[iel],
                dlmag[iel],
                predef_i[:, :, el.inodes],
                lflags,
                ddlmag,
                mdload,
                pnewdt,
            )
            if lflags[2] in (STIFF_AND_RHS, STIFF_ONLY, MASS_ONLY, MASS_AND_RHS):
                A[IX(eft, eft)] += A_e

            if lflags[2] in (STIFF_AND_RHS, RHS_ONLY, MASS_AND_RHS):
                rhs[eft] += rhs_e

    if lflags[2] in (STIFF_AND_RHS, RHS_ONLY):
        rhs[:] += Q

    return
