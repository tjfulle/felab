def apply_boundary_conditions(K, F, doftags, dofvals, u=None, du=None):
    """
    Apply boundary conditions to the global stiffness ``K`` and global
    force ``F``.

    Parameters
    ----------
    K : ndarray
        Global stiffness
    F : ndarray
        Global force

    Returns
    -------
    Kbc, Fbc : ndarray
        Boundary condition modified stiffness and force

    Notes
    -----
    Boundary conditions are applied in such a way that ``K`` remains
    symmetric by transferring columns associated with known degrees of
    freedom to ``F``. This method is intended to be called in the
    application code's ``Solve`` method.

    """
    if u is None:
        return _apply_bc_1(K, F, doftags, dofvals)
    return _apply_bc_2(K, F, doftags, dofvals, u, du)

def _apply_bc_1(K, F, doftags, dofvals):
    ubc = []
    Kbc,  Fbc = K.copy(), F.copy()
    for (i, I) in enumerate(doftags):
        ubc.append(dofvals[i])
        Fbc -= [K[k,I] * dofvals[i] for k in range(K.shape[0])]
        Kbc[I,:] = Kbc[:,I] = 0.
        Kbc[I,I] = 1.
    Fbc[doftags] = ubc
    return Kbc, Fbc

def _apply_bc_2(K, F, doftags, dofvals, u, du):
    ubc = []
    Kbc, Fbc = K.copy(), F.copy()
    for (i, I) in enumerate(doftags):
        u_cur = u[I] + du[I]
        ufac = dofvals[i] - u_cur
        ubc.append(ufac)
        Fbc -= [K[k,I] * ufac for k in range(K.shape[0])]
        Kbc[I,:] = Kbc[:,I] = 0.
        Kbc[I,I] = 1.
    Fbc[doftags] = ubc
    return Kbc, Fbc
