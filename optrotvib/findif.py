def wiberg_displacement_sizes(omegas, temp):
    c1 = 16.857
    c2 = 0.719384
    if temp < 0:
        raise ValueError("Temperature must be in Kelvin, therefore it can't be negative")
    if temp > 0:
        size_A0 = [np.sqrt((c1 / v) * (np.cosh(c2 * v / temp) / np.sinh(c2 * v / temp))) for v in omegas]
    if temp == 0:
        raise ValueError("Wiberg zero point displacement sizes have not been implemented yet")
    return [x / physconst['bohr2angstroms'] for x in size_A0]


def mort_displacment_sizes(S, avg_disp_per_atom=0.04):
    three_n, nmode = S.shape
    natom = three_n / 3
    d = 0.0
    for m in range(nmode):
        for n in range(natom):
            atom_sum = S[n, m]**2 + S[n + 1, m]**2 + S[n + 2, m]**2
            d += np.sqrt(atom_sum)

    delta = avg_disp_per_atom * (nmode * natom) / (np.sqrt(3) * d)
    return [delta] * nmode


def _get_TR_space(m, geom, space='TR', tol=LINEAR_A_TOL, verbose=False):
    """Form the idealized translation and rotation degrees of freedom from geometry `geom` and massed `m`.

    Parameters
    ----------
    m: array (natom)
        The masses of each atom
    geom: array (natom, 3)
        The positions of each atom
    space: {'T','R', 'TR'} default 'TR'
        Space to generate, either Translations, 'T', rotations 'R' or both 'TR' (default)
    tol: tolerance (to handle noisy linear geometries)
    verbose: {True, False} default False
        Set to True to enable extra printing
    """
    sqrtmmm = np.repeat(np.sqrt(m), 3)
    xxx = np.repeat(geom[:, 0], 3)
    yyy = np.repeat(geom[:, 1], 3)
    zzz = np.repeat(geom[:, 2], 3)

    z = np.zeros_like(m)
    i = np.ones_like(m)
    ux = np.ravel([i, z, z], order='F')
    uy = np.ravel([z, i, z], order='F')
    uz = np.ravel([z, z, i], order='F')

    T = [sqrtmmm * v for v in [ux, uy, uz]]
    R = [sqrtmmm * (yyy * uz - zzz * uy), sqrtmmm * (zzz * ux - xxx * uz), sqrtmmm * (xxx * uy - yyy * ux)]

    TRspace = []
    if 'T' in space:
        TRspace = TRspace + T
    if 'R' in space:
        TRspace = TRspace + R

    TRspace = np.vstack(TRspace)

    def orth(A, tol=tol):
        u, s, vh = np.linalg.svd(A, full_matrices=False)
        if verbose:
            print(s)
        M, N = A.shape
        eps = np.finfo(float).eps
        if not tol:
            tol = max(M, N) * np.amax(s) * eps
        num = np.sum(s > tol, dtype=int)
        Q = u[:, :num]
        return Q

    TRindep = orth(TRspace.T)
    TRindep = TRindep.T

    if verbose:
        print(TRindep.shape, '<--', TRspace.shape)
        print(np.linalg.norm(TRindep, axis=1))
        print('-' * 80)

    return TRindep


def _phase_cols_to_max_element(arr, tol=LINEAR_A_TOL, verbose=False):
    """Returns copy of 2D `arr` scaled such that, within cols, max(fabs)
    element is positive. If max(fabs) is pos/neg pair, scales so first
    element (within `tol`) is positive.

    """
    rephasing = []
    arr2 = np.zeros_like(arr)

    def phase_max_pos(arr):
        max_idx = 0
        max_val = -1
        for i, v in enumerate(arr):
            if abs(v) >= max_val:
                max_idx = i
                max_val = abs(v)
        # print("Max idx: {}".format(max_idx))
        # print("Value : {}".format(arr[max_idx]))
        phase = arr[max_idx] / abs(arr[max_idx])
        if phase < 0:
            rephase = True
        else:
            rephase = False
        arr *= phase

        return arr, rephase

    for v in range(arr.shape[1]):
        arr2[:, v], rep = phase_max_pos(arr[:, v].copy())
        if rep:
            rephasing.append(str(v))

    if rephasing and verbose:
        print('Negative modes rephased:', ', '.join(rephasing))

    return arr2


def get_vibinfo(hess, geom, mass):
    """Computes the vibinfo for this mass/geom/hess combo

    Parameters
    ----------
    hess : ndarray of float
        (3*natom, 3*natom) non-mass-weighted hessian in atomic units [Eh/a0/a0]
    geom : ndarray of float
        (natom, 3) geometry in [a0] at which the hessian was computed
    mass : ndarray of float
        (natom,) atomic masses [u].

    Returns
    -------
    dict
       Returns a dictionary of QCAspect objects (fields: lbl, units, data, comment). Contents summarized below.

    .. _`table:vibaspectinfo`:

    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | key           | description (lbl & comment)                | units     | data (real/imaginary modes)                          |
    +===============+============================================+===========+======================================================+
    | omega         | frequency                                  | cm^-1     | nd.array(ndof) complex (real/imag)                   |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | q             | normal mode, normalized mass-weighted      | a0 u^1/2  | ndarray(ndof, ndof) float                            |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | w             | normal mode, un-mass-weighted              | a0        | ndarray(ndof, ndof) float                            |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | x             | normal mode, normalized un-mass-weighted   | a0        | ndarray(ndof, ndof) float                            |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | mu            | reduced mass                               | u         | ndarray(ndof) float (+/+)                            |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | f             | force constant                             | Eh/a0     | ndarray(ndof) float (+/-)                            |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    """
    vibinfo = {}
    if (mass.shape[0] == geom.shape[0] == (hess.shape[0] // 3) == (hess.shape[1] // 3)) and (geom.shape[1] == 3):
        pass
    else:
        raise AttributeError(
            """Dimension mismatch amoung mass ({}), geometry ({}) and Hessian ({})""".format(
                mass.shape, geom.shape, hess.shape))

    mw = sum(mass)
    com = np.array([0.0, 0.0, 0.0])
    for i, m in enumerate(mass):
        com[0] += (geom[i, 0] * mass[i]) / mw
        com[1] += (geom[i, 1] * mass[i]) / mw
        com[2] += (geom[i, 2] * mass[i]) / mw
    if sum(com) > 1.0e-10:
        raise ValueError("Error, molecule must be at the center of mass")

    natom = len(mass)

    # get idealized rotations/translations
    TRspace = _get_TR_space(mass, geom, space='TR', tol=LINEAR_A_TOL)
    nrt = TRspace.shape[0]

    # form rot/tran projector
    P = np.identity(3 * natom)
    for irt in TRspace:
        P -= np.outer(irt, irt)

    # form mass weighting matrix
    sqrtmmm = np.repeat(np.sqrt(mass), 3)
    sqrtmmminv = np.divide(1.0, sqrtmmm)

    # mass weight the hessian
    nmwhess = hess.copy()
    mwhess = np.einsum('i,ij,j->ij', sqrtmmminv, nmwhess, sqrtmmminv)

    # project out translation/rotation
    mwhess_proj = np.dot(P.T, mwhess).dot(P)

    force_constant_au, qL = np.linalg.eigh(mwhess_proj)
    idx = np.argsort(force_constant_au)
    #idx = idx[:(3*natom]
    # sort
    force_constant_au = force_constant_au[idx]
    qL = qL[:, idx]
    # remove zeros
    idx = [i for i, v in enumerate(force_constant_au) if np.abs(v) > 1.0e-5]
    force_constant_au = force_constant_au[idx]
    qL = qL[:, idx]
    #qL = _phase_cols_to_max_element(qL)

    # populate vibinfo
    vibinfo['f'] = force_constant_au
    uconv_km = physconst['hartree2J'] / (physconst['bohr2m'] * physconst['bohr2m'] * physconst['amu2kg'])
    uconv_cm = 1.0 / (2.0 * np.pi * physconst['c'] * 100.00)
    # uconv_cm_1 = np.sqrt(physconst['na'] * physconst['hartree2J'] * 1.0e19) / (2 * np.pi * physconst['c'] *
    #         physconst['bohr2angstroms'])
    freq_cm_1 = np.sqrt(force_constant_au * uconv_km) * uconv_cm
    vibinfo['omega'] = freq_cm_1

    vibinfo['q'] = qL
    wL = np.einsum('i,ij->ij', sqrtmmminv, qL)
    vibinfo['w'] = wL

    reduced_mass = np.divide(1.0, np.linalg.norm(wL, axis=0)**2)
    vibinfo['mu'] = reduced_mass

    xL = np.sqrt(reduced_mass) * wL
    vibinfo['x'] = xL

    return vibinfo
