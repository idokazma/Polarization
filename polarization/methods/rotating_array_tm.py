"""Rotating-frame dipole-polarizability solver. Port of
``src/methods/RotatingArray_2D_TM.m``.

The MATLAB version keeps a local ``GRTM`` Green's helper that duplicates
``scalar_green`` but with a different prefactor. We preserve that helper
exactly (it is part of the locked-in numerical convention) rather than
reusing ``polarization.green.scalar``.
"""

from __future__ import annotations

from typing import Tuple

import numpy as np
from scipy.special import hankel1, jv


def _GRTM(params, ro: np.ndarray, rp: np.ndarray, k0: float, Rotan: float, nindex: float) -> complex:
    """Approximate rotating-medium 2D scalar Green's function in TM.

    ro, rp : 2-element row vectors (x, y).
    """
    mu0 = params.mu0
    r = float(np.linalg.norm(ro - rp))
    krn = k0 * r * nindex
    omega = k0 * params.c
    Gst = -0.25 * omega * mu0 * hankel1(0, krn)
    # z . (rp x ro)  -- MATLAB:  diff(ro .* flip(rp)) = ro(1)*rp(2) - ro(2)*rp(1)
    z_d_rp_cross_ro = ro[0] * rp[1] - ro[1] * rp[0]
    return Gst * np.exp(1j * (k0 * k0) * Rotan * z_d_rp_cross_ro)


def _cyl_pol_tm(params, aCyl: float, k0: float) -> complex:
    """Thin-cylinder electric polarizability under TM illumination."""
    omega_mu0 = k0 * params.c * params.mu0
    n1 = params.n_out
    n2 = params.n_in
    k0n1a = k0 * n1 * aCyl
    k0n2a = k0 * n2 * aCyl
    j0_2 = jv(0, k0n2a)
    j1_2 = jv(1, k0n2a)
    Dnum = -jv(1, k0n1a) * j0_2 * n1 + j1_2 * jv(0, k0n1a) * n2
    Den = -hankel1(1, k0n1a) * j0_2 * n1 + j1_2 * hankel1(0, k0n1a) * n2
    minus_b0_TM = Dnum / Den
    return 4.0 * minus_b0_TM / omega_mu0


def rotating_array_2d_tm(params, E_inc_z=None) -> Tuple[np.ndarray, np.ndarray, complex]:
    """Returns (PVector, Efields, alphaTM)."""
    p = params
    epsilon_r1 = p.er_out
    epsilon_r2 = p.er_in
    Rotan = p.OMEGA / p.omega

    sca_x = np.atleast_1d(np.asarray(p.sca_x, dtype=float))
    sca_y = np.atleast_1d(np.asarray(p.sca_y, dtype=float))
    rarray = np.column_stack((sca_x, sca_y))

    Npoints = rarray.shape[0]
    k0 = p.k0
    nb = p.n_out

    alphaTM = _cyl_pol_tm(p, p.radius, k0)
    InvalphaTM = 1.0 / alphaTM

    MM = InvalphaTM * np.eye(Npoints, dtype=complex)
    rhsV = np.zeros(Npoints, dtype=complex)
    src = np.array([[p.source_loc_x], [p.source_loc_y]], dtype=float)

    # Import here to avoid circular import
    from polarization.green.scalar import scalar_green

    for ni in range(Npoints):
        ri = rarray[ni]
        rhsV[ni] = scalar_green(src, ri.reshape(2, 1), p, False)[0]
        for nj in range(Npoints):
            if nj == ni:
                continue
            rj = rarray[nj]
            MM[ni, nj] = -_GRTM(p, ri, rj, k0, Rotan, nb)

    PVector = np.linalg.solve(MM, rhsV)
    # Eq. 22: E = I * i / (pi r^2 omega delta_eps)
    factor = 1j / (np.pi * (p.radius ** 2) * p.omega * (epsilon_r2 - epsilon_r1) * p.e0)
    Efields = PVector * factor

    return PVector, Efields, alphaTM
