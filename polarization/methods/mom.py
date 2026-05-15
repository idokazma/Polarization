"""Method of Moments solver. Port of the live path of ``src/methods/MoM.m``.

The MATLAB file is 374 lines of which ~half is commented-out scaffolding. We
port the active code path:

    [SOL, mean_E, mean_I, effective_radius, alpha_mom, mean_E0] = MoM(params)

The ``params.calc_full_sol`` branch (full-grid reconstruction) is implemented
but not yet exercised by any example.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np
from scipy.special import hankel1

from polarization.green.scalar import scalar_green_mom


class MoMResult(NamedTuple):
    SOL: object                 # -1 unless calc_full_sol, else 2D field
    mean_E: np.ndarray
    mean_I: np.ndarray
    effective_radius: float
    alpha_mom: complex
    mean_E0: np.ndarray


def mom(params) -> MoMResult:
    p = params

    X_mom, Y_mom = np.meshgrid(p.x_mom, p.y_mom)

    delta_eps = p.er_in - p.er_out
    loc = np.flatnonzero((np.sqrt(X_mom ** 2 + Y_mom ** 2).ravel(order="F")) < p.radius)
    Xm_flat = X_mom.ravel(order="F")
    Ym_flat = Y_mom.ravel(order="F")

    epsilon_bg = p.er_out
    lambda_ = p.lambda_
    k0 = 2.0 * np.pi / lambda_
    beta = 1.0

    dx = abs(p.x_mom[1] - p.x_mom[0]) if p.x_mom.size > 1 else 0.0
    dy = abs(p.y_mom[1] - p.y_mom[0]) if p.y_mom.size > 1 else 0.0
    V_i = dx * dy
    if len(loc) == 0:
        raise RuntimeError("MoM: no grid cells inside the scatterer; reduce wid_n/len_n or increase radius")
    effective_radius = float(np.sqrt(V_i * len(loc) / np.pi))

    R_i_eff = float(np.sqrt(V_i / np.pi))
    gamma = R_i_eff / k0 * hankel1(1, R_i_eff * k0) + 1j * 2.0 / (np.pi * k0 ** 2)

    # Build scatterer point cloud across all centres
    sca_x_arr = np.atleast_1d(np.asarray(p.sca_x, dtype=float))
    sca_y_arr = np.atleast_1d(np.asarray(p.sca_y, dtype=float))
    n_centres = sca_x_arr.size
    scatterers_x = np.concatenate([Xm_flat[loc] + sca_x_arr[pp] for pp in range(n_centres)])
    scatterers_y = np.concatenate([Ym_flat[loc] + sca_y_arr[pp] for pp in range(n_centres)])

    src_x = p.source_loc_x
    src_y = p.source_loc_y

    if not p.is_plane_wave:
        source_xy = np.array([[src_x], [src_y]], dtype=float)
        test_xy = np.vstack((scatterers_x, scatterers_y))
        E_bg_st = scalar_green_mom(source_xy, test_xy, p.n_out, k0, p.c, p.OMEGA)
    else:
        E_bg_st = 1.0 * np.exp(-1j * k0 * scatterers_x)

    # Solve linear system
    N = scatterers_x.size
    G_mat_dy = np.empty((N, N), dtype=complex)
    for ii in range(N):
        src_xy = np.array([[scatterers_x[ii]], [scatterers_y[ii]]], dtype=float)
        test_xy = np.vstack((scatterers_x, scatterers_y))
        G_mat_dy[:, ii] = scalar_green_mom(src_xy, test_xy, p.n_out, k0, p.c, p.OMEGA)
    G_mat_dy = -G_mat_dy * V_i * k0 ** 2 * delta_eps

    for ii in range(N):
        G_mat_dy[ii, ii] = 1.0 - (1j * np.pi / 2.0 * beta * gamma) * k0 ** 2 * delta_eps

    solution_mat_st = np.linalg.solve(G_mat_dy, E_bg_st)

    L = len(loc)
    mean_E = np.zeros(n_centres, dtype=complex)
    mean_E0 = np.zeros(n_centres, dtype=complex)
    mean_I = np.zeros(n_centres, dtype=complex)

    for ll in range(n_centres):
        sl = slice(ll * L, (ll + 1) * L)
        mean_E[ll] = solution_mat_st[sl].mean()
        mean_E0[ll] = E_bg_st[sl].mean()
        mean_I[ll] = mean_E[ll] * (-1j * p.omega * p.e0 * (p.er_in - p.er_out) * V_i * L)

    # The MATLAB code overwrites alpha_mom with the same mean over all
    # scatterer cells; preserve that behaviour.
    alpha_mom = mean_I / np.mean(E_bg_st)
    # Then collapses to a scalar via implicit broadcast in plot_fil_mom_pol;
    # leave as array for general n_centres.

    print(f"num of scatterers {scatterers_x.size}")

    SOL = -1
    if p.calc_full_sol == 1:
        # Optional full-grid reconstruction. Not exercised by any of the
        # bundled examples, but we provide a faithful port.
        X = p.X
        Y = p.Y
        if not p.is_plane_wave:
            grid_src = np.array([[src_x], [src_y]], dtype=float)
            grid_test = np.vstack((X.ravel(order="F"), Y.ravel(order="F")))
            SOL = scalar_green_mom(grid_src, grid_test, p.n_out, k0, p.c, p.OMEGA).reshape(
                X.shape, order="F"
            )
        else:
            SOL = np.exp(-1j * k0 * X)
        for i in range(N):
            src_xy = np.array([[scatterers_x[i]], [scatterers_y[i]]], dtype=float)
            test_xy = np.vstack((X.ravel(order="F"), Y.ravel(order="F")))
            add = solution_mat_st[i] * scalar_green_mom(
                src_xy, test_xy, p.n_out, k0, p.c, p.OMEGA
            ) * V_i * k0 ** 2 * delta_eps
            SOL = SOL + add.reshape(X.shape, order="F")

    return MoMResult(
        SOL=SOL,
        mean_E=mean_E,
        mean_I=mean_I,
        effective_radius=effective_radius,
        alpha_mom=alpha_mom,
        mean_E0=mean_E0,
    )
