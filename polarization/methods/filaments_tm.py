"""Filament method for the TM polarization across one or more cylinders.

Port of the live code path of ``src/methods/filaments_TM_multiple.m``.
Symbolic-toolbox / ``hit_plane != 0`` / ``calc_full_sol`` branches are
omitted because no example exercises them.

The MATLAB function returns:

    [full_fields, mean_E_fil, mean_I_fil, alpha_Fil, output_mat]

with ``output_mat`` a 1x10 cell array of per-scatterer scalars used by the
alpha-tensor presentation plots. We return a NamedTuple to keep the layout
explicit and the call site readable.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np

from polarization.green.scalar import scalar_green
from polarization.green.dyadic import dyadic_green


class FilamentResult(NamedTuple):
    full_fields: tuple                  # (E, Hx, Hy) sentinels unless calc_full_sol
    mean_E: np.ndarray
    mean_I: np.ndarray
    alpha: complex
    output_mat: np.ndarray              # shape (10,) of per-scatterer scalars (or arrays)


def filaments_tm_multiple(params) -> FilamentResult:
    p = params

    sca_x = np.atleast_1d(np.asarray(p.sca_x, dtype=float))
    sca_y = np.atleast_1d(np.asarray(p.sca_y, dtype=float))
    n_sca = sca_x.size
    R = p.radius

    # Build per-scatterer filaments + test points.
    # filaments_total columns: x, y, scatterer-index (0-based), inside-flag (0=outside, 1=inside)
    # test_point_total columns: x, y, scatterer-index, nhat_x, nhat_y, nhat_z
    filaments_blocks = []
    test_blocks = []

    theta = np.linspace(0.0, 2.0 * np.pi, p.N_filaments, endpoint=False)
    theta_test = np.linspace(0.0, 2.0 * np.pi, p.N_testpoints, endpoint=False)

    for k in range(n_sca):
        out_xy = np.vstack((
            p.R_out * R * np.cos(theta) + sca_x[k],
            p.R_out * R * np.sin(theta) + sca_y[k],
        ))
        in_xy = np.vstack((
            p.R_in * R * np.cos(theta) + sca_x[k],
            p.R_in * R * np.sin(theta) + sca_y[k],
        ))
        # inside filaments first (matches MATLAB stacking order)
        n_in = in_xy.shape[1]
        n_out = out_xy.shape[1]
        block_in = np.column_stack((
            in_xy[0], in_xy[1],
            np.full(n_in, k), np.ones(n_in),
        ))
        block_out = np.column_stack((
            out_xy[0], out_xy[1],
            np.full(n_out, k), np.zeros(n_out),
        ))
        filaments_blocks.append(block_in)
        filaments_blocks.append(block_out)

        # test points are on the cylinder surface, centred at the cylinder
        test_xy_local = np.vstack((R * np.cos(theta_test), R * np.sin(theta_test)))
        # normal vectors at each test point (unit radial)
        nhat_x = test_xy_local[0] / np.linalg.norm(test_xy_local, axis=0)
        nhat_y = test_xy_local[1] / np.linalg.norm(test_xy_local, axis=0)
        nhat_z = np.zeros_like(nhat_x)
        test_xy = test_xy_local + np.array([[sca_x[k]], [sca_y[k]]])
        block_test = np.column_stack((
            test_xy[0], test_xy[1],
            np.full(theta_test.size, k),
            nhat_x, nhat_y, nhat_z,
        ))
        test_blocks.append(block_test)

    filaments_total = np.vstack(filaments_blocks)
    test_point_total = np.vstack(test_blocks)
    Nt = test_point_total.shape[0]
    Nf = filaments_total.shape[0]

    # Incident fields at the test points
    test_xy_only = test_point_total[:, :2].T  # (2, Nt)

    if p.is_plane_wave:
        # Inline the closed forms used in generate_plane_wave so we don't
        # depend on a symbolic toolbox.
        if p.tm:
            E_inc_z = np.exp(-1j * p.k0 * test_point_total[:, 0])
            H_inc_x = np.zeros(Nt, dtype=complex)
            H_inc_y = (1j / (p.omega * (p.mu0 * p.mr_out))) * (
                -1j * p.k0 * np.exp(-1j * p.k0 * test_point_total[:, 0])
            )
        else:
            raise NotImplementedError("Filaments TE plane-wave path not ported")
    else:
        src = np.array([[p.source_loc_x], [p.source_loc_y]], dtype=float)
        E_inc_z = scalar_green(src, test_xy_only, p, False)
        H_inc_x, H_inc_y = dyadic_green(src, test_xy_only, p, False)

    # Assemble Mx, MHz
    Mx = np.zeros((Nt, Nf), dtype=complex)
    MHz = np.zeros((Nt, Nf), dtype=complex)
    Vex = np.zeros(Nt, dtype=complex)
    Vhz = np.zeros(Nt, dtype=complex)

    fil_sca_idx = filaments_total[:, 2].astype(int)
    fil_inside_flag = filaments_total[:, 3].astype(int)
    test_sca_idx = test_point_total[:, 2].astype(int)

    for i in range(Nt):
        nhat = test_point_total[i, 3:6]
        test_xy = test_point_total[i, :2].reshape(2, 1)

        # case 1: same scatterer, filament outside -> radiate inwards
        mask = (fil_sca_idx == test_sca_idx[i]) & (fil_inside_flag == 0)
        if np.any(mask):
            f_xy = filaments_total[mask, :2].T
            Mx[i, mask] = -scalar_green(f_xy, test_xy, p, True)
            Gx, Gy = dyadic_green(f_xy, test_xy, p, True)
            # cross(nhat, [Gx, Gy, 0]) z-component:  nhat_x*Gy - nhat_y*Gx
            MHz[i, mask] = -(nhat[0] * Gy - nhat[1] * Gx)

        # case 2: same scatterer, filament inside -> radiate outwards
        mask = (fil_sca_idx == test_sca_idx[i]) & (fil_inside_flag == 1)
        if np.any(mask):
            f_xy = filaments_total[mask, :2].T
            Mx[i, mask] = scalar_green(f_xy, test_xy, p, False)
            Gx, Gy = dyadic_green(f_xy, test_xy, p, False)
            MHz[i, mask] = (nhat[0] * Gy - nhat[1] * Gx)

        # case 3: other scatterer, filament inside
        mask = (fil_sca_idx != test_sca_idx[i]) & (fil_inside_flag == 1)
        if np.any(mask):
            f_xy = filaments_total[mask, :2].T
            Mx[i, mask] = scalar_green(f_xy, test_xy, p, False)
            Gx, Gy = dyadic_green(f_xy, test_xy, p, False)
            MHz[i, mask] = (nhat[0] * Gy - nhat[1] * Gx)

        Vex[i] = -E_inc_z[i]
        Vhz[i] = -(nhat[0] * H_inc_y[i] - nhat[1] * H_inc_x[i])

    M_total = np.vstack((Mx, MHz))
    V_total = np.concatenate((Vex, Vhz))
    I_solution, *_ = np.linalg.lstsq(M_total, V_total, rcond=None)

    # Compute field inside each scatterer on the grid points falling inside
    X = p.X
    Y = p.Y
    grid_inside = (X ** 2 + Y ** 2) <= 0.999999 * R * R
    n_in_pts = int(grid_inside.sum())
    E_sol_final = np.zeros((n_sca, n_in_pts), dtype=complex)
    Hx_sol_final = np.zeros((n_sca, n_in_pts), dtype=complex)
    Hy_sol_final = np.zeros((n_sca, n_in_pts), dtype=complex)

    X_in_template = X[grid_inside]
    Y_in_template = Y[grid_inside]

    # Inside fields come from filaments outside the scatterer (flag=0)
    for tt in range(Nf):
        if fil_inside_flag[tt] != 0:
            continue
        s = fil_sca_idx[tt]
        in_xy = np.vstack((X_in_template + sca_x[s], Y_in_template + sca_y[s]))
        src_xy = filaments_total[tt, :2].reshape(2, 1)
        E_sol_final[s] += I_solution[tt] * scalar_green(src_xy, in_xy, p, True)
        Hx, Hy = dyadic_green(src_xy, in_xy, p, True)
        Hx_sol_final[s] += I_solution[tt] * Hx
        Hy_sol_final[s] += I_solution[tt] * Hy

    mean_E_fil = E_sol_final.mean(axis=1)
    mean_I_fil_E_only = mean_E_fil * (
        -1j * p.omega * p.e0 * (p.er_in - p.er_out) * np.pi * R * R
    )

    # Compute alpha_Fil based on the incident field at the inside grid points
    # of the FIRST scatterer (the MATLAB code reuses the last in_loc_line).
    if n_sca >= 1:
        last_s = n_sca - 1
        in_xy = np.vstack((X_in_template + sca_x[last_s], Y_in_template + sca_y[last_s]))
        src_xy = np.array([[p.source_loc_x], [p.source_loc_y]], dtype=float)
        hit_field_inside = scalar_green(src_xy, in_xy, p, False)
    else:
        hit_field_inside = np.array([1.0 + 0j])
    alpha_Fil = np.mean(mean_I_fil_E_only) / np.mean(hit_field_inside)

    # Compute the incident fields inside each scatterer (used by output_mat)
    E_inc_z_inside = np.zeros((n_sca, n_in_pts), dtype=complex)
    H_inc_x_inside = np.zeros((n_sca, n_in_pts), dtype=complex)
    H_inc_y_inside = np.zeros((n_sca, n_in_pts), dtype=complex)
    src_xy = np.array([[p.source_loc_x], [p.source_loc_y]], dtype=float)
    for mm in range(n_sca):
        in_xy = np.vstack((X_in_template + sca_x[mm], Y_in_template + sca_y[mm]))
        E_inc_z_inside[mm] = scalar_green(src_xy, in_xy, p, False)
        Hx, Hy = dyadic_green(src_xy, in_xy, p, False)
        H_inc_x_inside[mm] = Hx
        H_inc_y_inside[mm] = Hy

    mean_Hx_inc = H_inc_x_inside.mean(axis=1)
    mean_Hy_inc = H_inc_y_inside.mean(axis=1)
    mean_Ez_inc = E_inc_z_inside.mean(axis=1)
    mean_Hx_fil = Hx_sol_final.mean(axis=1)
    mean_Hy_fil = Hy_sol_final.mean(axis=1)

    # H_current_add and the magnetic-current correction.
    # In MATLAB this branch only fires for single-scatterer arrays (the path
    # used by the alpha-tensor regression in TM_alpha_presentation_fields).
    if n_sca == 1:
        in_x = X_in_template + sca_x[0]
        in_y = Y_in_template + sca_y[0]
        final_current_H = -Hx_sol_final[0] * in_x - Hy_sol_final[0] * in_y
        rotation_current_H = -H_inc_x_inside[0] * in_x - H_inc_y_inside[0] * in_y
        H_current_add = final_current_H - rotation_current_H
    else:
        H_current_add = np.zeros(1, dtype=complex)

    mean_I_fil = (
        np.mean(E_sol_final, axis=1)
        * (-1j * p.omega * p.e0 * (p.er_in - p.er_out) * np.pi * R * R)
        + 1j * p.omega / (p.c ** 2) * p.OMEGA * np.mean(H_current_add) * np.pi * R * R
    )

    mean_Im_fil_x = (
        -1j * p.omega * p.OMEGA / (p.c ** 2) * sca_x * (np.mean(E_sol_final, axis=1) - mean_Ez_inc) * np.pi * R * R
    )
    mean_Im_fil_y = (
        -1j * p.omega * p.OMEGA / (p.c ** 2) * sca_y * (np.mean(E_sol_final, axis=1) - mean_Ez_inc) * np.pi * R * R
    )

    output_mat = np.array(
        [
            mean_I_fil,
            mean_I_fil_E_only,
            mean_Hx_inc,
            mean_Hy_inc,
            mean_Ez_inc,
            mean_Im_fil_x,
            mean_Im_fil_y,
            mean_E_fil,
            mean_Hx_fil,
            mean_Hy_fil,
        ],
        dtype=object,
    )

    full_fields = (-1, -1, -1)
    return FilamentResult(
        full_fields=full_fields,
        mean_E=mean_E_fil,
        mean_I=mean_I_fil,
        alpha=alpha_Fil,
        output_mat=output_mat,
    )
