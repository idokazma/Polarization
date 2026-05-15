"""Mie series for a single dielectric cylinder under TM plane-wave excitation.

Port of ``src/methods/Mie_Series_TM.m``. Valid only for the single-scatterer
case with the cylinder at the origin and plane-wave excitation; the runner
gates on these conditions before calling it.
"""

from __future__ import annotations

from typing import Tuple

import numpy as np
from scipy.special import hankel1, jv


def mie_series_tm(params, E_inc_z: np.ndarray) -> Tuple[Tuple[np.ndarray, np.ndarray, np.ndarray], complex]:
    p = params

    # cart2pol -> theta, rho
    theta = np.arctan2(p.Y, p.X)
    rho = np.hypot(p.X, p.Y)

    # Field at the origin
    centre_mask = (p.X == 0.0) & (p.Y == 0.0)
    E0 = E_inc_z[centre_mask]
    if E0.size == 0:
        # Fallback if the grid doesn't sample the origin exactly
        ix = np.argmin(np.abs(p.x))
        iy = np.argmin(np.abs(p.y))
        E0 = np.array([E_inc_z[iy, ix]])
    E_0 = E0[0]

    n_1 = p.n_out
    n_2 = p.n_in
    k0 = p.k0
    radius = p.radius

    etha_1 = np.sqrt(p.mr_out / p.er_out)
    etha_2 = np.sqrt(p.mr_in / p.er_in)

    m = np.arange(p.max_m + 1)

    # b_TM
    b_up_A = etha_2 * 0.5 * (jv(m - 1, k0 * n_1 * radius) - jv(m + 1, k0 * n_1 * radius)) * jv(m, k0 * n_2 * radius)
    b_up_B = etha_1 * 0.5 * (jv(m - 1, k0 * n_2 * radius) - jv(m + 1, k0 * n_2 * radius)) * jv(m, k0 * n_1 * radius)
    b_dn_C = etha_2 * 0.5 * (hankel1(m - 1, k0 * n_1 * radius) - hankel1(m + 1, k0 * n_1 * radius)) * jv(m, k0 * n_2 * radius)
    b_dn_D = etha_1 * 0.5 * (jv(m - 1, k0 * n_2 * radius) - jv(m + 1, k0 * n_2 * radius)) * hankel1(m, k0 * n_1 * radius)

    b_TM = ((-1.0) * (-1j) ** m) * (b_up_A - b_up_B) / (b_dn_C - b_dn_D)
    c_TM = b_TM * hankel1(m, k0 * n_1 * radius) / jv(m, k0 * n_2 * radius) \
           + ((-1j) ** m) * jv(m, k0 * n_1 * radius) / jv(m, k0 * n_2 * radius)

    sigma = np.full_like(m, 2.0, dtype=float)
    sigma[0] = 1.0

    out_mask = rho > radius
    in_mask = rho <= radius

    E_mid = np.zeros_like(p.X, dtype=complex)
    H_rho_mid = np.zeros_like(p.X, dtype=complex)
    H_theta_mid = np.zeros_like(p.X, dtype=complex)

    for i, mi in enumerate(m):
        # outside
        if np.any(out_mask):
            E_mid[out_mask] += sigma[i] * b_TM[i] * hankel1(mi, k0 * n_1 * rho[out_mask]) * np.cos(mi * theta[out_mask])
            H_rho_mid[out_mask] += sigma[i] * b_TM[i] * hankel1(mi, k0 * n_1 * rho[out_mask]) \
                * (-mi * np.sin(mi * theta[out_mask])) / rho[out_mask]
            H_theta_mid[out_mask] += sigma[i] * b_TM[i] \
                * 0.5 * k0 * n_1 * (hankel1(mi - 1, k0 * n_1 * rho[out_mask]) - hankel1(mi + 1, k0 * n_1 * rho[out_mask])) \
                * np.cos(mi * theta[out_mask])
        # inside
        if np.any(in_mask):
            E_mid[in_mask] += sigma[i] * c_TM[i] * jv(mi, k0 * n_2 * rho[in_mask]) * np.cos(mi * theta[in_mask])
            # MATLAB uses rho/(rho) safely because at rho=0 sin(mi*theta)=0; we guard against div-by-zero
            rho_in = rho[in_mask]
            safe = rho_in.copy()
            safe[safe == 0.0] = 1.0
            H_rho_mid[in_mask] += sigma[i] * c_TM[i] * jv(mi, k0 * n_2 * rho_in) \
                * (-mi * np.sin(mi * theta[in_mask])) / safe
            H_theta_mid[in_mask] += sigma[i] * c_TM[i] \
                * 0.5 * k0 * n_2 * (jv(mi - 1, k0 * n_2 * rho_in) - jv(mi + 1, k0 * n_2 * rho_in)) \
                * np.cos(mi * theta[in_mask])

    H_rho_mid[in_mask] = H_rho_mid[in_mask] / (1j * p.omega * p.mr_in * p.mu0) * E_0
    H_rho_mid[out_mask] = H_rho_mid[out_mask] / (1j * p.omega * p.mr_out * p.mu0) * E_0
    H_theta_mid[in_mask] = H_theta_mid[in_mask] / (1j * p.omega * p.mr_in * p.mu0) * E_0
    H_theta_mid[out_mask] = H_theta_mid[out_mask] / (1j * p.omega * p.mr_out * p.mu0) * E_0
    H_theta_mid = -H_theta_mid

    H_x_mid = np.cos(theta) * H_rho_mid - np.sin(theta) * H_theta_mid
    H_y_mid = np.sin(theta) * H_rho_mid + np.cos(theta) * H_theta_mid

    H_x_final = H_x_mid.copy()
    H_x_final[out_mask] = H_x_final[out_mask] + p.H_inc_x[out_mask]
    H_y_final = H_y_mid.copy()
    H_y_final[out_mask] = H_y_final[out_mask] + p.H_inc_y[out_mask]

    E_final = np.zeros_like(p.X, dtype=complex)
    E_final[out_mask] = E_0 * E_mid[out_mask] + E_inc_z[out_mask]
    E_final[in_mask] = E_0 * E_mid[in_mask]

    full_fields = (E_final, H_x_final, H_y_final)
    mean_E_mie = E_final[in_mask].mean() if np.any(in_mask) else 0.0 + 0.0j
    return full_fields, mean_E_mie
