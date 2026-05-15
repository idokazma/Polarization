"""Plane-wave incident field. Port of ``src/sources/generate_plane_wave.m``.

The MATLAB source builds ``Esym``/``Hsym`` via the Symbolic Math Toolbox and
``curl``; the closed-form expressions immediately above are identical, so we
drop the symbolic blocks entirely (they are only read by code paths that are
gated with ``if (false)`` or commented out in the MATLAB).
"""

from __future__ import annotations

import numpy as np

from polarization.config.parameters import SimParams


def generate_plane_wave(params: SimParams) -> SimParams:
    p = params
    if p.tm:
        p.E_inc_z = np.exp(-1j * p.k0 * p.X)
        p.E_inc_x = np.zeros_like(p.E_inc_z)
        p.E_inc_y = np.zeros_like(p.E_inc_z)

        p.H_inc_y = (1j / (p.omega * (p.mu0 * p.mr_out))) * (
            -1j * p.k0 * np.exp(-1j * p.k0 * p.X)
        )
        p.H_inc_x = np.zeros_like(p.H_inc_y)
        p.H_inc_z = np.zeros_like(p.H_inc_y)

    if p.te:
        if p.plane_wave_direction == 1:
            p.H_inc_z = np.exp(-1j * p.k0 * p.X)
            p.H_inc_x = np.zeros_like(p.H_inc_z)
            p.H_inc_y = np.zeros_like(p.H_inc_z)

            p.E_inc_y = (-1j / (p.omega * (p.e0 * p.er_out))) * (
                -1j * p.k0 * np.exp(-1j * p.k0 * p.X)
            )
            p.E_inc_x = np.zeros_like(p.E_inc_y)
            p.E_inc_z = np.zeros_like(p.E_inc_y)
        else:
            p.H_inc_z = np.exp(-1j * p.k0 * p.Y)
            p.H_inc_x = np.zeros_like(p.H_inc_z)
            p.H_inc_y = np.zeros_like(p.H_inc_z)

            p.E_inc_x = (1j / (p.omega * (p.e0 * p.er_out))) * (
                -1j * p.k0 * np.exp(-1j * p.k0 * p.Y)
            )
            p.E_inc_y = np.zeros_like(p.E_inc_x)
            p.E_inc_z = np.zeros_like(p.E_inc_x)

    return p
