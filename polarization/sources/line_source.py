"""Line-current incident field. Port of ``src/sources/generate_source_wave.m``."""

from __future__ import annotations

import numpy as np

from polarization.config.parameters import SimParams
from polarization.green.scalar import scalar_green
from polarization.green.dyadic import dyadic_green


def generate_source_wave(params: SimParams) -> SimParams:
    p = params
    src = np.array([[p.source_loc_x], [p.source_loc_y]], dtype=float)
    test = np.vstack((p.X.ravel(order="F"), p.Y.ravel(order="F")))

    if p.tm:
        E_inc_z = scalar_green(src, test, p, False)
        p.E_inc_z = E_inc_z.reshape(p.X.shape, order="F")
        p.E_inc_x = np.zeros_like(p.E_inc_z)
        p.E_inc_y = np.zeros_like(p.E_inc_z)

        Hx, Hy = dyadic_green(src, test, p, False)
        p.H_inc_x = Hx.reshape(p.X.shape, order="F")
        p.H_inc_y = Hy.reshape(p.X.shape, order="F")
        p.H_inc_z = np.zeros_like(p.H_inc_y)

    if p.te:
        H_inc_z = scalar_green(src, test, p, False)
        p.H_inc_z = H_inc_z.reshape(p.X.shape, order="F")
        p.H_inc_x = np.zeros_like(p.H_inc_z)
        p.H_inc_y = np.zeros_like(p.H_inc_z)

        Ex, Ey = dyadic_green(src, test, p, False)
        scale = 1j / (p.omega * (p.e0 * p.er_out))
        p.E_inc_x = (scale * Ex).reshape(p.X.shape, order="F")
        p.E_inc_y = (scale * Ey).reshape(p.X.shape, order="F")
        p.E_inc_z = np.zeros_like(p.E_inc_y)

    return p
