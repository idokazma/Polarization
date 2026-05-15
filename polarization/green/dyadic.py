"""Dyadic rotating-frame Green's function. Port of ``src/green/dyadic_green.m``."""

from __future__ import annotations

import numpy as np
from scipy.special import hankel1


def _coords(arr):
    arr = np.atleast_2d(np.asarray(arr))
    return arr[0, :], arr[1, :]


def dyadic_green(source, test, params, inside: bool):
    """Return (Gx, Gy) row vectors of magnetic dyadic Green's components."""
    if inside:
        n = params.n_in
        mu = params.mu0 * params.mr_in
    else:
        n = params.n_out
        mu = params.mu0 * params.mr_out

    c = params.c
    OMEGA = params.OMEGA
    k0 = params.k0
    omega = params.omega

    xs, ys = _coords(source)
    xt, yt = _coords(test)

    rot = np.exp(1j * k0 * OMEGA / c * (xs * yt - ys * xt))
    rr = np.sqrt((xt - xs) ** 2 + (yt - ys) ** 2)

    fact = 1.0 / (1j * omega * mu)
    fact_n = 1j * omega * OMEGA / (c * c)

    Gst = 1j / 4.0 * hankel1(0, k0 * n * rr)
    Gst_d = k0 * n * 1j / 4.0 * hankel1(1, k0 * n * rr)

    Gx_1 = ((ys - yt) / rr) * Gst_d
    Gx_2 = -(xt - xs) * Gst * fact_n

    Gy_1 = -((xs - xt) / rr) * Gst_d
    Gy_2 = -(yt - ys) * Gst * fact_n

    Gx = (Gx_1 + Gx_2) * rot * fact
    Gy = (Gy_1 + Gy_2) * rot * fact
    return Gx, Gy
