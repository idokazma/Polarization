"""Scalar rotating-frame Green's function.

Port of ``src/green/scalar_green.m`` and ``src/green/scalar_green_mom.m``.
The two variants share their math; ``scalar_green_mom`` accepts the medium
parameters explicitly so MoM can be called without a full ``SimParams``.

Convention (preserved from MATLAB):
    source, test : (2, N) arrays of (x, y) coordinates.
        Broadcasting follows MATLAB: when both are (2, N) of the same N, the
        function returns an N-vector of pointwise Green's values; when
        ``source`` is (2, 1) and ``test`` is (2, N) (or vice versa), it
        broadcasts.
"""

from __future__ import annotations

import numpy as np
from scipy.special import hankel1


def _coords(arr):
    """Return (x_row, y_row) row vectors from a (2, N) array."""
    arr = np.atleast_2d(np.asarray(arr))
    return arr[0, :], arr[1, :]


def scalar_green(source, test, params, inside: bool) -> np.ndarray:
    """Rotating 2D scalar Green's function for the params struct."""
    if inside:
        n = params.n_in
    else:
        n = params.n_out
    return scalar_green_mom(source, test, n, params.k0, params.c, params.OMEGA)


def scalar_green_mom(source, test, n, k0, c, OMEGA, _unused=None) -> np.ndarray:
    """Variant taking explicit (n, k0, c, OMEGA)."""
    xs, ys = _coords(source)
    xt, yt = _coords(test)

    rot = np.exp(1j * k0 * OMEGA / c * (xs * yt - ys * xt))
    rr = np.sqrt((xt - xs) ** 2 + (yt - ys) ** 2)

    Gst = 1j / 4.0 * hankel1(0, k0 * n * rr)
    return Gst * rot
