"""Tests for the Green's function modules."""

import numpy as np
import pytest

from polarization.config.parameters import SimParams, generate_parameters
from polarization.green.scalar import scalar_green, scalar_green_mom
from polarization.green.dyadic import dyadic_green


def _make_params():
    p = SimParams(len=1e-6, wid=1e-6, len_n=21, wid_n=21, OMEGA=0.0)
    return generate_parameters(p)


def test_scalar_green_real_positive_distance_matches_hankel():
    p = _make_params()
    src = np.array([[0.0], [0.0]])
    test = np.array([[1e-7], [0.0]])
    G = scalar_green(src, test, p, False)
    # Without rotation, the Green's function is 1j/4 * H_0^{(1)}(k0 n r)
    from scipy.special import hankel1
    r = 1e-7
    expected = 1j / 4 * hankel1(0, p.k0 * p.n_out * r)
    assert np.allclose(G, expected, rtol=1e-12)


def test_scalar_green_rotation_unitary():
    """For zero distance to swap of source/test, the rotation factor inverts."""
    p = _make_params()
    p.OMEGA = 1e9
    src = np.array([[1e-7], [0.0]])
    test = np.array([[0.0], [1e-7]])
    G_fwd = scalar_green(src, test, p, False)
    G_rev = scalar_green(test, src, p, False)
    # |G| is the same, phase differs by the rotation factor sign
    assert np.allclose(np.abs(G_fwd), np.abs(G_rev), rtol=1e-12)


def test_scalar_green_and_scalar_green_mom_agree():
    p = _make_params()
    src = np.array([[0.0], [0.0]])
    test = np.array([[2e-7, 3e-7], [1e-7, -1e-7]])
    G1 = scalar_green(src, test, p, False)
    G2 = scalar_green_mom(src, test, p.n_out, p.k0, p.c, p.OMEGA)
    assert np.allclose(G1, G2, rtol=1e-14)


def test_dyadic_green_returns_two_arrays():
    p = _make_params()
    src = np.array([[0.0], [0.0]])
    test = np.array([[2e-7, 3e-7], [1e-7, -1e-7]])
    Gx, Gy = dyadic_green(src, test, p, False)
    assert Gx.shape == (2,)
    assert Gy.shape == (2,)
    assert np.all(np.isfinite(Gx))
    assert np.all(np.isfinite(Gy))
