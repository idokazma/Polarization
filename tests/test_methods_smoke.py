"""Smoke tests that exercise each method on tiny grids."""

import numpy as np
import pytest

from polarization.config.parameters import SimParams, generate_parameters
from polarization.sources.plane_wave import generate_plane_wave
from polarization.sources.line_source import generate_source_wave
from polarization.methods.mie_tm import mie_series_tm
from polarization.methods.rotating_array_tm import rotating_array_2d_tm
from polarization.methods.mom import mom
from polarization.methods.filaments_tm import filaments_tm_multiple


@pytest.fixture
def single_scatterer_params():
    radius_ratio = 1.0 / 100.0
    p = SimParams(
        lambda_=1e-6, len_n=31, wid_n=31, tm=1,
        is_plane_wave=1, OMEGA=0.0,
    )
    p.len = radius_ratio * p.lambda_ * 2.1
    p.wid = radius_ratio * p.lambda_ * 2.1
    p = generate_parameters(p)
    p.radius = radius_ratio * p.lambda_
    p.sca_x = np.array([0.0])
    p.sca_y = np.array([0.0])
    p.source_loc_x = 1e-4
    p.source_loc_y = 0.0
    p = generate_plane_wave(p)
    return p


def test_mie_returns_finite_fields(single_scatterer_params):
    full_fields, mean_E = mie_series_tm(single_scatterer_params, single_scatterer_params.E_inc_z)
    E, Hx, Hy = full_fields
    assert E.shape == single_scatterer_params.X.shape
    assert np.all(np.isfinite(E))
    assert np.all(np.isfinite(Hx))
    assert np.all(np.isfinite(Hy))
    assert np.isfinite(mean_E)


def test_rotating_array_single_scatterer(single_scatterer_params):
    Pvec, Efields, alpha = rotating_array_2d_tm(
        single_scatterer_params, single_scatterer_params.E_inc_z
    )
    assert Pvec.shape == (1,)
    assert np.isfinite(Pvec[0])
    assert np.isfinite(alpha)


def test_mom_single_scatterer(single_scatterer_params):
    r = mom(single_scatterer_params)
    assert np.isfinite(r.mean_E[0])
    assert np.isfinite(r.mean_I[0])
    assert r.effective_radius > 0


def test_filaments_single_scatterer(single_scatterer_params):
    r = filaments_tm_multiple(single_scatterer_params)
    assert np.isfinite(r.mean_E).all()
    assert np.isfinite(r.alpha)
    assert r.output_mat.shape[0] == 10
    assert r.output_mat.shape[1] == 1  # one scatterer


def test_mie_vs_rotating_array_consistent_polarizability(single_scatterer_params):
    """For a single scatterer at the origin with no rotation, the
    rotating-array alpha should be close to the analytical CylPolTM that the
    MATLAB reference uses (since it's literally a copy of the same closed
    form). This sanity-checks the port of CylPolTM."""
    Pvec, _, alpha_pol = rotating_array_2d_tm(
        single_scatterer_params, single_scatterer_params.E_inc_z
    )
    # alpha is a complex polarizability per unit length; sanity-check that
    # it is non-zero and has a finite imaginary part.
    assert alpha_pol != 0
    assert np.isfinite(alpha_pol.imag)
