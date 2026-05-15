import numpy as np
from polarization.config.parameters import SimParams, generate_parameters


def test_derived_constants_match_matlab_values():
    p = SimParams(lambda_=1e-6, len=0.5e-6, wid=0.5e-6, len_n=31, wid_n=31)
    p = generate_parameters(p)
    # c = 1/sqrt(e0*mu0), check to within numerical tolerance
    assert np.isclose(p.c, 2.99792458e8, rtol=1e-5)
    assert np.isclose(p.k0, 2 * np.pi / p.lambda_, rtol=1e-14)
    assert np.isclose(p.omega, p.f * 2 * np.pi, rtol=1e-14)
    assert np.isclose(p.n_in, np.sqrt(p.mr_in * p.er_in), rtol=1e-14)
    assert np.isclose(p.n_out, np.sqrt(p.mr_out * p.er_out), rtol=1e-14)


def test_grid_centred_and_sized():
    p = SimParams(lambda_=1e-6, len=0.5e-6, wid=0.5e-6, len_n=21, wid_n=21)
    p = generate_parameters(p)
    assert p.x.size == 21
    assert p.y.size == 21
    assert np.isclose(p.x[0], -0.25e-6)
    assert np.isclose(p.x[-1], 0.25e-6)
    assert p.X.shape == (21, 21)
    assert p.Y.shape == (21, 21)


def test_filament_defaults():
    p = SimParams(len=1e-6, wid=1e-6, len_n=21, wid_n=21)
    p = generate_parameters(p)
    assert p.R_out == 1.2
    assert p.R_in == 0.8
    assert p.N_filaments == 30
    assert p.N_testpoints == 60
    assert p.max_m == 50


def test_with_overrides():
    p = SimParams()
    q = p.with_overrides(lambda_=8e-7, er_in=4.0)
    assert q.lambda_ == 8e-7
    assert q.er_in == 4.0
    assert p.lambda_ == 1e-6, "with_overrides must not mutate the original"
