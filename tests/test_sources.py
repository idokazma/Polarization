import numpy as np

from polarization.config.parameters import SimParams, generate_parameters
from polarization.sources.plane_wave import generate_plane_wave
from polarization.sources.line_source import generate_source_wave


def _params(grid=21):
    p = SimParams(lambda_=1e-6, len=1e-6, wid=1e-6, len_n=grid, wid_n=grid, tm=1, te=0)
    return generate_parameters(p)


def test_plane_wave_tm_E_amplitude_is_unity():
    p = _params()
    p = generate_plane_wave(p)
    assert np.allclose(np.abs(p.E_inc_z), 1.0, rtol=1e-14)


def test_plane_wave_h_y_magnitude_matches_intrinsic_admittance():
    p = _params()
    p = generate_plane_wave(p)
    # |H_y| = k0 / (omega mu0 mr_out)
    expected = p.k0 / (p.omega * p.mu0 * p.mr_out)
    assert np.allclose(np.abs(p.H_inc_y), expected, rtol=1e-12)


def test_line_source_E_finite_off_origin():
    p = _params()
    p.source_loc_x = 5e-6
    p.source_loc_y = 5e-6
    p = generate_source_wave(p)
    assert p.E_inc_z.shape == (21, 21)
    assert np.all(np.isfinite(p.E_inc_z))
    assert np.all(np.isfinite(p.H_inc_x))
    assert np.all(np.isfinite(p.H_inc_y))
