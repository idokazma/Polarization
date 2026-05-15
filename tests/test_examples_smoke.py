"""End-to-end smoke tests at minimal grid sizes."""

import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless

from polarization.examples.single_scatterer import single_scatterer
from polarization.examples.multiple_scatterers import multiple_scatterers


def test_single_scatterer_smoke():
    results = single_scatterer(
        overrides=dict(
            len_n=21, wid_n=21,
            omega_factors=np.array([-1e-5, 0.0, 1e-5]),
            shift_factors=np.array([1.0]),
            n_sources=3,
        ),
        plot=False,
    )
    # Shape is [ni=1][nj=3][nt=3]
    assert len(results) == 1
    assert len(results[0]) == 3
    assert len(results[0][0]) == 3
    r = results[0][0][0]
    assert r.pol.Pvec.size == 1
    assert r.filaments.mean_E.size == 1


def test_multiple_scatterers_smoke():
    results = multiple_scatterers(
        overrides=dict(
            len_n=51, wid_n=51,
            n_scatterers=3,
            omega_factors=np.array([-1e-5, 0.0, 1e-5]),
        ),
        plot=False,
    )
    assert len(results) == 1
    assert len(results[0]) == 3
    r = results[0][1]  # OMEGA=0 entry
    assert r.pol.Pvec.size == 3
    assert r.filaments.mean_E.size == 3
