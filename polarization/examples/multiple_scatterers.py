"""Vogel-spiral multi-scatterer comparison. Port of ``examples/multiple_scatterers.m``."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np

from polarization.config.parameters import SimParams, generate_parameters
from polarization.geometry.vogel import golden_angle_array
from polarization.sources.plane_wave import generate_plane_wave
from polarization.sources.line_source import generate_source_wave
from polarization.runners.eval_tm import eval_tm_results
from polarization.plotting.compare import plot_fil_mom_pol


@dataclass
class MultipleScattersConfig:
    eps_vec: float = 11.4
    radius_ratio: float = 2.0 / 100.0
    lambda_: float = 1e-6
    Iz: float = 1.0
    is_plane_wave: int = 0
    plane_wave_direction: int = 1
    len_n: int = 200
    wid_n: int = 200
    len: float = 0.5e-6
    wid: float = 0.5e-6
    n_scatterers: int = 15
    omega_factors: Sequence[float] = field(
        default_factory=lambda: 1e-1 * np.arange(-5e-4, 5e-4 + 1e-6, 5e-5)
    )
    shift_factors: Sequence[float] = field(default_factory=lambda: np.array([0.0]))


def multiple_scatterers(overrides: Optional[dict] = None, plot: bool = True):
    cfg = MultipleScattersConfig()
    if overrides:
        for k, v in overrides.items():
            if not hasattr(cfg, k):
                raise AttributeError(f"unknown multiple_scatterers override {k!r}")
            setattr(cfg, k, v)

    params = SimParams(
        te=0,
        tm=1,
        lambda_=cfg.lambda_,
        Iz=cfg.Iz,
        OMEGA=0.0,
        is_plane_wave=cfg.is_plane_wave,
        calc_full_sol=0,
        plane=0,
        hit_plane=0,
        len_n=cfg.len_n,
        wid_n=cfg.wid_n,
        len=cfg.len,
        wid=cfg.wid,
        plane_wave_direction=cfg.plane_wave_direction,
    )

    VogelArrayXY = golden_angle_array(max(cfg.n_scatterers, 100))
    params.sca_x = params.lambda_ * VogelArrayXY[0, : cfg.n_scatterers]
    params.sca_y = params.lambda_ * VogelArrayXY[1, : cfg.n_scatterers]
    print(np.max(np.hypot(params.sca_x, params.sca_y)))

    params = generate_parameters(params)
    omega_factors = np.asarray(cfg.omega_factors, dtype=float)
    params.OMEGA_vec = params.omega * omega_factors
    params.shift_vec = np.asarray(cfg.shift_factors, dtype=float) * params.lambda_

    print(f"OMEGA*rho_c/c        = {np.max(params.OMEGA_vec) * np.max(params.shift_vec) / 3e8:f}")
    print(
        f"OMEGA*max(rho_c)/c   = "
        f"{np.max(params.sca_x ** 2 + params.sca_y ** 2) * np.max(params.OMEGA_vec) / 3e8:f}"
    )

    ni = len(params.shift_vec)
    nj = len(params.OMEGA_vec)
    results: list = [[None for _ in range(nj)] for _ in range(ni)]

    for i in range(ni):
        params.er_in = cfg.eps_vec
        params.n_in = np.sqrt(params.mr_in * params.er_in)
        for j in range(nj):
            params.radius = cfg.radius_ratio * params.lambda_
            params.OMEGA = params.OMEGA_vec[j]
            params.sca_x = params.lambda_ * VogelArrayXY[0, : cfg.n_scatterers]
            params.sca_y = params.lambda_ * VogelArrayXY[1, : cfg.n_scatterers]
            params.source_loc_x = 0.0
            params.source_loc_y = 0.0

            if params.is_plane_wave:
                params = generate_plane_wave(params)
            else:
                params = generate_source_wave(params)

            results[i][j] = eval_tm_results(params)

    if plot:
        plot_fil_mom_pol(results, show=True)
    return results
