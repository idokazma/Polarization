"""Single rotating scatterer comparison. Port of ``examples/single_scatterer.m``.

Sweeps over rotation rates, lateral shifts, and source positions; runs MoM,
Filaments, and Polarizability theory for each configuration; hands the result
cube to the alpha-tensor presentation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np

from polarization.config.parameters import SimParams, generate_parameters
from polarization.sources.plane_wave import generate_plane_wave
from polarization.sources.line_source import generate_source_wave
from polarization.runners.eval_tm import eval_tm_results
from polarization.plotting.alpha_fields import tm_alpha_presentation_fields


@dataclass
class SingleScattererConfig:
    eps_vec: float = 11.4
    radius_ratio: float = 1.0 / 100.0
    lambda_: float = 1e-6
    Iz: float = 1.0
    is_plane_wave: int = 0
    plane_wave_direction: int = 1
    len_n: int = 61
    wid_n: int = 61
    grid_extent_factor: float = 2.1
    omega_factors: Sequence[float] = field(
        default_factory=lambda: np.arange(-3e-5, 3e-5 + 1e-6, 1e-5)
    )
    shift_factors: Sequence[float] = field(
        default_factory=lambda: np.array([1, 2, 4, 10, 20, 50, 100, 200])
    )
    source_distance_lambda: float = 200.0
    n_sources: int = 9


def single_scatterer(overrides: Optional[dict] = None, plot: bool = True):
    cfg = SingleScattererConfig()
    if overrides:
        for k, v in overrides.items():
            if not hasattr(cfg, k):
                raise AttributeError(f"unknown single_scatterer override {k!r}")
            setattr(cfg, k, v)

    params = SimParams(
        te=0,
        tm=1,
        lambda_=cfg.lambda_,
        Iz=cfg.Iz,
        OMEGA=0.0,
        multiscatterer=0,
        is_plane_wave=cfg.is_plane_wave,
        calc_full_sol=0,
        plane=0,
        hit_plane=0,
        len_n=cfg.len_n,
        wid_n=cfg.wid_n,
        len=cfg.radius_ratio * cfg.lambda_ * cfg.grid_extent_factor,
        wid=cfg.radius_ratio * cfg.lambda_ * cfg.grid_extent_factor,
        plane_wave_direction=cfg.plane_wave_direction,
    )
    params = generate_parameters(params)
    omega_factors = np.asarray(cfg.omega_factors, dtype=float)
    params.OMEGA_vec = params.omega * omega_factors
    params.shift_vec = np.asarray(cfg.shift_factors, dtype=float) * params.lambda_

    phase = 2.0 * np.pi * np.linspace(0.0, 1.0, cfg.n_sources + 1)[:-1]
    sources = cfg.source_distance_lambda * params.lambda_ * np.exp(1j * phase)

    ni = len(params.shift_vec)
    nj = len(params.OMEGA_vec)
    nt = len(sources)
    results: list = [[[None for _ in range(nt)] for _ in range(nj)] for _ in range(ni)]

    counter = 0
    total = ni * nj * nt
    for t in range(nt):
        for i in range(ni):
            params.er_in = cfg.eps_vec
            params.n_in = np.sqrt(params.mr_in * params.er_in)
            for j in range(nj):
                params.radius = cfg.radius_ratio * params.lambda_
                params.OMEGA = params.OMEGA_vec[j]
                params.sca_x = np.array([float(params.shift_vec[i])])
                params.sca_y = np.array([0.0])
                params.source_loc_x = float(np.real(sources[t]) + params.sca_x[0])
                params.source_loc_y = float(np.imag(sources[t]) + params.sca_y[0])

                counter += 1
                print(
                    f"Calculating... Overall: {100*counter/total:.1f}% "
                    f"Shift: {100*(i+1)/ni:.1f}%, "
                    f"Omega: {100*(j+1)/nj:.1f}%, "
                    f"Source: {100*(t+1)/nt:.1f}%."
                )

                if params.is_plane_wave:
                    params = generate_plane_wave(params)
                else:
                    params = generate_source_wave(params)

                results[i][j][t] = eval_tm_results(params)

    if plot:
        tm_alpha_presentation_fields(results, show=True)
    return results
