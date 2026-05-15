"""Spectral sweep over a rotating block-crystal array.

Port of ``examples/block_crystals.m``. Requires ``data/arraypoints.mat`` (a
MATLAB .mat file containing an ``ArrayPoints`` variable shaped (N, 2)).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

from polarization.config.parameters import SimParams, generate_parameters
from polarization.green.scalar import scalar_green
from polarization.methods.rotating_array_tm import rotating_array_2d_tm
from polarization.sources.plane_wave import generate_plane_wave
from polarization.sources.line_source import generate_source_wave


@dataclass
class BlockCrystalsConfig:
    eps_vec: float = 11.4
    radius_ratio: float = 1.0 / 10.0
    Iz: float = 1.0
    is_plane_wave: int = 0
    plane_wave_direction: int = 1
    len_n: int = 81
    wid_n: int = 81
    grid_extent_factor: float = 2.1
    array_scale: float = 0.75
    wavelengths: Sequence[float] = field(
        default_factory=lambda: np.arange(0.2, 3.8 + 1e-6, 0.005) * 1e-6
    )
    observed_points: np.ndarray = field(
        default_factory=lambda: np.array([[-3.0, 0.0, 3.0], [5.0, 5.0, 5.0]]) * 1e-6
    )
    source_distance: float = 3e-6
    OMEGA: float = 0.0


def block_crystals(overrides: Optional[dict] = None, plot: bool = True, data_path: Optional[Path] = None):
    cfg = BlockCrystalsConfig()
    if overrides:
        for k, v in overrides.items():
            if not hasattr(cfg, k):
                raise AttributeError(f"unknown block_crystals override {k!r}")
            setattr(cfg, k, v)

    if data_path is None:
        data_path = Path(__file__).resolve().parents[2] / "data" / "arraypoints.mat"
    data_path = Path(data_path)
    if not data_path.is_file():
        raise FileNotFoundError(
            f"Required data file not found: {data_path}\n"
            "See data/README.md for how to obtain or generate this file."
        )

    loaded = loadmat(str(data_path))
    ArrayPoints = loaded["ArrayPoints"]

    params = SimParams(
        te=0,
        tm=1,
        Iz=cfg.Iz,
        multiscatterer=1,
        is_plane_wave=cfg.is_plane_wave,
        calc_full_sol=0,
        plane=0,
        hit_plane=0,
        len_n=cfg.len_n,
        wid_n=cfg.wid_n,
        len=cfg.radius_ratio * 1e-6 * cfg.grid_extent_factor,
        wid=cfg.radius_ratio * 1e-6 * cfg.grid_extent_factor,
        plane_wave_direction=cfg.plane_wave_direction,
    )
    params.sca_x = ArrayPoints[:, 0] * 1e-6 * cfg.array_scale
    params.sca_y = ArrayPoints[:, 1] * 1e-6 * cfg.array_scale
    params.OMEGA_vec = np.array([cfg.OMEGA])
    params.shift_vec = np.array([0.0])
    params.OMEGA = cfg.OMEGA

    # Single source at angle 2*pi*3/4 (MATLAB sources(end-1) of 5 evenly-spaced angles)
    phase = 2.0 * np.pi * np.linspace(0.0, 1.0, 5)[-2]
    source = cfg.source_distance * np.exp(1j * phase)
    sources = np.array([source])

    observed_point = np.asarray(cfg.observed_points)
    n_obs = observed_point.shape[1]
    nw = len(cfg.wavelengths)
    E_sol = np.zeros((n_obs, nw), dtype=complex)
    E_hit = np.zeros((n_obs, nw), dtype=complex)

    for w, lam in enumerate(cfg.wavelengths):
        params.lambda_ = float(lam)
        params = generate_parameters(params)
        for t, s in enumerate(sources):
            params.source_loc_x = float(np.real(s))
            params.source_loc_y = float(np.imag(s))

            params.er_in = cfg.eps_vec
            params.n_in = np.sqrt(params.mr_in * params.er_in)
            params.radius = cfg.radius_ratio * 1e-6
            params.OMEGA = cfg.OMEGA

            print(f"Wavelength {w+1}/{nw}")

            if params.is_plane_wave:
                params = generate_plane_wave(params)
            else:
                params = generate_source_wave(params)

            Pvec, _, alpha_pol = rotating_array_2d_tm(params, params.E_inc_z)
            print(np.abs(alpha_pol))

            src_xy = np.array([[params.source_loc_x], [params.source_loc_y]], dtype=float)
            for pp in range(n_obs):
                test = observed_point[:, pp:pp + 1]
                E_hit[pp, w] = scalar_green(src_xy, test, params, False)[0]
                for tt in range(len(Pvec)):
                    sx = np.array([[params.sca_x[tt]], [params.sca_y[tt]]], dtype=float)
                    E_sol[pp, w] = E_sol[pp, w] - params.omega * params.mu0 / 4.0 / (1j / 4.0) * Pvec[tt] \
                        * scalar_green(sx, test, params, False)[0]

    fig = None
    if plot:
        fig, ax = plt.subplots()
        ax.plot(cfg.wavelengths, 10.0 * np.log10(np.abs(E_hit + E_sol) / np.abs(E_hit)).T)
        ax.set_xlabel(r"$\lambda$ [m]")
        ax.set_ylabel(r"$|E_{total}|/|E_{hit}|$ [dB]")
        ax.grid(True, which="both")
        plt.show()
    return E_sol, E_hit, fig
