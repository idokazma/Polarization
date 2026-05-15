"""Simulation parameter container.

Direct port of ``src/config/generate_parameters.m``. The MATLAB ``params``
struct is replaced by a mutable ``SimParams`` dataclass; ``generate_parameters``
fills in the derived fields (grid, physical constants, refractive indices,
filament defaults) exactly as the MATLAB version does.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np


@dataclass
class SimParams:
    # ---- Geometry / grid (caller sets) ----
    lambda_: float = 1e-6           # wavelength [m]
    len: float = 0.5e-6             # grid extent along y [m]
    wid: float = 0.5e-6             # grid extent along x [m]
    len_n: int = 200                # grid samples along y
    wid_n: int = 200                # grid samples along x

    # ---- Polarisation flags ----
    tm: int = 1
    te: int = 0

    # ---- Media (caller may override before generate_parameters) ----
    mr_in: float = 1.0
    mr_out: float = 1.0
    er_in: float = 11.4
    er_out: float = 1.0

    # ---- Scatterer geometry ----
    radius: float = 0.0
    sca_x: np.ndarray = field(default_factory=lambda: np.zeros(1))
    sca_y: np.ndarray = field(default_factory=lambda: np.zeros(1))

    # ---- Source ----
    Iz: float = 1.0
    is_plane_wave: int = 0
    plane_wave_direction: int = 1
    source_loc_x: float = 0.0
    source_loc_y: float = 0.0

    # ---- Rotation ----
    OMEGA: float = 0.0

    # ---- Sweeps (filled by example scripts) ----
    OMEGA_vec: np.ndarray = field(default_factory=lambda: np.zeros(1))
    shift_vec: np.ndarray = field(default_factory=lambda: np.zeros(1))

    # ---- Solver flags ----
    multiscatterer: int = 0
    calc_full_sol: int = 0
    plane: int = 0
    hit_plane: int = 0

    # ---- Derived (filled by generate_parameters) ----
    x: np.ndarray = field(default_factory=lambda: np.zeros(0))
    y: np.ndarray = field(default_factory=lambda: np.zeros(0))
    X: np.ndarray = field(default_factory=lambda: np.zeros(0))
    Y: np.ndarray = field(default_factory=lambda: np.zeros(0))
    x_mom: np.ndarray = field(default_factory=lambda: np.zeros(0))
    y_mom: np.ndarray = field(default_factory=lambda: np.zeros(0))
    e0: float = 0.0
    mu0: float = 0.0
    c: float = 0.0
    f: float = 0.0
    k0: float = 0.0
    omega: float = 0.0
    n_in: float = 0.0
    n_out: float = 0.0
    R_in: float = 0.8
    R_out: float = 1.2
    N_filaments: int = 30
    N_testpoints: int = 60
    max_m: int = 50

    # ---- Incident fields (filled by sources/*) ----
    E_inc_x: Optional[np.ndarray] = None
    E_inc_y: Optional[np.ndarray] = None
    E_inc_z: Optional[np.ndarray] = None
    H_inc_x: Optional[np.ndarray] = None
    H_inc_y: Optional[np.ndarray] = None
    H_inc_z: Optional[np.ndarray] = None

    # Used by single-scatterer Mie comparison
    direction: int = 0

    # ---- Helpers --------------------------------------------------------

    def with_overrides(self, **overrides) -> "SimParams":
        """Return a copy with overrides applied (analogue of merge_defaults)."""
        from dataclasses import replace
        valid = {k for k in self.__dataclass_fields__}
        clean = {}
        for k, v in overrides.items():
            if k not in valid:
                raise AttributeError(f"SimParams has no field {k!r}")
            clean[k] = v
        return replace(self, **clean)


def generate_parameters(params: SimParams) -> SimParams:
    """Fill the derived simulation constants, X/Y grids, refractive indices,
    and filament defaults. Mirrors ``src/config/generate_parameters.m``.
    """
    p = params
    p.y = np.linspace(0.0, p.len, p.len_n) - p.len / 2.0
    p.x = np.linspace(0.0, p.wid, p.wid_n) - p.wid / 2.0
    p.X, p.Y = np.meshgrid(p.x, p.y)            # matches MATLAB meshgrid layout

    p.y_mom = p.y.copy()
    p.x_mom = p.x.copy()

    # Physical constants
    p.e0 = 8.8541878128e-12
    p.mu0 = 1.25663706212e-6
    p.c = 1.0 / np.sqrt(p.e0 * p.mu0)
    p.f = p.c / p.lambda_
    p.k0 = 2.0 * np.pi / p.lambda_
    p.omega = p.f * 2.0 * np.pi

    p.n_in = np.sqrt(p.mr_in * p.er_in)
    p.n_out = np.sqrt(p.mr_out * p.er_out)

    # Filament method defaults
    p.R_out = 1.2
    p.R_in = 0.8
    p.N_filaments = 30
    p.N_testpoints = 60
    p.max_m = 50

    # Coerce scatterer coordinates to numpy arrays so downstream code can
    # treat single and multiple scatterer cases uniformly.
    p.sca_x = np.atleast_1d(np.asarray(p.sca_x, dtype=float))
    p.sca_y = np.atleast_1d(np.asarray(p.sca_y, dtype=float))

    return p
