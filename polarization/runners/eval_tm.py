"""Per-configuration orchestrator. Port of ``src/runners/eval_TM_results.m``."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from polarization.methods.mom import mom, MoMResult
from polarization.methods.rotating_array_tm import rotating_array_2d_tm
from polarization.methods.mie_tm import mie_series_tm
from polarization.methods.filaments_tm import filaments_tm_multiple, FilamentResult


@dataclass
class MoMOut:
    ok: bool = False
    E_sol: object = None
    mean_E: Optional[np.ndarray] = None
    mean_I: Optional[np.ndarray] = None
    alpha: Optional[np.ndarray] = None
    mean_E0: Optional[np.ndarray] = None


@dataclass
class PolOut:
    Pvec: np.ndarray = field(default_factory=lambda: np.zeros(0))
    E_sol: np.ndarray = field(default_factory=lambda: np.zeros(0))
    alpha: complex = 0j


@dataclass
class FilOut:
    E_sol: object = None
    mean_E: Optional[np.ndarray] = None
    mean_I: Optional[np.ndarray] = None
    alpha: complex = 0j
    field_current_mat: Optional[np.ndarray] = None


@dataclass
class MieOut:
    E_sol: object = None
    mean_E: complex = 0j


@dataclass
class TMResult:
    params: object = None
    mom: MoMOut = field(default_factory=MoMOut)
    pol: PolOut = field(default_factory=PolOut)
    filaments: FilOut = field(default_factory=FilOut)
    mie: Optional[MieOut] = None


def eval_tm_results(params) -> TMResult:
    """Run MoM, RotatingArray polarizability, Filaments (and Mie when valid)."""
    r = TMResult()

    # MoM
    try:
        mom_r: MoMResult = mom(params)
        r.mom = MoMOut(
            ok=True,
            E_sol=mom_r.SOL,
            mean_E=mom_r.mean_E,
            mean_I=mom_r.mean_I,
            alpha=mom_r.alpha_mom,
            mean_E0=mom_r.mean_E0,
        )
        params.radius = mom_r.effective_radius
    except Exception as exc:  # pragma: no cover - matches MATLAB warning behaviour
        import warnings
        warnings.warn(f"MoM failed: {exc}", RuntimeWarning)
        r.mom = MoMOut(ok=False)

    # Polarizability theory
    Pvec, E_sol_pol, alpha_pol = rotating_array_2d_tm(params, params.E_inc_z)
    r.pol = PolOut(Pvec=Pvec, E_sol=E_sol_pol, alpha=alpha_pol)

    # Mie (only single scatterer at origin with plane wave)
    sca_x = np.atleast_1d(params.sca_x)
    sca_y = np.atleast_1d(params.sca_y)
    if (
        params.is_plane_wave == 1
        and sca_x.size == 1
        and float(sca_x[0]) == 0.0
        and float(sca_y[0]) == 0.0
    ):
        full_fields, mean_E_mie = mie_series_tm(params, params.E_inc_z)
        r.mie = MieOut(E_sol=full_fields, mean_E=np.mean(mean_E_mie))
        print("done MIE")

    # Filaments
    fil: FilamentResult = filaments_tm_multiple(params)
    r.filaments = FilOut(
        E_sol=fil.full_fields,
        mean_E=fil.mean_E,
        mean_I=fil.mean_I,
        alpha=fil.alpha,
        field_current_mat=fil.output_mat,
    )

    r.params = deepcopy(params)
    return r
