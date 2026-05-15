"""Alpha-tensor regression and presentation plots.

Port of the live code path of ``src/plotting/TM_alpha_presentation_fields.m``.
The MATLAB script is 656 lines of which the bulk is duplicate / commented-out
scratch plotting. We keep the same regression formulae and produce one figure
per alpha component (angle + magnitude in a 3x3 grid).

The single-scatterer scenario produces a results cube indexed [i][j][t]:
    i = shift index (rho_c sweep)
    j = OMEGA index
    t = source position index
The MATLAB regression collapses ``t`` into the source sweep used by linsolve,
so the port preserves that.
"""

from __future__ import annotations

from typing import List, Sequence

import numpy as np
import matplotlib.pyplot as plt


def _extract_field_currents(results_ijt) -> tuple:
    """Pull the per-scenario quantities out of the filament output_mat.

    Returns arrays indexed by (rho, omega, source) for each filament-output
    component. Each result has output_mat[k] of shape (n_sca,) -- for the
    single-scatterer scenario n_sca == 1, so we squeeze."""
    ni = len(results_ijt)
    nj = len(results_ijt[0])
    nt = len(results_ijt[0][0])
    shape = (ni, nj, nt)

    Iz = np.zeros(shape, dtype=complex)
    Iez = np.zeros(shape, dtype=complex)
    Hx = np.zeros(shape, dtype=complex)
    Hy = np.zeros(shape, dtype=complex)
    Ez = np.zeros(shape, dtype=complex)
    Im_x = np.zeros(shape, dtype=complex)
    Im_y = np.zeros(shape, dtype=complex)
    Ez_sol = np.zeros(shape, dtype=complex)
    Hx_sol = np.zeros(shape, dtype=complex)
    Hy_sol = np.zeros(shape, dtype=complex)
    Ez_inc_mom = np.zeros(shape, dtype=complex)
    Ez_sol_mom = np.zeros(shape, dtype=complex)

    for i in range(ni):
        for j in range(nj):
            for t in range(nt):
                r = results_ijt[i][j][t]
                fc = r.filaments.field_current_mat
                Iz[i, j, t] = np.atleast_1d(fc[0])[0]
                Iez[i, j, t] = np.atleast_1d(fc[1])[0]
                Hx[i, j, t] = np.atleast_1d(fc[2])[0]
                Hy[i, j, t] = np.atleast_1d(fc[3])[0]
                Ez[i, j, t] = np.atleast_1d(fc[4])[0]
                Im_x[i, j, t] = np.atleast_1d(fc[5])[0]
                Im_y[i, j, t] = np.atleast_1d(fc[6])[0]
                Ez_sol[i, j, t] = np.atleast_1d(fc[7])[0]
                Hx_sol[i, j, t] = np.atleast_1d(fc[8])[0]
                Hy_sol[i, j, t] = np.atleast_1d(fc[9])[0]
                if r.mom.mean_E0 is not None:
                    Ez_inc_mom[i, j, t] = np.atleast_1d(r.mom.mean_E0)[0]
                if r.mom.mean_E is not None:
                    Ez_sol_mom[i, j, t] = np.atleast_1d(r.mom.mean_E)[0]
    return Iz, Iez, Hx, Hy, Ez, Im_x, Im_y, Ez_sol, Hx_sol, Hy_sol, Ez_inc_mom, Ez_sol_mom


def tm_alpha_presentation_fields(results, show: bool = True, savepath: str | None = None):
    """Regress the alpha tensor from the filament/MoM outputs and plot it.

    results : 3D list [i][j][t] of TMResult.
    """
    params = results[0][0][0].params
    Iz, Iez, Hx, Hy, Ez, Im_x, Im_y, Ez_sol, Hx_sol, Hy_sol, Ez_inc_mom, Ez_sol_mom = (
        _extract_field_currents(results)
    )

    ni, nj, nt = Ez.shape
    # Allocate 3x3 alpha tensors of shape (ni, nj)
    keys = [
        "ezez", "ezhx", "ezhy",
        "hxez", "hxhx", "hxhy",
        "hyez", "hyhx", "hyhy",
    ]
    alpha = {k: np.zeros((ni, nj), dtype=complex) for k in keys}
    alpha_mom = {k: np.zeros((ni, nj), dtype=complex) for k in ("ezez", "ezhx", "ezhy")}

    for rho_i in range(ni):
        for j in range(nj):
            # Linear regression over the source dimension t.
            # Rows are sources, columns are [Hx, Hy, Ez].
            A = np.column_stack((Hx[rho_i, j, :], Hy[rho_i, j, :], Ez[rho_i, j, :]))
            for label, target in (
                ("ez", Ez_sol[rho_i, j, :]),
                ("hx", Hx_sol[rho_i, j, :]),
                ("hy", Hy_sol[rho_i, j, :]),
            ):
                coef, *_ = np.linalg.lstsq(A, target, rcond=None)
                alpha[f"{label}hx"][rho_i, j] = coef[0]
                alpha[f"{label}hy"][rho_i, j] = coef[1]
                alpha[f"{label}ez"][rho_i, j] = coef[2]

            # MoM regression uses (Hx, Hy, Ez_inc_mom) -> Ez_sol_mom
            A_mom = np.column_stack(
                (Hx[rho_i, j, :], Hy[rho_i, j, :], Ez_inc_mom[rho_i, j, :])
            )
            coef, *_ = np.linalg.lstsq(A_mom, Ez_sol_mom[rho_i, j, :], rcond=None)
            alpha_mom["ezhx"][rho_i, j] = coef[0]
            alpha_mom["ezhy"][rho_i, j] = coef[1]
            alpha_mom["ezez"][rho_i, j] = coef[2]

    # Fix OMEGA=0 nan/zero by averaging neighbouring phases
    omega_vec = params.OMEGA_vec
    for j, omv in enumerate(omega_vec):
        if omv != 0.0:
            continue
        if 0 < j < nj - 1:
            for k in alpha:
                magnitude = np.abs(alpha[k][:, j])
                phase = 0.5 * (np.angle(alpha[k][:, j - 1]) + np.angle(alpha[k][:, j + 1]))
                alpha[k][:, j] = magnitude * np.exp(1j * phase)

    omega_norm = omega_vec / params.omega
    shift_vec = np.atleast_1d(params.shift_vec)

    # 3x3 angle figure
    fig_ang, axs = plt.subplots(3, 3, figsize=(12, 10))
    labels = {
        "ezez": r"$\sigma^{ee}_{zz}$",
        "ezhx": r"$\sigma^{em}_{z\rho}$",
        "ezhy": r"$\sigma^{em}_{z\theta}$",
        "hxez": r"$\sigma^{me}_{\rho z}$",
        "hxhx": r"$\sigma^{mm}_{\rho\rho}$",
        "hxhy": r"$\sigma^{mm}_{\rho\theta}$",
        "hyez": r"$\sigma^{me}_{\theta z}$",
        "hyhx": r"$\sigma^{mm}_{\theta\rho}$",
        "hyhy": r"$\sigma^{mm}_{\theta\theta}$",
    }
    layout = [
        ("ezez", (0, 0)), ("ezhx", (0, 1)), ("ezhy", (0, 2)),
        ("hxez", (1, 0)), ("hxhx", (1, 1)), ("hxhy", (1, 2)),
        ("hyez", (2, 0)), ("hyhx", (2, 1)), ("hyhy", (2, 2)),
    ]
    for key, (r, c) in layout:
        ax = axs[r, c]
        for rho_i in range(ni):
            ax.plot(omega_norm, np.angle(alpha[key][rho_i]), "x-",
                    label=f"{shift_vec[rho_i]/params.lambda_:.0f}" if c == 2 and r == 0 else None)
        ax.set_title(rf"$\angle$ {labels[key]}")
        ax.set_xlabel(r"$\bar{\Omega}$")
        ax.grid(True)
    fig_ang.suptitle("Alpha tensor — angle")
    fig_ang.tight_layout()

    # 3x3 magnitude figure
    fig_mag, axs2 = plt.subplots(3, 3, figsize=(12, 10))
    for key, (r, c) in layout:
        ax = axs2[r, c]
        for rho_i in range(ni):
            ax.plot(omega_norm, np.abs(alpha[key][rho_i]), "x-")
        ax.set_title(rf"$|{labels[key].strip('$')}|$")
        ax.set_xlabel(r"$\bar{\Omega}$")
        ax.grid(True)
    fig_mag.suptitle("Alpha tensor — magnitude")
    fig_mag.tight_layout()

    figs = [fig_ang, fig_mag]

    if savepath is not None:
        from pathlib import Path
        outdir = Path(savepath)
        outdir.mkdir(parents=True, exist_ok=True)
        fig_ang.savefig(outdir / "alpha_angle.png", dpi=120, bbox_inches="tight")
        fig_mag.savefig(outdir / "alpha_magnitude.png", dpi=120, bbox_inches="tight")

    if show:
        plt.show()
    return figs, alpha, alpha_mom
