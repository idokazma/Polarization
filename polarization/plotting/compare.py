"""Comparison plots across the OMEGA sweep. Port of ``plot_fil_mom_pol.m``."""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt


def plot_fil_mom_pol(results, show: bool = True, savepath: str | None = None):
    """Compare MoM / Filaments / Polarizability results across the OMEGA sweep.

    results : 2D list-of-lists [i][j] of TMResult; j is the OMEGA dimension.
    """
    ni = len(results)
    nj = len(results[0])
    params_final = results[0][0].params

    # Use mean_E from each method (one row per scatterer centre)
    n_centres = np.atleast_1d(params_final.sca_x).size
    E_MoM = np.zeros((n_centres, nj), dtype=complex)
    E_fil = np.zeros((n_centres, nj), dtype=complex)
    E_pol = np.zeros((n_centres, nj), dtype=complex)
    OMEGA_vals = np.zeros(nj)

    for i in range(ni):
        for j in range(nj):
            r = results[i][j]
            if r.mom.mean_E is not None:
                E_MoM[:, j] = np.atleast_1d(r.mom.mean_E)
            E_fil[:, j] = np.atleast_1d(r.filaments.mean_E)
            E_pol[:, j] = np.atleast_1d(r.pol.E_sol).ravel()[:n_centres] \
                if np.atleast_1d(r.pol.E_sol).size >= n_centres else np.atleast_1d(r.pol.E_sol)
            OMEGA_vals[j] = r.params.OMEGA_vec[j] / r.params.omega

    ttl = (rf"$\lambda$ = {params_final.lambda_*1e6:g} $\mu$m, "
           rf"Radius = {params_final.radius/params_final.lambda_:.4f}$\lambda$, "
           rf"$\epsilon$ = {params_final.er_in:g}")

    def panel(op, label):
        fig, ax = plt.subplots()
        for ll in range(E_MoM.shape[0]):
            ax.plot(OMEGA_vals, op(E_MoM[ll]), "x-", label="MoM" if ll == 0 else None)
            ax.plot(OMEGA_vals, op(E_fil[ll]), "s-", label="Fil" if ll == 0 else None)
            ax.plot(OMEGA_vals, op(E_pol[ll]), "d-", label="Pol" if ll == 0 else None)
        ax.set_title(f"{label}(E) vs $\\Omega/\\omega$, {ttl}")
        ax.set_xlabel(r"$\Omega/\omega$")
        ax.grid(True, which="both")
        ax.legend()
        return fig

    figs = [
        panel(np.abs, "abs"),
        panel(np.real, "real"),
        panel(np.imag, "imag"),
    ]

    # Scatterer array figure
    fig_arr, ax_arr = plt.subplots()
    sx = np.atleast_1d(params_final.sca_x) * 1e6
    sy = np.atleast_1d(params_final.sca_y) * 1e6
    sc = ax_arr.scatter(sx, sy, c=np.abs(E_MoM[:, 0]), s=200, cmap="jet")
    plt.colorbar(sc, ax=ax_arr)
    ax_arr.set_title("Scatterers Array")
    ax_arr.set_xlabel(r"x [$\mu$m]")
    ax_arr.set_ylabel(r"y [$\mu$m]")
    ax_arr.set_aspect("equal")
    ax_arr.grid(True)
    figs.append(fig_arr)

    if savepath is not None:
        from pathlib import Path
        outdir = Path(savepath)
        outdir.mkdir(parents=True, exist_ok=True)
        for k, f in enumerate(figs):
            f.savefig(outdir / f"compare_{k:02d}.png", dpi=120, bbox_inches="tight")
    if show:
        plt.show()
    return figs
