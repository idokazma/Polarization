"""Vogel (golden-angle) spiral generator. Port of ``src/geometry/GA_generator.m``."""

from __future__ import annotations

import numpy as np


def golden_angle_array(N: int = 100) -> np.ndarray:
    """Return a 2xN array of (x, y) coordinates on a Vogel spiral."""
    golden_ratio_sq = 1.5 + 0.5 * np.sqrt(5.0)
    alpha_GA = 2.0 * np.pi / golden_ratio_sq
    min_dist = 10.0
    a = (1.0 / 1.6) * min_dist

    ip = np.arange(1, N + 1, dtype=float)
    rp = np.sqrt(ip) * a
    alphap = alpha_GA * ip
    xy = np.vstack((rp * np.cos(alphap), rp * np.sin(alphap)))
    return xy
