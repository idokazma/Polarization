import numpy as np
from polarization.geometry.vogel import golden_angle_array


def test_golden_angle_shape():
    xy = golden_angle_array(100)
    assert xy.shape == (2, 100)


def test_golden_angle_first_point_at_radius_a():
    xy = golden_angle_array(1)
    expected_r = (1.0 / 1.6) * 10.0  # min_dist * (1/1.6)
    assert np.isclose(np.hypot(xy[0, 0], xy[1, 0]), expected_r)


def test_golden_angle_radii_monotonic():
    xy = golden_angle_array(50)
    r = np.hypot(xy[0], xy[1])
    # Radii must increase with sqrt(index)
    assert np.all(np.diff(r) > 0)
