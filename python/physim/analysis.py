"""Diagnostics computed from recorded trajectories."""

import numpy as np


def relative_energy_error(traj):
    """``(E(t) - E(0)) / |E(0)|`` for each recorded frame."""
    e = np.asarray(traj.energy)
    return (e - e[0]) / abs(e[0])
