import numpy as np
import pytest

from vec.boards import Board
from vec.evaluate import Stage

N_GENES = 60
GENES = [f"g{i}" for i in range(N_GENES)]


def simulate_stage(n_cells, means_by_type, props, rng, coords_scale=None, coords_shift=None):
    """Log-normalised Poisson counts from per-type mean profiles, plus an optional 3D blob per cell type."""
    types = rng.choice(len(props), size=n_cells, p=props)
    counts = rng.poisson(means_by_type[types])
    lib = counts.sum(1, keepdims=True).clip(min=1)
    X = np.log1p(1e4 * counts / lib).astype(np.float32)
    coords = None
    if coords_scale is not None:
        centres = coords_shift if coords_shift is not None else np.zeros((len(props), 3))
        coords = (centres[types] + rng.normal(size=(n_cells, 3)) * coords_scale).astype(np.float32)
    return Stage(X=X, coords=coords, celltype=np.array([f"type{t}" for t in types]))


@pytest.fixture(scope="session")
def rng():
    return np.random.default_rng(0)


@pytest.fixture(scope="session")
def profiles(rng):
    base = rng.gamma(2.0, 2.0, size=(3, N_GENES))
    later = base.copy()
    later[:, :8] *= 4.0        # genes going up
    later[:, 8:16] *= 0.2      # genes going down
    return base, later


def make_board(task, needs_coords):
    return Board(key=f"{task}:synthetic", task=task, genes=GENES, needs_coords=needs_coords, min_cells=10,
                 max_cells=100000)


@pytest.fixture(scope="session")
def boards():
    return {t: make_board(t, needs_coords=(t != "T1")) for t in ("T1", "T2", "T3")}
