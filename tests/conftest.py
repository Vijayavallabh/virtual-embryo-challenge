import os
import subprocess

import numpy as np
import pytest

from vec import provenance
from vec.boards import Board
from vec.evaluate import Stage
from vec.io import write_submission

N_GENES = 60
GENES = [f"g{i}" for i in range(N_GENES)]


def git(repo, *args):
    """git in a test's own repository, with no global or system configuration; returns its standard output."""
    env = {**os.environ, "GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1"}
    return subprocess.run(["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@example.org",
                           "-c", "commit.gpgsign=false", *args], check=True, capture_output=True, text=True,
                          env=env).stdout.strip()


def make_repo(path):
    """A git repository at `path` with one commit (method.py) and a clean tree."""
    path.mkdir(parents=True, exist_ok=True)
    git(path, "init", "-q")
    (path / "method.py").write_text("x = 1\n")
    git(path, "add", "method.py")
    git(path, "commit", "-qm", "init")
    return path


def build(adata, path, repo):
    """write_submission as one build in a new process: the build's start is recorded for `repo` (as the first import of
    vec.provenance records it for this package's repository), then the file is written and stamped against `repo`."""
    provenance._record_start(repo)
    return write_submission(adata, path, repo=repo)


@pytest.fixture
def own_repo(tmp_path_factory, monkeypatch):
    """A clean git repository in place of this package's (the default a builder's stamps are written against), with
    this test's build started there, so that no test stamps against the real checkout."""
    repo = make_repo(tmp_path_factory.mktemp("repo"))
    monkeypatch.setattr(provenance, "REPO_ROOT", repo)
    provenance._record_start(repo)
    return repo


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
