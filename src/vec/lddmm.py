"""Large-deformation diffeomorphic registration (LDDMM) by geodesic shooting, in plain torch float64, for row E53 of
experiments/harness/TRIES.md (as specified; branch campaign/moon-embryo-geodesic): the geodesic overlay form of
`vec.shape_interp` and its rotated-start placebo.

The deformation is a Gaussian-kernel velocity field carried by control points q with momenta p. Its geodesic is the
Hamiltonian flow of H(q, p) = 1/2 sum_ij p_i . p_j K(q_i, q_j), integrated with N_STEPS classical Runge-Kutta (RK4)
steps on [0, 1] from (q0, p0). The trajectory is stored at the 2 N_STEPS + 1 nodes 0, h/2, h, ..., 1 (h = 1 / N_STEPS;
a half node holds the mean of the RK4 step's second and third stage states), and the time-dependent field v_t(x) =
sum_j K(x, q_j(t)) p_j(t) reads the state at time t by linear interpolation between nodes. A point is carried from time
a to time b (either direction) by RK4 on dx/dt = v_t(x) with sub-steps cut at the step grid, so a full step reads the
stored node states exactly. The same stored trajectory carries points forward and back.

The fit: from p0 = 0, L-BFGS (strong Wolfe line search) on
    E(p0) = 1/2 <p0, K(q0, q0) p0> + GAMMA * data(phi_1(X), Y) / data(X, Y),
where X and Y are SUBSAMPLE cells of the earlier and the later bracket (stream [seed, SUBSAMPLE_STREAM]), q0 are
N_CONTROL control points chosen by k-means++ seeding (D^2 sampling) on X (stream [seed, CONTROL_STREAM]), and data is
the squared MMD between the two unit-mass empirical measures under the kernel sum over the position scales DATA_SCALES
of exp(-|x - y|^2 / (2 s^2)) exp(-|f - g|^2 / (2 b^2)), with f, g the cells' expression features and b the median
pairwise feature distance of the pooled subsamples. Every constant is fixed by the row; nothing is swept.

Inputs are in the earlier bracket's centred RMS-1 frame (`vec.t2_build.interpolate_clouds`), so SIGMA_V and the data
scales are in RMS units as given. torch is imported only when a fit runs (the gpu dependency group).
"""

from __future__ import annotations

import hashlib
import os
from dataclasses import dataclass, field

import numpy as np

SUBSAMPLE = 4000
N_CONTROL = 1500
SIGMA_V = 0.2
DATA_SCALES = (0.05, 0.1, 0.2, 0.4)
GAMMA = 10.0
N_STEPS = 10
LBFGS_ITERS = 150
SUBSAMPLE_STREAM = 51
CONTROL_STREAM = 52
CHUNK = 4096  # points carried per batch (bounds a Jacobian batch's autograd memory near 6 GB on the GPU)
CACHE_SIZE = 4
DEVICE_ENV = "VEC_LDDMM_DEVICE"


def _torch():
    import torch

    return torch


def device_name() -> str:
    """The fit's device: VEC_LDDMM_DEVICE if set, else cuda when torch sees a GPU, else cpu."""
    torch = _torch()
    want = os.environ.get(DEVICE_ENV)
    if want:
        return want
    return "cuda" if torch.cuda.is_available() else "cpu"


# ------------------------------------------------------------------------------------------------ numpy helpers
def subsample_rows(n: int, m: int, rng: np.random.Generator) -> np.ndarray:
    """min(m, n) distinct rows, sorted."""
    return np.sort(rng.choice(n, min(m, n), replace=False))


def kmeanspp(X, k: int, rng: np.random.Generator) -> np.ndarray:
    """k-means++ seeding (D^2 sampling) on the rows of X: k distinct row indices, in the order chosen. All rows when
    k >= len(X). A step where every remaining squared distance is zero draws uniformly among unchosen rows."""
    X = np.asarray(X, np.float64)
    n = len(X)
    if k >= n:
        return np.arange(n)
    chosen = np.empty(k, dtype=np.int64)
    chosen[0] = rng.integers(n)
    d2 = ((X - X[chosen[0]]) ** 2).sum(1)
    taken = np.zeros(n, dtype=bool)
    taken[chosen[0]] = True
    for i in range(1, k):
        w = np.where(taken, 0.0, d2)
        tot = float(w.sum())
        if tot > 0:
            j = int(rng.choice(n, p=w / tot))
        else:
            j = int(rng.choice(np.flatnonzero(~taken)))
        chosen[i] = j
        taken[j] = True
        d2 = np.minimum(d2, ((X - X[j]) ** 2).sum(1))
    return chosen


# ------------------------------------------------------------------------------------------------ torch kernels
def _sqdist(x, y):
    d = (x * x).sum(1)[:, None] + (y * y).sum(1)[None, :] - 2.0 * x @ y.T
    return d.clamp_min(0.0)


def _gauss(x, y, sigma: float):
    return (-_sqdist(x, y) / (2.0 * sigma * sigma)).exp()


def hamiltonian_rhs(q, p, sigma: float = SIGMA_V):
    """(dq/dt, dp/dt) of the geodesic equations for the Gaussian kernel: dq = K p, dp_i = sum_j (p_i . p_j) K_ij
    (q_i - q_j) / sigma^2."""
    K = _gauss(q, q, sigma)
    dq = K @ p
    W = K * (p @ p.T) / (sigma * sigma)
    dp = W.sum(1)[:, None] * q - W @ q
    return dq, dp


def shoot(q0, p0, sigma: float = SIGMA_V, n_steps: int = N_STEPS):
    """The geodesic from (q0, p0) by RK4: lists of q and p at the 2 n_steps + 1 nodes (times k h / 2)."""
    h = 1.0 / n_steps
    qs, ps = [q0], [p0]
    q, p = q0, p0
    for _ in range(n_steps):
        k1q, k1p = hamiltonian_rhs(q, p, sigma)
        q2, p2 = q + 0.5 * h * k1q, p + 0.5 * h * k1p
        k2q, k2p = hamiltonian_rhs(q2, p2, sigma)
        q3, p3 = q + 0.5 * h * k2q, p + 0.5 * h * k2p
        k3q, k3p = hamiltonian_rhs(q3, p3, sigma)
        q4, p4 = q + h * k3q, p + h * k3p
        k4q, k4p = hamiltonian_rhs(q4, p4, sigma)
        qs.append(0.5 * (q2 + q3))
        ps.append(0.5 * (p2 + p3))
        q = q + h / 6.0 * (k1q + 2 * k2q + 2 * k3q + k4q)
        p = p + h / 6.0 * (k1p + 2 * k2p + 2 * k3p + k4p)
        qs.append(q)
        ps.append(p)
    return qs, ps


def _state(qs, ps, tau: float, n_steps: int):
    """The control state at time tau, linear between the stored nodes."""
    half = 0.5 / n_steps
    pos = tau / half
    i = int(np.floor(pos + 1e-9))
    i = min(max(i, 0), len(qs) - 1)
    f = pos - i
    if abs(f) < 1e-9 or i == len(qs) - 1:
        return qs[i], ps[i]
    return (1 - f) * qs[i] + f * qs[i + 1], (1 - f) * ps[i] + f * ps[i + 1]


def _cuts(a: float, b: float, n_steps: int) -> list[float]:
    """a, the step grid's times strictly between a and b (in the direction of travel), b."""
    h = 1.0 / n_steps
    lo, hi = min(a, b), max(a, b)
    inner = [k * h for k in range(1, n_steps) if lo + 1e-12 < k * h < hi - 1e-12]
    if b < a:
        inner = inner[::-1]
    return [a, *inner, b]


def carry(x, qs, ps, a: float, b: float, sigma: float = SIGMA_V, n_steps: int = N_STEPS):
    """Points x (torch, n x 3) carried from time a to time b along the stored field (RK4, sub-steps at the grid)."""
    if a == b:
        return x
    cuts = _cuts(a, b, n_steps)
    for u, w in zip(cuts[:-1], cuts[1:], strict=True):
        d = w - u
        s0, sm, s1 = _state(qs, ps, u, n_steps), _state(qs, ps, 0.5 * (u + w), n_steps), _state(qs, ps, w, n_steps)

        def v(y, s):
            return _gauss(y, s[0], sigma) @ s[1]

        k1 = v(x, s0)
        k2 = v(x + 0.5 * d * k1, sm)
        k3 = v(x + 0.5 * d * k2, sm)
        k4 = v(x + d * k3, s1)
        x = x + d / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
    return x


def feature_bandwidth(F):
    """The median pairwise distance (i < j) of the rows of F (torch)."""
    torch = _torch()
    D = _sqdist(F, F).sqrt()
    mask = torch.ones(D.shape, dtype=torch.bool, device=D.device).triu(1)
    return D[mask].median()


class _Data:
    """The data term against a fixed target: sum over DATA_SCALES of the product-kernel squared MMD (unit masses)."""

    def __init__(self, X, Y, FX, FY, scales=DATA_SCALES):
        torch = _torch()
        self.Y, self.scales = Y, tuple(scales)
        F = torch.cat([FX, FY])
        b = feature_bandwidth(F)
        if not torch.isfinite(b) or float(b) <= 0:
            raise ValueError("the expression features have no spread (median pairwise distance 0)")
        self.bandwidth = float(b)
        self.kxx = _gauss(FX, FX, self.bandwidth)
        self.kxy = _gauss(FX, FY, self.bandwidth)
        dyy = _sqdist(Y, Y)
        kyy = _gauss(FY, FY, self.bandwidth)
        self.const = sum(((-dyy / (2 * s * s)).exp() * kyy).mean() for s in self.scales)

    def __call__(self, Xm):
        dxx, dxy = _sqdist(Xm, Xm), _sqdist(Xm, self.Y)
        out = self.const
        for s in self.scales:
            out = (out + ((-dxx / (2 * s * s)).exp() * self.kxx).mean()
                   - 2 * ((-dxy / (2 * s * s)).exp() * self.kxy).mean())
        return out


@dataclass
class GeodesicFit:
    """A fitted geodesic: control points q0, momenta p0 (numpy, float64), and the fit's flat facts (`record`)."""

    q0: np.ndarray
    p0: np.ndarray
    record: dict
    sigma: float = SIGMA_V
    n_steps: int = N_STEPS
    device: str = "cpu"
    _traj: tuple | None = field(default=None, repr=False)
    _memo: dict = field(default_factory=dict, repr=False)

    def _remember(self, kind: str, X: np.ndarray, a: float, b: float, compute) -> np.ndarray:
        """The same request (rows, times, kind) answers with the same array, so the two scorer builds of one build seed
        read bit-identical positions whatever the device's reduction order."""
        k = (kind, _key(X, seed=0), float(a), float(b))
        if k not in self._memo:
            self._memo[k] = compute()
        return self._memo[k].copy()

    def _trajectory(self):
        if self._traj is None:
            torch = _torch()
            q0 = torch.as_tensor(self.q0, dtype=torch.float64, device=self.device)
            p0 = torch.as_tensor(self.p0, dtype=torch.float64, device=self.device)
            with torch.no_grad():
                self._traj = shoot(q0, p0, self.sigma, self.n_steps)
        return self._traj

    def carry(self, X, a: float, b: float) -> np.ndarray:
        """Rows of X (numpy, n x 3) carried from time a to time b."""
        torch = _torch()
        X = np.asarray(X, np.float64)

        def compute():
            qs, ps = self._trajectory()
            out = np.empty_like(X)
            with torch.no_grad():
                for i in range(0, len(X), CHUNK):
                    x = torch.as_tensor(X[i:i + CHUNK], dtype=torch.float64, device=self.device)
                    out[i:i + CHUNK] = carry(x, qs, ps, a, b, self.sigma, self.n_steps).cpu().numpy()
            return out

        return self._remember("carry", X, a, b, compute)

    def jacobian_det(self, X, a: float, b: float) -> np.ndarray:
        """det D(map from a to b) at each row of X (autograd, one backward pass per output axis)."""
        torch = _torch()
        X = np.asarray(X, np.float64)

        def compute():
            qs, ps = self._trajectory()
            out = np.empty(len(X))
            for i in range(0, len(X), CHUNK):
                x = torch.as_tensor(X[i:i + CHUNK], dtype=torch.float64, device=self.device).requires_grad_(True)
                y = carry(x, qs, ps, a, b, self.sigma, self.n_steps)
                rows = [torch.autograd.grad(y[:, j].sum(), x, retain_graph=j < 2)[0] for j in range(3)]
                J = torch.stack(rows, dim=1)  # n x 3 (output) x 3 (input)
                out[i:i + CHUNK] = torch.linalg.det(J).detach().cpu().numpy()
            return out

        return self._remember("jacobian", X, a, b, compute)


def fit(Zb, Zt, Fb, Ft, seed: int, *, n_sub: int = SUBSAMPLE, n_control: int = N_CONTROL, sigma: float = SIGMA_V,
        scales=DATA_SCALES, gamma: float = GAMMA, n_steps: int = N_STEPS, iters: int = LBFGS_ITERS,
        device: str | None = None) -> GeodesicFit:
    """The registration from the earlier cloud Zb (with features Fb) to the target cloud Zt (Ft), at a build seed.
    The keyword arguments exist for the synthetic tests; the row's run uses the defaults."""
    torch = _torch()
    Zb, Zt = np.asarray(Zb, np.float64), np.asarray(Zt, np.float64)
    Fb, Ft = np.asarray(Fb, np.float64), np.asarray(Ft, np.float64)
    if len(Fb) != len(Zb) or len(Ft) != len(Zt) or Fb.shape[1] != Ft.shape[1]:
        raise ValueError("one feature row per cell and the same features in both brackets")
    rng = np.random.default_rng([int(seed), SUBSAMPLE_STREAM])
    ib = subsample_rows(len(Zb), n_sub, rng)
    it = subsample_rows(len(Zt), n_sub, rng)
    ic = kmeanspp(Zb[ib], n_control, np.random.default_rng([int(seed), CONTROL_STREAM]))
    dev = device or device_name()

    def T(a):
        return torch.as_tensor(np.ascontiguousarray(a), dtype=torch.float64, device=dev)

    X, Y, FX, FY = T(Zb[ib]), T(Zt[it]), T(Fb[ib]), T(Ft[it])
    q0 = X[torch.as_tensor(ic, device=dev)].clone()
    data = _Data(X, Y, FX, FY, scales)
    with torch.no_grad():
        d0 = float(data(X))
    if not np.isfinite(d0) or d0 <= 0:
        raise ValueError("the data term at the rigid start is not positive")
    K0 = _gauss(q0, q0, sigma)
    p = torch.zeros_like(q0, requires_grad=True)
    opt = torch.optim.LBFGS([p], lr=1.0, max_iter=iters, line_search_fn="strong_wolfe")

    def energy():
        qs, ps = shoot(q0, p, sigma, n_steps)
        Xm = carry(X, qs, ps, 0.0, 1.0, sigma, n_steps)
        kin = 0.5 * (p * (K0 @ p)).sum()
        return kin, data(Xm)

    def closure():
        opt.zero_grad()
        kin, d = energy()
        e = kin + gamma * d / d0
        e.backward()
        return e

    opt.step(closure)
    with torch.no_grad():
        kin, d1 = energy()
    p0 = p.detach()
    norms = p0.norm(dim=1)
    rec = {"sigma_v": float(sigma), "gamma": float(gamma), "n_steps": int(n_steps), "lbfgs_max_iter": int(iters),
           "lbfgs_n_iter": int(opt.state[opt._params[0]].get("n_iter", 0)),
           "feature_bandwidth": data.bandwidth, "data_rigid": d0, "data_final": float(d1),
           "data_ratio": float(d1) / d0, "kinetic": float(kin), "energy": float(kin) + gamma * float(d1) / d0,
           "p0_norm_mean": float(norms.mean()), "p0_norm_max": float(norms.max()), "device": str(dev)}
    return GeodesicFit(q0=q0.cpu().numpy(), p0=p0.cpu().numpy(), record=rec, sigma=float(sigma), n_steps=int(n_steps),
                       device=str(dev))


# ------------------------------------------------------------------------------------------------ cache
_CACHE: dict[str, GeodesicFit] = {}


def _key(*arrays, seed: int) -> str:
    h = hashlib.sha256()
    for a in arrays:
        a = np.ascontiguousarray(np.asarray(a, np.float64))
        h.update(np.asarray(a.shape, np.int64).tobytes())
        h.update(a.tobytes())
    h.update(np.int64(seed).tobytes())
    h.update(repr((SUBSAMPLE, N_CONTROL, SIGMA_V, DATA_SCALES, GAMMA, N_STEPS, LBFGS_ITERS, SUBSAMPLE_STREAM,
                   CONTROL_STREAM)).encode())
    return h.hexdigest()


def fit_cached(Zb, Zt, Fb, Ft, seed: int) -> GeodesicFit:
    """`fit` with the row's constants, kept for the last CACHE_SIZE input sets of this process, so the same build
    under two scorer builds fits once (the same inputs give the same fit)."""
    k = _key(Zb, Zt, Fb, Ft, seed=seed)
    if k not in _CACHE:
        if len(_CACHE) >= CACHE_SIZE:
            _CACHE.pop(next(iter(_CACHE)))
        _CACHE[k] = fit(Zb, Zt, Fb, Ft, seed, n_sub=SUBSAMPLE, n_control=N_CONTROL, sigma=SIGMA_V, scales=DATA_SCALES,
                        gamma=GAMMA, n_steps=N_STEPS, iters=LBFGS_ITERS)
    return _CACHE[k]
