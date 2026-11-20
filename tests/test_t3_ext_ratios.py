"""Log ratios, guards, the minimum-cells rule, shrinkage, chimeras and the frozen factor form of experiments/t3_ext,
on hand-made synthetic numbers only (no competition file, no outside data)."""

from __future__ import annotations

import importlib.util
import math
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / "experiments" / "t3_ext"
T3 = ROOT / "experiments" / "t3_decomp"
for p in (EXT, T3):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
import t3x_common as cm  # noqa: E402
import t3x_tables as tb  # noqa: E402


def _load(name):
    spec = importlib.util.spec_from_file_location(f"t3_ext_{name}", EXT / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    return m


pr = _load("process")


def test_constants_are_the_fixed_ones():
    assert (cm.CONFIDENT, cm.MIN_CONF_SHARE, cm.SHARE_RATIO_MAX, cm.SPEARMAN_MARGIN, cm.MIN_CELLS, cm.MIN_UMI,
            cm.PSEUDOCOUNT, cm.N_HALF, cm.CLIP, cm.STAGE_MAX) == (0.6, 0.5, 2.0, 0.05, 100, 20, 10.0, 200,
                                                                  (0.25, 4.0), 9.5)


def test_log_ratio_on_hand_checked_numbers():
    x_ko = np.array([100.0, 50.0, 10.0, 0.0, 30.0])
    x_ct = np.array([50.0, 50.0, 30.0, 0.0, 30.0])
    present = np.array([True, True, True, True, False])
    L = pr.log_ratio(x_ko, x_ct, present)
    # N over present genes only: 160 and 130
    assert L[0] == pytest.approx(math.log(110 / 160) - math.log(60 / 130))  # 0.398497
    assert L[0] == pytest.approx(0.398497, abs=1e-6)
    assert L[1] == pytest.approx(math.log(130 / 160))                          # -0.207639
    assert L[2] == 0.0 and L[3] == 0.0                                         # under 20 UMIs in an arm
    assert np.isnan(L[4])                                                      # absent gene
    big = pr.log_ratio(np.array([5000.0, 100.0]), np.array([100.0, 5000.0]), np.array([True, True]))
    assert big.tolist() == pytest.approx([math.log(4), math.log(0.25)])       # clipped to log [0.25, 4]


def test_shrinkage_on_hand_checked_numbers():
    # n_h = 2 / (1/100 + 1/400) = 160; cell term 160/360 = 0.4444; m = 2 -> 0.6667; min 0.4444
    assert cm.shrink_k([100, 400], 2, 1.0) == pytest.approx((160 / 360) * (2 / 3))   # 0.296296
    assert cm.shrink_k([100, 400], 2, 0.5) == pytest.approx(0.148148, abs=1e-6)
    # many cells, one library per arm: the cell term is capped by m/(m+1) = 0.5
    assert cm.shrink_k([1000, 1000], 1, 1.0) == pytest.approx(0.25)
    assert cm.shrink_k([1000, 1000], 4, 1.0) == pytest.approx(0.8 * 0.8)   # min(0.8333, 0.8) * 0.8
    assert cm.shrink_k([0, 1000], 2, 1.0) == 0.0


def test_design_factor():
    f = cm.design_factor
    assert f({"design": "germline", "control_relation": "littermate"}) == 1.0
    assert f({"design": "conditional", "control_relation": "same_capture"}) == 1.0
    assert f({"design": "conditional", "control_relation": "heterozygous"}) == 0.5
    assert f({"design": "germline", "control_relation": "same_series"}) == 0.5
    assert f({"design": "f0_crispr", "control_relation": "littermate"}) == 0.5
    assert f({"design": "chimera", "control_relation": "same_capture", "chimera": {"second_difference": True}}) == 1.0
    assert f({"design": "chimera", "control_relation": "same_capture"}) == 0.5


# ---------------------------------------------------------------------------------------------------------------------
# a synthetic dataset: three reference groups with distinct profiles

P = 40
RNG = np.random.default_rng(7)
GROUPS = ["CM", "Endoderm", "Neural"]
MU = {g: RNG.dirichlet(np.full(P, 0.5)) for g in GROUPS}
REF = {"genes": [f"g{i}" for i in range(P)], "groups": GROUPS, "shares": {"CM": 0.2, "Endoderm": 0.3, "Neural": 0.5},
       "mu": MU, "file_sha256": "x"}
FOLD = np.ones(P)
FOLD[:5] = 3.0


def _cells(arm, g, n, sample, stage=8.5, fold=None, conf=0.9, depth=3000, rng=None):
    rng = rng or np.random.default_rng(zlib.crc32(f"{arm}{g}{n}{sample}{stage}".encode()))
    p = MU[g] * (fold if fold is not None else 1.0)
    C = rng.multinomial(depth, p / p.sum(), size=n)
    obs = pd.DataFrame({"arm": arm, "stage": stage, "sample": sample, "call": g, "conf": conf},
                       index=[f"{sample}:{g}:{arm}:{i}" for i in range(n)])
    return C, obs


def _stack(parts):
    return np.vstack([c for c, _ in parts]), pd.concat([o for _, o in parts])


def test_unit_ratio_passes_a_clean_group_and_equals_k_times_the_hand_ratio():
    parts = [_cells("ko", "CM", 150, "K1", fold=FOLD), _cells("control", "CM", 150, "C1"),
             _cells("ko", "Neural", 150, "K1"), _cells("control", "Neural", 150, "C1")]
    C, obs = _stack(parts)
    present = np.ones(P, bool)
    Lg, qc = pr.unit_log_ratios(C, present, obs, REF, d=1.0)
    assert qc["status"] == "ok" and set(Lg) == {"CM", "Neural"}
    x_ko, x_ct = parts[0][0].sum(0), parts[1][0].sum(0)
    k = cm.shrink_k([150, 150], 1, 1.0)   # n_h 150 -> 150/350 = 0.4286 < 0.5; times 0.5
    assert k == pytest.approx((150 / 350) * 0.5)
    assert Lg["CM"] == pytest.approx(k * pr.log_ratio(x_ko, x_ct, present))
    enough = np.minimum(x_ko, x_ct) >= cm.MIN_UMI
    assert np.all(Lg["CM"][:5][enough[:5]] > 0) and np.all(Lg["CM"][~enough] == 0)
    assert qc["stages"]["E8.5"]["groups"]["CM"]["m"] == 1


def test_minimum_cells_rule():
    parts = [_cells("ko", "CM", 99, "K1"), _cells("control", "CM", 99, "C1"),
             _cells("ko", "Neural", 100, "K1"), _cells("control", "Neural", 100, "C1")]
    C, obs = _stack(parts)
    Lg, qc = pr.unit_log_ratios(C, np.ones(P, bool), obs, REF, d=1.0)
    assert set(Lg) == {"Neural"}
    assert "fewer than 100" in qc["stages"]["E8.5"]["groups"]["CM"]["skipped"]


def test_relabelling_rule_drops_a_group_whose_call_share_moves_more_than_2x():
    parts = [_cells("ko", "CM", 400, "K1"), _cells("ko", "Neural", 200, "K1"),
             _cells("control", "CM", 150, "C1"), _cells("control", "Neural", 450, "C1")]
    C, obs = _stack(parts)
    Lg, qc = pr.unit_log_ratios(C, np.ones(P, bool), obs, REF, d=1.0)
    assert "CM" not in Lg     # share 0.67 vs 0.25
    assert "relabelling" in qc["stages"]["E8.5"]["groups"]["CM"]["skipped"]


def test_a_stage_with_few_confident_cells_drops_and_no_group_means_left_out():
    parts = [_cells("ko", "CM", 150, "K1", conf=0.55), _cells("control", "CM", 150, "C1")]
    C, obs = _stack(parts)
    Lg, qc = pr.unit_log_ratios(C, np.ones(P, bool), obs, REF, d=1.0)
    assert Lg == {} and qc["status"] == "left_out"
    assert "fewer than half" in qc["stages"]["E8.5"]["skipped"]


def test_profile_guard_drops_a_group_whose_control_profile_is_another_groups():
    ko, _ = _cells("ko", "Neural", 150, "K1")
    ct, _ = _cells("control", "Neural", 150, "C1")
    obs = pd.concat([_cells("ko", "CM", 150, "K1")[1], _cells("control", "CM", 150, "C1")[1]])
    Lg, qc = pr.unit_log_ratios(np.vstack([ko, ct]), np.ones(P, bool), obs, REF, d=1.0)
    assert "CM" not in Lg
    assert qc["stages"]["E8.5"]["groups"]["CM"]["profile_guard"]["best"] == "Neural"


def test_stages_average_with_equal_weight_and_a_late_stage_refuses():
    parts = [_cells("ko", "CM", 150, "K1", fold=FOLD), _cells("control", "CM", 150, "C1"),
             _cells("ko", "CM", 300, "K2", stage=9.5, fold=FOLD), _cells("control", "CM", 300, "C2", stage=9.5)]
    C, obs = _stack(parts)
    present = np.ones(P, bool)
    Lg, _ = pr.unit_log_ratios(C, present, obs, REF, d=1.0)
    a = cm.shrink_k([150, 150], 1, 1.0) * pr.log_ratio(parts[0][0].sum(0), parts[1][0].sum(0), present)
    b = cm.shrink_k([300, 300], 1, 1.0) * pr.log_ratio(parts[2][0].sum(0), parts[3][0].sum(0), present)
    assert Lg["CM"] == pytest.approx((a + b) / 2)
    late = obs.copy()
    late["stage"] = 10.5
    with pytest.raises(cm.Refused):
        pr.unit_log_ratios(C, present, late, REF, d=1.0)


def test_chimeras_need_the_second_difference_and_take_it():
    present = np.ones(P, bool)
    parts = [_cells("ko_mut", "CM", 150, "P1", fold=FOLD), _cells("ko_host", "CM", 150, "P1"),
             _cells("wt_mut", "CM", 150, "P2"), _cells("wt_host", "CM", 150, "P2")]
    C, obs = _stack(parts)
    with pytest.raises(cm.Refused, match="second difference"):
        pr.unit_log_ratios(C[:300], present, obs.iloc[:300], REF, d=1.0, chimera=True)
    with pytest.raises(cm.Refused):
        pr.unit_log_ratios(C, present, obs, REF, d=1.0, chimera=False)
    Lg, qc = pr.unit_log_ratios(C, present, obs, REF, d=1.0, chimera=True)
    x = {a: p[0].sum(0) for a, p in zip(pr.CHIMERA_ARMS, parts, strict=True)}
    raw = pr.log_ratio(x["ko_mut"], x["ko_host"], present) - pr.log_ratio(x["wt_mut"], x["wt_host"], present)
    raw = np.where(np.min(np.vstack(list(x.values())), 0) < cm.MIN_UMI, 0.0, np.clip(raw, *cm.LOG_CLIP))
    assert Lg["CM"] == pytest.approx(cm.shrink_k([150] * 4, 1, 1.0) * raw)


# ---------------------------------------------------------------------------------------------------------------------
# the frozen factor's form


def test_pooled_finish_equals_the_frozen_gene_factor_on_a_hand_made_case():
    import t3_decompose as t3d

    rng = np.random.default_rng(3)
    mu_wt = rng.dirichlet(np.full(P, 2.0), size=3)
    ratio = np.exp(rng.normal(0, 0.4, P))
    mu_ko = mu_wt * ratio
    mu_ko /= mu_ko.sum(1, keepdims=True)
    # one group passing (mask), weight its share: the frozen ratio of weighted means equals the one-group ratio
    w, mask = np.array([0.2, 0.3, 0.5]), np.array([1.0, 0.0, 0.0])
    frozen, _ = t3d.gene_factor(mu_ko, mu_wt, w, mask, eps=0.0)
    ref = {"genes": REF["genes"], "groups": ["A", "B", "C"], "shares": dict(zip("ABC", w, strict=True)),
           "mu": dict(zip("ABC", mu_wt, strict=True))}
    Lg = {"A": np.log(mu_ko[0] / mu_wt[0])}
    ours, used = tb.finish_pooled(Lg, ref)
    assert used == ["A"]
    assert ours == pytest.approx(frozen, rel=1e-9)     # same form; the frozen eps (1e-5) is the only difference
    assert ours == pytest.approx(t3d.gene_factor(mu_ko, mu_wt, w, mask)[0], rel=0.05)
    # same ratio in two groups: still the frozen value (the rescale uses the passing groups' weighted profile)
    mu_ko2 = mu_wt * ratio
    frozen2, _ = t3d.gene_factor(mu_ko2, mu_wt, w, np.array([1.0, 1.0, 0.0]), eps=0.0)
    ours2, _ = tb.finish_pooled({"A": np.log(ratio), "B": np.log(ratio)}, ref)
    assert ours2 == pytest.approx(frozen2, rel=1e-9)


def test_finish_order_rescale_clip_then_absent_exactly_one():
    den = np.array([0.5, 0.3, 0.2])
    L = np.array([np.log(8.0), 0.0, np.nan])
    phi = tb.finish(L, den)
    # rescale on (8, 1, 1): den.phi / sum den = 4.0 + 0.3 + 0.2 = 4.5 -> (1.778, 0.222, 0.222); clip -> 0.25
    assert phi[0] == pytest.approx(8 / 4.5)
    assert phi[1] == pytest.approx(0.25)
    assert phi[2] == 1.0
    assert tb.finish(np.full(3, np.nan), den).tolist() == [1.0, 1.0, 1.0]


def test_a_written_pooled_table_is_accepted_by_the_frozen_load_factor(tmp_path):
    import build_submission as bs

    genes = REF["genes"]
    phi, _ = tb.finish_pooled({"CM": np.where(np.arange(P) < 5, 0.5, 0.0)}, REF)
    tb.write_table(tmp_path / "t", "ext_cardiac", genes, phi, {"CM": tb.finish(np.zeros(P), MU["CM"])}, {"x": 1})
    got, rec = bs.load_factor(tmp_path / "t", "ext_cardiac", genes, {})
    assert got.tolist() == pytest.approx(phi.tolist(), rel=0, abs=0)
    assert rec["recorded_in"].endswith("factors.json")
    with pytest.raises(bs.Refused):
        bs.load_factor(tmp_path / "t", "ext_cardiac", genes[::-1], {})
    pg = pd.read_csv(tmp_path / "t" / "ext_cardiac_per_group.csv")
    assert list(pg.columns) == ["group", "gene", "factor"] and set(pg["group"]) == {"CM"}
    assert pg["gene"].tolist() == genes
