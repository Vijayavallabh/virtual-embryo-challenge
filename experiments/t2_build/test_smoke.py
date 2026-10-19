"""The smoke test's shareable summary, on a synthetic summary (no competition data, no board file).

The synthetic summary has every field experiments/t2_build/smoke.py writes. Its count, sample, guard and check records
come from vec.t2_build itself (n_cells, sample_cells, cloud_guard, build_submission) on made-up clouds and boards with
planted cell limits; made-up emulator results go through smoke._summarise and the records through smoke._checks. The
shareable pair that smoke._write makes of it must keep only the allow-listed fields, and no planted limit, submitted
count or simple derivation of them may appear in its JSON or its text.

Run: uv run pytest -q experiments/t2_build/test_smoke.py
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pytest

from vec import t2_build as tb
from vec.boards import GROUPS, RANKED_METRICS, Board
from vec.evaluate import SCORER_BUILDS, _handedness_summary

sys.path.insert(0, str(Path(__file__).parent))
import smoke  # noqa: E402

GENES = [f"g{i}" for i in range(20)]
SEEDS = (0, 1, 2)
DIGEST, COMMIT = "abcdef" * 10 + "abcd", "deadbeef" * 5  # no digits: every digit in the output is a value
# Planted limits: no limit, count or ratio derived from them equals another number in the summary, so any match in the
# shareable pair is a leak. Each job runs one regime and stands for no real board. Job a: the count (max_cells) exceeds
# the pool, so rows repeat and are jittered. Job b: the count is an integer override (the rule then names the count) at
# most half the pool, so the block control draws the same count.
#        board key        min_cells truth max   pool  dropped override
PLANTED = {
    "a": ("T2:synthetic:a", 1013, 1427, 2749, 2011, 517, None),
    "b": ("T2:synthetic:b", 1013, 1129, 1877, 3907, 211, 1609),
}
TOKEN = re.compile(r"\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")
# The decimals specified for the shareable's guard ratios and scores, written out here independently of smoke.py.
DECIMALS = {"cover_ratio": 2, "rms_ratio": 3, "expected_rms_ratio": 3, "score_mean": 2, "score_sd": 2, "skills_mean": 3}


def cloud(n, seed):
    """A filled box cut into 15-unit z sections, rows stored section by section (as in the released files)."""
    rng = np.random.default_rng(seed)
    C = rng.uniform(-1.0, 1.0, size=(n, 3)) * np.array([120.0, 80.0, 60.0])
    C[:, 2] = np.round(C[:, 2] / 15.0) * 15.0
    return C[np.argsort(-C[:, 2], kind="stable")]


def expression(n, seed, scaled):
    """log1p of MERFISH-like counts (a median library of about 150), scaled to 1e4 per cell or not."""
    rng = np.random.default_rng(seed)
    counts = rng.poisson(rng.gamma(0.5, 15.0, size=(1, len(GENES))), size=(n, len(GENES))).astype(np.float64)
    counts[counts.sum(1) == 0, 0] = 1.0
    return np.log1p(1e4 * counts / counts.sum(1, keepdims=True) if scaled else counts)


def scored(board, n_pred, build, offset):
    """smoke._summarise of a made-up score_board_seeds result (its keys, made-up values) plus the timing smoke adds."""
    per = []
    for s in SEEDS:
        skills = {m: 0.7031 + offset + 0.0019 * s + 0.0117 * j for j, m in enumerate(RANKED_METRICS["T2"])}
        hand = {"dice_proper": 0.7113, "dice_mirror": 0.7051, "dice_gap": 0.0062, "parity_pred": 1, "parity_truth": -1,
                "prefix_mode": "mirror", "dice_prefix": 0.7051, "metric_seed": s, "local_parity": -1,
                "local_parity_n_true": board.truth_cells, "local_parity_seed": 0, "panel_dice_matches_build": True}
        per.append({"seed": s, "score": 50.0 + 20.0 * skills["d2_shape"], "skills": skills,
                    "metrics": {m: 1.0 - v for m, v in skills.items()},
                    "group_skill": {g: skills["d2_shape"] for g, _, _ in GROUPS["T2"]}, "anchors": {},
                    "extras": {"build": build, "handedness": hand, "truth_cells": 313, "ref_cells": 307,
                               "pred_cells": n_pred,
                               "inverted_anchors": ["occupancy_dice"] if build == "prefix" and s == 1 else []}})
    scores = np.array([p["score"] for p in per])
    res = {"board": board.key, "score_mean": float(scores.mean()), "score_sd": float(scores.std(ddof=1)),
           "per_seed": per, "build": build, "handedness": _handedness_summary([p["extras"]["handedness"] for p in per])}
    return smoke._summarise(res) | {"seconds": 1.5}


def synthetic_summary(score: bool) -> tuple[dict, dict]:
    """A summary with every field smoke.main writes (as _run_job and main assemble it), and per job the planted values
    and derivations the shareable pair must not hold: {name: value}."""
    summary = {"generated_utc": "2026-10-19T12:00:00+00:00", "git_commit": COMMIT, "host": "synthetic",
               "omp_num_threads": "8", "seeds": list(SEEDS), "sample_seed": 0, "anchors": "computed (stand-in targets)",
               "jobs": {}, "checks": {}}
    secrets = {}
    for i, (name, (key, min_cells, truth, max_cells, n_pool, dropped, override)) in enumerate(PLANTED.items()):
        board = Board(key=key, task="T2", genes=GENES, needs_coords=True, min_cells=min_cells, max_cells=max_cells,
                      truth_cells=truth)
        P = cloud(n_pool, seed=i)
        X = expression(n_pool, seed=i, scaled=name == "a")  # job b's cells were never scaled: a warning
        n, count = tb.n_cells(board, n_pool, override=override)
        s = tb.sample_cells(P, n, seed=0)
        _, rep = tb.build_submission(board, X[s.rows], s.coords, genes=GENES)
        guard = {"copy": tb.cloud_guard(s.coords, P, seed=0).as_dict()}
        m = min(n, n_pool // 2)
        start = (n_pool - m) // 2
        guard["block_control"] = tb.cloud_guard(P[start:start + m], P, seed=0).as_dict()
        target_rms = 1.2837 * tb.rms_radius(P)  # more decimals than the shareable's RMS ratios keep: rounding shows
        scaled = tb.rescale_to_rms(s.coords, target_rms)
        guard["oracle_size_control"] = tb.cloud_guard(scaled, P, seed=0, expected_rms_ratio=tb.rms_radius(scaled)
                                                      / tb.rms_radius(P)).as_dict()
        n_pred = {"copy": n, "block_control": m, "oracle_size_control": n}
        r = {
            "board": key, "board_spec": {"min_cells": min_cells, "max_cells": max_cells, "truth_cells": truth,
                                         "n_genes": len(GENES)},
            "inputs": {k: {"file": f"T2/{k}_{name}.h5ad", "sha256": DIGEST} for k in ("source", "target", "ref")},
            "pool": {"n_input": n_pool + dropped, "n_dropped_duplicates": dropped, "n_pool": n_pool,
                     "rms_radius": tb.rms_radius(P)},
            "count": count, "sample": s.record,
            "check_submission": {"ok": rep.ok, "errors": rep.errors, "warnings": rep.warnings, "n_cells": rep.n_cells},
            "block_control": {"n": m, "first_row": start},
            "oracle_size_control": {"target_rms_radius": target_rms, "scale": target_rms / tb.rms_radius(s.coords)},
            "guard": guard, "files_mb": {cond: round(k / 97.0, 1) for cond, k in n_pred.items()},
            "scores": ({cond: {b: scored(board, n_pred[cond], b, 0.031 * j + 0.0007 * i) for b in SCORER_BUILDS}
                        for j, cond in enumerate(guard)} if score else {}),
            "build_seconds": 3.2,
        }
        summary["checks"][name] = smoke._checks(board, n, r, score)
        r["seconds"] = 10.5
        summary["jobs"][name] = r
        planted = {"n": n, "max_cells": max_cells, "truth_cells": truth, "min_cells": min_cells,
                   "|n - pool|": abs(n - n_pool), "max_cells - truth_cells": max_cells - truth, "n - truth_cells":
                   n - truth, "n / pool": n / n_pool, "n / truth_cells": n / truth, "truth_cells / max_cells":
                   truth / max_cells, "truth_cells / pool": truth / n_pool}
        if n != max_cells:
            planted["n / max_cells"] = n / max_cells
        if m == n:  # the block control drew the submitted count
            planted |= {"block count": m, "block first row": start}
        if s.record["n_jittered"]:  # the jitter of the repeated rows replays to their number
            planted |= {k: s.record[k] for k in ("max_shift", "mean_shift_jittered")}
        secrets[name] = planted
    summary["all_checks_pass"] = all(all(c.values()) for c in summary["checks"].values())
    return summary, secrets


def hits(text: str, secrets: dict) -> list[tuple[str, str]]:
    """(token, 'job: secret') for every number in `text` that equals a planted value: integers exactly, ratios and
    shifts at full precision or rounded to 2, 3 or 4 decimals."""
    out = []
    for tok in TOKEN.findall(text):
        t = float(tok)
        for job, planted in secrets.items():
            for name, v in planted.items():
                if isinstance(v, int | np.integer):
                    match = t == v
                else:
                    match = abs(t - v) <= 1e-9 * abs(v) or any(t == round(v, k) for k in (2, 3, 4))
                if match:
                    out.append((tok, f"{job}: {name}"))
    return out


def leaf_paths(o, path=""):
    """Dotted paths of a JSON object's leaves; list items share the suffix []."""
    if isinstance(o, dict):
        for k, v in o.items():
            yield from leaf_paths(v, f"{path}.{k}" if path else k)
    elif isinstance(o, list):
        for v in o:
            yield from leaf_paths(v, f"{path}[]")
    else:
        yield path


def allowed_paths(summary: dict) -> set[str]:
    """The allow-list as specified for the shareable summary, written out here independently of smoke.py."""
    paths = {"generated_utc", "git_commit", "seeds[]", "shareable", "all_checks_pass"}
    for name, r in summary["jobs"].items():
        paths |= {f"checks.{name}.{c}" for c in summary["checks"][name]}
        p = f"jobs.{name}."
        paths |= {p + k for k in ("board", "count_rule", "pool.n_pool", "pool.n_dropped_duplicates",
                                  "check_submission.ok", "check_submission.warnings[]")}
        paths |= {f"{p}inputs.{k}.{f}" for k in r["inputs"] for f in ("file", "sha256")}
        for cond, g in r["guard"].items():
            paths |= {f"{p}guard.{cond}.{f}" for f in ("ok", "cover_ratio", "rms_ratio")}
            if g["expected_rms_source"] == "argument":  # passed in (the oracle's), not a reference statistic
                paths.add(f"{p}guard.{cond}.expected_rms_ratio")
        for cond, by_build in r["scores"].items():
            for b, v in by_build.items():
                q = f"{p}scores.{cond}.{b}."
                paths |= {q + f for f in ("score_mean", "score_sd", "inverted_anchors[]", "inverted_anchor_seeds[]")}
                paths |= {f"{q}skills_mean.{m}" for m in v["skills_mean"]}
    return paths


@pytest.fixture(scope="module", params=[True, False], ids=["scored", "no-score"])
def run(request, tmp_path_factory):
    """The synthetic summary, its planted values and the folder holding the four files smoke._write makes of it."""
    summary, secrets = synthetic_summary(score=request.param)
    out = tmp_path_factory.mktemp("smoke")
    smoke._write(out, summary, ",".join(str(s) for s in SEEDS))
    return summary, secrets, out


def test_synthetic_summary_carries_the_leaks_the_shareable_must_drop(run):
    summary, secrets, out = run
    full = json.loads((out / "smoke_summary.json").read_text())
    assert full["jobs"]["b"]["count"]["rule"] == "override 1609"
    assert summary["jobs"]["a"]["sample"]["n_jittered"] == 2749 - 2011
    g = full["jobs"]["b"]["guard"]  # one seed and one count: the reference draws coincide
    assert g["copy"]["ref_rms_ratio"] == g["block_control"]["ref_rms_ratio"]
    assert g["copy"]["expected_rms_ratio"] == g["block_control"]["expected_rms_ratio"]
    # the scan finds the planted values where they are: the full pair, which stays on the DGX
    for f, names in (("smoke_summary.json", ("n", "max_cells", "truth_cells", "min_cells", "n / pool")),
                     ("smoke_summary.txt", ("n", "truth_cells"))):
        found = {h for _, h in hits((out / f).read_text(), secrets)}
        assert {f"{job}: {k}" for job in PLANTED for k in names} <= found, f


def test_shareable_keeps_only_the_allow_list(run):
    summary, _, out = run
    shared = json.loads((out / "smoke_summary_shareable.json").read_text())
    paths, allowed = set(leaf_paths(shared)), allowed_paths(summary)
    assert paths <= allowed, sorted(paths - allowed)
    assert {p for p in allowed if not p.endswith("[]")} <= paths, sorted({p for p in allowed - paths})
    assert shared["jobs"]["a"]["count_rule"] == "max_cells"
    assert shared["jobs"]["b"]["count_rule"] == "override"  # the override's argument is the count
    assert any("implied library" in w for w in shared["jobs"]["b"]["check_submission"]["warnings"])


def test_shareable_rounds_guard_ratios_and_scores_to_the_specified_decimals(run):
    summary, _, out = run
    shared = json.loads((out / "smoke_summary_shareable.json").read_text())
    full = json.loads((out / "smoke_summary.json").read_text())
    leaves = []  # (field, shareable value, full summary value)
    for name, r in shared["jobs"].items():
        for cond, g in r["guard"].items():
            leaves += [(f, x, full["jobs"][name]["guard"][cond][f]) for f, x in g.items() if f in DECIMALS]
        for cond, by_build in r["scores"].items():
            for b, v in by_build.items():
                fv = full["jobs"][name]["scores"][cond][b]
                leaves += [(f, v[f], fv[f]) for f in ("score_mean", "score_sd")]
                leaves += [("skills_mean", x, fv["skills_mean"][m]) for m, x in v["skills_mean"].items()]
    assert [(f, x, v) for f, x, v in leaves if x != round(v, DECIMALS[f])] == []
    # every field is there, with made-up values that a dropped, wider or narrower round() would change
    fields = list(DECIMALS) if summary["jobs"]["a"]["scores"] else ["cover_ratio", "rms_ratio", "expected_rms_ratio"]
    for field in fields:
        k, vs = DECIMALS[field], [v for f, _, v in leaves if f == field]
        assert any(round(v, k + 1) != round(v, k) for v in vs), field  # unrounded or wider
        assert any(round(v, k - 1) != round(v, k) for v in vs), field  # narrower


def test_shareable_pair_holds_no_planted_limit_or_derivation(run):
    _, secrets, out = run
    for f in ("smoke_summary_shareable.json", "smoke_summary_shareable.txt"):
        assert hits((out / f).read_text(), secrets) == [], f


def test_written_pair_is_what_shareable_and_lines_give(run):
    summary, _, out = run
    seeds = ",".join(str(s) for s in SEEDS)
    shared = smoke._shareable(summary)
    lines = smoke._lines(shared, seeds)
    assert (out / "smoke_summary_shareable.json").read_text() == json.dumps(shared, indent=2, default=float) + "\n"
    assert (out / "smoke_summary_shareable.txt").read_text() == "\n".join(lines) + "\n"
    full_txt = (out / "smoke_summary.txt").read_text()
    assert full_txt.startswith("\n".join(lines) + "\n") and "counts (full report only)" in full_txt
    assert json.loads((out / "smoke_summary.json").read_text()) == json.loads(json.dumps(summary, default=float))
    assert sum(line.startswith("  score ") for line in lines) == (12 if summary["jobs"]["a"]["scores"] else 0)
