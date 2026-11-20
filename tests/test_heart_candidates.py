"""Row R16's heart candidate list (experiments/harness/TRIES.md, as proposed), pinned against the builder as it stands.

The best file so far is C, the full trim (row R10: `build_test.py heart --entry core_trim --expression copy --z-mode
raw`, then halve.py). The list is K1, the half trim (form B: `--entry s2_trim`, with the size field the builder ties
to it), and K2, the full-scope copy (form A: `--entry s1`, set heart2's command). Both are commands the builder has
on main; this file adds no builder option. On synthetic stages, a synthetic test board whose limit is never shown and
the size stand-in of tests/test_heart_trend_build.py: each candidate differs from C only in its named change (the
scope option with its size field), with the same expression step (the copy draw, no momentum or trend), seed, count,
raw size mode and halving; the left-off items stay refused; the DGX launcher runs exactly these commands. Nothing is
read from data/ and nothing is scored."""

from __future__ import annotations

import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_build_heart_core_trim import CORE_TRIM  # noqa: E402
from test_build_heart_iso_trim import BEFORE, ISO_TRIM, _build, _ctx  # noqa: E402
from test_heart_trend_build import SECRET_MAX, STAGES, install, sha, synthetic_stage  # noqa: E402

FULL_COPY = ["--entry", "s1", "--expression", "copy", "--z-mode", "raw"]  # set heart2's command
# name -> (argv, output stem as the launcher writes it, the scope's trim fraction or None, the size step's hypothesis)
LIST = {
    "C": (CORE_TRIM, "T2heart_s2core", 0.0, "core"),
    "K1": (ISO_TRIM, "T2heart_s2trimcopy", 0.5, "H2"),
    "K2": (FULL_COPY, "T2heart_s2iso", None, "H1"),
}
# What a candidate's record may differ in from C's: the scope option and the size it brings, nothing else
SCOPE_KEYS = {"entry", "trim", "scope_option", "size"}
LAUNCHER = REPO / "experiments" / "builds" / "heart_candidates.sh"


@pytest.fixture(scope="module")
def bt():
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import build_test

    return build_test


@pytest.fixture
def heart_env(bt, monkeypatch, own_repo):
    reads: list[str] = []
    install(bt, monkeypatch.setattr, reads)
    stand_in = bt.t2_build.size_or_hold
    calls: list[dict] = []

    def size(coords, size_json, hypothesis=None, **kw):
        calls.append({"hypothesis": hypothesis, **kw})
        return stand_in(coords, size_json, hypothesis, **kw)

    monkeypatch.setattr(bt.t2_build, "size_or_hold", size)
    return reads, calls


def _pool(frac):
    from vec import t2_build
    from vec.io import COORD_KEY

    ref = synthetic_stage(*STAGES["E10.5.h5ad"])
    if frac is None:
        return ref
    keep = t2_build.radial_trim_mask(np.asarray(ref.obsm[COORD_KEY], dtype=float)[:, :3],
                                     np.asarray(ref.obs["core"], dtype=bool), frac)
    return ref[keep].copy()


# ------------------------------------------------------------------------------------------------ the list
def test_every_command_on_the_list_is_one_the_builder_has(bt):
    """No option is added: C, K1 and K2 parse to the copy draw in the raw mode at the default seed, with no trend,
    depth, herding or cfm, and differ only in --entry."""
    parsed = {k: bt.build_parser().parse_args(["heart", *v[0], "--out", "runs/x/h.h5ad"]) for k, v in LIST.items()}
    for a in parsed.values():
        assert a.expression == "copy" and a.z_mode == "raw" and a.seed == bt.harness.BUILD_SEED_OFFSET + 1
        assert not (a.trend or a.depth or a.herd or a.cfm)
    differ = {k for k in vars(parsed["C"]) if len({str(getattr(a, k)) for a in parsed.values()}) > 1}
    assert differ == {"entry"}
    assert bt.HEART_CORE_TRIM == "heart " + " ".join(CORE_TRIM) and bt.HEART_ISO_TRIM == "heart " + " ".join(ISO_TRIM)
    assert bt.HEART_SCOPE_OPTIONS["s2_trim"]["size_source"] == "z_modes.raw.scope.H2.pred_rms_full"
    assert bt.HEART_SCOPE_OPTIONS["core_trim"]["size_source"] == "z_modes.raw.pred_rms_core"


def test_each_candidate_differs_from_the_best_file_only_in_its_scope_option(bt, heart_env, tmp_path, capsys):
    """Each file is the copy draw (vec.harness.draw_copy, same seed and count rule) from its own pool with the drawn
    cells' released expression; the size step gets its own hypothesis in the raw mode against the full E10.5 cloud;
    the row count is the count rule's for all three (it does not read the pool); every other record field is C's."""
    from vec import harness
    from vec.io import COORD_KEY, dense

    reads, calls = heart_env
    built = {}
    for name, (argv, stem, frac, hyp) in LIST.items():
        reads.clear()
        calls.clear()
        out, rec, text = _build(bt, tmp_path, stem, *argv)
        start = synthetic_stage(*STAGES["E9.5.h5ad"])
        _, _, X, coords = harness.draw_copy(_ctx(_pool(frac), start))
        got = ad.read_h5ad(out)
        np.testing.assert_array_equal(dense(got.X).astype(np.float32), X)
        np.testing.assert_array_equal(np.asarray(got.obsm[COORD_KEY], dtype=np.float32)[:, :3],
                                      (np.asarray(coords, dtype=np.float32).astype(float) * 1.25).astype(np.float32))
        assert reads == ["E10.5.h5ad", "E9.5.h5ad"]  # no earlier stage: no trend
        assert len(calls) == 1 and calls[0]["hypothesis"] == hyp and calls[0]["z_mode"] == "raw"
        assert calls[0]["reference_cloud"].n_obs == synthetic_stage(*STAGES["E10.5.h5ad"]).n_obs
        assert got.n_obs == SECRET_MAX and rec["n_cells"] == "max_cells"
        assert str(SECRET_MAX) not in text + capsys.readouterr().out
        built[name] = (out, rec)
    c = built["C"][1]
    for name in ("K1", "K2"):
        rec = built[name][1]
        assert set(rec) - SCOPE_KEYS == set(c) - SCOPE_KEYS
        assert {k: rec[k] for k in set(rec) - SCOPE_KEYS} == {k: c[k] for k in set(c) - SCOPE_KEYS}
        assert rec["expression"] == "copy" and rec["z_mode"] == "raw" and "trend" not in rec
    assert built["K1"][1]["scope_option"]["size_source_field"] == "z_modes.raw.scope.H2.pred_rms_full"
    assert built["K2"][1]["trim"] is None and "scope_option" not in built["K2"][1]


def test_the_halves_keep_the_same_rows_for_every_candidate(bt, heart_env, tmp_path):
    """halve.py with its default seed (1000001), as the uploaded files: the same row positions in all three files."""
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    import halve

    rows = {}
    for name, (argv, stem, _, _) in LIST.items():
        out, _, _ = _build(bt, tmp_path, stem, *argv)
        half = out.parent / f"{stem}_half.h5ad"
        digest = halve.halve(str(out), str(half))
        full, got = ad.read_h5ad(out), ad.read_h5ad(half)
        idx = halve.half_rows(full.n_obs, halve.DEFAULT_SEED)
        assert digest == sha(half) and list(got.obs_names) == list(full.obs_names[idx])
        rows[name] = idx
    np.testing.assert_array_equal(rows["C"], rows["K1"])
    np.testing.assert_array_equal(rows["C"], rows["K2"])


def test_k2_is_set_heart2s_command_byte_for_byte(bt, heart_env, tmp_path):
    out, _, _ = _build(bt, tmp_path, "k2", *FULL_COPY)
    assert sha(out) == BEFORE["s1_copy_raw"][0]


# ------------------------------------------------------------------------------------------------ left off the list
@pytest.mark.parametrize("argv", [
    [*CORE_TRIM, "--trend"], [*ISO_TRIM, "--trend"], [*FULL_COPY, "--trend"],  # K3, the trend step, at any scope
    [*CORE_TRIM, "--depth"], [*ISO_TRIM, "--herd"], [*FULL_COPY, "--cfm"],  # depth, herding, cfm
    ["--entry", "core_trim", "--expression", "copy", "--z-mode", "recorded"],  # another size mode for the full trim
    ["--entry", "core_trim", "--expression", "copy"],
])
def test_the_left_off_items_stay_refused_before_reading_anything(bt, heart_env, tmp_path, argv):
    with pytest.raises(SystemExit) as e:
        bt.main(["heart", *argv, "--out", str(tmp_path / "h.h5ad")])
    assert "refused" in str(e.value)
    assert heart_env[0] == [] and not (tmp_path / "h.h5ad").exists()


@pytest.mark.parametrize("entry", ["quarter_trim", "trim_0.25", "H0"])
def test_no_other_scope_value_parses(bt, entry):
    with pytest.raises(SystemExit):
        bt.build_parser().parse_args(["heart", "--entry", entry, "--expression", "copy", "--z-mode", "raw",
                                      "--out", "runs/x/h.h5ad"])


# ------------------------------------------------------------------------------------------------ the launcher
def test_the_launcher_builds_exactly_the_list():
    sh = LAUNCHER.read_text()
    args = dict(re.findall(r'^(C|K1|K2)_ARGS="([^"]+)"$', sh, flags=re.M))
    assert {k: v.split() for k, v in args.items()} == {k: v[0] for k, v in LIST.items()}
    assert "build C T2heart_s2core $C_ARGS" in sh
    assert "build \"$STEP\" T2heart_s2trimcopy $K1_ARGS" in sh and "for STEP in K1a K1b" in sh
    assert "build \"$STEP\" T2heart_s2iso $K2_ARGS" in sh and "for STEP in K2a K2b" in sh
    assert "HEART3_HALF=04cde6361c04840ee47fff32eb9a20f4781bd649387cc540ed9fdea71e65fa35" in sh
    assert "HEART2_HALF=dc436d9225bb28e7997d2ec51e95254cd125e569d44c059596298d0a7b6af5c2" in sh
    assert "experiments/builds/halve.py" in sh and "vec check" in sh and "BOARD=T2:heart:test_extrap" in sh
    assert "OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16" in sh
    assert "git status --porcelain" in sh and "S2_RULES.sha256" in sh and "'^| R16 '" in sh
    assert "PATH=$HOME/.local/bin:$PATH" in sh
    # every builder, halving and check call sends its output to /dev/null
    for line in sh.splitlines():
        if re.search(r"build_test\.py heart|halve\.py \"|vec check \"", line):
            assert "> /dev/null 2>&1" in line, line
    listed = re.findall(r"tests/test_[a-z0-9_]+\.py", sh)
    assert "tests/test_heart_candidates.py" in listed and all((REPO / f).exists() for f in listed)
