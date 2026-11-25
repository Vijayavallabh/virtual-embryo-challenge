"""`build_test.py t1` under rows R54 and R58 of experiments/harness/TRIES.md (the step-change round's T1 lines), each
one change to T1_s2rcMSB's command (row R26's: `t1 --g-mom --base lib3 --feng-groups recomputed --feng-fraction-mult
1.5 --direction-momentum all`):
  T1_s2rcMSB_LM16  ... --late-mix e15_17 --captain-ruling R54 (library 3's own chamber mix moved to row N64's shares)
  T1_s2rcMSBm20    ... --direction-momentum-damping 2.0 --captain-ruling R58 (MS's own damping 2.0)

Held here on synthetic stages and directions only (the fixtures of tests/test_t1_r51_build.py), with temporary ledgers
and a temporary shares file in place of the repository's (no test reads the repository's ledger or a real shares file):
(1) nothing that builds today changes: the builder at main 9a1b641 (read with `git show`) and this one route every
argument set of a grid alike unless the set gives --late-mix or --captain-ruling R58, and the earlier lines' tests keep
pinning every earlier command's bytes; (2) each new file is taken only with its row (rows R54's and R58's first cells
as the ledger holds them pass the builder's `_captain_order_names`), and every other use is refused before anything but
the ledger (and, for R54, its shares file) is read; (3) R58's file is MSB's command with MS's damping 2.0: with that
damping put back at B0's it is MSB's bytes (the checker's m20msb); (4) R54's file draws library 3 by quotas at the
shares file's moved mix: with the shares at library 3's own values its quotas are MSB's (identity on quotas, not bytes
alone) and so are its bytes, and with the moved mix put back at the pool's own it is MSB's bytes (the checker's
lm16own); every step after the draw is MSB's; only the shares file committed at HEAD with the sha256 row N64 names is
taken; (5) both build the same bytes twice. No released file, outside-atlas file, board or score is read."""

from __future__ import annotations

import hashlib
import importlib.util
import itertools
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for _p in ("experiments/external", "experiments/t1_m1", "experiments/t1_direction", "experiments/t1_hygiene",
           "experiments/builds", "experiments/t1_scope"):
    sys.path.insert(0, str(ROOT / _p))

import build_test as _bt  # noqa: E402  (the module itself, for the helpers)
import feng_relabel as fr  # noqa: E402
import late_mix  # noqa: E402
import test_t1_candidates_build as tcb  # noqa: E402
import test_t1_feng_groups_none  # noqa: E402, F401  (imported before any fixture, as tests/test_feng_relabel.py does)
import test_t1_r51_build as r51t  # noqa: E402
from conftest import git_env  # noqa: E402
from test_feng_relabel import _recomputed_env  # noqa: E402
from test_hidden_counts import t1_env  # noqa: E402, F401  (the synthetic build_t1 fixture)
from test_t1_candidates_build import RC, _build, _spy  # noqa: E402
from test_t1_msb_build import MSB, _three_libraries  # noqa: E402

from vec import harness  # noqa: E402
from vec.core import quotas  # noqa: E402

bt = tcb.bt  # the builder module's fixture

SEED = harness.BUILD_SEED_OFFSET + 1
PAIR = ["--feng-fraction-mult", "1.5", "--direction-momentum", "all"]
LM16 = [*PAIR, "--late-mix", "e15_17", "--captain-ruling", "R54"]
M20 = [*PAIR, "--direction-momentum-damping", "2.0", "--captain-ruling", "R58"]
BASE_CMD = "t1 --g-mom --base lib3 --feng-groups recomputed"
# Rows R54's and R58's first cells as the ledger holds them on main (9a1b641), verbatim
FIRST_R54 = ("R54 (the captain's order of 2026-11-24, ruling, no run: \"Yes, T1_s2rcMSB_LM16, my named exception, "
             "bound to T72: T1's composition and depth and the 17 November request rule may change for this file, and "
             "the strict form is set aside for Feng's pooled E15.5–E17.5 counts, at row R54's shares.\"; T1's late-mix "
             "line, the step-change round's T1 CLOSED design; the captain's named exception on T1's composition and "
             "depth, closed after the 4 November exposure, on the 17 November rule for requests on them, and on the "
             "strict form for Feng's permitted CD1 cells at E15.5 to E17.5 (pooled library GSM5795773), for one file, "
             "the banked command with library 3's chamber mix at the shares that row N64's step commits "
             "(`shares.json`, its sha256 written there), built only after row T72 reads PASS; the T1 builder is to "
             "take `--late-mix e15_17` only with `--captain-ruling R54`; decision log 2026-11-24 (www); not a try)")
FIRST_R58 = ("R58 (the captain's order of 2026-11-24, ruling, no run: \"Yes, T1_s2rcMSBm20 under R22, bound to its "
             "gate G0.\"; T1's line 11, the step-change round's T1 PORTFOLIO design; the T1 builder is to take "
             "`--feng-fraction-mult 1.5 --direction-momentum all --direction-momentum-damping 2.0` on "
             "`t1 --g-mom --base lib3 --feng-groups recomputed` only with `--captain-ruling R58`; decision log "
             "2026-11-24 (www); not a try)")
ROW_R54 = f"| {FIRST_R54} | T1 | late mix |"
ROW_R58 = f"| {FIRST_R58} | T1 | line 11 |"
# Numbers a refusal may show: the options' values, lib3, the rows' IDs, the files' names and the anchor's stages
NUMBERS_ALLOWED = {"0.5", "1.0", "1.5", "2.0", "1", "2", "3", "15", "17", "16", "20", "13", "23", "26", "51", "54",
                   "58", "64", "10", "5", "256"}


def _n64(sha: str | None) -> str:
    named = "" if sha is None else f"`shares.json sha256 {sha}`"
    return (f"| N64 (step 0 of row R54's late-mix line, target-free, not a try: `experiments/t1_scope/late_mix.py`) "
            f"| T1 | step 0 | written: {named} |")


def _ledger(r54: str | None = ROW_R54, r58: str | None = ROW_R58, n64: str | None = None) -> str:
    return r51t._ledger_text() + "".join(f"{row}\n" for row in (r54, r58, n64) if row)


@pytest.fixture
def ledger(bt, tmp_path, monkeypatch):
    """A temporary ledger holding row R51's (tests/test_t1_r51_build.py's) and rows R54 and R58 as main holds their
    first cells; the builder reads it in place of the repository's."""
    p = tmp_path / "TRIES.md"
    p.write_text(_ledger())
    monkeypatch.setattr(bt, "TRIES_FILE", p)
    return p


@pytest.fixture
def ms4(bt, t1_env, monkeypatch):  # noqa: F811
    _recomputed_env(bt, t1_env, monkeypatch)
    _three_libraries(t1_env)
    return _spy(bt, monkeypatch)


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _shares_record(env, cm: float | None = None, blood: float | None = None, own: dict | None = None) -> dict:
    """A shares file as row N64's step writes it, on the fixture's library 3 (its own mix by the step's own call),
    moved to `cm` and `blood` (default: library 3's own values)."""
    own = late_mix.own_shares(env.ref) if own is None else own
    heart = list(_bt.t1_scope.load_spec().heart)
    cm = _bt.t1_scope.share_of(own, "cm", heart) if cm is None else cm
    blood = float(own["blood"]) if blood is None else blood
    return {"row": "N64", "late_mix": "e15_17", "file": "T1_s2rcMSB_LM16", "own_shares": own, "cm_share_of_heart": cm,
            "blood_share_of_chamber": blood, "shares": late_mix.moved_mix(own, cm, blood)}


def _commit_shares(monkeypatch, tmp_path, rec: dict, tree: bytes | None = None) -> tuple[Path, str]:
    """`rec` written as the shares file, standing as committed at HEAD (feng_relabel.committed_at_head answers its
    bytes for its path; the tree's copy is `tree` when given); returns its path and sha256."""
    raw = (json.dumps(rec, indent=1, sort_keys=True) + "\n").encode()
    path = tmp_path / "n64_late_mix" / "shares.json"
    path.parent.mkdir(exist_ok=True)
    path.write_bytes(raw if tree is None else tree)
    monkeypatch.setattr(_bt, "LATE_MIX_SHARES", path)
    inner = fr.committed_at_head
    monkeypatch.setattr(fr, "committed_at_head", lambda p, inner=inner: raw if Path(p) == path else inner(p))
    return path, hashlib.sha256(raw).hexdigest()


@pytest.fixture
def lm(bt, ms4, t1_env, ledger, tmp_path, monkeypatch):  # noqa: F811
    """Row R54's setup at shares moved to cardiomyocytes 0.45 of heart cells and blood 0.10 of chamber cells: the
    shares file committed and named in row N64. Returns (record, sha256, spy)."""
    rec = _shares_record(t1_env, 0.45, 0.10)
    _, sha = _commit_shares(monkeypatch, tmp_path, rec)
    ledger.write_text(_ledger(n64=_n64(sha)))
    return rec, sha, ms4


def _n(env) -> int:
    return harness.count_for(env.board, env.ref.n_obs)


def _quotas(shares: dict, n: int) -> dict:
    return {c: q for c, q in quotas(shares, n).items() if q > 0}


def _no_read(*a, **k):
    pytest.fail("nothing is read before the refusal")


def _refused(monkeypatch, tmp_path, argv: list[str], match: str) -> str:
    """Refused before any stage, board, union, Feng file or step is read (the ledger and the shares file only), with
    no number in the message but the allowed ones and what the command line said."""
    for obj, name in ((_bt, "load_stage"), (_bt, "load_boards"), (_bt.hygiene, "load_union"),
                      (_bt.feng_dirs, "feng_directions"), (_bt, "_t1_hm_step"), (_bt, "_t1_ms_step")):
        monkeypatch.setattr(obj, name, _no_read)
    name = "T1_dry.h5ad" if "--board" in argv else "T1.h5ad"
    with pytest.raises(SystemExit, match=match) as e:
        _bt.main([*argv, "--out", str(tmp_path / name)])
    msg = str(e.value)
    given = set(re.findall(r"\d+(?:\.\d+)?", " ".join(argv)))
    assert set(re.findall(r"\d+(?:\.\d+)?", msg)) <= NUMBERS_ALLOWED | given, msg
    return msg


def _calls(monkeypatch) -> dict:
    return r51t._calls(_bt, monkeypatch)


def _parse(argv: list[str]):
    return _bt.build_parser().parse_args([*argv, "--out", "x.h5ad"])


# ------------------------------------------------------------------ the rows' first cells and the constants


def test_the_rows_first_cells_as_the_ledger_holds_them_pass_the_builders_check(tmp_path):
    p = tmp_path / "TRIES.md"
    p.write_text(f"# ledger\n{ROW_R54}\n{ROW_R58}\n")
    assert _bt._captain_order_names("R54", _bt.T1_R54_OPTION, p)
    assert _bt._captain_order_names("R58", _bt.T1_R58_OPTION, p)
    assert _bt.T1_R54_OPTION == "`--late-mix e15_17`" and _bt.LATE_MIX_CHOICES == ("e15_17",)
    assert _bt.T1_R58_OPTION == "`--feng-fraction-mult 1.5 --direction-momentum all --direction-momentum-damping 2.0`"
    assert _bt.T1_R54_COMMAND == f"{BASE_CMD} {' '.join(LM16)}" and _bt.T1_R58_COMMAND == f"{BASE_CMD} {' '.join(M20)}"
    for bad in (ROW_R54.replace("`--late-mix e15_17`", "--late-mix e15_17"),
                ROW_R54.replace("`--captain-ruling R54`", "--captain-ruling R54"),
                ROW_R54.replace("ruling, no run", "a try"), ROW_R54.replace("R54 (the captain's order", "R54 (a draft"),
                f"{ROW_R54}\n{ROW_R54}"):
        p.write_text(f"# ledger\n{bad}\n")
        assert not _bt._captain_order_names("R54", _bt.T1_R54_OPTION, p), bad[:40]
    p.write_text(f"# ledger\n{ROW_R58.replace('damping 2.0`', 'damping 1.0`')}\n")
    assert not _bt._captain_order_names("R58", _bt.T1_R58_OPTION, p)


def test_the_options_are_off_by_default(bt):
    a = _parse(["t1", *RC])
    assert a.late_mix is None and a.direction_momentum_damping is None and a.captain_ruling is None
    assert bt._t1_r54(a) is None and bt._t1_r58(a) is None


def test_the_named_sha_is_read_from_row_n64_alone(tmp_path):
    p = tmp_path / "TRIES.md"
    sha, other = "a" * 64, "b" * 64
    for text, want in ((_n64(sha), sha), (_n64(None), None), (f"{_n64(sha)}\n{_n64(sha)}", None),
                       (_n64(sha).replace("written:", f"written: `shares.json sha256 {other}` and"), None),
                       (_n64(sha).replace("written:", f"again `shares.json sha256 {sha}`,"), sha),
                       (_n64(sha).replace("| N64 (", "| N640 ("), None),
                       (_n64(sha).replace(f"sha256 {sha}", f"sha256 {sha[:-1]}"), None)):
        p.write_text(f"# ledger\n{text}\n")
        assert _bt._late_mix_named_sha(p) == want, text[:30]


# ------------------------------------------------------------------ row R58's file


def test_m20_is_taken_with_its_row_and_routes_nowhere_else(bt, ledger):
    a = _parse(["t1", *RC, *M20])
    want = {"direction_momentum": "all", "captain_ruling": "R58", "file": "T1_s2rcMSBm20", "base_file": "T1_s2rcMSB",
            "option": bt.T1_R58_OPTION, "feng_fraction_mult": 1.5, "g_mom_damping": 0.5, "held_momentum": None,
            "direction_momentum_damping": 2.0, "order": bt.T1_R58_ORDER, "command": bt.T1_R58_COMMAND}
    assert bt._t1_r58(a) == want and bt._r58_ordered(a) and bt._t1_r54(a) is None
    assert bt._bet_refusals(a) is None and bt._feng_recomputed_refusals(a) is None and bt._t1_pooled(a) is None
    for other in (bt._ruled_fraction_ordered, bt._ruled_damping_ordered, bt._ms_ordered, bt._msb_ordered,
                  bt._pooled_ordered, bt._hm_ordered, bt._ms25_ordered, bt._r51_ordered, bt._r54_ordered):
        assert not other(a), other.__name__


def test_m20_builds_the_same_bytes_twice_and_is_neither_msb_nor_m10(bt, ms4, ledger, tmp_path):
    f0, _ = _build(bt, tmp_path, "T1_s2rcMSB", MSB)
    fa, ra = _build(bt, tmp_path, "m20a", M20)
    fb, rb = _build(bt, tmp_path, "m20b", M20)
    f10, _ = _build(bt, tmp_path, "m10", r51t.ARGS["T1_s2rcMSBm10"])
    assert fa.read_bytes() == fb.read_bytes() and r51t._record_sha(ra) == r51t._record_sha(rb)
    assert len({_sha(f0), _sha(fa), _sha(f10)}) == 3


def test_m20s_steps_are_msbs_but_ms_damping(bt, ms4, ledger, tmp_path, monkeypatch):
    seen = _calls(monkeypatch)
    _build(bt, tmp_path, "T1_s2rcMSB", MSB)
    _build(bt, tmp_path, "T1_s2rcMSBm20", M20)
    assert seen["ms"][0] == {"fraction": 1.5 * 0.4, "seed": SEED, "kw": {}}
    assert seen["ms"][1] == {"fraction": 1.5 * 0.4, "seed": SEED,
                             "kw": {"damping": 2.0, "option_note": bt.T1_R58_NOTE}}
    assert seen["g_mom"][0] == seen["g_mom"][1] and seen["g_mom"][1]["damping"] == 0.5
    assert seen["hm"] == []


def test_m20s_record_is_msbs_plus_its_candidate_field(bt, ms4, ledger, tmp_path):
    _, r0 = _build(bt, tmp_path, "T1_s2rcMSB", MSB)
    _, r1 = _build(bt, tmp_path, "T1_s2rcMSBm20", M20)
    skip = {"t1_candidate", "g_mom"}
    assert {k: v for k, v in r1.items() if k not in skip} == {k: v for k, v in r0.items() if k not in skip}
    g0, g1 = dict(r0["g_mom"]), dict(r1["g_mom"])
    g0.pop("median_depth_scale"), g1.pop("median_depth_scale")
    assert g0 == g1
    cand = r1["t1_candidate"]
    assert cand["file"] == "T1_s2rcMSBm20" and cand["captain_ruling"] == "R58" and cand["applied_fraction"] == 1.5 * 0.4
    m0, m1 = dict(r0["t1_candidate"]["momentum"]), dict(cand["momentum"])
    assert m1.pop("damping") == 2.0 and m0.pop("damping") == 0.5
    assert m1.pop("damping_option") == bt.T1_R58_NOTE and m1 == m0
    assert "R58" not in json.dumps(r0) and "late_mix" not in r1


def test_m20_with_ms_damping_at_b0s_is_msbs_bytes_by_the_checker(bt, ms4, ledger, tmp_path):
    chk = _checker()
    f0, _ = _build(bt, tmp_path, "T1_s2rcMSB", MSB)
    out = tmp_path / "m20_msb.h5ad"
    before = (bt._t1_r58, bt._t1_ms_step, bt._t1_g_mom)
    assert chk.main(["m20msb", "t1", *RC, *M20, "--out", str(out)]) == 0
    assert out.read_bytes() == f0.read_bytes()
    assert (bt._t1_r58, bt._t1_ms_step, bt._t1_g_mom) == before
    assert chk.main(["m20msb", "t1", *RC, *MSB, "--out", str(tmp_path / "x.h5ad")]) == 1  # MSB's route, not R58's


@pytest.mark.parametrize("ruling", [None, "R51", "R26", "R23", "R54", "R44"])
def test_damping_2_without_r58_is_refused_as_before(bt, ledger, monkeypatch, tmp_path, ruling):
    tail = [] if ruling is None else ["--captain-ruling", ruling]
    _refused(monkeypatch, tmp_path, ["t1", *RC, *PAIR, "--direction-momentum-damping", "2.0", *tail],
             "--direction-momentum-damping 2.0 is refused: it takes 1.0 only, and only in row R51's file")


@pytest.mark.parametrize("argv, match", [
    ([*PAIR, "--captain-ruling", "R58"], "--captain-ruling R58 orders one file, T1_s2rcMSBm20"),
    (["--feng-fraction-mult", "2.0", "--direction-momentum", "all", "--direction-momentum-damping", "2.0",
      "--captain-ruling", "R58"], "orders one file"),
    ([*PAIR, "--direction-momentum-damping", "1.0", "--captain-ruling", "R58"], "orders one file"),
    ([*PAIR, "--direction-momentum-damping", "2", "--captain-ruling", "R58"], "orders one file"),
    ([*PAIR, "--direction-momentum-damping", "2.00", "--captain-ruling", "R58"], "orders one file"),
    (["--feng-fraction-mult", "1.5", "--direction-momentum-damping", "2.0", "--captain-ruling", "R58"],
     "orders one file"),
    (["--captain-ruling", "R58"], "orders one file"),
    ([*M20, "--g-mom-damping", "1.0"], "refused with --g-mom-damping: it is T1_s2rcMSB's command with one change"),
    ([*M20, "--held-momentum", "0.5"], "refused with --held-momentum"),
    ([*M20, "--released-bet", "share_trend", "--bet-strength", "0.5"], "refused with --released-bet, --bet-strength"),
    ([*M20, "--late-mix", "e15_17"], "--late-mix e15_17 is row R54's file T1_s2rcMSB_LM16 alone and needs "
                                     "--captain-ruling R54, not R58"),
    ([*M20, "--no-hold-unbracketed"], "refused with --no-hold-unbracketed"),
    ([*M20, "--seed", "1000002"], "changes T1_s2rcMSB only: it needs the command"),
    ([*M20, "--board", "T1:val"], "changes T1_s2rcMSB only"),
    ([*M20, "--cm-share", "0.6"], "--feng-groups recomputed builds one file"),
    ([*M20, "--herd"], "--feng-groups recomputed builds one file"),
    ([*M20, "--base", "pooled"], "--feng-groups recomputed builds one file"),
])
def test_every_other_use_of_r58_is_refused(bt, ledger, monkeypatch, tmp_path, argv, match):
    _refused(monkeypatch, tmp_path, ["t1", *RC, *argv], match)


@pytest.mark.parametrize("cmd", [["--g-mom", "--base", "lib3"], [], ["--g-mom", "--base", "lib3", "--feng-groups",
                                                                       "none"]])
def test_m20_on_another_command_is_refused(bt, ledger, monkeypatch, tmp_path, cmd):
    _refused(monkeypatch, tmp_path, ["t1", *cmd, *M20], "row R58's file T1_s2rcMSBm20 changes T1_s2rcMSB only")


@pytest.mark.parametrize("row", [None, ROW_R58.replace("`--captain-ruling R58`", "--captain-ruling R58"),
                                 ROW_R58.replace("ruling, no run", "a try"), f"{ROW_R58}\n{ROW_R58}",
                                 ROW_R58.replace("damping 2.0`", "damping 2.0 x`")])
def test_m20_without_its_row_is_refused(bt, ledger, monkeypatch, tmp_path, row):
    ledger.write_text(_ledger(r58=row))
    _refused(monkeypatch, tmp_path, ["t1", *RC, *M20], "no single row R58 whose first cell is the captain's order")


# ------------------------------------------------------------------ row R54's file


def test_the_steps_own_shares_are_the_builders_lib3_mix(bt, ms4, t1_env, ledger, tmp_path, monkeypatch):  # noqa: F811
    """late_mix.own_shares (row N64's step) takes library 3's mix exactly as the builder's lib3 path passes it to
    build_submission."""
    got = []
    inner = bt.t1_scope.build_submission
    monkeypatch.setattr(bt.t1_scope, "build_submission", lambda *a, **k: (got.append(dict(a[2])), inner(*a, **k))[1])
    _build(bt, tmp_path, "T1_s2rcMSB", MSB)
    assert got == [late_mix.own_shares(t1_env.ref)]


def test_lm16_at_library_3s_own_values_has_msbs_quotas_and_bytes(bt, ms4, t1_env, ledger, tmp_path, monkeypatch):  # noqa: F811
    """Identity on quotas: the shares file at library 3's own two shares gives the same quotas as the own mix
    (with_share rescales floats, so the quotas are compared, then the bytes)."""
    rec = _shares_record(t1_env)
    _, sha = _commit_shares(monkeypatch, tmp_path, rec)
    ledger.write_text(_ledger(n64=_n64(sha)))
    n = _n(t1_env)
    assert _quotas(rec["shares"], n) == _quotas(rec["own_shares"], n)
    assert max(abs(rec["shares"].get(c, 0.0) - v) for c, v in rec["own_shares"].items()) < 1e-12
    f0, _ = _build(bt, tmp_path, "T1_s2rcMSB", MSB)
    f1, r1 = _build(bt, tmp_path, "T1_s2rcMSB_LM16_own", LM16)
    assert f1.read_bytes() == f0.read_bytes()
    assert r1["late_mix"]["median_library_over_msb_draw"] == 1.0


def test_lm16_draws_library_3_at_the_moved_mix_by_quotas(bt, lm, t1_env, tmp_path):  # noqa: F811
    rec, sha, spy = lm
    f0, r0 = _build(bt, tmp_path, "T1_s2rcMSB", MSB)
    fa, ra = _build(bt, tmp_path, "LM16a", LM16)
    fb, rb = _build(bt, tmp_path, "LM16b", LM16)
    assert fa.read_bytes() == fb.read_bytes() and r51t._record_sha(ra) == r51t._record_sha(rb)
    assert fa.read_bytes() != f0.read_bytes()
    drawn = spy[-1]["before"]  # the population LM16's post_cells gets: its draw
    cats, counts = np.unique(drawn.obs["category"].astype(str).to_numpy(), return_counts=True)
    assert dict(zip(cats.tolist(), counts.tolist(), strict=True)) == _quotas(rec["shares"], _n(t1_env))
    assert ra["shares"] == pytest.approx(rec["shares"], abs=1e-12)
    assert ra["cm_share_of_heart"] == pytest.approx(0.45, abs=1e-12)
    assert ra["blood_share_of_chamber"] == pytest.approx(0.10, abs=1e-12)
    late = ra["late_mix"]
    assert late["shares_sha256"] == sha and late["row"] == "N64" and late["anchor"] == "e15_17"
    assert late["own_shares"] == rec["own_shares"] and late["cm_share_of_heart"] == 0.45
    assert np.isfinite(late["median_library_over_msb_draw"]) and late["median_library_over_msb_draw"] > 0
    cand = ra["t1_candidate"]
    assert cand["file"] == "T1_s2rcMSB_LM16" and cand["captain_ruling"] == "R54" and cand["shares_sha256"] == sha
    assert "own_shares" not in cand and "late_shares" not in cand and cand["command"] == bt.T1_R54_COMMAND
    skip = {"t1_candidate", "g_mom", "late_mix", "shares", "cm_share_of_heart", "blood_share_of_chamber", "realised"}
    assert {k: v for k, v in ra.items() if k not in skip} == {k: v for k, v in r0.items() if k not in skip}
    assert set(ra) - set(r0) == {"late_mix"}


def test_lm16s_steps_after_the_draw_are_msbs(bt, lm, tmp_path, monkeypatch):
    seen = _calls(monkeypatch)
    _build(bt, tmp_path, "T1_s2rcMSB", MSB)
    _build(bt, tmp_path, "T1_s2rcMSB_LM16", LM16)
    assert seen["ms"][0] == seen["ms"][1] == {"fraction": 1.5 * 0.4, "seed": SEED, "kw": {}}
    assert seen["g_mom"][0] == seen["g_mom"][1] and seen["hm"] == []
    rec, _, spy = lm
    assert spy[0]["fraction"] is None and spy[1]["fraction"] is None  # the atlas step inside MS's step, as MSB's


def test_lm16_with_the_moved_mix_put_back_is_msbs_bytes_by_the_checker(bt, lm, tmp_path):
    chk = _checker()
    f0, _ = _build(bt, tmp_path, "T1_s2rcMSB", MSB)
    out = tmp_path / "lm16_own.h5ad"
    before = (bt._t1_r54, bt._t1_late_mix_shares, bt._t1_ms_step, bt._t1_g_mom)
    assert chk.main(["lm16own", "t1", *RC, *LM16, "--out", str(out)]) == 0
    assert out.read_bytes() == f0.read_bytes()
    assert (bt._t1_r54, bt._t1_late_mix_shares, bt._t1_ms_step, bt._t1_g_mom) == before
    assert chk.main(["lm16own", "t1", *RC, *MSB, "--out", str(tmp_path / "x.h5ad")]) == 1  # MSB's route, not R54's
    for argv in ([], ["lm16own"], ["lm16own", "heart"], ["x", "t1"]):
        with pytest.raises(SystemExit, match="usage"):
            chk.main(argv)


def test_lm16_refuses_a_pool_mix_other_than_the_recorded_one(bt, ms4, t1_env, ledger, tmp_path, monkeypatch):  # noqa: F811
    rec = _shares_record(t1_env, 0.45, 0.10)
    own = dict(rec["own_shares"])
    own["cm"], own["blood"] = own["cm"] + 1e-9, own["blood"] - 1e-9
    rec = {**rec, "own_shares": own, "shares": late_mix.moved_mix(own, 0.45, 0.10)}
    _, sha = _commit_shares(monkeypatch, tmp_path, rec)
    ledger.write_text(_ledger(n64=_n64(sha)))
    with pytest.raises(SystemExit, match="the pool's own chamber mix is not the one row N64's step recorded"):
        _build(bt, tmp_path, "T1_s2rcMSB_LM16", LM16)


def test_lm16_refuses_moved_shares_other_than_the_recorded_ones(bt, ms4, t1_env, ledger, tmp_path, monkeypatch):  # noqa: F811
    rec = _shares_record(t1_env, 0.45, 0.10)
    rec = {**rec, "shares": late_mix.moved_mix(rec["own_shares"], 0.46, 0.10)}
    _, sha = _commit_shares(monkeypatch, tmp_path, rec)
    ledger.write_text(_ledger(n64=_n64(sha)))
    with pytest.raises(SystemExit, match="the moved shares are not row N64's recorded shares"):
        _build(bt, tmp_path, "T1_s2rcMSB_LM16", LM16)


@pytest.mark.parametrize("argv, match", [
    ([*PAIR, "--late-mix", "e15_17"], "--late-mix e15_17 is row R54's file T1_s2rcMSB_LM16 alone and needs "
                                      "--captain-ruling R54"),
    ([*PAIR, "--late-mix", "e15_17", "--captain-ruling", "R26"], "needs --captain-ruling R54, not R26"),
    ([*PAIR, "--late-mix", "e15_17", "--captain-ruling", "R51"], "needs --captain-ruling R54, not R51"),
    (["--feng-fraction-mult", "2.0", "--direction-momentum", "all", "--late-mix", "e15_17", "--captain-ruling", "R54"],
     "moves T1_s2rcMSB's chamber mix only: it needs the command"),
    (["--feng-fraction-mult", "1.5", "--late-mix", "e15_17", "--captain-ruling", "R54"], "moves T1_s2rcMSB's chamber"),
    (["--direction-momentum", "all", "--late-mix", "e15_17", "--captain-ruling", "R54"], "moves T1_s2rcMSB's chamber"),
    (["--late-mix", "e15_17", "--captain-ruling", "R54"], "moves T1_s2rcMSB's chamber mix only"),
    ([*LM16, "--g-mom-damping", "1.0"], "refused with --g-mom-damping: row R54's file is T1_s2rcMSB's command with its "
                                        "chamber mix moved"),
    ([*LM16, "--held-momentum", "0.5"], "refused with --held-momentum"),
    ([*LM16, "--direction-momentum-damping", "1.0"], "refused with --direction-momentum-damping"),
    ([*LM16, "--cm-share", "0.6"], "refused with --cm-share"),
    ([*LM16, "--blood-share", "0.1"], "refused with --blood-share"),
    ([*LM16, "--lib3-pair"], "refused with --lib3-pair"),
    ([*LM16, "--herd"], "refused with --herd"),
    ([*LM16, "--strong-null"], "refused with --strong-null"),
    ([*LM16, "--released-bet", "share_trend", "--bet-strength", "0.5"], "refused with --released-bet, --bet-strength"),
    ([*LM16, "--no-hold-unbracketed"], "R1 and R2 stand as in T1_s2rcMSB"),
    ([*LM16, "--no-feng-hygiene"], "R1 and R2 stand as in T1_s2rcMSB"),
    ([*LM16, "--feng-none-ruling", "x"], "R1 and R2 stand as in T1_s2rcMSB"),
    ([*LM16, "--seed", "1000002"], "moves T1_s2rcMSB's chamber mix only"),
    ([*LM16, "--board", "T1:val"], "moves T1_s2rcMSB's chamber mix only"),
])
def test_every_other_use_of_the_late_mix_is_refused(bt, lm, monkeypatch, tmp_path, argv, match):
    _refused(monkeypatch, tmp_path, ["t1", *RC, *argv], match)


@pytest.mark.parametrize("cmd", [["--g-mom", "--base", "lib3"], ["--g-mom", "--base", "pooled", "--feng-groups",
                                                                   "recomputed"],
                                 ["--base", "lib3", "--feng-groups", "recomputed"],
                                 ["--g-mom", "--base", "mid", "--feng-groups", "recomputed"]])
def test_lm16_on_another_command_is_refused(bt, lm, monkeypatch, tmp_path, cmd):
    _refused(monkeypatch, tmp_path, ["t1", *cmd, *LM16], "moves T1_s2rcMSB's chamber mix only")


@pytest.mark.parametrize("r54", [None, ROW_R54.replace("`--late-mix e15_17`", "--late-mix e15_17"),
                                 ROW_R54.replace("ruling, no run", "a try"), f"{ROW_R54}\n{ROW_R54}"])
def test_lm16_without_its_row_is_refused(bt, lm, ledger, monkeypatch, tmp_path, r54):
    _, sha, _ = lm
    ledger.write_text(_ledger(r54=r54, n64=_n64(sha)))
    _refused(monkeypatch, tmp_path, ["t1", *RC, *LM16], "no single row R54 whose first cell is the captain's ruling")


def test_lm16_takes_only_the_shares_file_row_n64_names(bt, ms4, t1_env, ledger, tmp_path, monkeypatch):  # noqa: F811
    rec = _shares_record(t1_env, 0.45, 0.10)
    path, sha = _commit_shares(monkeypatch, tmp_path, rec)
    argv = ["t1", *RC, *LM16]
    for n64 in (None, _n64(None), _n64("c" * 64).replace("written:", f"or `shares.json sha256 {sha}`")):
        ledger.write_text(_ledger(n64=n64))
        _refused(monkeypatch, tmp_path, argv, "has no single row N64 naming exactly one")
    ledger.write_text(_ledger(n64=_n64("d" * 64)))
    _refused(monkeypatch, tmp_path, argv, "the committed shares file's sha256 is not the one row N64 names")
    ledger.write_text(_ledger(n64=_n64(sha)))
    path.write_bytes(path.read_bytes() + b" ")  # the tree's copy differs from HEAD's
    _refused(monkeypatch, tmp_path, argv, "the tree's copy of the shares file shares.json is not the one committed")
    path.unlink()
    _refused(monkeypatch, tmp_path, argv, "the tree's copy of the shares file")
    monkeypatch.setattr(fr, "committed_at_head", lambda p: None)
    _refused(monkeypatch, tmp_path, argv, "is not committed at HEAD")


@pytest.mark.parametrize("change", [{"row": "N65"}, {"late_mix": "e14"}, {"cm_share_of_heart": 1.0},
                                    {"blood_share_of_chamber": 0.0}, {"shares": {"cm": 0.5}}, {"own_shares": None}])
def test_lm16_refuses_a_shares_file_that_is_not_the_steps(bt, ms4, t1_env, ledger, tmp_path, monkeypatch, change):  # noqa: F811
    rec = {**_shares_record(t1_env, 0.45, 0.10), **change}
    _, sha = _commit_shares(monkeypatch, tmp_path, rec)
    ledger.write_text(_ledger(n64=_n64(sha)))
    _refused(monkeypatch, tmp_path, ["t1", *RC, *LM16], "the shares file is not row N64's step output for e15_17")


def test_the_cm_and_blood_share_options_stay_refused_on_lib3(bt, lm, monkeypatch, tmp_path):
    """Row R6's refusal is untouched: --cm-share and --blood-share on --base lib3 (without --late-mix) as before."""
    _refused(monkeypatch, tmp_path, ["t1", *RC, "--cm-share", "0.6"], "--feng-groups recomputed builds one file")
    with pytest.raises(SystemExit, match="--base lib3 is refused with --cm-share or --blood-share"):
        bt.main(["t1", "--base", "lib3", "--cm-share", "0.6", "--out", str(tmp_path / "T1.h5ad")])


# ------------------------------------------------------------------ the old builder against this one


OLD_COMMIT = "9a1b641"  # main as merged into this branch, without these options
OLD_SHA = "b4476d41ad295edf5e11601ddcb6365a58f19002d895efb66450cb35949473d7"  # its experiments/builds/build_test.py


class _Reached(Exception):
    """The command passed every refusal and reached its first read (the argument names which)."""


def _reach(where: str):
    def read(*a, **k):
        raise _Reached(where)

    return read


def _old_builder(tmp_path: Path):
    """The builder at OLD_COMMIT, read with `git show` (skipped where that commit is not in the checkout)."""
    env = {**git_env(), "GIT_OPTIONAL_LOCKS": "0"}
    res = subprocess.run(["git", "-C", str(ROOT), "show", f"{OLD_COMMIT}:experiments/builds/build_test.py"],
                         capture_output=True, env=env)
    if res.returncode != 0:
        pytest.skip(f"commit {OLD_COMMIT} is not available in this checkout")
    assert hashlib.sha256(res.stdout).hexdigest() == OLD_SHA
    path = tmp_path / "old" / "build_test_before_lm_m20.py"
    path.parent.mkdir()
    path.write_bytes(res.stdout)
    spec = importlib.util.spec_from_file_location("build_test_before_lm_m20", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


CMDS = [RC, ["--g-mom", "--base", "lib3"], ["--g-mom", "--base", "pooled", "--feng-groups", "recomputed"], []]
MULTS = [None, "1.5", "2.0"]
DMS = [None, "all"]
DMD = [None, "1.0", "2.0"]
LATE = [None, "e15_17"]
RULINGS = [None, "R54", "R58", "R51", "R26", "R23", "R22"]
EXTRAS = [[], ["--seed", "1000002"], ["--cm-share", "0.6"], ["--g-mom-damping", "1.0"], ["--held-momentum", "0.5"]]
GRID_SETS = 4 * 3 * 2 * 3 * 2 * 7 * 5


def _grid() -> list[list[str]]:
    out = []
    for cmd, mult, dm, dmd, late, ruling, extra in itertools.product(CMDS, MULTS, DMS, DMD, LATE, RULINGS, EXTRAS):
        argv = ["t1", *cmd]
        for flag, value in (("--feng-fraction-mult", mult), ("--direction-momentum", dm),
                            ("--direction-momentum-damping", dmd), ("--late-mix", late), ("--captain-ruling", ruling)):
            argv += [] if value is None else [flag, value]
        out.append([*argv, *extra])
    return out


def _kind(argv: list[str], named: set[str]) -> str:
    """"new54" / "new58" for exactly row R54's / R58's command under a ledger naming the row, "touch" for any other set
    giving --late-mix or --captain-ruling R58, "other" for every other set."""
    tail = argv[len(["t1", *RC]):] if argv[:len(RC) + 1] == ["t1", *RC] else None
    if tail == LM16 and "R54" in named:
        return "new54"
    if tail == M20 and "R58" in named:
        return "new58"
    if "--late-mix" in argv or ("--captain-ruling" in argv and argv[argv.index("--captain-ruling") + 1] == "R58"):
        return "touch"
    return "other"


def test_the_old_builder_and_this_one_route_every_other_argument_set_alike(bt, tmp_path, monkeypatch):
    """Side by side on the grid, under three synthetic ledgers (rows R54 and R58 as main holds them with row N64
    naming the shares file; without them; with R58 alone): every set that gives neither --late-mix nor --captain-ruling
    R58 is refused with the same message, or reaches the same first read with the same candidate, MS, pooled and row
    R51 values, in both builders, and this one's R54 and R58 routes return None for it; R54's and R58's commands under
    their rows are taken by this one with their values (the old builder refuses each); every other set giving --late-mix
    or --captain-ruling R58 is refused by both. Nothing is built: every read is stopped where it would start."""
    old = _old_builder(tmp_path)
    assert not hasattr(old, "_t1_r54") and hasattr(bt, "_t1_r54") and hasattr(old, "_t1_r51")
    ledger = tmp_path / "TRIES.md"
    record = tmp_path / "relabel" / "MANIFEST.json"  # never written: the relabelling's record is the stand-in below
    shares = {"row": "N64", "late_mix": "e15_17", "own_shares": {"cm": 0.5, "blood": 0.5}, "cm_share_of_heart": 0.6,
              "blood_share_of_chamber": 0.4, "shares": {"cm": 0.6, "blood": 0.4}}
    raw = (json.dumps(shares, sort_keys=True) + "\n").encode()
    path = tmp_path / "n64_late_mix" / "shares.json"
    path.parent.mkdir()
    path.write_bytes(raw)
    sha = hashlib.sha256(raw).hexdigest()
    monkeypatch.setattr(bt, "LATE_MIX_SHARES", path)
    monkeypatch.setattr(fr, "RECORD", record)
    monkeypatch.setattr(fr, "committed_at_head", lambda p: json.dumps({"route": "recomputed"}).encode()
                        if Path(p) == record else raw if Path(p) == path else None)
    monkeypatch.setattr(fr, "later_ruling", lambda *a, **k: "no later ruling (synthetic)")
    monkeypatch.setattr(bt.hygiene, "load_union", _reach("union"))
    monkeypatch.setattr(bt.strong_null, "check_approved", _reach("strong_null"))
    monkeypatch.setattr(bt.feng_dirs, "feng_directions", _reach("feng"))
    spies, parsers = {}, {}
    for mod in (old, bt):
        parsers[mod] = mod.build_parser()
        monkeypatch.setattr(mod, "build_parser", lambda p=parsers[mod]: p)
        monkeypatch.setattr(mod, "TRIES_FILE", ledger)
        monkeypatch.setattr(mod, "load_stage", _reach("stage"))
        monkeypatch.setattr(mod, "load_boards", _reach("board"))
        new = ("_t1_r54", "_t1_r58") if mod is bt else ()
        for name in ("_t1_ms", "_t1_candidate", "_t1_pooled", "_t1_r51", *new):
            inner = getattr(mod, name)
            monkeypatch.setattr(mod, name, lambda a, inner=inner, name=name: spies.setdefault(name, inner(a)))
    out = tmp_path / "out"
    out.mkdir()

    def route(mod, argv):
        spies.clear()
        name = "T1_dry.h5ad" if "--board" in argv else "T1.h5ad"
        try:
            mod.main([*argv, "--out", str(out / name)])
        except SystemExit as e:
            return ("refused", str(e.code)), dict(spies)
        except _Reached as e:
            return (("taken", e.args[0], spies.get("_t1_ms"), spies.get("_t1_candidate"), spies.get("_t1_pooled"),
                     spies.get("_t1_r51")), dict(spies))
        raise AssertionError(f"{argv} neither refused nor read")

    n64 = _n64(sha)
    ledgers = {"after": (_ledger(n64=n64), {"R54", "R58"}), "before": (_ledger(r54=None, r58=None), set()),
               "r58_only": (_ledger(r54=None, n64=n64), {"R58"})}
    grid = _grid()
    assert len({tuple(x) for x in grid}) == len(grid) == GRID_SETS
    count = {"other": 0, "touch": 0, "new54": 0, "new58": 0, "taken_alike": 0, "refused_alike": 0}
    with open(tmp_path / "stderr.txt", "w") as err, pytest.MonkeyPatch.context() as quiet:
        quiet.setattr(sys, "stderr", err)  # the old parser's usage errors on --late-mix go to a file
        for text, named in ledgers.values():
            ledger.write_text(text)
            for argv in grid:
                kind = _kind(argv, named)
                (o, _), (n, seen) = route(old, argv), route(bt, argv)
                count[kind] += 1
                if kind == "other":
                    assert n == o, argv
                    assert seen.get("_t1_r54") is None and seen.get("_t1_r58") is None, argv
                    count["taken_alike" if o[0] == "taken" else "refused_alike"] += 1
                    continue
                assert o[0] == "refused", argv  # the old builder took none of these
                if kind == "new54":
                    assert n == ("taken", "union", None, None, None, None), argv
                    assert seen["_t1_r54"]["file"] == "T1_s2rcMSB_LM16" and seen["_t1_r54"]["shares_sha256"] == sha
                elif kind == "new58":
                    assert n == ("taken", "union", None, None, None, None), argv
                    assert seen["_t1_r58"]["file"] == "T1_s2rcMSBm20"
                else:
                    assert n[0] == "refused", argv
    assert sum(count[k] for k in ("other", "touch", "new54", "new58")) == 3 * GRID_SETS
    assert count["taken_alike"] and count["refused_alike"] and count["touch"]
    assert count["new54"] == 1 and count["new58"] == 2  # R54's under "after"; R58's under "after" and "r58_only"


def _checker():
    import t1_lm_m20_check

    return t1_lm_m20_check


def test_the_checker_prints_nothing_and_writes_nothing_of_its_own():
    src = (ROOT / "experiments" / "builds" / "t1_lm_m20_check.py").read_text()
    body = src.split('"""', 2)[2]
    code = "\n".join(line.split("#")[0] for line in body.splitlines())
    for word in ("print(", "unlink", "rmtree", "os.remove", "write_text", "write_bytes", "upload", "submit",
                 ".build.json", "read_text"):
        assert word not in code, word
