"""The test-board build launcher of row R60's seven files (T1's line 12, each one change to T1_s2rcMS20's command) and
row R46's pair T1_s2rcB25MSH05.

experiments/builds/t1_l12_line.sh first rebuilds twenty-one known files (t1_r51_line.sh's nineteen, T1_s2rcMS20 and
T1_s2rcB25MS) and stops unless each gives its known bytes; then it builds each new command once on its own path with
its one setting at its base's value (experiments/builds/t1_l12_check.py base) and stops unless each gives its base's
bytes (T1_s2rcMS20's for row R60's, T1_s2rcB25MS's for row R46's); then it builds the eight files twice each from one
clean commit and compares each pair. These tests hold its command strings, hashes and step order to that, check that
it deletes nothing and reads no output of a test build but exit codes, sha256 and the stamps' commit and clean fields,
and run its gate's paths for real in throwaway repositories under a temporary home, with a stand-in `uv` that fails,
so that a call passing the gate stops at the tests step and builds nothing. Nothing here reads the repository's own
ledger or rules file for its state of today; build_test is imported only to check that it takes the commands the
launcher runs, against a temporary ledger."""

import json
import re
import subprocess
import sys
from pathlib import Path

import pytest
import test_t1_r51_launcher as r51l
from test_t1_l12_build import REAL_R46, REAL_R60
from test_t1_l67pb_launcher import _sums

REPO = Path(__file__).resolve().parents[1]
LAUNCHER = REPO / "experiments" / "builds" / "t1_l12_line.sh"
BASE = "--g-mom --base lib3 --feng-groups recomputed"
MS20_OPT = "--feng-fraction-mult 2.0 --direction-momentum all"
# step -> (file name, args variable, the args, hash variable, sha256, exit code): t1_r51_line.sh's nineteen, then
# T1_s2rcMS20 (row R51's file) and T1_s2rcB25MS (row R45's file)
KNOWN = dict(r51l.KNOWN)
KNOWN |= {
    "MS20": ("T1_s2rcMS20", "MS20", f"{BASE} {MS20_OPT} --captain-ruling R51", "MS20_SHA",
             "629d3a1cffa8151f0345a169ba97420e4ee3b43e379bc3bf7d33ab947dcb62ef", 39),
    "B25MS": ("T1_s2rcB25MS", "B25MS", f"{BASE} --feng-fraction-mult 2.5 --direction-momentum all --captain-ruling R45",
              "B25MS_SHA", "3127a2d9c92aa064e3f0fdf400447a074b7964ceacf2274da8a14ada7425f411", 50),
}
# file -> (args variable, step stem, its args)
NEW = {"T1_s2rcMS20d10": ("D10", "d10", "$L12 --g-mom-damping 1.0 --captain-ruling R60"),
       "T1_s2rcMS20H05": ("H05", "H05", "$L12 --held-momentum 0.5 --captain-ruling R60"),
       "T1_s2rcMS20d20": ("D20", "d20", "$L12 --g-mom-damping 2.0 --captain-ruling R60"),
       "T1_s2rcMS20m10": ("M10", "m10", "$L12 --direction-momentum-damping 1.0 --captain-ruling R60"),
       "T1_s2rcMS20d15": ("D15", "d15", "$L12 --g-mom-damping 1.5 --captain-ruling R60"),
       "T1_s2rcMS20H10": ("H10", "H10", "$L12 --held-momentum 1.0 --captain-ruling R60"),
       "T1_s2rcMS20m20": ("M20", "m20", "$L12 --direction-momentum-damping 2.0 --captain-ruling R60"),
       "T1_s2rcB25MSH05": ("PAIR", "pair", f"{BASE} --feng-fraction-mult 2.5 --direction-momentum all --held-momentum "
                                           "0.5 --captain-ruling R46")}
OPTS = {f: args.replace("$L12", f"{BASE} {MS20_OPT}").replace(f"{BASE} ", "").rsplit(" --captain-ruling", 1)[0]
        for f, (_, _, args) in NEW.items()}
# new step -> (file name, args variable, base hash variable, exit code when it is not the base's bytes)
BASECHECKS = {f"{stem}base": (f"{f}_base", key, "B25MS_SHA" if key == "PAIR" else "MS20_SHA", 51 + i)
              for i, (f, (key, stem, _)) in enumerate(NEW.items())}
PAIRS = {(f"{stem}a", f"{stem}b"): (f, key, 41 + i) for i, (f, (key, stem, _)) in enumerate(NEW.items())}
STEPS = ["gate", "tests", *KNOWN, *BASECHECKS, *(s for pair in PAIRS for s in pair), "done"]


def _sh() -> str:
    return LAUNCHER.read_text()


def _args(sh: str) -> dict[str, str]:
    return dict(re.findall(r'^([A-Z0-9]+)_ARGS="([^"]+)"$', sh, flags=re.M))


def test_bash_parses_the_launcher():
    out = subprocess.run(["bash", "-n", str(LAUNCHER)], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr


def test_the_launcher_reproduces_21_files_checks_the_bases_bytes_then_builds_each_new_file_twice():
    sh = _sh()
    assert f'L12="{BASE} {MS20_OPT}"' in sh
    want = {key: args for _, key, args, _, _, _ in KNOWN.values()} | {key: a for key, _, a in NEW.values()}
    assert _args(sh) == want
    for step, (name, key, _, var, sha, code) in KNOWN.items():
        assert f'identity {step} {name} "${var}" {code} ${key}_ARGS' in sh
        assert f"{var}={sha}" in sh and sh.count(sha) == 1
    for step, (name, key, var, code) in BASECHECKS.items():
        assert f'base_check {step} {name} "${var}" {code} ${key}_ARGS' in sh
    for (a, b), (name, key, code) in PAIRS.items():
        assert f'for STEP in {a} {b}; do begin; build "$STEP" {name} ${key}_ARGS || finish $?; step_ok; done' in sh
        assert f'cmp -s "$O/{a}/{name}.h5ad" "$O/{b}/{name}.h5ad" || finish {code}' in sh
    marks = ("STEP=gate", "STEP=tests", *(f"identity {k} " for k in KNOWN), *(f"base_check {k} " for k in BASECHECKS),
             *(f"for STEP in {a} {b}" for a, b in PAIRS), "STEP=done")
    order = [sh.index(s) for s in marks]
    assert order == sorted(order)
    listing = sh[sh.index("(cd \"$O\" && sha256sum"):sh.index('> "$S/SHA256SUMS"')]
    assert re.findall(r"([A-Za-z0-9]+)/\*\.h5ad", listing) == STEPS[2:-1]  # every built folder, once, in order
    assert 'got=$(sha256sum "$O/$1/$2.h5ad" | cut -d\' \' -f1)' in sh
    assert "OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16" in sh
    assert sh.index("NUM_THREADS=8") < sh.index("STEP=tests") < sh.index("NUM_THREADS=16") < sh.index("identity D ")
    assert "NUMPY_MADVISE_HUGEPAGE=0" in sh and 'TMPDIR="$S/tmp"' in sh
    first_uv = sh.index("uv run")
    assert sh.index("export PATH=$HOME/.local/bin:$PATH") < first_uv and ".light" not in sh
    assert sh.index("unset UV_PROJECT_ENVIRONMENT") < first_uv
    assert 'PY="uv run --frozen --no-sync python -I -B"' in sh
    assert '[ "$head" = "$COMMIT" ]' in sh and '[ "$st" = "$head True" ] || return 15' in sh
    assert "--seed" not in sh and "--board T1" not in sh and "vec check" in sh and "BOARD=T1:test" in sh
    assert "upload" not in sh.replace("Nothing here uploads", "").replace("the uploads are gated", "")
    assert "cd " not in sh.replace('(cd "$O" && sha256sum', "")
    assert "S2_RULES.md" not in sh.replace("docs/S2_RULES.sha256", "")


def test_the_base_builds_go_through_the_checker_and_every_build_is_checked_and_stamped():
    sh = _sh()
    base = sh[sh.index("build_base() {"):sh.index("same() {")]
    assert ('$PY experiments/builds/t1_l12_check.py base t1 "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 12'
            in base) and 'checked "$d" "$name"' in base
    build = sh[sh.index("build() {"):sh.index("build_base() {")]
    assert ('$PY experiments/builds/build_test.py t1 "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 11'
            in build) and 'checked "$d" "$name"' in build


def test_the_two_new_pins_are_the_committed_lists():
    """T1_s2rcMS20's and T1_s2rcB25MS's hashes as the committed run lists hold them, built twice with one sha256."""
    r51 = _sums("t1r51_0cdcc6f_r1.SHA256SUMS")
    hm = _sums("t1hm_20d5129_h1.SHA256SUMS")
    assert r51["MS20a/T1_s2rcMS20.h5ad"] == r51["MS20b/T1_s2rcMS20.h5ad"] == KNOWN["MS20"][4]
    assert hm["MS25a/T1_s2rcB25MS.h5ad"] == hm["MS25b/T1_s2rcB25MS.h5ad"] == KNOWN["B25MS"][4]
    r51_sh = (REPO / "experiments" / "builds" / "t1_r51_line.sh").read_text()
    assert f'MS20_ARGS="{KNOWN["MS20"][2]}"' in r51_sh
    hm_sh = (REPO / "experiments" / "builds" / "t1_hm_line.sh").read_text()
    assert f'MS25_ARGS="{KNOWN["B25MS"][2].replace("R45", "$ROW_M")}"' in hm_sh


def test_the_builder_takes_the_commands_the_launcher_runs(tmp_path, monkeypatch):
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    sys.path.insert(0, str(REPO / "experiments" / "external"))
    import build_test as bt
    import feng_relabel as fr
    from test_t1_hm_build import ROW_M

    ledger = tmp_path / "TRIES.md"
    ledger.write_text(r51l.ROWS_ALL + r51l.ROW_R51 + f"{ROW_M}\n| {REAL_R60} | T1 |\n"
                      f"| {REAL_R46} | T1 |\n")
    monkeypatch.setattr(bt, "TRIES_FILE", ledger)
    monkeypatch.setattr(fr, "committed_at_head", lambda p: json.dumps({"route": "recomputed"}).encode())
    args = {k: v.replace("$L12", f"{BASE} {MS20_OPT}") for k, v in _args(_sh()).items()}
    for _, key, *_ in KNOWN.values():
        a = bt.build_parser().parse_args(["t1", *args[key].split(), "--out", "x.h5ad"])
        assert bt._t1_l12(a) is None  # no earlier file takes row R60's or R46's route
    for f, (key, _, _) in NEW.items():
        a = bt.build_parser().parse_args(["t1", *args[key].split(), "--out", "x.h5ad"])
        assert a.board is None and a.seed == bt.harness.BUILD_SEED_OFFSET + 1
        got = bt._t1_r51(a)  # nothing read but the ledger
        assert got["file"] == f and got["option"] == f"`{OPTS[f]}`"
        assert got["captain_ruling"] == ("R46" if key == "PAIR" else "R60")
        assert bt._bet_refusals(a) is None and bt._feng_recomputed_refusals(a) is None and bt._t1_pooled(a) is None


def test_the_tests_step_runs_the_new_tests_and_t1_r51_line_shs():
    sh = _sh()
    step = sh[sh.index("STEP=tests"):sh.index("step_ok", sh.index("STEP=tests"))]
    listed = re.findall(r"tests/test_[a-z0-9_]+\.py", step)
    r51 = (REPO / "experiments" / "builds" / "t1_r51_line.sh").read_text()
    r51_step = r51[r51.index("STEP=tests"):r51.index("step_ok", r51.index("STEP=tests"))]
    assert {"tests/test_t1_l12_launcher.py", "tests/test_t1_l12_build.py",
            *re.findall(r"tests/test_[a-z0-9_]+\.py", r51_step)} == set(listed)
    assert len(listed) == len(set(listed)) and all((REPO / f).exists() for f in listed)
    assert "|| finish 10" in step and '--basetemp="$S/tmp/pytest"' in step


def test_the_launcher_deletes_nothing_and_opens_no_quote_in_a_usage_message():
    code = "\n".join(line.split("#")[0] for line in _sh().splitlines()).splitlines()
    assert not any(re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate|mv)\s", line) for line in code)
    assert not any(re.search(r"find .*-delete|git (clean|reset|checkout|stash)", line) for line in code)
    usages = re.findall(r"\$\{\d+:\?([^}]*)\}", _sh())
    assert len(usages) == 2
    for msg in usages:
        assert "'" not in msg and '"' not in msg, msg


def test_the_documented_exit_codes_are_the_ones_used():
    sh = _sh()
    head = sh[:sh.index("set -u")]
    used = {int(c) for c in re.findall(r"(?:finish|exit) (\d+)", sh)}
    used |= {code for *_, code in KNOWN.values()} | {code for *_, code in BASECHECKS.values()}
    used |= {code for *_, code in PAIRS.values()}
    for c in used:
        assert re.search(rf"(^|\D){c}(\D|$)", head), c


# ------------------------------------------------------------------ the gate, run for real


ROW_R60 = f"| {REAL_R60} | T1 | line 12 |\n"
ROW_R46 = f"| {REAL_R46} | T1 | the pair |\n"


def _repo(tmp_path: Path, tries: str = ROW_R60 + ROW_R46, **kw):
    return r51l._repo(tmp_path, tries=tries, launcher=_sh(), **kw)


def _scratch(tmp_path: Path, run: str) -> Path:
    return tmp_path / "home" / "v" / "scratch" / "t1l12" / run


def _stopped_at_the_gate(tmp_path: Path, run: str, code: int) -> None:
    s = _scratch(tmp_path, run)
    assert (s / "gate.exit").read_text().strip() == str(code)
    assert sorted(p.name for p in s.iterdir()) == ["gate.exit", "tmp"]
    assert not list((s / "tmp").iterdir())


def test_a_call_that_passes_the_gate_writes_the_head_and_stops_at_the_failing_tests_step(tmp_path):
    repo, env, head = _repo(tmp_path)
    out = r51l._run(repo, env, head, "r1")
    assert out.returncode == 10, out.stderr
    s = _scratch(tmp_path, "r1")
    assert (s / "gate.exit").read_text().strip() == "0" and (s / "tests.exit").read_text().strip() == "10"
    assert (s / "head").read_text().strip() == head
    assert sorted(p.name for p in s.iterdir()) == ["gate.exit", "head", "tests.exit", "tests.out", "tmp"]
    o = repo / "runs" / f"t1l12_{head[:7]}_r1"
    assert o.is_dir() and not list(o.iterdir())
    again = r51l._run(repo, env, head, "r1")
    assert again.returncode == 91 and "never reused" in again.stderr


def test_a_wrong_commit_a_changed_file_and_odd_names_are_refused(tmp_path):
    repo, env, head = _repo(tmp_path)
    out = r51l._run(repo, env, "0" * 40, "r1")
    assert out.returncode == 96 and "HEAD is not" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 96)
    assert r51l._run(repo, env, head[:12], "r2").returncode == 96
    for bad in ("a/b", "r 1", "$(id)"):
        assert r51l._run(repo, env, head, bad).returncode == 90
    assert r51l._run(repo, env, head).returncode != 0
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"# ledger\n{ROW_R60}{ROW_R46}more\n")
    out = r51l._run(repo, env, head, "r3")
    assert out.returncode == 98 and "tracked files differ" in out.stderr
    _stopped_at_the_gate(tmp_path, "r3", 98)


OPT_STRINGS = [o for f, o in OPTS.items() if f != "T1_s2rcB25MSH05"]


@pytest.mark.parametrize("tries, code, says", [
    ("", 95, "exactly one row R60"),
    (ROW_R46, 95, "exactly one row R60"),
    (ROW_R60 + ROW_R60 + ROW_R46, 95, "exactly one row R60"),
    (ROW_R60.replace("R60 (the captain's order", "R60 (a proposal") + ROW_R46, 95, "not the captain's ruling"),
    (ROW_R60.replace("ruling, no run", "a try") + ROW_R46, 95, "does not say ruling, no run"),
    *[(ROW_R60.replace(f"`{o}`", f"`{o} --herd`") + ROW_R46, 95, f"does not name `{o}`") for o in OPT_STRINGS],
    (ROW_R60.replace("`--captain-ruling R60`", "`--captain-ruling R51`") + ROW_R46, 95,
     "does not name `--captain-ruling R60`"),
    (ROW_R60, 94, "exactly one row R46"),
    (ROW_R60 + ROW_R46 + ROW_R46, 94, "exactly one row R46"),
    (ROW_R60 + ROW_R46.replace("R46 (the captain's", "R46 (a proposal of"), 94, "not the captain's ruling"),
    (ROW_R60 + ROW_R46.replace("ruling, no run", "a try"), 94, "does not say ruling, no run"),
    (ROW_R60 + ROW_R46.replace("T1_s2rcB25MSH05", "T1_s2rcB25MSH05x"), 94, "does not name T1_s2rcB25MSH05"),
    (ROW_R60 + ROW_R46.replace("T1_s2rcB25MSH05", "T1_s2rcB25MS"), 94, "does not name T1_s2rcB25MSH05"),
])
def test_each_gate_refusal_on_the_two_rows(tmp_path, tries, code, says):
    repo, env, head = _repo(tmp_path, tries=tries)
    out = r51l._run(repo, env, head, "r1")
    assert out.returncode == code, out.stderr
    assert says in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", code)
    assert not (repo / "runs").exists()


@pytest.mark.parametrize("kw, code", [({"good_sha": False}, 97), ({"relabel": False}, 92), ({"venv": None}, 88),
                                      ({"venv": "other"}, 88)])
def test_the_other_gate_refusals(tmp_path, kw, code):
    repo, env, head = _repo(tmp_path, **kw)
    out = r51l._run(repo, env, head, "r1")
    assert out.returncode == code, out.stderr
    _stopped_at_the_gate(tmp_path, "r1", code)
