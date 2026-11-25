"""experiments/builds/heart_trend_line_h.sh, the test-board builds of the heart's trend-form line (row R55, the
captain's order of 2026-11-24; the row fixed, no ROW argument), each form only when row H35's outcome cell, read by
code, says KEPT for it. Every call runs in a throwaway git repository under the test's own folder, with a temporary
home and a stand-in uv (and, where named, a stand-in sha256sum) first on PATH: no launcher runs with a real commit,
nothing reads a competition file and nothing is built. The uv stand-in is tests/test_heart_trend_line_f_launcher's;
the sha256sum stand-in gives the six identity steps' committed hashes. The gate's two inline readers are held to the
builder's R55 check and to the screen's own reader of row H35 on the same ledgers."""

from __future__ import annotations

import importlib.util
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_heart_organ_launchers import _calls, _env, _git_fn  # noqa: E402
from test_heart_trend_line_f_launcher import STUB_UV_F  # noqa: E402

LAUNCHER = REPO_ROOT / "experiments" / "builds" / "heart_trend_line_h.sh"
G_LAUNCHER = REPO_ROOT / "experiments" / "builds" / "heart_trend_line_g.sh"
G_LIST = REPO_ROOT / "experiments" / "builds" / "results" / "heartg_b2bbcf5_g1.SHA256SUMS"
FIRST = "H35 (the counted screen of row R55's trend forms: ten tries; G0 first, target-free, not a try)"
FORMS = ("category", "contrast", "whole", "portion")
TAGS = {"category": "tk", "contrast": "tx", "whole": "tw", "portion": "tp"}
# the six identity steps: step -> (file name, its sha256, its half's sha256), as heartg_b2bbcf5_g1.SHA256SUMS gives them
IDENTITY = {
    "D": ("T2heart_organ", "76e21cd2c754362c043ec32ca863a1e1b7f2cfde549e0366d6104354bfffe18b",
          "1b4342c151c19dfc827d4b7a310c6a0fa8908f1e369836ebffb22c27c93b4f00"),
    "OTS05": ("T2heart_organ_ts05", "c3b5fc9c7b47d0d0aafc4e5876bac7f06988cb2d628a0030b4b643806da8bceb",
              "f918e76731aac18eda3b46c8a35c18f26ae00a1a6a1db36d33e7515704cc43ea"),
    "OTS10": ("T2heart_organ_ts10", "efaaabe0b744d82dc40135ed9269a6baf8938c9cdd83f85fb3917c0b102e7a2b",
              "ee9e0e68d62e8d1448f1dbbbce1146dc5fb4877e94d77fc289fa82744acb06a9"),
    "OTS075": ("T2heart_organ_ts075", "109f7f26966c4a3cd0cd18807032ea8eb8f7786df833cd7b3064a21298b1cdbb",
               "1a66b5a5a9d5388d806f548d5c51570345f02e340a6ef2d4ce3490d1e2673f5d"),
    "OTS125": ("T2heart_organ_ts125", "79d0b813db758ce88bda22ed1f8c8b58a560214360bfbfc2795c7eee372419da",
               "adc3243665949f9802056216d86443e2ea0d45d07534013f99e26afb7e57ade0"),
    "OTS15": ("T2heart_organ_ts15", "ade351c1627d0eb5bdbc43cb613cf1b914053e19e4efca0dd943f939b1035d08",
              "f818931834d08a6906e010e8af9f077ba15be779aba8965cbaa880a9bdf1196c"),
}
STOP_CODES = {"D": (20, 21), "OTS05": (22, 23), "OTS10": (24, 25), "OTS075": (26, 27), "OTS125": (28, 29),
              "OTS15": (32, 33)}
ARGS = {"D": "--entry organ_scope --expression copy --z-mode raw",
        "OTS05": "--entry organ_scope --trend-strength 0.5 --z-mode raw --captain-ruling R29",
        "OTS10": "--entry organ_scope --trend-strength 1.0 --z-mode raw --captain-ruling R29",
        "OTS075": "--entry organ_scope --trend-strength 0.75 --z-mode raw --captain-ruling R49",
        "OTS125": "--entry organ_scope --trend-strength 1.25 --z-mode raw --captain-ruling R49",
        "OTS15": "--entry organ_scope --trend-strength 1.5 --z-mode raw --captain-ruling R52",
        **{f: f"--entry organ_scope --trend-strength 1.25 --trend-form {f} --z-mode raw --captain-ruling R55"
           for f in FORMS}}
R55_FIRST = ("R55 (the captain's order of 2026-11-24, ruling, no run: a synthetic stand-in; the heart builder is to "
             "take `--entry organ_scope` with `--trend-strength 1.25` and `--trend-form category`, "
             "`--trend-form contrast`, `--trend-form whole` or `--trend-form portion` only with "
             "`--captain-ruling R55`; not a try)")
R55_ROW = f"| {R55_FIRST} | heart | x |\n"


def outcome(words: dict) -> str:
    kept = any(w == "KEPT" for w in words.values())
    return ("OUTCOME H35: " + ("" if kept else "STOP: ") + "; ".join(f"{f} {words[f]}" for f in FORMS) + ".")


def h35_row(cell: str) -> str:
    return f"| {FIRST} | heart | x | x | x | x | x | x | x | x | {cell} |\n"


ALL_KEPT = h35_row(f'**Run.** In the read\'s fixed words: "{outcome({f: "KEPT" for f in FORMS})}" "A screen."')

# A stand-in sha256sum that gives the committed hashes for the identity steps' files and halves (each one replaceable
# through the environment: STUB_SHA_<STEP> and STUB_SHA_<STEP>_HALF) and the real hash otherwise.
STUB_SHA_HEAD = r"""#!/usr/bin/env bash
REAL=__REAL__
if [ "${1:-}" = -c ]; then exec "$REAL" "$@"; fi
for f in "$@"; do
  case "$(basename "$f")" in
"""
STUB_SHA_TAIL = r"""    *) "$REAL" "$f" ;;
  esac
done
"""


def _stub_sha() -> str:
    cases = "".join(f'    {name}.h5ad) echo "${{STUB_SHA_{step}:-{full}}}  $f" ;;\n'
                    f'    {name}_half.h5ad) echo "${{STUB_SHA_{step}_HALF:-{half}}}  $f" ;;\n'
                    for step, (name, full, half) in IDENTITY.items())
    return STUB_SHA_HEAD.replace("__REAL__", shutil.which("sha256sum")) + cases + STUB_SHA_TAIL


def _bin(tmp_path, sha_stub=False):
    b = tmp_path / ("bin_sha" if sha_stub else "bin")
    b.mkdir(exist_ok=True)
    (b / "uv").write_text(STUB_UV_F)
    (b / "uv").chmod(0o755)
    if sha_stub:
        (b / "sha256sum").write_text(_stub_sha())
        (b / "sha256sum").chmod(0o755)
    return b


def _ledger(h35: str = ALL_KEPT, r55: str = R55_ROW) -> str:
    return "| R29 (other) | heart | x |\n| R49 (other) | heart | x |\n| R52 (other) | heart | x |\n" + r55 + h35


RECORDS = Path("experiments") / "harness" / "results" / "h35_heart_trend_forms"
# files the screen read, as the launcher names them (tracked stand-ins in the throwaway repository)
SCREENED_FILES = ("src/a.py", "configs/c.yaml", "experiments/builds/build_test.py", "experiments/builds/halve.py",
                  "experiments/coords/size_procedure.py", "experiments/s2_hedge/payoff.py",
                  "experiments/harness/h25_heart_forms.py", "experiments/harness/h26_heart_scope_forms.py",
                  "experiments/harness/h27_heart_forecast_size.py", "experiments/harness/h28_heart_trend_strength.py",
                  "experiments/harness/h29b_heart_organ.py", "experiments/harness/h34_heart_organ_trend.py",
                  "experiments/harness/h35_heart_trend_forms.py")


def _records(repo, run_commit, words=None, **over):
    """The screen's records as the main session commits them with the outcome: G0's and the read's, both naming the
    run commit (G0's from a clean tree), the read keeping the forms `words` keeps (default: all four)."""
    words = words or {f: "KEPT" for f in FORMS}
    kept = any(w == "KEPT" for w in words.values())
    g0 = {"row": "H35", "axis": "h35-heart-trend-forms:g0", "commit": run_commit, "verdict": "PASS",
          "provenance": {"commit": run_commit, "tree_clean": True, "utc": "x"}, **over.get("g0", {})}
    read = {"row": "H35", "commit": run_commit, "void": False, "outcome": "KEPT" if kept else "STOP",
            "words": words, **over.get("read", {})}
    (repo / RECORDS).mkdir(parents=True, exist_ok=True)
    (repo / RECORDS / "g0.json").write_text(json.dumps(g0))
    (repo / RECORDS / "reading.json").write_text(json.dumps(read))


@pytest.fixture()
def repo(tmp_path):
    """A throwaway repository: the rules file and its pin, stand-ins for the files the screen read and a ledger with
    row R55 as the captain's order for the four files and row H35 open, committed as the screen's run commit; on top,
    one commit with row H35's outcome cell keeping all four and the screen's records naming the run commit."""
    repo = tmp_path / "repo"
    for d in ("docs", "experiments/harness", "src"):
        (repo / d).mkdir(parents=True)
    (tmp_path / "home").mkdir()
    (repo / "docs" / "S2_RULES.md").write_text("rules\n")
    sha = subprocess.run(["sha256sum", "docs/S2_RULES.md"], cwd=repo, capture_output=True, text=True).stdout
    (repo / "docs" / "S2_RULES.sha256").write_text(sha)
    for f in SCREENED_FILES:
        (repo / f).parent.mkdir(parents=True, exist_ok=True)
        (repo / f).write_text("x = 1\n")
    (repo / "experiments" / "harness" / "TRIES.md").write_text(_ledger(h35_row("pending")))
    (repo / ".gitignore").write_text("runs/\n")
    git = _git_fn(repo, tmp_path)
    git("init", "-q")
    git("add", ".")
    git("commit", "-q", "-m", "the screen's run commit")
    run_commit = git("rev-parse", "HEAD")
    (repo / "experiments" / "harness" / "TRIES.md").write_text(_ledger())
    _records(repo, run_commit)
    git("add", ".")
    git("commit", "-q", "-m", "the outcome and the records")
    repo_run[repo] = run_commit
    return repo, git, git("rev-parse", "HEAD")


repo_run: dict = {}  # repository -> its screen's run commit


def _relaunch(repo, git, text, words=None):
    """Commit a new ledger (and, with `words`, the read's record keeping those forms) and return the new HEAD."""
    (repo / "experiments" / "harness" / "TRIES.md").write_text(text)
    if words is not None:
        _records(repo, repo_run[repo], words)
    git("commit", "-q", "-am", "ledger")
    return git("rev-parse", "HEAD")


def _launch(tmp_path, *args, sha_stub=True, extra=None):
    r = subprocess.run(["bash", str(LAUNCHER), *args], cwd=tmp_path / "repo", capture_output=True, text=True,
                       env=_env(tmp_path, _bin(tmp_path, sha_stub), extra), timeout=180)
    return r.returncode, r.stdout, r.stderr


def _scratch(tmp_path, run):
    return tmp_path / "home" / "v" / "scratch" / "hearth" / run


def _exit(tmp_path, run, step):
    p = _scratch(tmp_path, run) / f"{step}.exit"
    return p.read_text().strip() if p.exists() else None


def _steps(forms):
    return ("gate", "tests", *IDENTITY, *(f"{TAGS[f]}{x}" for f in forms for x in "ab"), "done")


def test_the_stand_in_sha256sum_gives_the_committed_hashes(tmp_path):
    b = _bin(tmp_path, sha_stub=True)
    for step, (name, full, half) in IDENTITY.items():
        for f, digest in ((f"x/{name}.h5ad", full), (f"x/{name}_half.h5ad", half)):
            out = subprocess.run([str(b / "sha256sum"), f], capture_output=True, text=True).stdout
            assert out == f"{digest}  {f}\n", step


def test_the_launcher_end_to_end_with_stand_ins(tmp_path, repo):
    repo, git, head = repo
    code, out, err = _launch(tmp_path, head, "r1")
    assert code == 0, err
    assert out == "" and err == ""
    for step in _steps(FORMS):
        assert _exit(tmp_path, "r1", step) == "0", step
    S = _scratch(tmp_path, "r1")
    rel = Path("runs") / f"hearth_{head[:7]}_r1"
    assert (S / "head").read_text().strip() == head and (S / "forms").read_text().split() == list(FORMS)
    assert (S / "screened").read_text().strip() == repo_run[repo] != head
    assert (S / "stamps").read_text().splitlines() == [
        *(f"{step}/{name} {head} True" for step, (name, _, _) in IDENTITY.items()),
        *(f"{TAGS[f]}{x}/T2heart_organ_ts125_{TAGS[f]} {head} True" for f in FORMS for x in "ab")]
    sums = (S / "SHA256SUMS").read_text()
    assert sums == (repo / rel / "SHA256SUMS").read_text()
    names = [ln.split()[1] for ln in sums.splitlines()]
    assert names == [*(f"{step}/{name}{h}.h5ad" for step, (name, _, _) in IDENTITY.items() for h in ("", "_half")),
                     *(f"{TAGS[f]}{x}/T2heart_organ_ts125_{TAGS[f]}{h}.h5ad" for f in FORMS for x in "ab"
                       for h in ("", "_half"))]
    digests = {ln.split()[1]: ln.split()[0] for ln in sums.splitlines()}
    for step, (name, full, half) in IDENTITY.items():
        assert digests[f"{step}/{name}.h5ad"] == full and digests[f"{step}/{name}_half.h5ad"] == half, step
    for t in TAGS.values():
        for h in ("", "_half"):
            assert digests[f"{t}a/T2heart_organ_ts125_{t}{h}.h5ad"] == digests[f"{t}b/T2heart_organ_ts125_{t}{h}.h5ad"]
    calls = _calls(tmp_path)
    py = "run --frozen --no-sync python -I -B experiments/builds/build_test.py heart"
    assert [c for c in calls if "build_test.py" in c] == [
        *(f"{py} {ARGS[step]} --out {rel}/{step}/{name}.h5ad" for step, (name, _, _) in IDENTITY.items()),
        *(f"{py} {ARGS[f]} --out {rel}/{TAGS[f]}{x}/T2heart_organ_ts125_{TAGS[f]}.h5ad" for f in FORMS for x in "ab")]
    halves = [c for c in calls if "experiments/builds/halve.py" in c]
    assert len(halves) == 6 + 8 and all("--seed" not in c for c in halves)  # halve.py's default seed, 1000001
    checks = [c for c in calls if " vec check " in c]
    assert len(checks) == 2 * (6 + 8) and all(c.endswith("--board T2:heart:test_extrap") for c in checks)
    tests = [c for c in calls if " -m pytest " in c]
    assert len(tests) == 1
    for f in ("test_heart_trend_line_h.py", "test_heart_trend_line_h_launcher.py", "test_h35_heart_trend_forms.py",
              "test_run_h35_launcher.py", "test_heart_trend_line_g.py", "test_heart_trend_line_f.py",
              "test_heart_trend_line_d.py", "test_build_heart_organ_scope.py", "test_h28_heart_trend_strength.py",
              "test_halve.py", "test_hidden_counts.py", "test_s2_rules.py"):
        assert f"tests/{f}" in tests[0], f
    assert _launch(tmp_path, head, "r1")[0] == 91  # a run folder is never reused


@pytest.mark.parametrize("words, built", [
    ({"category": "DROP", "contrast": "KEPT", "whole": "IDLE", "portion": "KEPT"}, ("contrast", "portion")),
    ({"category": "KEPT", "contrast": "IDLE", "whole": "IDLE", "portion": "IDLE"}, ("category",)),
    ({"category": "IDLE", "contrast": "DROP", "whole": "KEPT", "portion": "DROP"}, ("whole",)),
])
def test_only_the_forms_row_h35_keeps_are_built(tmp_path, repo, words, built):
    repo, git, _ = repo
    head = _relaunch(repo, git, _ledger(h35_row(f'"{outcome(words)}"')), words)
    code, _, err = _launch(tmp_path, head, "k1")
    assert code == 0, err
    for step in _steps(built):
        assert _exit(tmp_path, "k1", step) == "0", step
    assert (_scratch(tmp_path, "k1") / "forms").read_text().split() == list(built)
    built_args = [c for c in _calls(tmp_path) if "--trend-form" in c]
    assert [re.search(r"--trend-form (\w+)", c).group(1) for c in built_args] == [f for f in built for _ in "ab"]
    for f in set(FORMS) - set(built):
        assert _exit(tmp_path, "k1", f"{TAGS[f]}a") is None


@pytest.mark.parametrize("cell, code", [
    ("OUTCOME H35: VOID: v3: no gate record; nothing is built.", 87),
    (outcome({"category": "IDLE", "contrast": "DROP", "whole": "IDLE", "portion": "IDLE"}), 86),
    ("pending", 92), ("", 92),
    ("OUTCOME H35: category KEPT; contrast KEPT; whole KEPT; portion KEPT", 92),  # no full stop
    ("OUTCOME H35: STOP: category KEPT; contrast IDLE; whole IDLE; portion IDLE.", 92),
    ("OUTCOME H35: category IDLE; contrast IDLE; whole IDLE; portion IDLE.", 92),  # none kept without STOP
    ("OUTCOME H35: category KEPT; contrast KEPT; whole KEPT; portion KEPT. OUTCOME H35: VOID: x.", 92),
    ("OUTCOME H35: category **KEPT**; contrast KEPT; whole KEPT; portion KEPT.", 92),
])
def test_nothing_is_built_unless_row_h35s_outcome_cell_keeps_a_form(tmp_path, repo, cell, code):
    repo, git, _ = repo
    head = _relaunch(repo, git, _ledger(h35_row(cell)))
    got, _, err = _launch(tmp_path, head, "n1")
    assert got == code and "refused: row H35" in err
    assert _exit(tmp_path, "n1", "gate") == str(code) and _exit(tmp_path, "n1", "tests") is None
    assert not any("build_test.py" in c or " -m pytest " in c for c in _calls(tmp_path))
    assert not (repo / "runs").exists()


@pytest.mark.parametrize("ledger", [
    _ledger(""),  # no row H35
    _ledger(ALL_KEPT + ALL_KEPT),  # two
    _ledger(ALL_KEPT.replace(FIRST, FIRST.replace("ten tries", "eight tries"))),
    _ledger(ALL_KEPT.replace("| heart | x |", "| heart |", 1)),  # ten cells
    _ledger(ALL_KEPT.replace("| H35 (", "| H350 (")),
])
def test_a_missing_or_malformed_row_h35_builds_nothing(tmp_path, repo, ledger):
    repo, git, _ = repo
    head = _relaunch(repo, git, ledger)
    got, _, err = _launch(tmp_path, head, "m1")
    assert got == 92 and "missing or malformed" in err
    assert not any("build_test.py" in c for c in _calls(tmp_path))


def test_row_r55_must_name_each_kept_forms_file(tmp_path, repo):
    repo, git, _ = repo
    words = {"category": "IDLE", "contrast": "KEPT", "whole": "IDLE", "portion": "KEPT"}
    no_portion = R55_ROW.replace(" or `--trend-form portion`", "")
    head = _relaunch(repo, git, _ledger(h35_row(f'"{outcome(words)}"'), no_portion))
    got, _, err = _launch(tmp_path, head, "q1")
    assert got == 95 and "naming portion" in err and _exit(tmp_path, "q1", "gate") == "95"
    # the same row orders a set that does not need portion
    words = {"category": "KEPT", "contrast": "KEPT", "whole": "IDLE", "portion": "DROP"}
    head = _relaunch(repo, git, _ledger(h35_row(f'"{outcome(words)}"'), no_portion), words)
    assert _launch(tmp_path, head, "q2")[0] == 0
    for i, bad in enumerate(("", R55_ROW + R55_ROW, R55_ROW.replace("| R55 (", "| R56 ("),
                             R55_ROW.replace("`--captain-ruling R55`", "`--captain-ruling R49`"),
                             R55_ROW.replace("ruling, no run", "awaiting a ruling"),
                             R55_ROW.replace("R55 (the captain's order", "R55 (the captain's named exception"))):
        head = _relaunch(repo, git, _ledger(ALL_KEPT, bad))
        assert _launch(tmp_path, head, f"bad{i}")[0] == 95, i


@pytest.mark.parametrize("change", ["no_g0", "no_read", "g0_not_json", "read_not_an_object", "other_commit_in_read",
                                    "unclean_g0", "provenance_other_commit", "void_read", "stop_read", "other_kept",
                                    "short_commit", "not_in_history"])
def test_the_screens_records_must_name_one_clean_run_commit_that_keeps_these_forms(tmp_path, repo, change):
    repo, git, _ = repo
    run = repo_run[repo]
    rec = repo / RECORDS
    if change == "no_g0":
        (rec / "g0.json").unlink()
    elif change == "no_read":
        (rec / "reading.json").unlink()
    elif change == "g0_not_json":
        (rec / "g0.json").write_text("{not json")
    elif change == "read_not_an_object":
        (rec / "reading.json").write_text("[]")
    elif change == "other_commit_in_read":
        _records(repo, run, read={"commit": "d" * 40})
    elif change == "unclean_g0":
        _records(repo, run, g0={"provenance": {"commit": run, "tree_clean": False}})
    elif change == "provenance_other_commit":
        _records(repo, run, g0={"provenance": {"commit": "d" * 40, "tree_clean": True}})
    elif change == "void_read":
        _records(repo, run, read={"void": True})
    elif change == "stop_read":
        _records(repo, run, {f: "IDLE" for f in FORMS})
    elif change == "other_kept":  # the read keeps one form, row H35's outcome cell all four
        _records(repo, run, {"category": "KEPT", "contrast": "IDLE", "whole": "IDLE", "portion": "DROP"})
    elif change == "short_commit":
        _records(repo, run[:12], g0={"provenance": {"commit": run[:12], "tree_clean": True}})
    else:  # a well-formed commit that is not in the history
        _records(repo, "e" * 40)
    git("add", "-A")
    git("commit", "-q", "-m", f"records: {change}")
    got, _, err = _launch(tmp_path, git("rev-parse", "HEAD"), "s1")
    assert got == 84 and "refused: the screen's" in err, (change, err)
    assert _exit(tmp_path, "s1", "gate") == "84" and not (_scratch(tmp_path, "s1") / "screened").exists()
    assert not any("build_test.py" in c or " -m pytest " in c for c in _calls(tmp_path))


@pytest.mark.parametrize("path", [*SCREENED_FILES, "src/vec/new_module.py", "configs/new.yaml"])
def test_a_file_the_screen_read_that_changed_since_its_run_commit_refuses_the_builds(tmp_path, repo, path):
    repo, git, _ = repo
    (repo / path).parent.mkdir(parents=True, exist_ok=True)
    (repo / path).write_text("x = 2\n")
    git("add", "-A")
    git("commit", "-q", "-m", "a change after the screen")
    got, _, err = _launch(tmp_path, git("rev-parse", "HEAD"), "c1")
    assert got == 85 and f"differs between its run commit {repo_run[repo]}" in err, (path, err)
    assert not any("build_test.py" in c or " -m pytest " in c for c in _calls(tmp_path))


def test_other_changes_since_the_screens_run_commit_build(tmp_path, repo):
    """The ledger, the screen's records, other lines' files, tests and documents may change after the screen's run
    commit: the builds go ahead."""
    repo, git, _ = repo
    for path in ("docs/notes.md", "tests/test_x.py", "experiments/harness/h36_heart_run_profile.py",
                 "experiments/builds/heart_xpc_build.sh", str(RECORDS / "gate.json"), str(RECORDS / "reading.txt")):
        (repo / path).parent.mkdir(parents=True, exist_ok=True)
        (repo / path).write_text("another line's file\n")
    ledger = repo / "experiments" / "harness" / "TRIES.md"
    ledger.write_text(ledger.read_text() + "| R56 (the captain's order of another line) | heart | x |\n")
    git("add", "-A")
    git("commit", "-q", "-m", "other work after the screen")
    code, _, err = _launch(tmp_path, git("rev-parse", "HEAD"), "o1")
    assert code == 0, err


@pytest.mark.parametrize("step", list(IDENTITY))
@pytest.mark.parametrize("half", [False, True])
def test_the_launcher_stops_unless_each_identity_step_gives_its_committed_bytes(tmp_path, repo, step, half):
    _, _, head = repo
    var = f"STUB_SHA_{step}" + ("_HALF" if half else "")
    code = STOP_CODES[step][int(half)]
    assert _launch(tmp_path, head, "s", extra={var: "0" * 64})[0] == code
    assert _exit(tmp_path, "s", step) == str(code)
    later = _steps(FORMS)[_steps(FORMS).index(step) + 1:]
    assert all(_exit(tmp_path, "s", s) is None for s in later)
    assert not any("--trend-form" in c for c in _calls(tmp_path))  # no form is built


def test_the_launcher_stops_at_file_d_on_its_real_hash(tmp_path, repo):
    _, _, head = repo
    assert _launch(tmp_path, head, "r2", sha_stub=False)[0] == 20
    assert _exit(tmp_path, "r2", "D") == "20" and _exit(tmp_path, "r2", "OTS05") is None
    assert not any("--trend-strength" in c for c in _calls(tmp_path))


def test_the_launcher_fails_on_a_build_a_check_a_stamp_the_tests_or_two_different_files(tmp_path, repo):
    _, _, head = repo
    assert _launch(tmp_path, head, "b1", extra={"STUB_BUILD_EXIT": "1"})[0] == 11
    assert _exit(tmp_path, "b1", "D") == "11"
    assert _launch(tmp_path, head, "b2", extra={"STUB_CHECK_EXIT": "1"})[0] == 13
    assert _launch(tmp_path, head, "b3", extra={"STUB_DIRTY": "1"})[0] == 15
    assert _launch(tmp_path, head, "b4", extra={"STUB_PYTEST_EXIT": "1"})[0] == 10
    assert _exit(tmp_path, "b4", "tests") == "10" and _exit(tmp_path, "b4", "D") is None
    # one form's two builds differ: done fails, every step before it passed
    assert _launch(tmp_path, head, "b5", extra={"STUB_NONDET": "1", "STUB_NONDET_ON": "--trend-form whole "})[0] == 31
    assert _exit(tmp_path, "b5", "done") == "31"
    assert all(_exit(tmp_path, "b5", s) == "0" for s in _steps(FORMS)[:-1])


def test_the_launcher_gate_refusals(tmp_path, repo):
    repo, git, head = repo
    assert _launch(tmp_path, head[:12], "g0")[0] == 90
    for odd in ("a/b", "it's", "$(id)"):
        assert _launch(tmp_path, head, odd)[0] == 90
    for i, extra in enumerate(("R55", "category", "")):  # the row is fixed and the forms come from row H35
        assert _launch(tmp_path, head, f"x{i}", extra)[0] == 90, extra
        assert not _scratch(tmp_path, f"x{i}").exists(), extra
    for args in ((), (head,), (head, "")):
        code = _launch(tmp_path, *args)[0]
        assert code != 0 and not (tmp_path / "home" / "v" / "scratch" / "hearth").exists(), args
    assert _launch(tmp_path, "d" * 40, "g1")[0] == 96 and _exit(tmp_path, "g1", "gate") == "96"
    (repo / "src" / "a.py").write_text("x = 2\n")
    assert _launch(tmp_path, head, "g2")[0] == 98
    git("checkout", "-q", "--", "src/a.py")
    (repo / "runs").mkdir()
    (repo / "runs" / f"hearth_{head[:7]}_g3").mkdir()
    assert _launch(tmp_path, head, "g3")[0] == 93
    (repo / "docs" / "S2_RULES.md").write_text("rules, changed\n")
    git("commit", "-q", "-am", "rules changed without a new pin")
    assert _launch(tmp_path, git("rev-parse", "HEAD"), "pin")[0] == 97
    assert not any("build_test.py" in c or " -m pytest " in c for c in _calls(tmp_path))


# ================================================================================================ the inline readers
def _inline(sh: str, start: str, end: str) -> str:
    i = sh.index(start) + len("python3 -I -B -c '")
    return sh[i:sh.index(end, i)]


def _r55_code() -> str:
    start = ("python3 -I -B -c 'import re\nimport sys\nrows = [ln for ln in open(sys.argv[1], encoding=\"utf-8\")"
             ".read().splitlines() if ln.startswith(\"| R55 \")]")
    return _inline(LAUNCHER.read_text(), start, "' experiments/harness/TRIES.md \"$form\"")


def test_the_gates_r55_check_is_the_builders(tmp_path, monkeypatch):
    """The gate's inline R55 check and the builder's _heart_tf_ordered take and refuse the same ledgers, form by
    form."""
    sys.path.insert(0, str(REPO_ROOT / "experiments" / "builds"))
    import build_test as bt

    code = _r55_code()
    ledgers = {"ok": _ledger(), "only": R55_ROW, "none": "| R27 (x) | T1 | x |\n", "two": _ledger() + R55_ROW,
               "other": _ledger().replace("| R55 (", "| R56 ("), "longer": _ledger().replace("| R55 (", "| R550 ("),
               "ruled": _ledger().replace("`--captain-ruling R55`", "`--captain-ruling R49`"),
               "wider": _ledger().replace("`--trend-strength 1.25`", "`--trend-strength 1.250`"),
               "inside": _ledger().replace("`--trend-form whole`", "`--trend-form whole --z-mode raw`"),
               "no_portion": _ledger().replace(" or `--trend-form portion`", ""),
               "kind": _ledger().replace("R55 (the captain's order", "R55 (the captain's named exception"),
               "unpiped": _ledger().replace("| R55 (", "|R55 ("),
               "outside": _ledger().replace(", `--trend-form contrast`", "").replace(
                   "| heart | x |\n| " + FIRST, "| `--trend-form contrast` | x |\n| " + FIRST),
               **{f"drop{i}": _ledger().replace(s, "x") for i, s in enumerate((
                   "`--entry organ_scope`", "`--trend-strength 1.25`", "`--captain-ruling R55`", "ruling, no run",
                   "the captain's order"))}}
    seen = set()
    for name, text in ledgers.items():
        path = tmp_path / f"{name}.md"
        path.write_text(text)
        monkeypatch.setattr(bt, "TRIES_FILE", path)
        for form in FORMS:
            gate = subprocess.run(["python3", "-I", "-B", "-c", code, str(path), form]).returncode == 0
            assert gate == bt._heart_tf_ordered(form), (name, form)
            seen.add(gate)
    assert seen == {True, False}


def test_the_gates_r55_check_takes_the_committed_ledger():
    """The inline R55 check passes on experiments/harness/TRIES.md as committed (read by code only) for each form."""
    for form in FORMS:
        r = subprocess.run(["python3", "-I", "-B", "-c", _r55_code(), str(REPO_ROOT / "experiments" / "harness" /
                                                                            "TRIES.md"), form])
        assert r.returncode == 0, form


def test_the_gates_h35_reader_is_the_screens(tmp_path):
    """The gate's inline reader of row H35 and the screen's own (h35_heart_trend_forms.screen_outcome) allow the same
    forms on the same ledgers (tests/test_h35_heart_trend_forms.py holds it on every word set as well)."""
    spec = importlib.util.spec_from_file_location("h35_for_line_h_launcher", REPO_ROOT / "experiments" / "harness" /
                                                  "h35_heart_trend_forms.py")
    h35 = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = h35
    spec.loader.exec_module(h35)
    assert h35.FIRST_CELL == FIRST and h35.FORMS == FORMS
    code = _inline(LAUNCHER.read_text(), "python3 -I -B -c '\nimport re", "' experiments/harness/TRIES.md)")
    cells = [outcome({"category": a, "contrast": b, "whole": "IDLE", "portion": "DROP"})
             for a in ("KEPT", "IDLE") for b in ("KEPT", "DROP")]
    for i, text in enumerate([_ledger(h35_row(c)) for c in cells] + [_ledger(""), _ledger(h35_row("pending"))]):
        path = tmp_path / f"h{i}.md"
        path.write_text(text)
        r = subprocess.run(["python3", "-I", "-B", "-c", code, str(path)], capture_output=True, text=True)
        want = h35.screen_outcome(text)
        assert (r.returncode == 0) == bool(want and want["allowed"]), i
        assert r.stdout.split() == (want["allowed"] if want else []), i


# ================================================================================================ as text
def _body(sh: str) -> str:
    return "\n".join(ln for ln in sh.splitlines() if not ln.lstrip().startswith("#"))


def test_the_launcher_deletes_nothing_parses_and_keeps_the_g_launchers_conventions():
    sh, g = LAUNCHER.read_text(), G_LAUNCHER.read_text()
    body = _body(sh)
    assert not re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|truncate|shred|mv)\s", body)
    assert not re.search(r"\bscp\b|\brsync\b|\bcurl\b|\bwget\b|\bssh\b|vec submit|\bupload", body)
    for msg in re.findall(r"\$\{[A-Za-z0-9_]+:\?([^}]*)\}", sh):
        assert "'" not in msg and '"' not in msg, msg
    assert len(re.findall(r"\$\{[A-Za-z0-9_]+:\?([^}]*)\}", sh)) == 2 and "${3" not in sh and "$ROW" not in sh
    assert "python -I -B" in sh and 'UV_PROJECT_ENVIRONMENT="$PWD/.venv"' in sh and "NUMPY_MADVISE_HUGEPAGE=0" in sh
    assert 'TMPDIR="$S/tmp"' in sh and "OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8" in sh
    assert "OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16" in sh and "set -o pipefail" in sh
    assert "/mnt/" not in sh and "vec-campaign/" not in body
    assert "S=$HOME/v/scratch/hearth/$RUN" in sh and "O=runs/hearth_${head:0:7}_$RUN" in sh
    assert body.index("export PATH=$HOME/.local/bin:$PATH") < body.index("uv run")
    assert body.index('[ "$#" = 2 ]') < body.index('mkdir -p "$HOME')
    assert subprocess.run(["bash", "-n", str(LAUNCHER)], capture_output=True).returncode == 0
    # the shared parts are heart_trend_line_g.sh's, line for line
    for start, end in (("finish() {", "step_ok() {"), ("stamp() {", "identity() {")):
        mine = sh[sh.index(start):sh.index(end, sh.index(start))]
        assert mine in g, start
    for line in ('[ "$head" = "$COMMIT" ] || { echo "refused: HEAD is not $COMMIT" >&2; finish 96; }',
                 "sha256sum -c docs/S2_RULES.sha256 > /dev/null 2>&1 || finish 97",
                 '$PY experiments/builds/build_test.py heart "$@" --out "$d/$name.h5ad" > /dev/null 2>&1 || return 11',
                 '  [ "$st" = "$head True" ] || return 15'):
        assert line in sh and line in g, line


def test_the_screened_files_are_the_screens_and_the_builds():
    """SCREEN_FILES names src, configs, the builder (whose context G0 takes) and halve.py, and the screen's scripts as
    they load each other, each one present; the test fixture's stand-ins are those paths."""
    sh = LAUNCHER.read_text()
    listed = [p for m in re.finditer(r'^SCREEN_FILES="(?:\$SCREEN_FILES )?([^"]*)"$', sh, re.M)
              for p in m.group(1).split()]
    scripts = [f"experiments/harness/{s}.py" for s in ("h25_heart_forms", "h26_heart_scope_forms",
                                                        "h27_heart_forecast_size", "h28_heart_trend_strength",
                                                        "h29b_heart_organ", "h34_heart_organ_trend",
                                                        "h35_heart_trend_forms")]
    assert listed == ["src", "configs", "experiments/builds/build_test.py", "experiments/builds/halve.py",
                      "experiments/coords/size_procedure.py", "experiments/s2_hedge/payoff.py", *scripts]
    assert all((REPO_ROOT / p).exists() for p in listed)
    assert [p for p in listed if p.startswith("experiments/")] == list(SCREENED_FILES[2:])
    loads = {"h35_heart_trend_forms": ("h34_heart_organ_trend.py", "experiments\" / \"builds\" / \"build_test.py"),
             "h34_heart_organ_trend": ("h28_heart_trend_strength.py", "h29b_heart_organ.py"),
             "h29b_heart_organ": ("h27_heart_forecast_size.py",),
             "h28_heart_trend_strength": ("h27_heart_forecast_size.py",),
             "h27_heart_forecast_size": ("h26_heart_scope_forms.py",),
             "h26_heart_scope_forms": ("import h25_heart_forms", "\"s2_hedge\" / \"payoff.py\""),
             "h25_heart_forms": ("import halve", "import size_procedure")}
    for script, names in loads.items():
        text = (REPO_ROOT / "experiments" / "harness" / f"{script}.py").read_text()
        for name in names:
            assert name in text, (script, name)


def test_the_launcher_reads_only_exit_sha_and_stamp():
    sh = LAUNCHER.read_text()
    for line in (ln for ln in sh.splitlines() if not ln.lstrip().startswith("#")):
        if line.lstrip().startswith(("$PY experiments/builds/", "uv run --frozen --no-sync vec check")):
            assert "> /dev/null 2>&1" in line, line
    body = _body(sh)
    assert 'print(s["commit"], s["clean"])' in body and '"size"' not in body
    assert not re.search(r"(^|[\s;&|(])(stat|du|wc|ls|cat|head|tail|less|h5ls)\s", body)
    for step in IDENTITY:
        assert f'{step}_ARGS="{ARGS[step]}"' in sh, step
    for f, t in TAGS.items():
        assert f'{t.upper()}_ARGS="{ARGS[f]}"' in sh, f
    assert "BOARD=T2:heart:test_extrap" in sh


def test_the_identity_hashes_are_the_committed_list():
    """Every identity step's hashes are the ones experiments/builds/results/heartg_b2bbcf5_g1.SHA256SUMS records (the
    last build of all six, OTS15's two builds one sha256); the launcher holds exactly these."""
    sh = LAUNCHER.read_text()
    listed: dict[str, set] = {}
    for ln in G_LIST.read_text().splitlines():
        digest, name = ln.split()
        listed.setdefault(name.split("/", 1)[1], set()).add(digest)
    for step, (name, full, half) in IDENTITY.items():
        assert listed[f"{name}.h5ad"] == {full} and listed[f"{name}_half.h5ad"] == {half}, step
        assert re.search(rf"^{step}_SHA={full}\b", sh, re.M), step
        assert re.search(rf"^{step}_HALF_SHA={half}\b", sh, re.M), step


def test_the_launchers_commands_are_the_builders():
    sys.path.insert(0, str(REPO_ROOT / "experiments" / "builds"))
    import build_test as bt

    sh = LAUNCHER.read_text()
    assert f'D_ARGS="{bt.HEART_ORGAN_SCOPE[len("heart "):]}"' in sh
    for s, var in ((0.5, "OTS05_ARGS"), (1.0, "OTS10_ARGS")):
        assert f'{var}="{bt.HEART_OTS_COMMANDS[s][len("heart "):]}"' in sh
    for s, var in ((0.75, "OTS075_ARGS"), (1.25, "OTS125_ARGS")):
        assert f'{var}="{bt.HEART_FURTHER_COMMANDS[s].format(row=bt.HEART_FURTHER_ROW)[len("heart "):]}"' in sh
    assert f'OTS15_ARGS="{bt.HEART_OTS15_COMMAND[len("heart "):]}"' in sh
    for f in FORMS:
        assert f'{bt.HEART_TF_TAGS[f].upper()}_ARGS="{bt.HEART_TF_COMMANDS[f][len("heart "):]}"' in sh
        assert "T2heart_organ_ts125_$tag" in sh and bt.HEART_TF_TAGS[f] == TAGS[f]
