"""The test-board build launcher of T1 lines 7 and 6 of row R22 and row R37's pooled base (experiments/harness/TRIES.md;
the captain's orders of 2026-11-24, ruling rows R39 and R40, and his named exception of 2026-11-24, row R37).

experiments/builds/t1_l67pb_line.sh first rebuilds fifteen known files (the eleven of t1_b2d_line.sh, T1_s2rcB30,
T1_s2rcCm25, T1_s2rcB2d15 and T1_s2rcB2d20) and stops unless each gives its known bytes; then it builds T1_s2rcCm30,
T1_s2rcB25d15, T1_s2rcB25d20 and T1_s2rc_pooled twice each from one clean commit and compares each pair. These tests
hold its command strings, hashes and step order to that, check the fifteen hashes against the committed build records,
check that it deletes nothing and reads no output of a test build, and run its gate's paths for real in throwaway
repositories under a temporary home, with a stand-in `uv` that fails, so that a call passing the gate stops at the tests
step and builds nothing. Nothing here reads the repository's own ledger or rules file for its state of today."""

import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
LAUNCHER = REPO / "experiments" / "builds" / "t1_l67pb_line.sh"
RESULTS = REPO / "experiments" / "builds" / "results"
BASE = "--g-mom --base lib3 --feng-groups recomputed"
# step -> (file name, args variable, args after BASE, hash variable, sha256, exit code)
KNOWN = {
    "D": ("T1_s2rc", "BASE", "", "R13_SHA", "ad8f24c71c90ca17bce5f8876d503471ab3faa5b17a0878addb00aff689fca4f", 20),
    "A": ("T1_s2rcA", "A", " --g-mom-damping 1.0", "A_SHA",
          "76865826f5a3741c7585f81c4241048ca1a95b5791ebb3c8feb1f2b536678ce4", 21),
    "B": ("T1_s2rcB", "B", " --feng-fraction-mult 1.5", "B_SHA",
          "c7240cb6a7f7ee93d037f437832efbaad004f1fbe608d80072f8efd5d20a0062", 22),
    "C": ("T1_s2rcC", "C", " --g-mom-damping 1.0 --feng-fraction-mult 1.5", "C_SHA",
          "31477ef7df5d5374a74a8249c570217020c593916ab4d3780049b7b3f87b0620", 23),
    "X2": ("T1_s2rcB2", "X2", " --feng-fraction-mult 2.0 --captain-ruling R20", "X2_SHA",
           "52897fe375b8a81afffa788592047c30361824148b00a8176514385060518837", 24),
    "X25": ("T1_s2rcB25", "X25", " --feng-fraction-mult 2.5 --captain-ruling R22", "X25_SHA",
            "ccb7c457c0e8ac8ac668bcea17bc61ed22f414f58ab5a63e946cb2f08c4f4e77", 25),
    "A15": ("T1_s2rcA15", "A15", " --g-mom-damping 1.5 --captain-ruling R22", "A15_SHA",
            "e2334137b54d98b04e9dece90213987e61fe9aa15a3fb8d85e6c262d8fb5ba75", 26),
    "A20": ("T1_s2rcA20", "A20", " --g-mom-damping 2.0 --captain-ruling R22", "A20_SHA",
            "372713686782be273a52427140eab652517e35b679f54fc2a969484dee730218", 27),
    "Cm20": ("T1_s2rcCm20", "CM20", " --g-mom-damping 1.0 --feng-fraction-mult 2.0 --captain-ruling R22", "CM20_SHA",
             "acce5f5f7097269bffebbfd583833e4c9c1c5ae77f1d5294dde6d13e197f863e", 28),
    "MS": ("T1_s2rc_MS", "MS", " --direction-momentum all --captain-ruling R23", "MS_SHA",
           "6103d10535665634ee0cad3390590a5df132d1805682fb1d309bb19581fdaa5e", 29),
    "MSB": ("T1_s2rcMSB", "MSB", " --feng-fraction-mult 1.5 --direction-momentum all --captain-ruling R26", "MSB_SHA",
            "78bb17eea5721e590cfd2c9f707e9cd6f2f496e081acede4f697dc4035871da2", 30),
    "X30": ("T1_s2rcB30", "X30", " --feng-fraction-mult 3.0 --captain-ruling R22", "X30_SHA",
            "b639de29bd13d7005450ca53e7cd3dacb2e362a55beca8545937858e6582a024", 31),
    "Cm25": ("T1_s2rcCm25", "CM25", " --g-mom-damping 1.0 --feng-fraction-mult 2.5 --captain-ruling R22", "CM25_SHA",
             "c5fb96a64075b6fbfb2f3066604d622ea8b89ac3a5e4f8a99d0fea7d93011dfd", 32),
    "B2d15": ("T1_s2rcB2d15", "B2D15", " --g-mom-damping 1.5 --feng-fraction-mult 2.0 --captain-ruling R22",
              "B2D15_SHA", "fe4e7403d6257cc265a02d0a13d40cdb60a97e9fdc0647d6b5c7e318aedaed6c", 33),
    "B2d20": ("T1_s2rcB2d20", "B2D20", " --g-mom-damping 2.0 --feng-fraction-mult 2.0 --captain-ruling R22",
              "B2D20_SHA", "bd12c66312ba32ab938b4ee71bd5c79e58a154d6fc56434012e6dad907e668dd", 34),
}
# new step pair -> (file name, args variable, exit code of the pair's comparison)
PAIRS = {("Cm30a", "Cm30b"): ("T1_s2rcCm30", "CM30", 41), ("B25d15a", "B25d15b"): ("T1_s2rcB25d15", "B25D15", 42),
         ("B25d20a", "B25d20b"): ("T1_s2rcB25d20", "B25D20", 43),
         ("pooleda", "pooledb"): ("T1_s2rc_pooled", "POOLED", 44)}
STEPS = ["gate", "tests", "depth", *KNOWN, *(s for pair in PAIRS for s in pair), "done"]
OPT_CM30 = "`--g-mom-damping 1.0 --feng-fraction-mult 3.0`"
OPT_B25D15 = "`--g-mom-damping 1.5 --feng-fraction-mult 2.5`"
OPT_B25D20 = "`--g-mom-damping 2.0 --feng-fraction-mult 2.5`"
OPT_POOLED = "`--base pooled` with `--feng-groups recomputed`"
NEW = {"CM30": f"{BASE} --g-mom-damping 1.0 --feng-fraction-mult 3.0 --captain-ruling R22",
       "B25D15": f"{BASE} --g-mom-damping 1.5 --feng-fraction-mult 2.5 --captain-ruling R22",
       "B25D20": f"{BASE} --g-mom-damping 2.0 --feng-fraction-mult 2.5 --captain-ruling R22",
       "POOLED": "--g-mom --base pooled --feng-groups recomputed --captain-ruling R37"}


def _sh() -> str:
    return LAUNCHER.read_text()


def _args(sh: str) -> dict[str, str]:
    return dict(re.findall(r'^([A-Z0-9]+)_ARGS="([^"]+)"$', sh, flags=re.M))


# ------------------------------------------------------------------ what the launcher runs


def test_bash_parses_the_launcher():
    out = subprocess.run(["bash", "-n", str(LAUNCHER)], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr


def test_the_launcher_reproduces_fifteen_files_then_builds_each_new_file_twice():
    sh = _sh()
    want = {key: f"{BASE}{tail}" for _, key, tail, _, _, _ in KNOWN.values()} | NEW
    assert _args(sh) == want
    for step, (name, key, _, var, sha, code) in KNOWN.items():
        assert f"build {step} {name} ${key}_ARGS || finish $?" in sh
        assert f'same {step} {name} "${var}" {code}' in sh
        assert f"{var}={sha}" in sh and sh.count(sha) == 1
    assert 'got=$(sha256sum "$O/$1/$2.h5ad" | cut -d\' \' -f1)' in sh
    assert '[ "$got" = "$3" ] || {' in sh and 'finish "$4"; }' in sh
    for (a, b), (name, key, code) in PAIRS.items():
        assert f'for STEP in {a} {b}; do begin; build "$STEP" {name} ${key}_ARGS || finish $?; step_ok; done' in sh
        assert f'cmp -s "$O/{a}/{name}.h5ad" "$O/{b}/{name}.h5ad" || finish {code}' in sh
    order = [sh.index(s) for s in ("STEP=gate", "STEP=tests", "STEP=depth", *(f"STEP={k}\n" for k in KNOWN),
                                   *(f"for STEP in {a} {b}" for a, b in PAIRS), "STEP=done")]
    assert order == sorted(order)
    listing = sh[sh.index("(cd \"$O\" && sha256sum"):sh.index('> "$S/SHA256SUMS"')]
    assert re.findall(r"([A-Za-z0-9]+)/\*\.h5ad", listing) == STEPS[3:-1]  # every built folder, once, in order
    assert 'cp "$S/SHA256SUMS" "$O/SHA256SUMS"' in sh
    assert "vec check" in sh and "BOARD=T1:test" in sh
    assert "--seed" not in sh and "--board T1" not in sh  # seed and board: the builder's defaults, as the references
    assert "OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16" in sh  # the builds, as the references
    assert "OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8" in sh  # the tests and the depth record
    assert sh.index("NUM_THREADS=8") < sh.index("STEP=tests") < sh.index("STEP=depth") < sh.index("NUM_THREADS=16")
    assert sh.index("NUM_THREADS=16") < sh.index("STEP=D\n")
    assert "NUMPY_MADVISE_HUGEPAGE=0" in sh and 'TMPDIR="$S/tmp"' in sh
    assert sh.index("PATH=$HOME/.local/bin:$PATH") < sh.index("uv run")  # PATH set before the first uv call
    assert 'PY="uv run --frozen --no-sync python -I -B"' in sh
    assert '[ "$head" = "$COMMIT" ]' in sh  # the full hash, never a prefix
    assert '[ "$st" = "$head True" ] || return 15' in sh  # a build counts only with a stamp of this commit, clean
    assert 'echo "$(basename "$d")/$name $st" >> "$S/stamps"' in sh
    assert "upload" not in sh.replace("Nothing here uploads anything", "")
    assert "/v/vec-campaign" not in sh.replace("(~/v/vec-campaign-d, never the reference checkout)", "")
    assert "cd " not in sh.replace('(cd "$O" && sha256sum', "")  # runs from the folder it is called in
    assert "S2_RULES.md" not in sh.replace("docs/S2_RULES.sha256", "")  # the rules file is checked by its hash only


def test_the_launcher_reads_nothing_a_test_builds_log_size_or_shape_would_reveal():
    sh = _sh()
    for line in sh.splitlines():  # every builder and check call: output to /dev/null, read by its exit code only
        if re.search(r'build_test\.py t1|vec check "', line):
            assert "> /dev/null 2>&1" in line, line
    code = "\n".join(line.split("#")[0] for line in sh.splitlines())
    for word in (".build.json", "du ", "stat ", "wc ", "ls ", "h5ls", "read_h5ad", "n_obs", "n_cells", "head -",
                 "tail ", "cat "):
        assert word not in code, word
    assert "tests.out" in sh  # the tests' log: synthetic stages only
    assert 'print(s["commit"], s["clean"])' in sh and "stamp.json" in sh  # of the stamp, its commit and clean only
    assert "curl" not in sh and "wget" not in sh and "git fetch" not in sh and "git pull" not in sh  # fetches nothing


def test_every_step_writes_running_first_then_its_code():
    sh = _sh()
    assert 'begin() { echo running > "$S/$STEP.exit"; }' in sh
    assert 'step_ok() { echo 0 > "$S/$STEP.exit"; }\nbegin\n\nhead=$(git rev-parse HEAD)\n' in sh
    for step in ("tests", "depth", *KNOWN, "done"):
        at = sh.index(f"STEP={step}\n")
        assert sh[at:].split("\n", 2)[1] == "begin", step  # the line after STEP=<step>
    for a, b in PAIRS:
        assert f"for STEP in {a} {b}; do begin; build" in sh
    assert sh.index("exit 91") < sh.index("STEP=gate") and sh.index("exit 90") < sh.index("STEP=gate")
    gate = sh[sh.index("STEP=gate"):sh.index("STEP=tests")]
    assert gate.index('echo "$head" > "$S/head"') > gate.index("finish 93")  # the head file after every check


def test_the_depth_step_runs_the_record_script_alone_and_never_stops_the_run():
    sh = _sh()
    step = sh[sh.index("STEP=depth\n"):sh.index("export OMP_NUM_THREADS=16")]
    assert step == ('STEP=depth\nbegin\n$PY experiments/builds/t1_pool_depth_record.py > "$S/pool_depth.json" '
                    '2> "$S/pool_depth.err"\ncode=$?\n[ "$code" = 0 ] && cp "$S/pool_depth.json" '
                    '"$O/pool_depth.json"\necho "$code" > "$S/$STEP.exit"\n\n')
    assert "finish" not in step and "exit " not in step  # no gate: the run goes on whatever its code


def test_the_documented_exit_codes_are_the_ones_used():
    sh = _sh()
    doc = sh[sh.index("# Exit codes:"):sh.index("set -u")]
    used = {int(c) for c in re.findall(r"(?:finish|exit|return) (\d+)", sh)}
    used |= {code for *_, code in KNOWN.values()}
    documented = {int(c) for c in re.findall(r"(?<![\d.])(\d\d)(?= )", doc)}
    assert used == documented, (used, documented)


def _sums(name: str) -> dict[str, str]:
    """A committed SHA256SUMS list of experiments/builds/results: path -> sha256."""
    words = (RESULTS / name).read_text().split()
    return dict(zip(words[1::2], words[0::2], strict=True))


def test_the_fifteen_known_hashes_are_the_recorded_ones():
    """The fifteen identity hashes are the ones the earlier runs' committed lists hold (no build here)."""
    sha = {step: v[4] for step, v in KNOWN.items()}
    cand = _sums("t1cand_1420c70.SHA256SUMS")
    assert (cand["T1_s2rcA.h5ad"], cand["T1_s2rcB.h5ad"], cand["T1_s2rcC.h5ad"]) == (sha["A"], sha["B"], sha["C"])
    x2 = _sums("t1x2_3d18841.SHA256SUMS")
    assert x2["D/T1_s2rc.h5ad"] == sha["D"] and x2["X2a/T1_s2rcB2.h5ad"] == x2["X2b/T1_s2rcB2.h5ad"] == sha["X2"]
    line1 = _sums("t1line_e3cc0e4_r1.SHA256SUMS")
    assert line1["X25a/T1_s2rcB25.h5ad"] == line1["X25b/T1_s2rcB25.h5ad"] == sha["X25"]
    assert line1["X30a/T1_s2rcB30.h5ad"] == line1["X30b/T1_s2rcB30.h5ad"] == sha["X30"]
    damp = _sums("t1damp_9a10157_r1.SHA256SUMS")
    assert damp["X15a/T1_s2rcA15.h5ad"] == damp["X15b/T1_s2rcA15.h5ad"] == sha["A15"]
    assert damp["X20a/T1_s2rcA20.h5ad"] == damp["X20b/T1_s2rcA20.h5ad"] == sha["A20"]
    cline = _sums("t1cline_2a1c248_c1.SHA256SUMS")
    assert cline["Cm20a/T1_s2rcCm20.h5ad"] == cline["Cm20b/T1_s2rcCm20.h5ad"] == sha["Cm20"]
    assert cline["Cm25a/T1_s2rcCm25.h5ad"] == cline["Cm25b/T1_s2rcCm25.h5ad"] == sha["Cm25"]
    ms = _sums("t1ms_6bc36d9_r1.SHA256SUMS")
    assert ms["MSa/T1_s2rc_MS.h5ad"] == ms["MSb/T1_s2rc_MS.h5ad"] == sha["MS"]
    msc = _sums("t1msc_a3e879a_m1.SHA256SUMS")
    assert msc["MSBa/T1_s2rcMSB.h5ad"] == msc["MSBb/T1_s2rcMSB.h5ad"] == sha["MSB"]
    b2d = _sums("t1b2dline_7e7d9a6_b1.SHA256SUMS")
    assert b2d["B2d15a/T1_s2rcB2d15.h5ad"] == b2d["B2d15b/T1_s2rcB2d15.h5ad"] == sha["B2d15"]
    assert b2d["B2d20a/T1_s2rcB2d20.h5ad"] == b2d["B2d20b/T1_s2rcB2d20.h5ad"] == sha["B2d20"]
    # the eleven of t1_b2d_line.sh's run, as that run rebuilt them
    for step in ("D", "A", "B", "C", "X2", "X25", "A15", "A20", "Cm20", "MS", "MSB"):
        assert b2d[f"{step}/{KNOWN[step][0]}.h5ad"] == sha[step], step


def test_the_builder_takes_the_commands_the_launcher_runs(tmp_path, monkeypatch):
    sys.path.insert(0, str(REPO / "experiments" / "builds"))
    sys.path.insert(0, str(REPO / "experiments" / "external"))
    import build_test as bt
    import feng_relabel as fr

    ledger = tmp_path / "TRIES.md"  # a temporary ledger with the ruling rows' orders, not the repository's
    ledger.write_text("# ledger\n| R20 (the captain's order of X; ruling, no run; `--feng-fraction-mult 2.0` "
                      "`--captain-ruling R20`) |\n"
                      "| R22 (the captain's order of X; ruling, no run; `--feng-fraction-mult 2.5` "
                      "`--feng-fraction-mult 3.0` `--g-mom-damping 1.5` `--g-mom-damping 2.0` "
                      "`--g-mom-damping 1.0 --feng-fraction-mult 2.0` `--g-mom-damping 1.0 --feng-fraction-mult 2.5` "
                      "`--g-mom-damping 1.5 --feng-fraction-mult 2.0` `--g-mom-damping 2.0 --feng-fraction-mult 2.0` "
                      f"{OPT_CM30}, {OPT_B25D15} and {OPT_B25D20}, each only with `--captain-ruling R22`) |\n"
                      "| R23 (the captain's order of X; ruling, no run; `--direction-momentum all` "
                      "`--captain-ruling R23`) |\n"
                      "| R26 (the captain's order of X; ruling, no run; `--feng-fraction-mult 1.5 --direction-momentum "
                      "all` `--captain-ruling R26`) |\n"
                      f"| R37 (the captain's named exception of X; ruling, no run; {OPT_POOLED} only with "
                      "`--captain-ruling R37`) |\n")
    monkeypatch.setattr(bt, "TRIES_FILE", ledger)
    monkeypatch.setattr(fr, "committed_at_head", lambda p: json.dumps({"route": "recomputed"}).encode())
    args = _args(_sh())
    want = {"BASE": None, "A": (1.0, 1.0, None), "B": (0.5, 1.5, None), "C": (1.0, 1.5, None),
            "X2": (0.5, 2.0, "R20"), "X25": (0.5, 2.5, "R22"), "A15": (1.5, 1.0, "R22"), "A20": (2.0, 1.0, "R22"),
            "CM20": (1.0, 2.0, "R22"), "MS": None, "MSB": (0.5, 1.5, None), "X30": (0.5, 3.0, "R22"),
            "CM25": (1.0, 2.5, "R22"), "B2D15": (1.5, 2.0, "R22"), "B2D20": (2.0, 2.0, "R22"),
            "CM30": (1.0, 3.0, "R22"), "B25D15": (1.5, 2.5, "R22"), "B25D20": (2.0, 2.5, "R22"), "POOLED": None}
    assert set(want) == set(args)
    for key, values in want.items():
        a = bt.build_parser().parse_args(["t1", *args[key].split(), "--out", "x.h5ad"])
        assert a.board is None and a.seed == bt.harness.BUILD_SEED_OFFSET + 1
        ms, cand, pooled = bt._t1_ms(a), bt._t1_candidate(a), bt._t1_pooled(a)  # the refusals pass: nothing is read
        assert bt._bet_refusals(a) is None and bt._feng_recomputed_refusals(a) is None  # but the ledger
        got = None if cand is None else (cand["g_mom_damping"], cand["feng_fraction_mult"], cand.get("captain_ruling"))
        assert got == values, key
        assert (ms is not None) == (key in ("MS", "MSB")) and (pooled is not None) == (key == "POOLED")
        if key == "CM30":
            assert "seventh T1 search line" in cand["order"]
        if key.startswith("B25D"):
            assert "sixth T1 search line" in cand["order"]
    a = bt.build_parser().parse_args(["t1", *args["POOLED"].split(), "--out", "x.h5ad"])
    assert a.base == "pooled" and bt._t1_pooled(a)["captain_ruling"] == "R37"


def test_the_gate_lines():
    sh = _sh()
    gate = sh[sh.index("STEP=gate"):sh.index("STEP=tests")]
    assert "git status --porcelain --untracked-files=no" in gate
    assert "sha256sum -c docs/S2_RULES.sha256 > /dev/null 2>&1 || finish 97" in gate
    assert "[ \"$(grep -c '^| R22 ' experiments/harness/TRIES.md)\" = 1 ]" in gate
    assert "first=$(grep '^| R22 ' experiments/harness/TRIES.md | cut -d'|' -f2)" in gate
    assert (f"for opt in '{OPT_CM30}' '{OPT_B25D15}' \\\n  '{OPT_B25D20}' '`--captain-ruling R22`'; do") in gate
    assert "[ \"$(grep -c '^| R37 ' experiments/harness/TRIES.md)\" = 1 ]" in gate
    assert "first37=$(grep '^| R37 ' experiments/harness/TRIES.md | cut -d'|' -f2)" in gate
    assert "grep -q \"^ R37 (the captain's named exception\"" in gate
    assert f"for opt in '{OPT_POOLED}' '`--captain-ruling R37`'; do" in gate
    assert "[ -d data/derived/feng_relabel_c ]" in gate
    assert "O=runs/t1l67pb_${head:0:7}_$RUN" in sh and "S=$HOME/v/scratch/t1l67pb/$RUN" in sh
    codes = [int(c) for c in re.findall(r"finish (\d+)", gate)]
    assert codes == [96, 98, 97, 95, 95, 94, 94, 94, 92, 93]


def test_the_tests_step_runs_the_builders_and_launchers_targeted_tests():
    sh = _sh()
    step = sh[sh.index("STEP=tests"):sh.index("step_ok", sh.index("STEP=tests"))]
    listed = re.findall(r"tests/test_[a-z0-9_]+\.py", step)
    assert {"tests/test_t1_l67pb_launcher.py", "tests/test_t1_l67pb_build.py", "tests/test_t1_b2d_line_launcher.py",
            "tests/test_t1_line5_build.py", "tests/test_t1_line3_build.py", "tests/test_t1_damping_line_build.py",
            "tests/test_t1_fraction_line_build.py", "tests/test_t1_fraction2_build.py", "tests/test_t1_ms_build.py",
            "tests/test_t1_msb_build.py", "tests/test_t1_candidates_build.py", "tests/test_hidden_counts.py",
            "tests/test_s2_rules.py"} <= set(listed)
    assert len(listed) == len(set(listed)) and all((REPO / f).exists() for f in listed)
    assert "|| finish 10" in step and '--basetemp="$S/tmp/pytest"' in step


def test_the_launcher_deletes_nothing_and_opens_no_quote_in_a_usage_message():
    code = [line.split("#")[0] for line in _sh().splitlines()]
    assert not any(re.search(r"(^|[\s;&|(])(rm|rmdir|unlink|shred|truncate|mv)\s", line) for line in code)
    assert not any(re.search(r"find .*-delete|git (clean|reset|checkout|stash)", line) for line in code)
    assert "rm " not in _sh()
    usages = re.findall(r"\$\{\d+:\?([^}]*)\}", _sh())
    assert len(usages) == 2
    for msg in usages:  # an apostrophe there opens a quote
        assert "'" not in msg and '"' not in msg, msg


# ------------------------------------------------------------------ the gate, run for real


ROW_R22 = (f"| R22 (the captain's order of a later day; ruling, no run; {OPT_CM30}, {OPT_B25D15} and {OPT_B25D20}, "
           "each only with `--captain-ruling R22`) | T1 | the lines |\n")
ROW_R37 = (f"| R37 (the captain's named exception of a later day; ruling, no run; {OPT_POOLED} only with "
           "`--captain-ruling R37`) | T1 | the pooled base |\n")


def _env(tmp_path: Path) -> dict:
    env = {k: v for k, v in os.environ.items() if k not in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE")}
    env["HOME"] = str(tmp_path / "home")
    env["GIT_CEILING_DIRECTORIES"] = str(tmp_path)
    return env


def _git(repo: Path, env: dict, *argv: str) -> str:
    return subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", "-c", "commit.gpgsign=false", *argv],
                          cwd=repo, env=env, capture_output=True, text=True, check=True).stdout.strip()


def _repo(tmp_path: Path, tries: str = ROW_R22 + ROW_R37, rules: str = "the amendment names row R22\n",
          relabel: bool = True, good_sha: bool = True) -> tuple[Path, dict, str]:
    """A throwaway repository holding the launcher and the files its gate reads, and a home whose ~/.local/bin holds a
    `uv` that fails (the launcher puts that folder first on PATH): returns the repository, the env and its HEAD."""
    env = _env(tmp_path)
    fake = tmp_path / "home" / ".local" / "bin"
    fake.mkdir(parents=True, exist_ok=True)
    (fake / "uv").write_text("#!/bin/sh\nexit 1\n")
    (fake / "uv").chmod(0o755)
    repo = tmp_path / "repo"
    (repo / "experiments" / "builds").mkdir(parents=True)
    (repo / "experiments" / "harness").mkdir(parents=True)
    (repo / "docs").mkdir()
    (repo / "experiments" / "builds" / "t1_l67pb_line.sh").write_text(_sh())
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"# ledger\n{tries}")
    (repo / "docs" / "S2_RULES.md").write_text(rules)
    sha = hashlib.sha256(rules.encode()).hexdigest() if good_sha else "0" * 64
    (repo / "docs" / "S2_RULES.sha256").write_text(f"{sha}  docs/S2_RULES.md\n")
    (repo / ".gitignore").write_text("data/\nruns/\n")
    if relabel:
        (repo / "data" / "derived" / "feng_relabel_c").mkdir(parents=True)
    _git(repo, env, "init", "-q")
    _git(repo, env, "add", "-A")
    _git(repo, env, "commit", "-q", "-m", "throwaway")
    return repo, env, _git(repo, env, "rev-parse", "HEAD")


def _run(repo: Path, env: dict, *argv: str) -> subprocess.CompletedProcess:
    return subprocess.run(["bash", "experiments/builds/t1_l67pb_line.sh", *argv], cwd=repo, env=env,
                          capture_output=True, text=True)


def _scratch(tmp_path: Path, run: str) -> Path:
    return tmp_path / "home" / "v" / "scratch" / "t1l67pb" / run


def _stopped_at_the_gate(tmp_path: Path, run: str, code: int) -> None:
    s = _scratch(tmp_path, run)
    assert (s / "gate.exit").read_text().strip() == str(code)
    assert sorted(p.name for p in s.iterdir()) == ["gate.exit", "tmp"]  # no head file, no later step's mark
    assert not list((s / "tmp").iterdir())


def test_a_call_that_passes_the_gate_writes_the_head_and_stops_at_the_failing_tests_step(tmp_path):
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, head, "r1")
    assert out.returncode == 10, out.stderr
    s = _scratch(tmp_path, "r1")
    assert (s / "gate.exit").read_text().strip() == "0" and (s / "tests.exit").read_text().strip() == "10"
    assert (s / "head").read_text().strip() == head
    assert sorted(p.name for p in s.iterdir()) == ["gate.exit", "head", "tests.exit", "tests.out", "tmp"]
    o = repo / "runs" / f"t1l67pb_{head[:7]}_r1"
    assert o.is_dir() and not list(o.iterdir())  # nothing built, no depth record
    before = {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()}
    again = _run(repo, env, head, "r1")
    assert again.returncode == 91 and "never reused" in again.stderr
    assert {p.name: p.read_bytes() for p in s.iterdir() if p.is_file()} == before
    other = _run(repo, env, head, "r2")
    assert other.returncode == 10 and (repo / "runs" / f"t1l67pb_{head[:7]}_r2").is_dir()


def test_a_wrong_commit_and_a_prefix_are_refused_and_a_run_name_is_never_reused(tmp_path):
    repo, env, head = _repo(tmp_path)
    out = _run(repo, env, "0" * 40, "r1")
    assert out.returncode == 96 and "HEAD is not" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 96)
    assert not (repo / "runs").exists()
    again = _run(repo, env, head, "r1")
    assert again.returncode == 91 and "never reused" in again.stderr
    _stopped_at_the_gate(tmp_path, "r1", 96)
    short = _run(repo, env, head[:12], "r2")
    assert short.returncode == 96
    _stopped_at_the_gate(tmp_path, "r2", 96)


def test_a_missing_or_odd_run_name_is_refused_before_any_folder_is_made(tmp_path):
    repo, env, head = _repo(tmp_path)
    assert _run(repo, env, head).returncode != 0
    assert _run(repo, env).returncode != 0
    odd = _run(repo, env, head, "a/b")
    assert odd.returncode == 90 and not (tmp_path / "home" / "v").exists()


def test_a_changed_tracked_file_is_refused(tmp_path):
    repo, env, head = _repo(tmp_path)
    (repo / "experiments" / "harness" / "TRIES.md").write_text(f"# ledger\n{ROW_R22}{ROW_R37}more\n")
    out = _run(repo, env, head, "r1")
    assert out.returncode == 98 and "tracked files differ" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 98)


R22_BARE = ("| R22 (the captain's order; `--g-mom-damping 1.0` `--g-mom-damping 1.5` `--g-mom-damping 2.0` "
            "`--feng-fraction-mult 2.5` `--feng-fraction-mult 3.0` `--captain-ruling R22`) |\n")


@pytest.mark.parametrize("kw, code, says", [
    ({"good_sha": False}, 97, None),
    ({"tries": ROW_R37}, 95, "exactly one row R22"),
    ({"tries": f"| R22b (a row with another ID; {OPT_CM30} {OPT_B25D15} {OPT_B25D20} `--captain-ruling R22`) |\n"
               + ROW_R37}, 95, "exactly one row R22"),
    ({"tries": ROW_R22 + ROW_R22 + ROW_R37}, 95, "exactly one row R22"),  # two rows with one ID
    ({"tries": R22_BARE + ROW_R37}, 95, f"first cell does not name {OPT_CM30}"),  # the bare options only
    ({"tries": f"| R22 (the captain's order; {OPT_CM30} and {OPT_B25D20}, with `--captain-ruling R22`) |\n" + ROW_R37},
     95, f"first cell does not name {OPT_B25D15}"),
    ({"tries": f"| R22 (the captain's order; {OPT_CM30} and {OPT_B25D15}, with `--captain-ruling R22`) |\n" + ROW_R37},
     95, f"first cell does not name {OPT_B25D20}"),
    ({"tries": f"| R22 (the captain's order; {OPT_CM30} {OPT_B25D15} {OPT_B25D20}) |\n" + ROW_R37}, 95,
     "first cell does not name `--captain-ruling R22`"),
    ({"tries": f"| R22 (the captain's order) | T1 | {OPT_CM30} {OPT_B25D15} {OPT_B25D20} `--captain-ruling R22` |\n"
               + ROW_R37}, 95, f"first cell does not name {OPT_CM30}"),  # named outside the first cell only
    ({"tries": f"| R22 (the captain's order) |\n| R39 ({OPT_CM30} {OPT_B25D15} {OPT_B25D20} `--captain-ruling R22`) |\n"
               + ROW_R37}, 95, f"first cell does not name {OPT_CM30}"),  # named in another row only
    ({"tries": ROW_R22}, 94, "exactly one row R37"),
    ({"tries": ROW_R22 + ROW_R37 + ROW_R37}, 94, "exactly one row R37"),
    ({"tries": ROW_R22 + f"| R37 (the captain's order of a later day; ruling, no run; {OPT_POOLED} only with "
                         "`--captain-ruling R37`) |\n"}, 94, "is not the captain's named exception"),
    ({"tries": ROW_R22 + "| R37 (the captain's named exception; `--base pooled` `--feng-groups recomputed` "
                         "`--captain-ruling R37`) |\n"}, 94, f"first cell does not name {OPT_POOLED}"),
    ({"tries": ROW_R22 + f"| R37 (the captain's named exception; {OPT_POOLED}) |\n"}, 94,
     "first cell does not name `--captain-ruling R37`"),
    ({"tries": ROW_R22 + f"| R37 (the captain's named exception) | T1 | {OPT_POOLED} `--captain-ruling R37` |\n"}, 94,
     f"first cell does not name {OPT_POOLED}"),  # named outside the first cell only
    ({"relabel": False}, 92, "feng_relabel_c is not in this checkout"),
])
def test_each_gate_refusal(tmp_path, kw, code, says):
    repo, env, head = _repo(tmp_path, **kw)
    out = _run(repo, env, head, "r1")
    assert out.returncode == code, out.stderr
    if says:
        assert says in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", code)
    assert not (repo / "runs" / f"t1l67pb_{head[:7]}_r1").exists()


# A stand-in `uv` that lets a call run past the tests step: the tests pass, the depth record prints a stand-in record
# (or fails, with DEPTH_FAILS set), a build writes a stand-in file whose stamp names HEAD and clean, `vec check` passes
# and the stamp reader is python itself. No stand-in file can have a known file's bytes, so the run stops at step D.
PASSING_UV = """#!/bin/sh
shift 3
case "$1" in
  vec) exit 0 ;;
  python)
    shift 3
    case "$1" in
      -m) exit 0 ;;
      -c) shift; exec python3 -I -B -c "$@" ;;
      experiments/builds/t1_pool_depth_record.py)
        [ -n "${DEPTH_FAILS:-}" ] && exit 1
        echo '{"pools": "stand-in"}'; exit 0 ;;
      experiments/builds/build_test.py)
        out=""; prev=""
        for a in "$@"; do [ "$prev" = "--out" ] && out=$a; prev=$a; done
        echo "stand-in $*" > "$out"
        printf '{"commit": "%s", "clean": true}' "$(git rev-parse HEAD)" > "$out.stamp.json"
        exit 0 ;;
    esac ;;
esac
exit 1
"""


@pytest.mark.parametrize("depth_fails", [False, True])
def test_a_call_past_the_tests_step_writes_the_depth_record_and_stops_at_the_first_unknown_bytes(tmp_path,
                                                                                               depth_fails):
    repo, env, head = _repo(tmp_path)
    (tmp_path / "home" / ".local" / "bin" / "uv").write_text(PASSING_UV)
    if depth_fails:
        env["DEPTH_FAILS"] = "1"
    out = _run(repo, env, head, "r1")
    assert out.returncode == 20, out.stderr
    assert "stopped: step D does not give the known bytes of T1_s2rc" in out.stderr
    s = _scratch(tmp_path, "r1")
    codes = {p.name: p.read_text().strip() for p in s.glob("*.exit")}
    assert codes == {"gate.exit": "0", "tests.exit": "0", "depth.exit": "1" if depth_fails else "0", "D.exit": "20"}
    assert (s / "stamps").read_text() == f"D/T1_s2rc {head} True\n"
    o = repo / "runs" / f"t1l67pb_{head[:7]}_r1"
    assert sorted(p.name for p in o.iterdir()) == (["D"] if depth_fails else ["D", "pool_depth.json"])
    assert sorted(p.name for p in (o / "D").iterdir()) == ["T1_s2rc.h5ad", "T1_s2rc.h5ad.stamp.json"]
    assert (o / "D" / "T1_s2rc.h5ad").read_text() == (f"stand-in experiments/builds/build_test.py t1 {BASE} --out "
                                                      f"runs/t1l67pb_{head[:7]}_r1/D/T1_s2rc.h5ad\n")
    if not depth_fails:
        assert (o / "pool_depth.json").read_text() == (s / "pool_depth.json").read_text() == '{"pools": "stand-in"}\n'


def test_an_existing_output_folder_is_refused_and_left_as_it_was(tmp_path):
    repo, env, head = _repo(tmp_path)
    o = repo / "runs" / f"t1l67pb_{head[:7]}_r1"
    o.mkdir(parents=True)
    (o / "keep").write_text("x")
    out = _run(repo, env, head, "r1")
    assert out.returncode == 93 and "output folder" in out.stderr
    _stopped_at_the_gate(tmp_path, "r1", 93)
    assert sorted(p.name for p in o.iterdir()) == ["keep"]
