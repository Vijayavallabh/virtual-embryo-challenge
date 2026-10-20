"""Build stamps (vec.provenance), as vec.io.write_submission writes them, against each test's own git repository."""

import errno
import hashlib
import json
import os
import socket
import sys
from datetime import timedelta

import numpy as np
import pytest
from conftest import GENES, build, git, make_repo
from typer.testing import CliRunner

from vec import cli, provenance
from vec import submit as vs
from vec.boards import REPO_ROOT, Board
from vec.io import make_submission, write_submission
from vec.provenance import STAMP_MAX_BYTES, STAMP_SUFFIX, read_stamp, stamp_path


@pytest.fixture
def repo(tmp_path):
    """A git repository with one commit and a clean tree."""
    return make_repo(tmp_path / "repo")


def _adata(seed=0, n=40):
    rng = np.random.default_rng(seed)
    counts = rng.poisson(2.0, size=(n, len(GENES)))
    return make_submission(np.log1p(1e4 * counts / counts.sum(1, keepdims=True).clip(min=1)), GENES)


def _stamp(path):
    return json.loads(stamp_path(path).read_text())


def test_write_submission_stamps_the_bytes_it_wrote_with_the_commit_and_a_clean_tree(repo, tmp_path):
    out = tmp_path / "preds" / "s1.h5ad"
    path = build(_adata(), out, repo)
    data, s, head = path.read_bytes(), _stamp(path), git(repo, "rev-parse", "HEAD")
    assert sorted(p.name for p in out.parent.iterdir()) == ["s1.h5ad", "s1.h5ad.stamp.json"]  # no temporary file left
    assert (s["version"], s["file"], s["sha256"], s["size"]) == (1, "s1.h5ad", hashlib.sha256(data).hexdigest(),
                                                                 len(data))
    assert (s["commit"], s["clean"], s["why"], s["unclean"], s["unclean_count"]) == (head, True, None, [], 0)
    assert s["start_commit"] == head and s["started_utc"] <= s["created_utc"]
    assert s["hostname"] == socket.gethostname() and s["argv"] == sys.argv
    assert all(vs.datetime.fromisoformat(s[k]).utcoffset() == timedelta(0) for k in ("started_utc", "created_utc"))
    assert read_stamp(stamp_path(path)) == s


def test_stamp_false_writes_no_stamp(repo, tmp_path, monkeypatch):
    monkeypatch.setattr(provenance, "git_state", lambda *a, **k: pytest.fail("the state was read"))
    path = write_submission(_adata(), tmp_path / "out" / "tmp.h5ad", stamp=False, repo=repo)
    assert [p.name for p in path.parent.iterdir()] == ["tmp.h5ad"]


def test_a_dirty_tree_is_recorded_with_its_paths(repo):
    head = git(repo, "rev-parse", "HEAD")
    (repo / "method.py").write_text("x = 2\n")  # an uncommitted change to a tracked file
    (repo / "notes.txt").write_text("an untracked file outside the output directories\n")
    for allowed in ("results/out.csv", "runs/r1/log.txt", "experiments/t9/results/a.json"):
        (repo / allowed).parent.mkdir(parents=True, exist_ok=True)
        (repo / allowed).write_text("outputs may stay untracked\n")
    out = repo / "preds" / "s1.h5ad"  # inside the repository, outside the output directories, and not ignored here
    for _ in range(2):  # a rebuild: the file and its earlier stamp, there when it starts, are its own, not unclean code
        s = _stamp(build(_adata(), out, repo))
        assert (s["commit"], s["clean"]) == (head, False)
        assert (sorted(s["unclean"]), s["unclean_count"]) == (["method.py", "notes.txt"], 2)
        assert s["why"] == f"uncommitted or untracked: {', '.join(s['unclean'])}"
    for i in range(25):
        (repo / f"n{i:02d}.txt").write_text("y\n")
    s = _stamp(build(_adata(), out, repo))
    assert len(s["unclean"]) == provenance.STAMP_MAX_PATHS == 20 and s["unclean_count"] == 27
    assert s["why"] == f"uncommitted or untracked: {', '.join(s['unclean'][:5])} ..."
    git(repo, "checkout", "-q", "--", "method.py")
    for p in repo.glob("n[0-9][0-9].txt"):
        p.unlink()
    (repo / "notes.txt").unlink()
    s = _stamp(build(_adata(), out, repo))
    assert (s["clean"], s["why"], s["unclean"], s["unclean_count"]) == (True, None, [], 0)


def test_the_stamp_compares_the_state_when_the_build_started(repo, tmp_path, monkeypatch):
    head, out = git(repo, "rev-parse", "HEAD"), tmp_path / "s1.h5ad"
    s = _stamp(write_submission(_adata(), out, repo=repo))  # no start recorded for this repository
    assert (s["clean"], s["why"], s["commit"], s["start_commit"], s["started_utc"], s["unclean"]) == (
        False, "no record of the repository's state when the build started", head, "unknown", None, [])
    with monkeypatch.context() as m:
        m.setenv("PATH", "")  # no git when the build starts
        assert provenance._record_start(repo).commit == "unknown"
    s = _stamp(write_submission(_adata(), out, repo=repo))
    assert (s["clean"], s["why"], s["commit"], s["start_commit"], s["unclean"]) == (
        False, "git could not tell the repository's state when the build started", head, "unknown", [])

    (repo / "notes.txt").write_text("y\n")  # untracked when the build starts, gone when the file is written
    provenance._record_start(repo)
    (repo / "notes.txt").unlink()
    (repo / "method.py").write_text("x = 2\n")  # and an uncommitted change only when it is written
    s = _stamp(write_submission(_adata(), out, repo=repo))
    assert (s["clean"], s["why"], s["unclean"], s["unclean_count"]) == (
        False, "uncommitted or untracked: method.py", ["method.py", "notes.txt"], 2)
    git(repo, "checkout", "-q", "--", "method.py")
    s = _stamp(write_submission(_adata(), out, repo=repo))  # the same build: it started with notes.txt there
    assert (s["clean"], s["why"], s["unclean"]) == (
        False, "uncommitted or untracked when the build started: notes.txt", ["notes.txt"])

    provenance._record_start(repo)
    git(repo, "commit", "-q", "--allow-empty", "-m", "a commit that lands while the build runs")
    later = git(repo, "rev-parse", "HEAD")
    s = _stamp(write_submission(_adata(), out, repo=repo))
    assert (s["clean"], s["why"], s["commit"], s["start_commit"], s["unclean"]) == (
        False, f"HEAD moved during the build, from {head[:12]} to {later[:12]}", later, head, [])
    s = _stamp(build(_adata(), out, repo))  # a new build, at the new HEAD
    assert (s["clean"], s["commit"], s["start_commit"]) == (True, later, later)


def test_head_moving_away_and_back_within_one_second_is_seen(repo, tmp_path, monkeypatch):
    """A reflog entry records its time in seconds, so a checkout away and back within the second of the newest entry
    repeats that entry exactly: the stamp also compares how many entries there are."""
    monkeypatch.setenv("GIT_COMMITTER_DATE", "1700000000 +0000")  # every reflog entry in the same second
    branch = git(repo, "rev-parse", "--abbrev-ref", "HEAD")
    git(repo, "checkout", "-qb", "side")
    git(repo, "checkout", "-q", branch)
    start = provenance._record_start(repo)
    git(repo, "checkout", "-q", "side")
    git(repo, "checkout", "-q", branch)
    now = provenance.repo_state(repo)
    assert (now.commit, now.reflog) == (start.commit, (start.reflog[0] + 2, start.reflog[1]))  # the same newest entry
    s = _stamp(write_submission(_adata(), tmp_path / "s1.h5ad", repo=repo))
    assert (s["clean"], s["why"]) == (
        False, f"HEAD moved during the build and back to {now.commit[:12]} (checkout: moving from side to {branch})")


def test_a_committed_ledger_with_new_rows_is_exempt_as_at_upload(repo, tmp_path):
    ledger = repo / provenance.LEDGER_PATH
    ledger.parent.mkdir()
    ledger.write_text("utc_time\n")
    git(repo, "add", provenance.LEDGER_PATH)
    git(repo, "commit", "-qm", "ledger")
    ledger.write_text("utc_time\n2026-11-20T00:00:00+00:00\n")  # a format check appended a row, before the build began
    assert _stamp(build(_adata(), tmp_path / "s1.h5ad", repo))["clean"]
    assert vs.LEDGER == REPO_ROOT / provenance.LEDGER_PATH  # the ledger vec.submit writes


def test_a_git_failure_never_fails_the_write(repo, tmp_path, monkeypatch):
    plain = tmp_path / "plain"
    plain.mkdir()
    monkeypatch.setenv("GIT_CEILING_DIRECTORIES", str(tmp_path))  # git never looks above it for a repository

    def boom(*args, **kwargs):
        raise RuntimeError("an unexpected failure")

    provenance._record_start(repo)
    for where, patch in ((plain, None), (repo, ("setenv", "PATH", "")), (repo, ("setattr", "git_state", boom))):
        with monkeypatch.context() as m:
            if patch and patch[0] == "setenv":
                m.setenv(*patch[1:])  # no git on the PATH
            elif patch:
                m.setattr(provenance, *patch[1:])
            path = write_submission(_adata(), tmp_path / "out" / "s1.h5ad", repo=where)
        s = _stamp(path)
        assert path.is_file() and s["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest(), where
        assert (s["commit"], s["clean"], s["unclean"], s["unclean_count"]) == ("unknown", False, None, None), patch
        assert s["why"] == "git cannot tell the repository's state when the file was written"


def test_a_failed_stamp_write_leaves_the_earlier_stamp_and_no_temporary_file(repo, tmp_path, monkeypatch):
    path = build(_adata(0), tmp_path / "s1.h5ad", repo)
    before = stamp_path(path).read_text()
    write_submission(_adata(1), path, stamp=False)  # rebuilt; stamping it fails at the rename

    def full(src, dst):
        raise OSError(errno.ENOSPC, "No space left on device")

    monkeypatch.setattr(os, "replace", full)
    with pytest.raises(OSError, match="No space"):
        provenance.write_stamp(path, repo)
    monkeypatch.undo()
    assert stamp_path(path).read_text() == before
    assert sorted(p.name for p in tmp_path.iterdir() if p.is_file()) == ["s1.h5ad", "s1.h5ad.stamp.json"]


def test_the_default_repository_is_the_packages(repo, tmp_path, monkeypatch):
    assert provenance.REPO_ROOT == REPO_ROOT == provenance.Path(provenance.__file__).resolve().parents[2]
    assert REPO_ROOT in provenance._STARTS  # its start: recorded when this process first imported vec.provenance
    monkeypatch.setattr(provenance, "REPO_ROOT", repo)  # stands in for this package's repository
    provenance._record_start()  # as at that import
    s = _stamp(write_submission(_adata(), tmp_path / "s1.h5ad"))
    assert s["commit"] == git(repo, "rev-parse", "HEAD") and s["clean"]


def test_a_stamp_sits_next_to_the_real_file_a_link_points_at(repo, tmp_path):
    target = tmp_path / "runs" / "s1.h5ad"
    build(_adata(0), target, repo)
    link = tmp_path / "latest.h5ad"
    link.symlink_to(target)
    build(_adata(1), link, repo)  # written through the link: the target holds the new bytes
    assert not stamp_path(link).exists() and _stamp(target)["file"] == "s1.h5ad"
    assert _stamp(target)["sha256"] == hashlib.sha256(target.read_bytes()).hexdigest()


def test_vec_baseline_stamps_its_file(own_repo, tmp_path, monkeypatch):
    board = Board(key="T1:x", task="T1", genes=GENES, needs_coords=False, min_cells=10, max_cells=1000)
    monkeypatch.setattr(cli, "get_board", lambda key: board)
    last, out = tmp_path / "last.h5ad", tmp_path / "pred.h5ad"
    _adata(3).write_h5ad(last)
    res = CliRunner().invoke(cli.app, ["baseline", "copy_last", "--board", "T1:x", "--last", str(last), "--out",
                                       str(out)])
    assert res.exit_code == 0, res.output
    s = read_stamp(stamp_path(out))
    assert s["sha256"] == hashlib.sha256(out.read_bytes()).hexdigest() and s["clean"]
    assert s["commit"] == git(own_repo, "rev-parse", "HEAD")


def test_vec_baseline_keeps_n_random_rows(own_repo, tmp_path, monkeypatch):
    import anndata as ad

    board = Board(key="T1:x", task="T1", genes=GENES, needs_coords=False, min_cells=10, max_cells=1000)
    monkeypatch.setattr(cli, "get_board", lambda key: board)
    last = tmp_path / "last.h5ad"
    _adata(3).write_h5ad(last)

    def run(*extra, out):
        return CliRunner().invoke(cli.app, ["baseline", "copy_last", "--board", "T1:x", "--last", str(last),
                                            *extra, "--out", str(out)])

    files = [tmp_path / f"p{i}.h5ad" for i in range(3)]
    for f, seed in zip(files, (0, 0, 1), strict=True):
        res = run("--n-cells", "12", "--seed", str(seed), out=f)
        assert res.exit_code == 0, res.output
    def dense(f):
        X = ad.read_h5ad(f).X
        return np.asarray(X.todense() if hasattr(X, "todense") else X)

    x = [dense(f) for f in files]
    assert all(m.shape == (12, len(GENES)) for m in x)
    assert np.array_equal(x[0], x[1]) and not np.array_equal(x[0], x[2])  # the seed fixes the rows
    assert read_stamp(stamp_path(files[0]))["sha256"] == hashlib.sha256(files[0].read_bytes()).hexdigest()
    for bad in ("41", "-1"):
        res = run("--n-cells", bad, out=tmp_path / "bad.h5ad")
        assert res.exit_code != 0 and not (tmp_path / "bad.h5ad").exists()


def test_read_stamp_refuses_what_is_not_a_stamp_of_this_version(repo, tmp_path):
    path = stamp_path(build(_adata(), tmp_path / "s1.h5ad", repo))
    good = json.loads(path.read_text())
    cases = [(b"{" + b" " * STAMP_MAX_BYTES + b"}", "larger than 64 KiB"), (b"not json", "not JSON"),
             (b"\xff\xfe\x00", "not JSON"), (b"[1]", "no version"),
             (json.dumps({**good, "version": True}), "no version"),
             (json.dumps({**good, "version": 2}), "version is 2"),
             (json.dumps({**good, "sha256": good["sha256"].upper()}), "malformed sha256"),
             (json.dumps({**good, "commit": "HEAD"}), "malformed commit"),
             (json.dumps({**good, "commit": good["commit"][:12]}), "malformed commit"),
             (json.dumps({**good, "clean": 1}), "malformed clean"),
             (json.dumps({**good, "why": ["moved"]}), "malformed why"),
             (json.dumps({**good, "size": -1, "unclean": "method.py", "unclean_count": "1"}),
              "malformed size, unclean, unclean_count")]
    for raw, message in cases:
        path.write_bytes(raw if isinstance(raw, bytes) else raw.encode())
        with pytest.raises(ValueError, match=message):
            read_stamp(path)
    path.write_text(json.dumps({**good, "commit": "unknown", "unclean": None, "unclean_count": None}))
    assert read_stamp(path)["commit"] == "unknown"
    with pytest.raises(FileNotFoundError):
        read_stamp(tmp_path / ("other.h5ad" + STAMP_SUFFIX))


def test_is_ancestor(repo):
    first = git(repo, "rev-parse", "HEAD")
    (repo / "method.py").write_text("x = 2\n")
    git(repo, "commit", "-qam", "second")
    second = git(repo, "rev-parse", "HEAD")
    assert provenance.is_ancestor(repo, first) is True and provenance.is_ancestor(repo, second) is True
    git(repo, "checkout", "-q", first)
    assert provenance.is_ancestor(repo, second) is False
    for name in ("f" * 40, "HEAD", "--all", second[:12]):  # a commit the repository lacks; not a full hex name
        assert provenance.is_ancestor(repo, name) is None


def test_the_clean_tree_rule_keeps_its_names_in_vec_submit():
    for name in ("OUTPUT_DIRS", "_output_path", "git_state", "_repo_paths", "unclean_paths"):
        assert getattr(vs, name) is getattr(provenance, name), name
