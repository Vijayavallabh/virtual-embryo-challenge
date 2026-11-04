"""Submission client against a local mock of the challenge API (127.0.0.1 only; no test reaches the internet)."""

import dataclasses
import errno
import hashlib
import io
import json
import os
import re
import shlex
import shutil
import signal
import socket
import subprocess
import sys
import tempfile
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pytest
from conftest import GENES, build, make_repo
from conftest import git as _git
from rich.console import Console
from typer.testing import CliRunner

from vec import cli
from vec import submit as vs
from vec.boards import load_boards
from vec.io import CheckReport, make_submission

TOKEN = "tok-FAKE-never-print-3f9c1a77"
T1 = "T1:test"
X001 = "X001"
REAL_DGX_RUN = vs.dgx_run  # the dgx fixture replaces it


def _parse_multipart(content_type: str, body: bytes):
    boundary = content_type.split("boundary=", 1)[1].encode()
    fields, filename, data = {}, None, b""
    for part in body.split(b"--" + boundary)[1:]:
        if part.startswith(b"--"):
            break
        head, _, payload = part[2:].partition(b"\r\n\r\n")
        payload = payload[:-2]
        disp = head.decode()
        name = re.search(r'; name="([^"]*)"', disp).group(1)
        fn = re.search(r'; filename="([^"]*)"', disp)
        if fn:
            filename, data = fn.group(1), payload
        else:
            fields[name] = payload.decode()
    return fields, filename, data


def _bundle(fields):
    if fields["task"] == "T2":
        return f"T2:{fields['setting']}:test_{fields['mode']}"
    if fields["task"] == "T3":
        return f"T3:{fields['target']}"
    return T1


class MockApi:
    """The API as far as the client uses it: idempotency (identical replay, 409 on a changed payload), quota, a
    request log, and scripted POST outcomes: "503" (nothing stored), "store_503" (stored, response lost), "409",
    "307" (a redirect, never followed), "reject" (stored as rejected), "echo_401" (an error detail that echoes the
    token), "echo_401_long" (the same with the token across the client's 300-character cut), "echo_ok" (a 200 format
    answer that echoes it), "echo_id" (stored under an id that echoes it), "markup" (stored with a status holding rich
    markup), "store_unreadable" (stored, answered 200 with an HTML page), "format_stored" (a format check stored as a
    submission), "format_unconfirmed" (a format answer saying format_only is false). With echo_status set, the status
    refresh (GET /challenge/submissions/{id}) echoes the token."""

    def __init__(self):
        self.requests, self.subs, self.keys, self.script = [], [], {}, []
        self.phase = {"phase": "p3", "label": "P3", "split": "test", "accepts_submissions": True}
        self.quota = None  # None: computed from subs; a dict body; or an int HTTP status
        self.subs_status = 200
        self.echo_status = False

    def official_per_board(self):
        """Assumed shape (not published): every current board with its phase-wide used count and limit."""
        out = {b: {"used": 0, "limit": 2} for b in (T1, "T2:heart:test_extrap")}
        for s in self.subs:
            if s["status"] != "rejected":
                out.setdefault(s["bundle_key"], {"used": 0, "limit": 2})["used"] += 1
        return out

    def posts(self, official=None):
        p = [r for r in self.requests if r["method"] == "POST"]
        return p if official is None else [r for r in p if ("format_only" not in r["fields"]) == official]


def _handler(api: MockApi):
    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *args):
            pass

        def _send(self, code, body):
            raw = body if isinstance(body, bytes) else json.dumps(body).encode()
            self.send_response(code)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(raw)))
            self.end_headers()
            self.wfile.write(raw)

        def _authed(self):
            return self.headers.get("Authorization") == f"Bearer {TOKEN}"

        def do_GET(self):
            api.requests.append({"method": "GET", "path": self.path, "auth": self.headers.get("Authorization")})
            if self.path == "/challenge/phase":
                return self._send(200, api.phase)
            if not self._authed():
                return self._send(401, {"detail": "sign in first"})
            if self.path == "/challenge/integration/quota":
                if isinstance(api.quota, int):
                    return self._send(api.quota, {"detail": "quota unavailable"})
                daily = {"T1": {"used": 0, "limit": 8, "uploads": 0, "upload_limit": 16}}
                return self._send(200, api.quota if api.quota is not None else
                                  {"daily": daily, "quota_credit": 0, "official_per_board": api.official_per_board()})
            if self.path == "/challenge/submissions":
                return self._send(api.subs_status, {"registered": True, "submissions": api.subs})
            if self.path.startswith("/challenge/submissions/"):
                if api.echo_status:
                    return self._send(200, {"id": self.path.rsplit("/", 1)[1],
                                            "status": f"scored {self.headers.get('Authorization')}"})
                hit = [s for s in api.subs if s["id"] == self.path.rsplit("/", 1)[1]]
                return self._send(200, hit[0]) if hit else self._send(404, {"detail": "no such submission"})
            return self._send(404, {"detail": "not found"})

        def do_POST(self):
            body = self.rfile.read(int(self.headers["Content-Length"]))
            fields, filename, data = _parse_multipart(self.headers["Content-Type"], body)
            rec = {"method": "POST", "path": self.path, "auth": self.headers.get("Authorization"),
                   "key": self.headers.get("Idempotency-Key"), "fields": fields, "filename": filename,
                   "sha256": hashlib.sha256(data).hexdigest()}
            api.requests.append(rec)
            action = api.script.pop(0) if api.script else ""
            auth = self.headers.get("Authorization")
            if action in ("echo_401", "echo_401_long"):
                text = "x" * 270 if action == "echo_401_long" else "rejected credentials"
                return self._send(401, {"detail": f"{text} {auth}"})
            if action == "echo_ok":  # a 200 answer whose text echoes the token
                return self._send(200, {"ok": False, "problems": [f"bad header {auth}"], "format_only": True})
            if self.path != "/challenge/submissions" or not self._authed():
                return self._send(401, {"detail": "sign in first"})
            if action in ("503", "409", "307"):
                return self._send(int(action), {"detail": {"503": "busy", "409": "accept the Rules first",
                                                           "307": "moved"}[action]})
            if fields.get("format_only") == "true" and action != "format_stored":
                return self._send(200, {"ok": True, "problems": [], "format_only": action != "format_unconfirmed"})
            digest, key = json.dumps([fields, filename, rec["sha256"]]), rec["key"]
            if key in api.keys:
                stored, resp = api.keys[key]
                return self._send(200, resp) if stored == digest else self._send(409, {"detail": "key reused"})
            sub = {"id": f"sub-{auth.split()[-1]}" if action == "echo_id" else f"sub{len(api.subs) + 1}",
                   "task": fields["task"], "bundle_key": _bundle(fields),
                   "status": {"reject": "rejected", "markup": "queued [/bold]"}.get(action, "queued"),
                   "sha256": rec["sha256"], "filename": filename}
            api.subs.append(sub)
            resp = {**sub, "phase": "p3", "bytes": len(data), "problems": [], "used_credit": False}
            if key:
                api.keys[key] = (digest, resp)
            if action == "store_503":
                return self._send(503, {"detail": "lost"})
            return self._send(200, b"<html>accepted</html>" if action == "store_unreadable" else resp)

        def do_PUT(self):
            api.requests.append({"method": self.command, "path": self.path, "fields": {}})
            self._send(405, {"detail": "not allowed"})

        do_DELETE = do_PATCH = do_PUT

    return Handler


@pytest.fixture(autouse=True)
def no_internet(monkeypatch):
    real = socket.create_connection

    def local_only(address, *args, **kwargs):
        if address[0] not in ("127.0.0.1", "localhost"):
            raise AssertionError(f"a test tried to reach {address[0]}")
        return real(address, *args, **kwargs)

    monkeypatch.setattr(socket, "create_connection", local_only)


@pytest.fixture
def api():
    mock = MockApi()
    server = ThreadingHTTPServer(("127.0.0.1", 0), _handler(mock))
    threading.Thread(target=server.serve_forever, kwargs={"poll_interval": 0.02}, daemon=True).start()
    mock.url = f"http://127.0.0.1:{server.server_address[1]}"
    yield mock
    server.shutdown()
    server.server_close()


@pytest.fixture
def env(tmp_path, monkeypatch, api):
    """Isolated boards, disclosure register, ledger and a clean git repo; the module defaults point at them."""
    boards = tmp_path / "boards"
    boards.mkdir()
    (boards / "genes.txt").write_text("\n".join(GENES) + "\n")
    spec = {"genes_file": "genes.txt", "needs_coords": False, "min_cells": 10, "max_cells": 500}
    index = {T1: {"key": T1, "task": "T1", "split": "test", **spec},
             "T1:val": {"key": "T1:val", "task": "T1", "split": "val", **spec},
             "T2:heart:test_extrap": {**spec, "key": "T2:heart:test_extrap", "task": "T2", "setting": "heart",
                                      "split": "test_extrap", "needs_coords": True},
             "T2:heart:test": {**spec, "key": "T2:heart:test", "task": "T2", "setting": "heart", "split": "test"}}
    (boards / "index.json").write_text(json.dumps(index))
    register = tmp_path / "EXTERNAL_DATA.md"
    register.write_text("| ID | Source |\n|---|---|\n| X001 | veckit (local scoring) |\n| X002 | Feng 2022 |\n")
    repo = make_repo(tmp_path / "repo")
    monkeypatch.setenv("VEC_TOKEN", TOKEN)
    monkeypatch.setattr(vs, "API_BASE", api.url)
    monkeypatch.setattr(vs, "BOARDS_DIR", boards)
    monkeypatch.setattr(vs, "EXTERNAL_DATA", register)
    monkeypatch.setattr(vs, "LEDGER", tmp_path / "submissions.csv")
    monkeypatch.setattr(vs, "REPO_ROOT", repo)
    # wider than any path (PATH_MAX is 4096): no printed line wraps and no table cell is cropped, wherever pytest puts
    # its tmp files
    monkeypatch.setattr(cli, "console", Console(emoji=False, highlight=False, width=10_000))
    return SimpleNamespace(tmp=tmp_path, boards=boards, repo=repo, ledger=tmp_path / "submissions.csv", api=api)


def _pred(env, name, seed, coords=False, where=None):
    """A valid prediction file in `where` (default env.tmp), built as in a new process and stamped by write_submission
    against the test's repository (clean, at its HEAD, unless a test has changed it)."""
    rng = np.random.default_rng(seed)
    counts = rng.poisson(2.0, size=(40, len(GENES)))
    X = np.log1p(1e4 * counts / counts.sum(1, keepdims=True).clip(min=1))
    xyz = rng.normal(size=(40, 3)) if coords else None
    return build(make_submission(X, GENES, xyz), (where or env.tmp) / name, env.repo)


def _official(path, **kw):
    return vs.submit(path, kw.pop("board", T1), official=True, confirm=vs.sha256_file(path)[:12],
                     disclose=kw.pop("disclose", X001), sleep=lambda s: None, **kw)


def _checked(env, name, seed):
    """A prediction file that has passed the server format check (the official guard needs one)."""
    path = _pred(env, name, seed)
    assert vs.submit(path, T1).server_status == "format_ok"
    return path


def test_default_sends_format_only_true_and_no_idempotency_key(env):
    f = _pred(env, "a.h5ad", 0)
    res = vs.submit(f, T1)
    (post,) = env.api.posts()
    assert post["fields"] == {"task": "T1", "format_only": "true"} and post["key"] is None
    assert post["filename"] == "a.h5ad" and post["sha256"] == vs.sha256_file(f)
    assert res.ok and res.server_status == "format_ok" and not env.api.subs
    gets = [(r["path"], r["auth"]) for r in env.api.requests if r["method"] == "GET"]
    assert gets == [("/challenge/phase", None)]  # a format check reads only the public phase, without the token

    heart = _pred(env, "heart.h5ad", 1, coords=True)
    vs.submit(heart, "T2:heart:test_extrap")
    assert env.api.posts()[-1]["fields"] == {"task": "T2", "setting": "heart", "mode": "extrap", "format_only": "true"}


def test_official_needs_both_flags_and_the_right_confirm(env):
    f = _checked(env, "a.h5ad", 0)
    sha = vs.sha256_file(f)
    with pytest.raises(vs.SubmitRefused, match="confirm"):
        vs.submit(f, T1, official=True, disclose=X001)
    with pytest.raises(vs.SubmitRefused, match="confirm") as e:
        vs.submit(f, T1, official=True, confirm="0" * 12, disclose=X001)
    assert e.value.plan.failed == ["confirm"]
    res = vs.submit(f, T1, confirm=sha[:12])  # --confirm without --official stays a format check
    assert res.server_status == "format_ok" and any("ignored" in w for w in res.plan.warnings)
    assert not env.api.posts(official=True)

    res = vs.submit(f, T1, official=True, confirm=sha[:12], disclose=X001)
    (post,) = env.api.posts(official=True)
    assert res.ok and res.server_id == "sub1" and res.server_status == "queued"
    assert post["fields"] == {"task": "T1"} and post["key"] == vs.idempotency_key(T1, sha) == f"vec.T1:test.{sha}"


def test_third_distinct_upload_refused(env):
    files = [_checked(env, f"p{i}.h5ad", i) for i in range(3)]
    assert _official(files[0]).ok and _official(files[1]).ok
    with pytest.raises(vs.SubmitRefused, match="budget") as e:
        _official(files[2])
    assert "2 of 2 spent" in next(g.detail for g in e.value.plan.guards if g.name == "budget")
    assert len(env.api.posts(official=True)) == 2


def test_server_quota_counts_uploads_the_ledger_never_saw(env):
    f = _checked(env, "a.h5ad", 0)
    env.api.subs += [{"id": f"web{i}", "task": "T1", "bundle_key": T1, "status": "scored", "sha256": str(i) * 64}
                     for i in range(2)]  # two uploads from the website
    with pytest.raises(vs.SubmitRefused, match="budget"):
        _official(f)
    env.api.subs.clear()
    env.api.quota = {"official_per_board": {T1: {"used": 2, "limit": 2}}}  # the quota endpoint alone says full
    with pytest.raises(vs.SubmitRefused, match="budget"):
        _official(f)
    assert not env.api.posts(official=True)


def test_ledger_counts_uploads_the_server_does_not_report_yet(env):
    f = _checked(env, "a.h5ad", 0)
    with vs.Ledger(env.ledger) as led:
        for sha in ("a" * 64, "b" * 64):
            plan = vs.Plan(f, T1, True, sha, "x.h5ad", {"task": "T1"}, vs.idempotency_key(T1, sha), [X001], "",
                           "abc")
            led.append(vs._row(plan, "queued", "old1", 200))
    with pytest.raises(vs.SubmitRefused, match="budget"):
        _official(f)


AMBIGUOUS_ENTRY = {T1: {"used": 0, "limit": 2}, "T1": {"used": 0, "limit": 2}}


@pytest.mark.parametrize("quota", [500, {"daily": {"T1": {"used": 0, "limit": 8}}}, {"official_per_board": 2},
                                   {"official_per_board": {T1: 1}}, {"official_per_board": {"T1:val": {"used": 0}}},
                                   {"official_per_board": AMBIGUOUS_ENTRY}])
def test_unknown_quota_refuses_official(env, quota):
    f = _checked(env, "a.h5ad", 0)
    env.api.quota = quota
    with pytest.raises(vs.SubmitRefused, match="quota") as e:
        _official(f)
    assert "budget" in e.value.plan.failed
    assert not env.api.posts(official=True)


def test_quota_shapes_read(env):
    board = load_boards(env.boards)[T1]
    assert vs.official_quota({"official_per_board": {T1: {"used": 1, "limit": 2}}}, board) == (1, 2)
    assert vs.official_quota({"official_per_board": {"T1": {"remaining": 2, "limit": 2}}}, board) == (0, 2)
    assert isinstance(vs.official_quota({"official_per_board": {T1: {"used": True, "limit": 2}}}, board), str)


def test_unreadable_submission_list_refuses_official(env):
    f = _checked(env, "a.h5ad", 0)
    env.api.subs_status = 500
    with pytest.raises(vs.SubmitRefused, match="submission list"):
        _official(f)


def test_retry_after_a_lost_response_is_identical_and_spends_nothing_twice(env):
    f = _checked(env, "a.h5ad", 0)
    env.api.script = ["store_503"]
    sleeps = []
    res = vs.submit(f, T1, official=True, confirm=vs.sha256_file(f)[:12], disclose=X001, sleep=sleeps.append)
    first, second = env.api.posts(official=True)
    assert sleeps == [2.0] and len(env.api.subs) == 1
    assert {k: first[k] for k in ("fields", "filename", "key", "sha256")} == \
        {k: second[k] for k in ("fields", "filename", "key", "sha256")}
    assert res.ok and res.server_id == "sub1"
    # the mock answers a changed payload under the same key with 409, so identity is what made the retry succeed
    changed = dataclasses.replace(res.plan, upload_name="other.h5ad")
    assert vs._post(env.api.url, TOKEN, changed, f, f.stat().st_size).status == 409


def test_a_retry_in_a_later_call_reuses_the_first_filename_fields_and_key(env):
    g = _checked(env, "g.h5ad", 5)
    env.api.script = ["503"] * 4
    res = _official(g)
    assert res.server_status == "http_error" and res.http_status == 503 and not res.ok
    renamed = env.tmp / "renamed.h5ad"
    shutil.copy(g, renamed)
    shutil.copy(vs.stamp_path(g), vs.stamp_path(renamed))  # the same bytes, so the same build stamp
    res = _official(renamed, model="a new name")
    posts = env.api.posts(official=True)
    assert res.ok and len(posts) == 5 and len(env.api.subs) == 1
    assert {p["filename"] for p in posts} == {"g.h5ad"} and len({p["key"] for p in posts}) == 1
    assert all(p["fields"] == {"task": "T1"} for p in posts)
    assert any("retry of an earlier official attempt" in w for w in res.plan.warnings)


def test_4xx_is_not_retried_automatically(env):
    f = _checked(env, "a.h5ad", 0)
    env.api.script = ["409"]
    res = _official(f)
    assert res.server_status == "http_error" and res.http_status == 409 and len(env.api.posts(official=True)) == 1
    assert _official(f).ok  # a manual retry after the cause is fixed goes through with the same request


def test_same_file_replays_the_earlier_record(env):
    f = _checked(env, "a.h5ad", 0)
    assert _official(f).server_id == "sub1"
    res = _official(f)
    assert res.replayed and res.server_id == "sub1" and res.server_status == "queued"
    assert len(env.api.posts(official=True)) == 1
    last = vs.read_ledger(env.ledger)[-1]
    assert "nothing was sent" in last["detail"] and last["mode"] == "replay"


def test_rejected_and_unstored_uploads_do_not_count(env):
    f1, f2, f3, f4 = (_checked(env, f"p{i}.h5ad", i) for i in range(4))
    env.api.script = ["reject"]
    assert _official(f1).server_status == "rejected"
    assert _official(f2).ok
    env.api.script = ["503"] * 4  # never stored, but the client cannot know: the outcome is ambiguous ...
    assert not _official(f3).ok
    with pytest.raises(vs.SubmitRefused, match="budget"):  # ... so it counts while the server does not list it
        _official(f4)
    vs.resolve(vs.sha256_file(f3)[:12], T1, "not-stored", "the portal lists p0 (rejected) and p1 only")
    assert _official(f4).ok  # a person's check settled it
    with pytest.raises(vs.SubmitRefused, match="budget"):
        _official(f3)


def test_dirty_tree_refuses_official_but_outputs_and_the_ledger_do_not(env, monkeypatch):
    f = _checked(env, "a.h5ad", 0)
    (env.repo / "method.py").write_text("x = 2\n")
    with pytest.raises(vs.SubmitRefused, match="clean tree") as e:
        _official(f)
    assert "method.py" in next(g.detail for g in e.value.plan.guards if g.name == "clean tree")
    assert not env.api.posts(official=True)

    (env.repo / "method.py").write_text("x = 1\n")
    (env.repo / "results").mkdir()
    (env.repo / "results" / "scratch.txt").write_text("untracked output\n")
    ledger = env.repo / "results" / "submissions.csv"  # a committed ledger that a format check then modifies
    monkeypatch.setattr(vs, "LEDGER", ledger)
    vs.submit(f, T1)
    _git(env.repo, "add", "results/submissions.csv")
    _git(env.repo, "commit", "-qm", "ledger")
    vs.submit(f, T1)
    assert _git(env.repo, "status", "--porcelain", "--untracked-files=no")  # the ledger is modified ...
    assert _official(f).ok  # ... and exempt


def test_other_guards(env):
    f = _checked(env, "a.h5ad", 0)
    with pytest.raises(vs.SubmitRefused, match="disclosure"):
        _official(f, disclose="")
    with pytest.raises(vs.SubmitRefused, match="disclosure"):
        vs.submit(f, T1, disclose="X001,X999")  # unknown ids refuse even a format check
    with pytest.raises(vs.SubmitRefused, match="server format check"):
        _official(_pred(env, "unchecked.h5ad", 9))
    with pytest.raises(vs.SubmitRefused, match="phase"):
        _official(f, board="T1:val")  # the server would route this upload to the test board
    env.api.phase = {**env.api.phase, "accepts_submissions": False}
    with pytest.raises(vs.SubmitRefused, match="phase"):
        _official(f)
    with pytest.raises(vs.SubmitRefused, match="board"):
        vs.submit(f, "T9:nope")
    assert not env.api.posts(official=True)


def test_board_fields_mirror_the_portal(env):
    def board(key, task, **kw):
        return vs.Board(key=key, task=task, genes=GENES, needs_coords=task != "T1", **kw)

    assert vs.board_fields(board("T1:val", "T1", split="val")) == {"task": "T1"}
    heart = board("T2:heart:val_extrap", "T2", setting="heart", split="val_extrap")
    assert vs.board_fields(heart) == {"task": "T2", "setting": "heart", "mode": "extrap"}
    assert vs.subtask_key(heart) == "T2/heart/extrap" and vs.board_split(heart) == "val"
    ko = board("T3:gata4", "T3", target_ko="gata4", split="val")
    assert vs.board_fields(ko) == {"task": "T3", "target": "gata4"} and vs.board_split(ko) == "val"
    for bad in (board("T3:gata4", "T3", target_ko="other"), board("T2:heart:test", "T2"), board("T4:x", "T4")):
        with pytest.raises(ValueError):
            vs.board_fields(bad)
    assert vs.board_split(board("T3:koA", "T3")) is None  # no published split: official uploads refuse

    f = _pred(env, "heart.h5ad", 3, coords=True)  # a board key without a mode is refused, not crashed on
    for kw in ({}, {"official": True, "confirm": vs.sha256_file(f)[:12], "disclose": X001}):
        with pytest.raises(vs.SubmitRefused, match="board") as e:
            vs.submit(f, "T2:heart:test", **kw)
        assert "cannot derive upload fields" in next(g.detail for g in e.value.plan.guards if g.name == "board")
    assert not env.api.requests


def test_local_check_refuses_before_sending(env):
    bad = _pred(env, "coords_missing.h5ad", 0)
    with pytest.raises(vs.SubmitRefused, match="local format check") as e:
        vs.submit(bad, "T2:heart:test_extrap")
    assert "spatial_3D" in next(g.detail for g in e.value.plan.guards if g.name == "local format check")
    renamed = env.tmp / "pred.txt"
    shutil.copy(bad, renamed)
    with pytest.raises(vs.SubmitRefused, match="h5ad"):
        vs.check_local(renamed, T1, env.boards)
    assert not env.api.posts()
    assert {r["server_status"] for r in vs.read_ledger(env.ledger)} == {"refused"}


def test_ledger_rows(env):
    f = _pred(env, "a.h5ad", 0)
    vs.submit(f, T1, disclose="X001,X002", note="first check")
    _official(f, disclose="x001", note="s1 for T1")
    rows = vs.read_ledger(env.ledger)
    assert list(rows[0]) == vs.LEDGER_FIELDS
    head = _git(env.repo, "rev-parse", "HEAD")
    fmt, sending, done = rows
    assert fmt["board"] == T1 and fmt["mode"] == "format" and fmt["server_status"] == "format_ok"
    assert fmt["file"].endswith("a.h5ad") and fmt["sha256"] == vs.sha256_file(f) and fmt["git_commit"] == head
    assert fmt["disclosure"] == "X001;X002" and fmt["note"] == "first check" and not fmt["idempotency_key"]
    assert sending["server_status"] == "sending" and sending["mode"] == "official"
    assert done["server_status"] == "queued" and done["server_id"] == "sub1" and done["http_status"] == "200"
    assert done["disclosure"] == "X001" and done["note"] == "s1 for T1" and done["git_commit"] == head
    assert json.loads(done["fields"]) == {"task": "T1"}
    assert done["idempotency_key"] == vs.idempotency_key(T1, fmt["sha256"])


def _output(res):
    text = res.output
    try:
        text += res.stderr
    except (ValueError, AttributeError):
        pass
    return text


def test_cli_and_token_never_in_output_or_ledger(env):
    runner = CliRunner()
    f = _pred(env, "p" * 200 + "/a.h5ad", 0)  # a long path, so that a wrapped line would break the hint below
    sha = vs.sha256_file(f)
    outs = []

    res = runner.invoke(cli.app, ["submit", str(f), "--board", T1, "--note", "cli check"])
    outs.append(_output(res))
    assert res.exit_code == 0, outs[-1]
    assert "Format check on T1:test (format_only=true)" in outs[-1] and "server format check passed" in outs[-1]
    assert f"vec submit {f} --board {T1} --official --confirm {sha[:12]} --disclose <ids from" in outs[-1]  # one line

    res = runner.invoke(cli.app, ["submit", str(f), "--board", T1, "--official"])
    outs.append(_output(res))
    assert res.exit_code == 1 and "FAIL confirm" in outs[-1] and "nothing was sent" in outs[-1]

    res = runner.invoke(cli.app, ["submit", str(f), "--board", T1, "--official", "--confirm", sha[:12],
                                  "--disclose", "X001,X002"])
    outs.append(_output(res))
    assert res.exit_code == 0, outs[-1]
    assert "OFFICIAL upload to T1:test" in outs[-1] and "uploaded submission sub1 (queued)" in outs[-1]

    res = runner.invoke(cli.app, ["submit", "status", "--board", T1])
    outs.append(_output(res))
    assert res.exit_code == 0 and "1/2" in outs[-1] and sha[:12] in outs[-1]

    env.api.script = ["echo_401"]  # a server error detail that echoes the token is scrubbed before printing
    res = runner.invoke(cli.app, ["submit", str(_pred(env, "b.h5ad", 1)), "--board", T1])
    outs.append(_output(res))
    assert res.exit_code == 1 and "<token>" in outs[-1]

    ledger = env.ledger.read_text()
    assert "<token>" in ledger
    for text in [*outs, ledger]:
        assert TOKEN not in text
    assert all(r["auth"] in (None, f"Bearer {TOKEN}") for r in env.api.requests)


def test_status_without_a_token_shows_the_ledger(env, monkeypatch):
    f = _checked(env, "a.h5ad", 0)
    monkeypatch.delenv("VEC_TOKEN")
    rep = vs.status(T1, token_file=env.tmp / "no-token")
    assert any(n.startswith("auth:") for n in rep.notes)
    (board,) = rep.boards
    assert board.format_ok == 1 and board.budget.server_used is None and board.budget.reason
    with pytest.raises(vs.vec_data.AuthError):
        vs.submit(f, T1, token_file=env.tmp / "no-token")


def test_boards_update_prints_only_the_summary(env, monkeypatch):
    calls = []

    def fake_update():
        calls.append(1)
        return {"boards": [T1, "T1:val"], "reduced": [T1], "t1_composition_written": False,
                "t1_composition_unreleased_stages": 1}

    monkeypatch.setattr("vec.boards.update_boards", fake_update)
    res = CliRunner().invoke(cli.app, ["boards", "update"])
    assert res.exit_code == 0 and calls == [1]
    assert "boards (2): T1:test, T1:val" in res.output and "submission spec (1): T1:test" in res.output
    assert "not written" in res.output and "anchor" not in res.output.lower()

    for phase, said in (("closed", "once the final phase is open"), ("unreadable", "could not be read")):
        monkeypatch.setattr("vec.boards.update_boards", lambda p=phase: {**fake_update(), "final_phase": p})
        res = CliRunner().invoke(cli.app, ["boards", "update"])  # why E10.5 is held back
        assert res.exit_code == 0 and said in " ".join(res.output.split())

    monkeypatch.setattr(cli, "load_boards", lambda: load_boards(env.boards))
    res = CliRunner().invoke(cli.app, ["boards"])  # unchanged: lists the boards
    assert res.exit_code == 0 and "T2:heart:test_extrap" in res.output


def test_token_only_goes_to_https_or_loopback(env):
    f = _pred(env, "a.h5ad", 0)
    with pytest.raises(vs.SubmitRefused, match="https"):
        vs.submit(f, T1, api_base="http://example.org")
    assert not env.api.requests


def test_token_echoed_in_a_successful_answer_is_scrubbed(env, capsys):
    f = _pred(env, "a.h5ad", 0)
    env.api.script = ["echo_ok"]
    res = vs.submit(f, T1)
    assert res.server_status == "format_failed" and TOKEN not in res.detail and "<token>" in res.detail
    assert TOKEN not in env.ledger.read_text()


def test_lost_response_counts_against_the_budget_while_the_server_list_lags(env):
    a = _checked(env, "a.h5ad", 0)
    assert _official(a).ok
    lost = "b" * 64  # sent, response lost, and not listed by the server yet
    with vs.Ledger(env.ledger) as led:
        plan = vs.Plan(a, T1, True, lost, "b.h5ad", {"task": "T1"}, vs.idempotency_key(T1, lost), [X001], "", "abc")
        led.append(vs._row(plan, "no_response", "", None, "attempt 4 of 4: timed out"))
    c = _checked(env, "c.h5ad", 2)
    with pytest.raises(vs.SubmitRefused, match="budget"):
        _official(c)
    assert len(env.api.posts(official=True)) == 1


def test_untracked_files_block_official_except_outputs_and_ignored_files(env):
    f = _checked(env, "a.h5ad", 0)
    (env.repo / ".gitignore").write_text("data/\n")
    _git(env.repo, "add", ".gitignore")
    _git(env.repo, "commit", "-qm", "ignore data")
    blocking = ["scripts/new_step.py", "run.py", "notes.txt", "experiments/t9/conf.yml", "experiments/t9/nb.ipynb",
                "src/new method.py", "experiments/t9/results.py"]
    allowed = ["results/x.csv", "runs/r1/log.txt", "experiments/t9/results/out.json", "experiments/results/a.csv",
               "data/loader.py"]
    for p in blocking + allowed:
        (env.repo / p).parent.mkdir(parents=True, exist_ok=True)
        (env.repo / p).write_text("y = 1\n")
    assert sorted(vs.git_state(env.repo)[1]) == sorted(blocking)  # unquoted, as git status -z gives them
    with pytest.raises(vs.SubmitRefused, match="clean tree") as e:
        _official(f)
    assert "notes.txt" in next(g.detail for g in e.value.plan.guards if g.name == "clean tree")
    for p in blocking:
        (env.repo / p).unlink()
    assert _official(f).ok
    assert not [p for p in allowed if not (env.repo / p).exists()]


def test_paths_with_spaces_compare_exactly(env, monkeypatch):
    ledger = env.repo / "results" / "sub missions.csv"
    monkeypatch.setattr(vs, "LEDGER", ledger)
    f = _checked(env, "a.h5ad", 0)
    _git(env.repo, "add", "results/sub missions.csv")
    (env.repo / "my method.py").write_text("x = 1\n")
    _git(env.repo, "add", "my method.py")
    _git(env.repo, "commit", "-qm", "names with spaces")
    vs.submit(f, T1)  # the committed ledger is modified again ...
    assert vs.git_state(env.repo, ignore=[ledger])[1] == []  # ... and stays exempt
    (env.repo / "my method.py").write_text("x = 2\n")
    assert vs.git_state(env.repo, ignore=[ledger])[1] == ["my method.py"]


# ---------------------------------------------------------------------------------------------------------------------
# Build stamps: the code that built the file, not the code checked out when it is sent


def _guard(plan, name):
    return next(g for g in plan.guards if g.name == name)


def _refused_for_the_stamp(path):
    """An official upload of `path` refused by the guard "build stamp" alone; returns its detail."""
    with pytest.raises(vs.SubmitRefused, match="build stamp") as e:
        _official(path)
    assert e.value.plan.failed == ["build stamp"]
    return _guard(e.value.plan, "build stamp").detail


def test_an_official_upload_needs_the_files_build_stamp(env):
    f = _checked(env, "a.h5ad", 0)
    vs.stamp_path(f).unlink()
    res = vs.submit(f, T1)  # a format check never requires one
    assert res.server_status == "format_ok" and vs.read_ledger(env.ledger)[-1]["build_commit"] == ""
    assert res.plan.warnings[-1].startswith("an official upload of this file would fail build stamp: no build stamp "
                                            "next to the file (a.h5ad.stamp.json)")
    assert "no build stamp next to the file" in _refused_for_the_stamp(f)
    assert not env.api.posts(official=True) and vs.read_ledger(env.ledger)[-1]["detail"] == "failed: build stamp"


def test_a_file_rebuilt_without_a_new_stamp_is_refused(env):
    f = _checked(env, "a.h5ad", 0)
    old = json.loads(vs.stamp_path(f).read_text())
    shutil.copyfile(_bigger_pred(env, "rebuilt.h5ad", 7), f)  # new bytes at the same path; the old stamp stays
    sha = vs.sha256_file(f)
    res = vs.submit(f, T1)
    assert res.server_status == "format_ok" and vs.read_ledger(env.ledger)[-1]["build_commit"] == ""
    assert f"the build stamp is for other bytes (sha256 {old['sha256'][:12]}, the file's is {sha[:12]})" in \
        res.plan.warnings[-1]
    assert "rebuilt without a new stamp" in _refused_for_the_stamp(f)
    assert not env.api.posts(official=True)


def test_a_file_built_from_a_dirty_tree_is_refused_even_when_the_tree_is_clean_again(env):
    head = _git(env.repo, "rev-parse", "HEAD")
    (env.repo / "method.py").write_text("x = 2\n")
    f = _pred(env, "a.h5ad", 0)  # built with an uncommitted change
    (env.repo / "method.py").write_text("x = 1\n")  # undone before the upload: the clean tree guard passes
    assert vs.submit(f, T1).server_status == "format_ok"
    assert vs.read_ledger(env.ledger)[-1]["build_commit"] == f"{head}+dirty"
    assert "not clean (uncommitted or untracked: method.py)" in _refused_for_the_stamp(f)
    assert vs.read_ledger(env.ledger)[-1]["build_commit"] == f"{head}+dirty"  # the refusal row names it too
    stamp = json.loads(vs.stamp_path(f).read_text())
    vs.stamp_path(f).write_text(json.dumps({**stamp, "why": None}))  # a stamp that gives no reason: told from its paths
    assert "not clean (uncommitted or untracked: method.py)" in _refused_for_the_stamp(f)
    assert not env.api.posts(official=True)


def test_a_file_built_at_a_commit_that_is_not_an_ancestor_of_head_is_refused(env):
    main = _git(env.repo, "rev-parse", "--abbrev-ref", "HEAD")
    _git(env.repo, "checkout", "-qb", "side")
    (env.repo / "method.py").write_text("x = 3\n")
    _git(env.repo, "commit", "-qam", "an experiment on a branch that is not merged")
    side = _git(env.repo, "rev-parse", "HEAD")
    f = _checked(env, "a.h5ad", 0)  # built clean at the branch's commit
    _git(env.repo, "checkout", "-q", main)
    assert vs.read_ledger(env.ledger)[-1]["build_commit"] == side
    assert f"build commit {side[:12]} is neither the local HEAD nor an ancestor of it" in _refused_for_the_stamp(f)
    stamp = json.loads(vs.stamp_path(f).read_text())
    vs.stamp_path(f).write_text(json.dumps({**stamp, "commit": "e" * 40}))  # a commit this repository does not have
    assert "build commit eeeeeeeeeeee is not in the local repository" in _refused_for_the_stamp(f)
    vs.stamp_path(f).write_text(json.dumps(stamp))
    assert not env.api.posts(official=True)
    _git(env.repo, "merge", "-q", "side")  # merged: the build commit is now HEAD's ancestor (here, HEAD itself)
    assert _official(f).ok


def test_a_file_built_clean_at_an_older_commit_is_sent_and_recorded_under_it(env):
    built_at = _git(env.repo, "rev-parse", "HEAD")
    f = _checked(env, "a.h5ad", 0)
    (env.repo / "method.py").write_text("x = 2\n")
    _git(env.repo, "commit", "-qam", "a later fix, after the file was built")
    head = _git(env.repo, "rev-parse", "HEAD")
    res = _official(f)
    assert res.ok and _guard(res.plan, "build stamp").detail == \
        f"built from a clean tree at {built_at[:12]}, the local HEAD or an ancestor of it"
    row = vs.read_ledger(env.ledger)[-1]
    assert (row["mode"], row["server_status"], row["git_commit"], row["build_commit"]) == \
        ("official", "queued", head, built_at)
    link = env.tmp / "links" / "s1.h5ad"  # a link: the stamp is read next to the file it points at
    link.parent.mkdir()
    link.symlink_to(f)
    res = _official(link)
    assert res.replayed and res.plan.build_commit == built_at and _guard(res.plan, "build stamp").ok


def test_the_files_own_stamp_is_exempt_from_the_clean_tree_rule_but_no_other_stamp(env):
    (env.repo / ".gitignore").write_text("*.h5ad\n")  # as in this repository: prediction files are never tracked
    _git(env.repo, "add", ".gitignore")
    _git(env.repo, "commit", "-qm", "ignore predictions")
    preds = env.repo / "preds"  # inside the tree, outside results/ and runs/
    f = _pred(env, "s1.h5ad", 0, where=preds)
    assert vs.submit(f, T1).server_status == "format_ok"
    other = _pred(env, "s2.h5ad", 2)
    shutil.copy(vs.stamp_path(other), preds / ("s2.h5ad" + vs.STAMP_SUFFIX))  # another file's stamp lies there too
    with pytest.raises(vs.SubmitRefused, match="clean tree") as e:
        _official(f)
    assert e.value.plan.failed == ["clean tree"] and _guard(e.value.plan, "clean tree").detail == \
        "uncommitted or untracked: preds/s2.h5ad.stamp.json"
    (preds / ("s2.h5ad" + vs.STAMP_SUFFIX)).unlink()
    assert _official(f).ok  # its own stamp, untracked next to it, is exempt like the ledger


def test_cli_shows_the_build_commit_in_the_plan_and_in_status(env):
    runner = CliRunner()
    built_at = _git(env.repo, "rev-parse", "HEAD")
    f = _checked(env, "a.h5ad", 0)
    sha = vs.sha256_file(f)
    (env.repo / "method.py").write_text("x = 2\n")
    _git(env.repo, "commit", "-qam", "a later fix")
    head = _git(env.repo, "rev-parse", "HEAD")
    res = runner.invoke(cli.app, ["submit", str(f), "--board", T1, "--official", "--confirm", sha[:12], "--disclose",
                                  X001])
    out = _output(res)
    assert res.exit_code == 0, out
    assert f"  commit    {head}\n  build     {built_at}\n" in out
    assert re.search(rf"PASS build stamp +built from a clean tree at {built_at[:12]}", out)
    res = runner.invoke(cli.app, ["submit", "status", "--board", T1])
    assert res.exit_code == 0 and f"official {T1}: {sha[:12]} sub1 (queued), commit {head}, build {built_at}" in \
        _output(res)
    vs.stamp_path(f).unlink()
    res = runner.invoke(cli.app, ["submit", str(f), "--board", T1])  # a format check without a stamp: shown, not needed
    assert res.exit_code == 0 and "  build     -\n" in _output(res) and "would fail build stamp" in _output(res)


# A builder in a process of its own, as a build on the DGX runs: it imports vec.io (whose import of vec.provenance
# records the repository's state as the build's start) and the method, says it is ready, waits, then writes the file
BUILDER = r"""
import sys
import time
from pathlib import Path

import vec.boards

repo, ready, go, out = map(Path, sys.argv[1:5])
vec.boards.REPO_ROOT = repo  # the test's repository stands in for this package's; vec.provenance is not imported yet
import numpy as np
from vec.io import make_submission, write_submission  # the build starts here

sys.path.insert(0, str(repo))
import method  # the code that builds the file, as the tree holds it now

ready.touch()
while not go.exists():
    time.sleep(0.01)
genes = sys.argv[5].split(",")
counts = np.random.default_rng(0).poisson(2.0, size=(40, len(genes)))
a = make_submission(np.log1p(1e4 * counts / counts.sum(1, keepdims=True).clip(min=1)), genes)
a.uns["method_x"] = method.x  # what the code that ran computed
write_submission(a, out)
"""


def _build_in_a_new_process(env, name, during=lambda: None):
    """Build `name` in env.tmp with BUILDER in a new process, which stamps it against env.repo; `during` runs once the
    build has started and before it writes the file. Returns the file."""
    ready, go, err, out = (env.tmp / f"{name}{ext}" for ext in (".ready", ".go", ".err", ""))
    environ = {**{k: v for k, v in os.environ.items() if "TOKEN" not in k.upper()},
               "PYTHONDONTWRITEBYTECODE": "1"}  # no __pycache__ in the test's repository, which ignores nothing
    with err.open("wb") as stderr:
        proc = subprocess.Popen([sys.executable, "-c", BUILDER, str(env.repo), str(ready), str(go), str(out),
                                 ",".join(GENES)], cwd=env.tmp, env=environ, stderr=stderr)
        try:
            deadline = time.monotonic() + 120
            while not ready.exists():
                assert proc.poll() is None and time.monotonic() < deadline, err.read_text()
                time.sleep(0.01)
            during()
            go.touch()
            proc.wait(timeout=120)
        finally:
            if proc.poll() is None:
                proc.kill()
                proc.wait()
    assert proc.returncode == 0, err.read_text()
    return out


def test_a_file_whose_build_saw_a_new_commit_land_is_never_sent_officially(env):
    """A build starts at A; a fix is committed (a pull, a sync) before it writes the file. A's code built the file, so
    its stamp must not say that it was built clean at the new HEAD."""
    a = _git(env.repo, "rev-parse", "HEAD")

    def pull():
        (env.repo / "method.py").write_text("x = 2\n")
        _git(env.repo, "commit", "-qam", "a fix, synced while the build ran")

    f = _build_in_a_new_process(env, "s1.h5ad", during=pull)
    b = _git(env.repo, "rev-parse", "HEAD")
    s = vs.read_stamp(vs.stamp_path(f))
    assert ad.read_h5ad(f).uns["method_x"] == 1 and not s["clean"], s  # built by A's code
    assert (s["commit"], s["start_commit"], s["why"]) == (b, a, f"HEAD moved during the build, from {a[:12]} to "
                                                                f"{b[:12]}")
    res = vs.submit(f, T1)  # a format check records the build commit as not clean, and warns
    assert res.server_status == "format_ok" and vs.read_ledger(env.ledger)[-1]["build_commit"] == f"{b}+dirty"
    assert "would fail build stamp: the build was not clean (HEAD moved during the build" in res.plan.warnings[-1]
    assert f"(HEAD moved during the build, from {a[:12]} to {b[:12]})" in _refused_for_the_stamp(f)
    assert not env.api.posts(official=True)


@pytest.mark.parametrize("change", ["nothing", "stash", "checkout and back"])
def test_only_a_build_that_ran_on_one_clean_commit_is_clean(env, change):
    head, branch = _git(env.repo, "rev-parse", "HEAD"), _git(env.repo, "rev-parse", "--abbrev-ref", "HEAD")
    during = lambda: None  # noqa: E731
    if change == "stash":  # the build starts with an uncommitted change, stashed before it writes the file
        (env.repo / "method.py").write_text("x = 3\n")
        during = lambda: _git(env.repo, "stash", "-q")  # noqa: E731
    elif change == "checkout and back":  # HEAD leaves while the build runs, and comes back before it writes the file
        _git(env.repo, "checkout", "-qb", "side")
        (env.repo / "method.py").write_text("x = 4\n")
        _git(env.repo, "commit", "-qam", "an experiment")
        _git(env.repo, "checkout", "-q", branch)

        def during():
            _git(env.repo, "checkout", "-q", "side")
            _git(env.repo, "checkout", "-q", branch)

    f = _build_in_a_new_process(env, "s1.h5ad", during)
    s = vs.read_stamp(vs.stamp_path(f))
    assert s["clean"] is (change == "nothing"), s
    assert (s["commit"], s["start_commit"]) == (head, head) and _git(env.repo, "status", "--porcelain") == ""
    if change == "nothing":
        assert s["why"] is None and vs.submit(f, T1).server_status == "format_ok" and _official(f).ok
        assert vs.read_ledger(env.ledger)[-1]["build_commit"] == head
        return
    if change == "stash":
        assert ad.read_h5ad(f).uns["method_x"] == 3  # built by the uncommitted change
        assert (s["why"], s["unclean"], s["unclean_count"]) == (
            "uncommitted or untracked when the build started: method.py", ["method.py"], 1)
    else:
        assert s["why"].startswith(f"HEAD moved during the build and back to {head[:12]} (checkout: moving from side")
    assert vs.submit(f, T1).server_status == "format_ok"
    assert vs.read_ledger(env.ledger)[-1]["build_commit"] == f"{head}+dirty"
    assert f"the build was not clean ({s['why']})" in _refused_for_the_stamp(f)
    assert not env.api.posts(official=True)


# ---------------------------------------------------------------------------------------------------------------------
# What is sent: a private snapshot of the confirmed bytes


def _bigger_pred(env, name, seed):
    """Another valid prediction with more cells, so that its size differs too."""
    rng = np.random.default_rng(seed)
    counts = rng.poisson(2.0, size=(55, len(GENES)))
    X = np.log1p(1e4 * counts / counts.sum(1, keepdims=True).clip(min=1))
    return build(make_submission(X, GENES), env.tmp / name, env.repo)


def test_a_retry_sends_the_confirmed_bytes_even_if_the_file_is_rewritten(env):
    a = _checked(env, "a.h5ad", 0)
    sha = vs.sha256_file(a)
    other = _bigger_pred(env, "other.h5ad", 7)
    env.api.script = ["503"]  # attempt 1: nothing stored

    def rewrite(_seconds):  # an scp or a pipeline rewrites the prediction in place during the back-off
        shutil.copyfile(other, a)

    res = vs.submit(a, T1, official=True, confirm=sha[:12], disclose=X001, sleep=rewrite)
    assert vs.sha256_file(a) == vs.sha256_file(other)  # the file on disk did change
    assert [p["sha256"] for p in env.api.posts(official=True)] == [sha, sha]
    assert res.ok and [s["sha256"] for s in env.api.subs] == [sha]


def test_a_change_after_the_plan_is_refused_and_the_snapshot_is_private_and_removed(env, monkeypatch):
    tmp = env.tmp / "tmpdir"
    tmp.mkdir()
    monkeypatch.setattr(tempfile, "tempdir", str(tmp))
    a, b = _checked(env, "a.h5ad", 0), _checked(env, "b.h5ad", 1)
    other = _pred(env, "other.h5ad", 7)
    with pytest.raises(vs.SubmitRefused, match="file changed before upload"):
        _official(a, on_plan=lambda plan: shutil.copyfile(other, a))
    assert not env.api.posts(official=True) and vs.read_ledger(env.ledger)[-1]["detail"] == "file changed before upload"
    seen = []

    def look(_seconds):
        (snap_dir,) = tmp.iterdir()
        (snap,) = snap_dir.iterdir()
        seen.append((snap_dir.stat().st_mode & 0o777, snap.stat().st_mode & 0o777, vs.sha256_file(snap)))

    env.api.script = ["503"]
    assert vs.submit(b, T1, official=True, confirm=vs.sha256_file(b)[:12], disclose=X001, sleep=look).ok
    assert seen == [(0o700, 0o600, vs.sha256_file(b))] and not list(tmp.iterdir())


# ---------------------------------------------------------------------------------------------------------------------
# The token: redacted once, in _call, from everything the server sends back; refused if malformed


def test_the_token_is_redacted_from_ids_replays_guard_details_and_status(env):
    runner = CliRunner()
    f = _checked(env, "a.h5ad", 0)
    sha = vs.sha256_file(f)
    official = ["submit", str(f), "--board", T1, "--official", "--confirm", sha[:12], "--disclose", X001]
    env.api.script = ["echo_id"]  # the server's id echoes the token
    outs = [_output(runner.invoke(cli.app, official))]
    assert "uploaded submission sub-<token> (queued)" in outs[-1]
    env.api.echo_status = True  # so does the replay's status refresh
    outs.append(_output(runner.invoke(cli.app, official)))
    assert "already sent: submission sub-<token> (scored Bearer <token>)" in outs[-1]
    env.api.subs.append({"id": f"web-{TOKEN}", "task": "T1", "bundle_key": T1, "status": f"scored {TOKEN}",
                         "sha256": "1" * 64})
    env.api.quota = {"official_per_board": {f"T1:test {TOKEN}": {"used": 0, "limit": 2}}}  # an unknown shape
    with pytest.raises(vs.SubmitRefused, match="quota") as e:
        _official(_checked(env, "b.h5ad", 1))
    details = " ".join(g.detail for g in e.value.plan.guards)
    assert "keys: T1:test <token>" in details
    outs.append(_output(runner.invoke(cli.app, ["submit", "status"])))
    assert "web-<token>" in outs[-1] and "T1:test <token>" in outs[-1]
    for text in [*outs, details, env.ledger.read_text()]:
        assert TOKEN not in text


def test_a_token_across_the_cut_leaves_no_prefix(env):
    env.api.script = ["echo_401_long"]
    res = CliRunner().invoke(cli.app, ["submit", str(_pred(env, "a.h5ad", 0)), "--board", T1])
    assert res.exit_code == 1 and "<token>" in env.ledger.read_text()
    for text in (_output(res).replace("\n", ""), env.ledger.read_text()):
        assert TOKEN[:8] not in text


def test_a_malformed_token_is_refused_before_anything_is_sent(env, monkeypatch):
    broken = "tok-FAKE-part1\npart2-never-print"  # e.g. after `echo NEW >> ~/.config/vec/token`
    monkeypatch.setenv("VEC_TOKEN", broken)
    f = _pred(env, "a.h5ad", 0)
    with pytest.raises(vs.vec_data.AuthError):
        vs.submit(f, T1)
    res = CliRunner().invoke(cli.app, ["submit", str(f), "--board", T1])
    assert res.exit_code == 2 and "auth error" in _output(res)
    assert "part1" not in _output(res) and "part2-never-print" not in _output(res) and not env.api.requests
    reply = vs._get(env.api.url, "/challenge/integration/quota", broken)  # http.client quotes it in a ValueError
    assert reply.status == 0 and "ValueError" in reply.detail
    assert "part1" not in reply.detail and "part2-never-print" not in reply.detail


def test_server_text_cannot_inject_markup(env):
    f = _checked(env, "a.h5ad", 0)
    sha = vs.sha256_file(f)
    env.api.script = ["markup"]  # stored, with the status "queued [/bold]"
    res = CliRunner().invoke(cli.app, ["submit", str(f), "--board", T1, "--official", "--confirm", sha[:12],
                                       "--disclose", X001])
    assert res.exit_code == 0, _output(res)
    assert "uploaded submission sub1 (queued [/bold])" in _output(res)
    env.api.subs.append({"id": "[red]web1", "task": "T1", "bundle_key": T1, "status": "[/]", "sha256": "1" * 64})
    res = CliRunner().invoke(cli.app, ["submit", "status", "--board", T1])
    assert res.exit_code == 0, _output(res)
    assert "[red]web1 [/]" in _output(res)


# ---------------------------------------------------------------------------------------------------------------------
# Budget: ambiguous until settled, resolutions, format checks the server stored, slots by server id


def _hide(env, path):
    """Take the stored submission of `path` off the server's list and quota (they lag); returns it."""
    sha = vs.sha256_file(path)
    hidden = [s for s in env.api.subs if s["sha256"] == sha]
    env.api.subs[:] = [s for s in env.api.subs if s["sha256"] != sha]
    return hidden


def _row(sha, status, http="", sid="", mode="official", when="2000-01-01T00:00:00+00:00"):
    return {**dict.fromkeys(vs.LEDGER_FIELDS, ""), "utc_time": when, "board": T1, "mode": mode, "sha256": sha,
            "server_status": status, "http_status": http, "server_id": sid}


def test_a_4xx_after_an_ambiguous_attempt_does_not_free_the_slot(env):
    a, b, c = (_checked(env, f"p{i}.h5ad", i) for i in range(3))
    assert _official(a).ok
    env.api.script = ["store_503", "409"]  # attempt 1 stored, its answer lost; the retry answered 409
    assert _official(b).http_status == 409
    hidden = _hide(env, b)  # neither the list nor the quota shows b yet
    with pytest.raises(vs.SubmitRefused, match="budget"):
        _official(c)
    env.api.subs += hidden
    assert len(env.api.posts(official=True)) == 3  # a, and b's two attempts


@pytest.mark.parametrize("attempts", [
    [("sending", ""), ("http_error", "307")],
    [("sending", ""), ("http_error", "503"), ("sending", ""), ("http_error", "429")],
    [("sending", ""), ("no_response", ""), ("sending", ""), ("http_error", "409")],
    [("sending", ""), ("unconfirmed", "200")],
    [("sending", ""), ("queued", "200")],  # a 2xx with a status but no id
    [("sending", "")],  # a crash mid-upload
])
def test_an_ambiguous_upload_counts_until_listed_or_resolved_however_old(attempts):
    sha, ident = "b" * 64, {"task": "T1"}
    rows = [_row(sha, s, http) for s, http in attempts]
    budget = vs.board_budget(T1, rows, (0, 2), [], ident)
    assert budget.spent == 1 and budget.ambiguous == {sha}
    listed = [{"id": "sub7", "bundle_key": T1, "status": "rejected", "sha256": sha}]
    assert vs.board_budget(T1, rows, (0, 2), listed, ident).spent == 0  # a listing settles it
    resolved = [*rows, _row(sha, "not-stored", mode="resolve")]
    assert vs.board_budget(T1, resolved, (0, 2), [], ident).spent == 0  # so does a person
    assert vs.board_budget(T1, resolved, (0, 2), [{**listed[0], "status": "queued"}], ident).spent == 1  # listing wins
    assert vs.board_budget(T1, [*resolved, *rows], (0, 2), [], ident).spent == 1  # a new attempt counts again


def test_resolve_records_a_decision_locally(env, monkeypatch):
    a, b, c = (_checked(env, f"p{i}.h5ad", i) for i in range(3))
    assert _official(a).ok
    env.api.script = ["503"] * 4  # nothing stored, but the client cannot know
    assert not _official(b).ok
    with pytest.raises(vs.SubmitRefused, match="budget"):
        _official(c)
    sha_a, sha_b = vs.sha256_file(a), vs.sha256_file(b)
    monkeypatch.delenv("VEC_TOKEN")
    n, runner = len(env.api.requests), CliRunner()
    res = runner.invoke(cli.app, ["submit", "resolve", sha_b[:12], "--board", T1, "--as", "not-stored", "--note",
                                  "portal lists p0 only"])
    assert res.exit_code == 0, _output(res)
    assert len(env.api.requests) == n  # local: nothing sent, no token read
    row = vs.read_ledger(env.ledger)[-1]
    assert (row["mode"], row["sha256"], row["server_status"], row["note"]) == ("resolve", sha_b, "not-stored",
                                                                              "portal lists p0 only")
    for sha, decision in ((sha_b, "not-stored"), (sha_b, "acknowledged"), (sha_a, "not-stored"), ("f" * 12, "x")):
        res = runner.invoke(cli.app, ["submit", "resolve", sha[:12], "--board", T1, "--as", decision, "--note", "x"])
        assert res.exit_code == 1 and "nothing was recorded" in _output(res), (sha, decision)
    assert vs.read_ledger(env.ledger)[-1] == row
    monkeypatch.setenv("VEC_TOKEN", TOKEN)
    assert _official(c).ok  # b no longer counts


def test_a_format_check_stored_as_a_submission_counts_and_blocks_the_board(env):
    a, b = _pred(env, "a.h5ad", 0), _pred(env, "b.h5ad", 1)
    env.api.script = ["format_stored"]  # the server ignored format_only
    res = vs.submit(a, T1)
    assert res.server_status == "UNEXPECTED_SUBMISSION" and not res.ok and len(env.api.subs) == 1
    lagging = vs.board_budget(T1, vs.read_ledger(env.ledger), (0, 2), [], {"task": "T1"})
    assert lagging.spent == 1
    n = len(env.api.posts())
    for kw in ({}, {"official": True, "confirm": vs.sha256_file(b)[:12], "disclose": X001}):
        with pytest.raises(vs.SubmitRefused, match="unresolved upload"):
            vs.submit(b, T1, **kw)
    assert len(env.api.posts()) == n  # no POST of any kind until a person has looked
    with pytest.raises(vs.SubmitRefused, match="acknowledged"):  # the server gave it an id: it is stored
        vs.resolve(vs.sha256_file(a)[:12], T1, "not-stored", "portal")
    vs.resolve(vs.sha256_file(a)[:12], T1, "acknowledged", "portal shows the stray submission")
    assert vs.submit(b, T1).server_status == "format_ok"
    assert vs.board_budget(T1, vs.read_ledger(env.ledger), (0, 2), [], {"task": "T1"}).spent == 1  # still counted


def test_a_format_check_answered_format_unconfirmed_counts_as_ambiguous_until_resolved(env):
    a, b = _pred(env, "a.h5ad", 0), _pred(env, "b.h5ad", 1)
    env.api.script = ["format_unconfirmed"]
    assert vs.submit(a, T1).server_status == "format_unconfirmed"
    sha = vs.sha256_file(a)
    lagging = vs.board_budget(T1, vs.read_ledger(env.ledger), (0, 2), [], {"task": "T1"})
    assert lagging.spent == 1 and lagging.ambiguous == {sha}
    with pytest.raises(vs.SubmitRefused, match="unresolved upload"):
        vs.submit(b, T1)
    vs.resolve(sha[:12], T1, "not-stored", "portal shows no new submission")
    assert vs.board_budget(T1, vs.read_ledger(env.ledger), (0, 2), [], {"task": "T1"}).spent == 0
    assert vs.submit(b, T1).server_status == "format_ok"


def test_an_identical_retry_of_an_unanswered_upload_is_not_counted_twice(env):
    a, b, c = (_checked(env, f"p{i}.h5ad", i) for i in range(3))
    assert _official(a).ok
    env.api.script = ["503"] * 4  # b: nothing stored, and the client cannot know
    assert not _official(b).ok
    res = _official(b)  # the identical request under the same key: at most b's own slot
    assert res.ok and not res.replayed and len(env.api.subs) == 2
    with pytest.raises(vs.SubmitRefused, match="budget"):
        _official(c)


def test_an_official_2xx_without_an_id_is_ambiguous_and_never_replayed_from_the_ledger(env):
    f, g, h = (_checked(env, f"p{i}.h5ad", i) for i in range(3))
    env.api.script = ["store_unreadable"]  # stored, but the answer is an intermediary's HTML page
    first = _official(f)
    assert first.server_status == "unconfirmed" and not first.ok
    hidden = _hide(env, f)
    again = _official(f)  # not listed yet: the identical request goes again, and the key returns the original
    env.api.subs += hidden
    assert not again.replayed and again.ok and again.server_id == "sub1" and len(env.api.subs) == 1
    assert _official(g).ok
    with pytest.raises(vs.SubmitRefused, match="budget"):
        _official(h)


def test_server_entries_count_by_id():
    copies = [{"id": f"web{i}", "bundle_key": T1, "status": "queued", "sha256": "d" * 64} for i in range(2)]
    assert vs.board_budget(T1, [], (1, 2), copies, {"task": "T1"}).spent == 2  # two stored copies of one file


def test_quota_entries_that_disagree_or_carry_other_counts_are_unknown(env):
    board = load_boards(env.boards)[T1]
    for entry in ({"used": 1, "limit": 2, "remaining": 0}, {"used": 0, "limit": 2, "pending": 2},
                  {"used": 0, "limit": 2, "final": True}, {"used": 0, "limit": 2, "by_status": {"queued": 1}}):
        assert isinstance(vs.official_quota({"official_per_board": {T1: entry}}, board), str), entry
    assert vs.official_quota({"official_per_board": {T1: {"used": 1, "limit": 2, "remaining": 1}}}, board) == (1, 2)
    assert vs.official_quota({"official_per_board": {T1: {"used": 1, "limit": 2, "board": T1}}}, board) == (1, 2)


# ---------------------------------------------------------------------------------------------------------------------
# Strays (format checks the server stored): where they count, what they block; slots settled one by one; no replay
# of a format check's submission


def test_a_stray_counts_and_blocks_on_the_board_the_server_named_and_on_the_requested_one(env):
    env.api.phase = {**env.api.phase, "split": "val"}  # so a format check on T1:val goes out
    a, b = _pred(env, "a.h5ad", 0), _pred(env, "b.h5ad", 1)
    env.api.script = ["format_stored"]  # stored anyway, and the mock files every T1 upload under T1:test
    res = vs.submit(a, "T1:val")
    assert res.server_status == "UNEXPECTED_SUBMISSION" and env.api.subs[0]["bundle_key"] == T1
    rows, sha = vs.read_ledger(env.ledger), vs.sha256_file(a)
    assert [(r["board"], r["server_id"]) for r in rows if r["server_status"] == "UNEXPECTED_SUBMISSION"] == \
        [("T1:val", "sub1"), (T1, "sub1")]
    out = _output(CliRunner().invoke(cli.app, ["submit", "status"]))
    for board in ("T1:val", T1):  # the list lags: the ledger alone counts it, on both boards
        assert vs.board_budget(board, rows, (0, 2), [], {"task": "T1"}).spent == 1 and vs._unresolved(rows, board)
        assert f"vec submit resolve {sha[:12]} --board {board} --as acknowledged" in res.detail
        assert (f"unresolved {board}: a format check of {sha[:12]} came back as UNEXPECTED_SUBMISSION; nothing is "
                f"sent to any board until a person has checked the portal and run: vec submit resolve {sha[:12]} "
                f"--board {board} --as acknowledged --note TEXT") in out
    vs.resolve(sha[:12], T1, "acknowledged", "the portal lists it on T1:test")
    with pytest.raises(vs.SubmitRefused, match="unresolved upload"):  # T1:val still waits for a person
        vs.submit(b, "T1:val")
    vs.resolve(sha[:12], "T1:val", "acknowledged", "the same submission, recorded under the board the request named")
    assert vs.submit(b, "T1:val").server_status == "format_ok"


def test_format_checks_on_a_board_outside_the_phase_split_are_refused(env):
    f = _pred(env, "a.h5ad", 0)
    with pytest.raises(vs.SubmitRefused, match="phase") as e:  # the test phase: a stray would land on T1:test
        vs.submit(f, "T1:val")
    assert "split 'test'" in next(g.detail for g in e.value.plan.guards if g.name == "phase")
    env.api.phase = "down"  # answered, but not with an object
    with pytest.raises(vs.SubmitRefused, match="phase"):
        vs.submit(f, T1)
    assert not env.api.posts() and {r["server_status"] for r in vs.read_ledger(env.ledger)} == {"refused"}
    env.api.phase = {**MockApi().phase, "accepts_submissions": False}  # a format check only needs the split
    assert vs.submit(f, T1).server_status == "format_ok"


def test_a_stray_on_one_board_blocks_every_post_to_every_board(env):
    a, heart = _pred(env, "a.h5ad", 0), _pred(env, "heart.h5ad", 1, coords=True)
    env.api.script = ["format_stored"]
    assert vs.submit(a, T1).server_status == "UNEXPECTED_SUBMISSION"
    n, sha = len(env.api.posts()), vs.sha256_file(a)
    for kw in ({}, {"official": True, "confirm": vs.sha256_file(heart)[:12], "disclose": X001}):
        with pytest.raises(vs.SubmitRefused, match="unresolved upload") as e:
            vs.submit(heart, "T2:heart:test_extrap", **kw)
        detail = next(g.detail for g in e.value.plan.guards if g.name == "unresolved upload")
        assert f"{sha[:12]} on {T1} (UNEXPECTED_SUBMISSION)" in detail and "nothing is sent to any board" in detail
        assert f"vec submit resolve {sha[:12]} --board {T1} --as acknowledged --note TEXT" in detail
    assert len(env.api.posts()) == n
    vs.resolve(sha[:12], T1, "acknowledged", "the portal shows the stray on T1:test")
    assert vs.submit(heart, "T2:heart:test_extrap").server_status == "format_ok"


def test_each_listed_submission_settles_one_ledger_slot_not_the_whole_file():
    sha, ident = "b" * 64, {"task": "T1"}
    rows = [_row(sha, "sending"), _row(sha, "no_response"),  # stored as sub1, the answer lost
            _row(sha, "UNEXPECTED_SUBMISSION", "200", "sub2", mode="format"),  # and a format check stored as sub2
            _row(sha, "acknowledged", sid="sub2", mode="resolve")]
    entry = {"bundle_key": T1, "status": "queued", "sha256": sha}
    assert vs.board_budget(T1, rows, (0, 2), [], ident).spent == 2
    assert vs.board_budget(T1, rows, (1, 2), [{**entry, "id": "sub2"}], ident).spent == 2  # sub1 not listed yet
    both = [{**entry, "id": "sub1"}, {**entry, "id": "sub2"}]
    assert vs.board_budget(T1, rows, (2, 2), both, ident).spent == 2  # nothing counts twice
    rejected = {**entry, "id": "sub1", "status": "rejected"}  # the ambiguous attempt, listed as rejected
    assert vs.board_budget(T1, rows, (0, 2), [rejected], ident).spent == 1
    answered = [_row(sha, "queued", "200", "sub1")]  # an answer, while only a website copy of the file is listed
    assert vs.board_budget(T1, answered, (1, 2), [{**entry, "id": "web7"}], ident).spent == 2


def test_an_official_retry_never_replays_a_format_checks_submission(env):
    b, c = _checked(env, "b.h5ad", 1), _checked(env, "c.h5ad", 2)
    env.api.script = ["store_503", "503", "503", "503"]  # stored as sub1, every answer lost
    assert not _official(b).ok
    env.api.script = ["format_stored"]  # a later format check of b is stored too, as sub2
    assert vs.submit(b, T1).server_status == "UNEXPECTED_SUBMISSION"
    vs.resolve(vs.sha256_file(b)[:12], T1, "acknowledged", "the portal shows the stray format check")
    hidden = list(env.api.subs)
    env.api.subs.clear()  # the list and the quota lag
    res = _official(b)  # the identical request under b's key, which the server answers with sub1
    assert not res.replayed and res.ok and res.server_id == "sub1"
    with pytest.raises(vs.SubmitRefused, match="budget"):  # sub1 and sub2 both count
        _official(c)
    env.api.subs += hidden
    assert len(env.api.subs) == 2


def test_an_official_upload_of_a_file_the_board_holds_as_a_stray_is_refused(env):
    a = _checked(env, "a.h5ad", 0)
    env.api.script = ["format_stored"]  # a later format check of a is stored anyway
    assert vs.submit(a, T1).server_status == "UNEXPECTED_SUBMISSION"
    vs.resolve(vs.sha256_file(a)[:12], T1, "acknowledged", "the portal shows the stray")
    n = len(env.api.posts())
    with pytest.raises(vs.SubmitRefused, match="already held") as e:  # the same bytes would take a second slot
        _official(a)
    assert len(env.api.posts()) == n and vs.read_ledger(env.ledger)[-1]["detail"] == "failed: already held"
    guard = next(g for g in e.value.plan.guards if g.name == "already held")
    assert "already holds this file as submission sub" in guard.detail


def test_making_the_snapshot_directory_fails_with_a_refusal(env, monkeypatch):
    a = _pred(env, "a.h5ad", 0)
    monkeypatch.setattr(vs.tempfile, "tempdir", str(env.tmp / "no" / "such" / "dir"))
    with pytest.raises(vs.SubmitRefused, match="cannot make a temporary directory") as e:
        vs.submit(a, T1)
    assert "free some space" in str(e.value) and not env.api.posts()
    assert vs.read_ledger(env.ledger)[-1]["detail"] == "cannot make a temporary directory (FileNotFoundError)"


def test_a_replay_never_settles_an_ambiguous_upload(env):
    b = _checked(env, "b.h5ad", 1)
    env.api.script = ["store_503", "503", "503", "503"]  # stored as sub1, every answer lost
    assert not _official(b).ok
    res = _official(b)  # the server lists it: replayed, nothing is sent
    rows = vs.read_ledger(env.ledger)
    assert res.replayed and res.server_id == "sub1" and rows[-1]["mode"] == "replay"
    assert vs._ledger_files(rows, T1)[vs.sha256_file(b)].ambiguous  # only a listing, an answer or a person settles it
    assert vs.board_budget(T1, rows, (1, 2), env.api.subs, {"task": "T1"}).spent == 1  # listed: one slot
    assert vs.board_budget(T1, rows, (0, 2), [], {"task": "T1"}).spent == 1  # and still one if the list lags


# ---------------------------------------------------------------------------------------------------------------------
# A file that cannot be read; site text in the board commands


def test_a_file_that_cannot_be_read_is_not_blamed_on_the_disk(env):
    a, b = _pred(env, "a.h5ad", 0), _pred(env, "b.h5ad", 1)
    with pytest.raises(vs.SubmitRefused, match="cannot read") as e:  # gone between the plan and the snapshot
        vs.submit(a, T1, on_plan=lambda plan: plan.file.unlink())
    assert "FileNotFoundError" in str(e.value) and "free some space" not in str(e.value)
    assert vs.read_ledger(env.ledger)[-1]["detail"] == "cannot read the file (FileNotFoundError)"

    class Failing(io.RawIOBase):  # a disk error while the copy reads it
        def readable(self):
            return True

        def readinto(self, _buffer):
            raise OSError(errno.EIO, "Input/output error")

    class Flaky(type(b)):
        def open(self, *args, **kwargs):
            return Failing()

    with pytest.raises(vs.SubmitRefused, match="cannot read") as e:
        vs.submit(b, T1, on_plan=lambda plan: setattr(plan, "file", Flaky(plan.file)))
    assert "Input/output error" in str(e.value) and "free some space" not in str(e.value)
    assert vs.read_ledger(env.ledger)[-1]["detail"] == "cannot read the file (OSError)" and not env.api.posts()


def test_site_text_is_never_markup_in_the_board_commands(env, monkeypatch):
    runner = CliRunner()
    board = dataclasses.replace(load_boards(env.boards)[T1], key="T1:x[/]", floor_model="[red]copy_last")
    monkeypatch.setattr(cli, "load_boards", lambda: {board.key: board})
    res = runner.invoke(cli.app, ["boards"])
    assert res.exit_code == 0, _output(res)
    assert "T1:x[/]" in res.output and "[red]copy_last" in res.output
    # a Board with no floor: load_boards now gives a published null its task's floor, so only code builds this one
    blank = dataclasses.replace(board, floor_model=None)
    monkeypatch.setattr(cli, "load_boards", lambda: {blank.key: blank})
    res = runner.invoke(cli.app, ["boards"])
    assert res.exit_code == 0, _output(res)
    monkeypatch.setattr(cli, "load_boards", lambda: {board.key: board})
    summary = {"boards": ["T1:x[/]", "[bold]T2"], "reduced": ["[/]"], "t1_composition_written": False,
               "t1_composition_unreleased_stages": "[i]1"}
    monkeypatch.setattr("vec.boards.update_boards", lambda: summary)
    res = runner.invoke(cli.app, ["boards", "update"])
    assert res.exit_code == 0, _output(res)
    assert "boards (2): T1:x[/], [bold]T2" in res.output and "spec (1): [/]" in res.output
    assert "describes [i]1 unreleased" in res.output
    monkeypatch.setattr(cli, "get_board", lambda key: board)
    report = CheckReport(errors=["gene panel mismatch: 1 missing (e.g. ['G[/]'])"], warnings=["[red]odd"], n_cells=3)
    monkeypatch.setattr(cli, "check_submission", lambda pred, b: report)
    res = runner.invoke(cli.app, ["check", str(env.tmp / "x.h5ad"), "--board", board.key])
    assert res.exit_code == 1, _output(res)
    assert "error gene panel mismatch: 1 missing (e.g. ['G[/]'])" in res.output and "warning [red]odd" in res.output


# ---------------------------------------------------------------------------------------------------------------------
# Files built on the DGX: pulled into a private copy for one call, provenance for official uploads, no token there


@pytest.fixture
def dgx(env, monkeypatch):
    """The DGX for these tests: dgx_run runs each command with bash, as ssh would, in a temporary home whose DGX_REPO is
    a clone of the local repo (the same HEAD), and logs it. With `cat_bytes` set, a cat sends those bytes instead; with
    `stamp_bytes`, a cat of a build stamp does."""
    home = env.tmp / "dgx-home"
    repo = home / vs.DGX_REPO.removeprefix("~/")
    repo.parent.mkdir(parents=True)
    _git(repo.parent, "clone", "-q", str(env.repo), repo.name)
    (repo / "results").mkdir()
    fake = SimpleNamespace(home=home, repo=repo, commands=[], cat_bytes=None, stamp_bytes=None)

    def run(command, stderr):
        fake.commands.append(command)
        if fake.stamp_bytes is not None and command.startswith("cat -- ") and \
                shlex.split(command)[-1].endswith(vs.STAMP_SUFFIX):
            (home / "other").write_bytes(fake.stamp_bytes)
            command = "cat -- other"
        if fake.cat_bytes is not None and command.startswith("cat -- "):
            (home / "other").write_bytes(fake.cat_bytes)
            command = "cat -- other"
        return subprocess.Popen(["bash", "-c", command], cwd=home, env={"HOME": str(home), "PATH": os.environ["PATH"]},
                                stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=stderr)

    monkeypatch.setattr(vs, "dgx_run", run)
    return fake


def _on_dgx(env, dgx, rel, seed):
    """A valid prediction file at `rel` from the DGX repo root, with its build stamp next to it, as a builder there
    leaves them (built clean at the HEAD the DGX clone shares with the local repo)."""
    path = dgx.repo / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    built = _pred(env, f"dgx{seed}.h5ad", seed)
    shutil.move(built, path)
    shutil.move(vs.stamp_path(built), vs.stamp_path(path))
    return path


def _private_tmp(env, monkeypatch):
    """An empty directory for tempfile, to see what a call leaves behind."""
    tmp = env.tmp / "tmpdir"
    tmp.mkdir()
    monkeypatch.setattr(tempfile, "tempdir", str(tmp))
    return tmp


def test_a_format_check_from_the_dgx_sends_a_private_copy_and_leaves_nothing_behind(env, dgx, monkeypatch):
    tmp = _private_tmp(env, monkeypatch)
    rel = "results/t1/pred.h5ad"
    f = _on_dgx(env, dgx, rel, 0)
    sha, seen = vs.sha256_file(f), []

    def look(plan):
        seen.append((plan.file.parent.stat().st_mode & 0o777, plan.file.stat().st_mode & 0o777,
                     vs.sha256_file(plan.file)))

    res = vs.submit(rel, T1, from_dgx=True, on_plan=look)
    (post,) = env.api.posts()
    assert res.server_status == "format_ok" and seen == [(0o700, 0o600, sha)] and not list(tmp.iterdir())
    assert post["fields"] == {"task": "T1", "format_only": "true"} and post["filename"] == "pred.h5ad"
    assert post["sha256"] == sha
    report, cat, cat_stamp = dgx.commands
    assert report.startswith("bash -c ") and cat == f"cat -- {shlex.quote(str(f.resolve()))}"
    assert cat_stamp == f"cat -- {shlex.quote(str(f.resolve()) + vs.STAMP_SUFFIX)}"
    row = vs.read_ledger(env.ledger)[-1]
    assert (row["file"], row["upload_name"], row["sha256"]) == (rel, "pred.h5ad", sha)

    def boom(_plan):
        raise RuntimeError("interrupted")

    with pytest.raises(RuntimeError):  # errors too leave nothing behind
        vs.submit(rel, T1, from_dgx=True, on_plan=boom)
    assert not list(tmp.iterdir()) and len(env.api.posts()) == 1


def test_an_official_upload_from_the_dgx_needs_the_local_head_and_a_clean_dgx_tree(env, dgx, monkeypatch):
    tmp = _private_tmp(env, monkeypatch)
    rel, head = "results/pred.h5ad", _git(env.repo, "rev-parse", "HEAD")
    sha = vs.sha256_file(_on_dgx(env, dgx, rel, 0))
    assert vs.submit(rel, T1, from_dgx=True).server_status == "format_ok"

    def official():
        return vs.submit(rel, T1, official=True, confirm=sha[:12], disclose=X001, from_dgx=True, sleep=lambda s: None)

    def refused(match):
        with pytest.raises(vs.SubmitRefused, match="dgx provenance") as e:
            official()
        assert e.value.plan.failed == ["dgx provenance"] and not list(tmp.iterdir())
        assert match in next(g.detail for g in e.value.plan.guards if g.name == "dgx provenance")

    (dgx.repo / "method.py").write_text("x = 2\n")
    _git(dgx.repo, "commit", "-qam", "a newer commit on the DGX")
    refused(f"is not the local HEAD {head[:12]}")
    _git(dgx.repo, "reset", "-q", "--hard", head)
    (dgx.repo / "method.py").write_text("x = 3\n")
    refused("uncommitted or untracked on the DGX: method.py")
    _git(dgx.repo, "checkout", "-q", "--", "method.py")
    (dgx.repo / "notes.txt").write_text("an untracked file outside results/\n")
    refused("uncommitted or untracked on the DGX: notes.txt")
    (dgx.repo / "notes.txt").unlink()
    assert not env.api.posts(official=True)
    res = official()
    assert res.ok and res.server_id == "sub1" and len(env.api.posts(official=True)) == 1 and not list(tmp.iterdir())
    assert next(g for g in res.plan.guards if g.name == "dgx provenance").ok
    row = vs.read_ledger(env.ledger)[-1]
    assert (row["mode"], row["file"], row["origin"]) == ("official", rel, f"PrakashDGX_H2:{rel}@{head}")


def test_a_pull_that_does_not_match_the_report_is_refused_and_leaves_nothing_behind(env, dgx, monkeypatch):
    tmp = _private_tmp(env, monkeypatch)
    data = _on_dgx(env, dgx, "results/pred.h5ad", 0).read_bytes()
    for other in (data[:-1] + bytes([data[-1] ^ 1]), data[:100]):  # one byte changed (the same size); cut short
        dgx.cat_bytes = other
        with pytest.raises(vs.SubmitRefused, match="transfer mismatch"):
            vs.submit("results/pred.h5ad", T1, from_dgx=True)
        assert not list(tmp.iterdir())
    assert not env.api.requests and not env.ledger.exists()  # refused before any plan: nothing sent or recorded


def test_dgx_paths_with_spaces_quotes_or_shell_syntax_stay_one_argument(env, dgx):
    rel = "results/my pred's \"v2\" $(touch pwned) `touch pwned2`.h5ad"
    sha = vs.sha256_file(_on_dgx(env, dgx, rel, 0))
    res = vs.submit(rel, T1, from_dgx=True)
    assert res.server_status == "format_ok" and env.api.posts()[-1]["sha256"] == sha
    sent = env.api.posts()[-1]["filename"]
    assert sent == vs.upload_name(Path(rel)) == "my_pred_s__v2____touch_pwned___touch_pwned2_.h5ad"
    assert vs.read_ledger(env.ledger)[-1]["file"] == rel and not list(env.tmp.rglob("pwned*"))
    absolute = _on_dgx(env, dgx, "results/abs.h5ad", 1)
    assert vs.submit(str(absolute), T1, from_dgx=True).server_status == "format_ok"  # absolute, inside the repo
    assert vs.read_ledger(env.ledger)[-1]["file"] == "results/abs.h5ad"


def test_dgx_paths_outside_the_repo_or_with_a_newline_are_refused(env, dgx):
    outside = _on_dgx(env, dgx, "../outside.h5ad", 0).resolve()  # next to the DGX repo, not in it
    (dgx.repo / "results" / "link.h5ad").symlink_to(outside)
    for path in ("../outside.h5ad", "results/../../outside.h5ad", str(outside), "results/link.h5ad"):
        with pytest.raises(vs.SubmitRefused, match="outside the repository"):
            vs.submit(path, T1, from_dgx=True)
    with pytest.raises(vs.SubmitRefused, match="No such file"):
        vs.submit("results/missing.h5ad", T1, from_dgx=True)
    nl = _on_dgx(env, dgx, "results/new\nline.h5ad", 1)
    (dgx.repo / "results" / "nl.h5ad").symlink_to(nl)
    with pytest.raises(vs.SubmitRefused, match="holds a newline"):  # a link to such a name: refused after the report
        vs.submit("results/nl.h5ad", T1, from_dgx=True)
    n = len(dgx.commands)
    for path in ("results/new\nline.h5ad", "results/a\0b.h5ad", ""):
        with pytest.raises(vs.SubmitRefused, match="newline or NUL"):
            vs.submit(path, T1, from_dgx=True)
    assert len(dgx.commands) == n  # refused before any remote command
    assert not [c for c in dgx.commands if c.startswith("cat ")] and not env.api.requests


def test_the_token_never_reaches_a_remote_command_or_its_environment(env, dgx, monkeypatch):
    monkeypatch.setenv("GITHUB_TOKEN", "gh-FAKE-never-sent")
    monkeypatch.setenv("hf_token_path", "/somewhere")
    sha = vs.sha256_file(_on_dgx(env, dgx, "results/pred.h5ad", 0))
    assert vs.submit("results/pred.h5ad", T1, from_dgx=True).ok  # the token is loaded and used for this upload
    res = vs.submit("results/pred.h5ad", T1, official=True, confirm=sha[:12], disclose=X001, from_dgx=True)
    assert res.ok and len(dgx.commands) == 6 and not [c for c in dgx.commands if TOKEN in c or "gh-FAKE" in c]
    seen = {}  # the real runner: the exact ssh command line, and an environment without any TOKEN variable
    monkeypatch.setattr(subprocess, "Popen", lambda argv, **kw: seen.update(argv=argv, **kw))
    with tempfile.TemporaryFile() as err:
        REAL_DGX_RUN("cat -- pyproject.toml", err)
    assert seen["argv"] == ["ssh", "-o", "BatchMode=yes", "PrakashDGX_H2", "--", "cat -- pyproject.toml"]
    assert not [k for k in seen["env"] if "TOKEN" in k.upper()] and seen["env"]["PATH"] == os.environ["PATH"]
    assert TOKEN not in json.dumps(seen["env"]) and seen["stdin"] == subprocess.DEVNULL


def test_the_ledger_origin_column_names_local_files_and_the_dgx_state(env, dgx):
    vs.submit(_pred(env, "a.h5ad", 0), T1)
    rel, head = "results/pred.h5ad", _git(env.repo, "rev-parse", "HEAD")
    _on_dgx(env, dgx, rel, 1)
    vs.submit(rel, T1, from_dgx=True)
    (dgx.repo / "method.py").write_text("x = 9\n")  # recorded, but a format check does not require a clean DGX tree
    res = vs.submit(rel, T1, from_dgx=True)
    assert res.server_status == "format_ok"
    assert [w for w in res.plan.warnings if "would fail dgx provenance" in w and "method.py" in w]
    rows = vs.read_ledger(env.ledger)
    assert list(rows[0])[-2:] == ["origin", "build_commit"]
    assert [r["origin"] for r in rows] == ["local", f"PrakashDGX_H2:{rel}@{head}", f"PrakashDGX_H2:{rel}@{head}+dirty"]
    assert rows[0]["file"].endswith("a.h5ad") and [r["file"] for r in rows[1:]] == [rel, rel]
    assert {r["upload_name"] for r in rows[1:]} == {"pred.h5ad"} and rows[2]["git_commit"] == head


def test_cli_from_dgx(env, dgx):
    runner, rel = CliRunner(), "results/pred.h5ad"
    sha = vs.sha256_file(_on_dgx(env, dgx, rel, 0))
    res = runner.invoke(cli.app, ["submit", rel, "--board", T1, "--from-dgx"])
    out = _output(res)
    assert res.exit_code == 0, out
    assert f"file      PrakashDGX_H2:{rel}, pulled for this call (sent as pred.h5ad" in out
    assert f"origin    PrakashDGX_H2:{rel}@{_git(env.repo, 'rev-parse', 'HEAD')}" in out
    assert f"vec submit {rel} --board {T1} --from-dgx --official --confirm {sha[:12]} --disclose <ids" in out
    res = runner.invoke(cli.app, ["submit", "../outside.h5ad", "--board", T1, "--from-dgx"])
    assert res.exit_code == 1 and "nothing was sent" in _output(res)
    for args in (["submit", "--help"], ["submit", "upload", "--help"]):
        assert "--from-dgx" in runner.invoke(cli.app, args).output


def test_hangup_and_terminate_unwind_the_call_so_no_private_copy_remains(env, dgx, monkeypatch):
    """A closed terminal or a dropped session (SIGHUP), or SIGTERM, mid-call: the call unwinds and removes its copy."""
    tmp = _private_tmp(env, monkeypatch)
    rel = "results/pred.h5ad"
    _on_dgx(env, dgx, rel, 0)
    before = {s: signal.getsignal(s) for s in (signal.SIGHUP, signal.SIGTERM)}
    for sig in (signal.SIGHUP, signal.SIGTERM):
        def hang_up(plan, sig=sig):
            assert plan.file.is_file() and plan.file.parent.parent == tmp  # the pulled copy exists right now
            os.kill(os.getpid(), sig)
            time.sleep(5)  # the handler interrupts this at once

        with pytest.raises(SystemExit) as e:
            vs.submit(rel, T1, from_dgx=True, on_plan=hang_up)
        assert e.value.code == 128 + sig and not list(tmp.iterdir())
    assert {s: signal.getsignal(s) for s in before} == before and not env.api.posts()  # handlers restored, nothing sent


def test_copies_left_by_a_killed_call_are_removed_by_the_next_one(env, monkeypatch):
    tmp = _private_tmp(env, monkeypatch)
    now, old = time.time(), time.time() - (vs.STALE_COPY_HOURS + 1) * 3600
    for name, age in (("vec-dgx-old", old), ("vec-submit-old", old), ("vec-dgx-new", now), ("other-old", old)):
        (tmp / name).mkdir()
        (tmp / name / "copy").write_bytes(b"x")
        os.utime(tmp / name, (age, age))
    (tmp / "vec-dgx-link").symlink_to(tmp / "other-old")  # a link is never followed
    os.utime(tmp / "vec-dgx-link", (old, old), follow_symlinks=False)
    assert sorted(vs.remove_stale_copies(now)) == ["vec-dgx-old", "vec-submit-old"]
    assert sorted(p.name for p in tmp.iterdir()) == ["other-old", "vec-dgx-link", "vec-dgx-new"]
    assert (tmp / "other-old" / "copy").is_file()
    os.utime(tmp / "vec-dgx-new", (old, old))
    vs.submit(_pred(env, "a.h5ad", 0), T1)  # every call clears what killed calls left
    assert sorted(p.name for p in tmp.iterdir()) == ["other-old", "vec-dgx-link"]


def test_an_unreadable_dgx_status_is_marked_in_the_origin_and_refuses_official(env, dgx):
    rel = "results/pred.h5ad"
    sha = vs.sha256_file(_on_dgx(env, dgx, rel, 0))
    assert vs.submit(rel, T1, from_dgx=True).server_status == "format_ok"
    (dgx.repo / ".git" / "index").write_bytes(b"not an index")  # git status fails; git rev-parse HEAD still works
    res = vs.submit(rel, T1, from_dgx=True)
    head = _git(env.repo, "rev-parse", "HEAD")
    assert res.server_status == "format_ok" and vs.read_ledger(env.ledger)[-1]["origin"] == \
        f"PrakashDGX_H2:{rel}@{head}+unknown"
    with pytest.raises(vs.SubmitRefused, match="dgx provenance") as e:
        vs.submit(rel, T1, official=True, confirm=sha[:12], disclose=X001, from_dgx=True, sleep=lambda s: None)
    assert "cannot read the git state" in next(g.detail for g in e.value.plan.guards if g.name == "dgx provenance")


def test_a_linked_dgx_file_keeps_the_name_it_was_given(env, dgx):
    """As for a local file, the upload name and the format check use the name given, not the link target's."""
    target = _on_dgx(env, dgx, "runs/build/output", 0)  # a valid prediction under a name without .h5ad
    (dgx.repo / "results" / "latest.h5ad").symlink_to(target)
    res = vs.submit("results/latest.h5ad", T1, from_dgx=True)
    (post,) = env.api.posts()
    assert res.server_status == "format_ok" and post["filename"] == "latest.h5ad"
    row = vs.read_ledger(env.ledger)[-1]
    assert (row["file"], row["upload_name"]) == ("runs/build/output", "latest.h5ad")  # the ledger names the target


def test_the_build_stamp_next_to_a_links_target_comes_along_from_the_dgx(env, dgx, monkeypatch):
    tmp = _private_tmp(env, monkeypatch)
    target = _on_dgx(env, dgx, "runs/build/s1.h5ad", 0)
    (dgx.repo / "results" / "latest.h5ad").symlink_to(target)
    sha, head, seen = vs.sha256_file(target), _git(env.repo, "rev-parse", "HEAD"), []

    def look(plan):
        stamp, names = vs.stamp_path(plan.file), sorted(p.name for p in plan.file.parent.iterdir())
        seen.append((names, stamp.read_bytes(), stamp.stat().st_mode & 0o777))

    res = vs.submit("results/latest.h5ad", T1, from_dgx=True, on_plan=look)
    assert res.server_status == "format_ok" and not [w for w in res.plan.warnings if "build stamp" in w]
    assert seen == [(["latest.h5ad", "latest.h5ad" + vs.STAMP_SUFFIX], vs.stamp_path(target).read_bytes(), 0o600)]
    assert dgx.commands[-1] == f"cat -- {shlex.quote(str(target.resolve()) + vs.STAMP_SUFFIX)}"
    assert vs.read_ledger(env.ledger)[-1]["build_commit"] == head and not list(tmp.iterdir())
    res = vs.submit("results/latest.h5ad", T1, official=True, confirm=sha[:12], disclose=X001, from_dgx=True,
                    sleep=lambda s: None)
    assert res.ok and _guard(res.plan, "build stamp").ok and _guard(res.plan, "dgx provenance").ok
    row = vs.read_ledger(env.ledger)[-1]
    assert (row["mode"], row["origin"], row["build_commit"]) == ("official", f"PrakashDGX_H2:runs/build/s1.h5ad@{head}",
                                                                 head)
    assert not list(tmp.iterdir())


def test_a_dgx_file_without_a_usable_stamp_is_checked_but_never_sent_officially(env, dgx):
    rel = "results/pred.h5ad"
    f = _on_dgx(env, dgx, rel, 0)
    sha, stamp = vs.sha256_file(f), vs.stamp_path(f)
    good = stamp.read_bytes()
    stamp.unlink()
    oversized = good.rstrip(b"\n") + b" " * vs.STAMP_MAX_BYTES + b"\n"  # still valid JSON, but too large to pull
    for make, why in ((lambda: None, f"no build stamp next to {rel} on PrakashDGX_H2"),
                      (lambda: stamp.symlink_to(f"{f.name}.elsewhere.json"), "is a link or not a regular file"),
                      (lambda: stamp.mkdir(), "is a link or not a regular file"),
                      (lambda: stamp.write_bytes(oversized), f"is {len(oversized)} bytes, more than a build stamp's "
                                                             "64 KiB")):
        make()
        n = len(dgx.commands)
        res = vs.submit(rel, T1, from_dgx=True)  # a format check works without it, with a warning
        assert res.server_status == "format_ok" and len(dgx.commands) == n + 2, why  # the report and the file's cat
        assert "an official upload of this file would fail build stamp: " in res.plan.warnings[-1]
        assert why in res.plan.warnings[-1] and vs.read_ledger(env.ledger)[-1]["build_commit"] == ""
        with pytest.raises(vs.SubmitRefused, match="build stamp") as e:
            vs.submit(rel, T1, official=True, confirm=sha[:12], disclose=X001, from_dgx=True)
        assert e.value.plan.failed == ["build stamp"] and why in _guard(e.value.plan, "build stamp").detail
        if stamp.is_dir() and not stamp.is_symlink():
            stamp.rmdir()
        else:
            stamp.unlink(missing_ok=True)
    assert not env.api.posts(official=True)
    stamp.write_bytes(good)
    assert vs.submit(rel, T1, official=True, confirm=sha[:12], disclose=X001, from_dgx=True).ok


def test_a_stamp_pull_that_does_not_match_the_report_is_refused_and_leaves_nothing_behind(env, dgx, monkeypatch):
    tmp = _private_tmp(env, monkeypatch)
    data = vs.stamp_path(_on_dgx(env, dgx, "results/pred.h5ad", 0)).read_bytes()
    for other, n in ((data[:-2] + b"x\n", len(data)), (data + b" " * vs.STAMP_MAX_BYTES, vs.STAMP_MAX_BYTES + 1)):
        dgx.stamp_bytes = other  # one byte changed (the same size); more than any stamp, cut off after the cap
        with pytest.raises(vs.SubmitRefused, match=f"transfer mismatch: the pull of results/pred.h5ad.stamp.json gave "
                                                   f"{n} bytes"):
            vs.submit("results/pred.h5ad", T1, from_dgx=True)
        assert not list(tmp.iterdir())
    assert not env.api.requests and not env.ledger.exists()  # refused before any plan: nothing sent or recorded


def test_a_dgx_files_own_stamp_is_exempt_from_dgx_provenance_but_no_other_stamp(env, dgx):
    (env.repo / ".gitignore").write_text("*.h5ad\n")  # as in this repository: prediction files are never tracked
    _git(env.repo, "add", ".gitignore")
    _git(env.repo, "commit", "-qm", "ignore predictions")
    _git(dgx.repo, "pull", "-q")  # the DGX at the same HEAD
    rel = "preds/s1.h5ad"  # inside the DGX tree, outside results/ and runs/
    sha = vs.sha256_file(_on_dgx(env, dgx, rel, 0))
    _on_dgx(env, dgx, "preds/s2.h5ad", 1)  # another file and its stamp next to it
    assert vs.submit(rel, T1, from_dgx=True).server_status == "format_ok"

    def official():
        return vs.submit(rel, T1, official=True, confirm=sha[:12], disclose=X001, from_dgx=True, sleep=lambda s: None)

    with pytest.raises(vs.SubmitRefused, match="dgx provenance") as e:
        official()
    assert e.value.plan.failed == ["dgx provenance"]
    assert _guard(e.value.plan, "dgx provenance").detail == "uncommitted or untracked on the DGX: preds/s2.h5ad" + \
        vs.STAMP_SUFFIX
    vs.stamp_path(dgx.repo / "preds" / "s2.h5ad").unlink()
    assert official().ok  # its own stamp, untracked next to it, is exempt like the ledger
