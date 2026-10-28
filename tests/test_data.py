import json
import socket
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import pytest

from vec.data import AuthError, DataFile, load_token

TOKEN = "tok-FAKE-never-print-3f9c1a77"
BROKEN = "tok-FAKE-part1\npart2-never-print"  # e.g. after `echo NEW >> ~/.config/vec/token`


@pytest.fixture(autouse=True)
def no_internet(monkeypatch):
    real = socket.create_connection

    def local_only(address, *args, **kwargs):
        if address[0] not in ("127.0.0.1", "localhost"):
            raise AssertionError(f"a test tried to reach {address[0]}")
        return real(address, *args, **kwargs)

    monkeypatch.setattr(socket, "create_connection", local_only)


def test_token_file_must_be_private(tmp_path, monkeypatch):
    monkeypatch.delenv("VEC_TOKEN", raising=False)
    path = tmp_path / "token"
    path.write_text("abc\n")
    path.chmod(0o644)
    with pytest.raises(AuthError, match="readable by others"):
        load_token(path)
    path.chmod(0o600)
    assert load_token(path) == "abc"


def test_token_accepts_pasted_session_json(tmp_path, monkeypatch):
    monkeypatch.delenv("VEC_TOKEN", raising=False)
    path = tmp_path / "token"
    path.write_text(json.dumps({"token": "xyz", "account": {}, "expiresAt": 0}))
    path.chmod(0o600)
    assert load_token(path) == "xyz"


@pytest.mark.parametrize("key", ["../etc/passwd", "a/../../b", "x y.h5ad", "a;rm -rf.h5ad"])
def test_unsafe_keys_rejected(key):
    with pytest.raises(ValueError):
        _ = DataFile(key=key, task="T1", role="train", stage="E8.5", modality="scRNA", bytes=1).relpath


def test_tilde_never_reaches_a_quoted_remote_path():
    from vec.data import resolve_dest

    assert resolve_dest("~/v/proj/data/raw", remote=True) == "v/proj/data/raw"
    assert resolve_dest("/abs/path", remote=True) == "/abs/path"
    assert not resolve_dest("~/x", remote=False).startswith("~")


def test_data_module_never_posts():
    import inspect

    import vec.data

    src = inspect.getsource(vec.data)
    assert "method=\"POST\"" not in src and "data=" not in src.replace("metadata", "")


def test_safe_key_kept():
    assert DataFile(key="T1/E8.5_RNA.h5ad", task="T1", role="train", stage="E8.5", modality="scRNA",
                    bytes=1).relpath == "T1/E8.5_RNA.h5ad"


def test_api_get_never_follows_a_redirect_with_the_token(monkeypatch):
    import threading
    from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

    import vec.data as vd

    seen = []

    class Target(BaseHTTPRequestHandler):
        def log_message(self, *args):
            pass

        def do_GET(self):
            seen.append(self.headers.get("Authorization"))
            self.send_response(200)
            self.end_headers()
            self.wfile.write(b"{}")

    target = ThreadingHTTPServer(("127.0.0.1", 0), Target)

    class Redirect(Target):
        def do_GET(self):
            self.send_response(302)
            self.send_header("Location", f"http://localhost:{target.server_port}/stolen")
            self.end_headers()

    first = ThreadingHTTPServer(("127.0.0.1", 0), Redirect)
    for srv in (target, first):
        threading.Thread(target=srv.serve_forever, daemon=True).start()
    try:
        monkeypatch.setattr(vd, "API_BASE", f"http://127.0.0.1:{first.server_port}")
        with pytest.raises(RuntimeError, match="302"):
            vd.api_get("/challenge/data/manifest", "tok-FAKE")
        assert seen == []
    finally:
        target.shutdown()
        first.shutdown()


@pytest.mark.parametrize("raw", [BROKEN, "tok-FAKE part2-never-print", "tok-FAKE\tpart2-never-print",
                                 "tok-FAKE\x00part2-never-print", "\ufefftok-FAKE-part2-never-print",
                                 "tok-FAKE-\u00e9-part2-never-print", '{"token": "tok-FAKE part2-never-print"}',
                                 '{"token": 5}', '{"tok-FAKE-part2-never-print', "", " \n "])
def test_a_malformed_token_is_refused_without_echoing_it(tmp_path, monkeypatch, raw):
    monkeypatch.delenv("VEC_TOKEN", raising=False)
    path = tmp_path / "token"
    path.write_text(raw + "\n")
    path.chmod(0o600)
    with pytest.raises(AuthError) as e:
        load_token(path)
    assert "part1" not in str(e.value) and "part2-never-print" not in str(e.value)
    monkeypatch.setenv("VEC_TOKEN", BROKEN)
    with pytest.raises(AuthError) as e:
        load_token(path)
    assert "part1" not in str(e.value) and "part2-never-print" not in str(e.value)


def test_a_malformed_token_never_reaches_the_terminal(tmp_path, monkeypatch):
    from typer.testing import CliRunner

    import vec.data as vd
    from vec import cli

    monkeypatch.setattr(vd, "API_BASE", "http://127.0.0.1:9")  # loopback, whatever the code under test does
    monkeypatch.setenv("VEC_TOKEN", BROKEN)
    res = CliRunner().invoke(cli.app, ["data", "manifest", "--token-file", str(tmp_path / "unused")])
    assert res.exit_code == 2 and "auth error" in res.output
    assert "part1" not in res.output and "part2-never-print" not in res.output


def test_redact_reaches_keys_values_and_escaped_forms():
    from vec.data import redact

    body = {"detail": [f"bad {TOKEN}", {f"k {TOKEN}": 1}], "n": 2}
    assert redact(body, TOKEN) == {"detail": ["bad <token>", {"k <token>": 1}], "n": 2}
    assert redact(f"Invalid header value {repr(f'Bearer {BROKEN}'.encode())}", BROKEN) == \
        "Invalid header value b'Bearer <token>'"
    assert redact(body, "") == body


def test_api_get_redacts_the_token_from_server_errors(monkeypatch):
    import vec.data as vd

    answers = [(401, "invalid credentials "), (500, "x" * 170 + " ")]  # the second puts the token across the cut

    class Echo(BaseHTTPRequestHandler):
        def log_message(self, *args):
            pass

        def do_GET(self):
            code, text = answers.pop(0)
            raw = json.dumps({"detail": text + self.headers.get("Authorization")}).encode()
            self.send_response(code)
            self.send_header("Content-Length", str(len(raw)))
            self.end_headers()
            self.wfile.write(raw)

    srv = ThreadingHTTPServer(("127.0.0.1", 0), Echo)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    try:
        monkeypatch.setattr(vd, "API_BASE", f"http://127.0.0.1:{srv.server_port}")
        with pytest.raises(AuthError) as e:
            vd.api_get("/challenge/data/manifest", TOKEN)
        assert "<token>" in str(e.value) and TOKEN not in str(e.value)
        with pytest.raises(RuntimeError) as e:
            vd.api_get("/challenge/data/manifest", TOKEN)
        assert TOKEN[:8] not in str(e.value)
    finally:
        srv.shutdown()


def test_api_get_turns_bad_header_values_and_urls_into_redacted_errors(monkeypatch):
    import vec.data as vd

    monkeypatch.setattr(vd, "API_BASE", "http://127.0.0.1:9")
    with pytest.raises(RuntimeError) as e:  # http.client quotes the whole header value in its ValueError
        vd.api_get("/challenge/data/manifest", BROKEN)
    assert "part1" not in str(e.value) and "part2-never-print" not in str(e.value) and e.value.__suppress_context__
    monkeypatch.setattr(vd, "API_BASE", "http://127.0.0.1:notaport")
    with pytest.raises(RuntimeError, match="InvalidURL"):
        vd.api_get("/challenge/data/manifest", TOKEN)


def _wide_cli(monkeypatch):
    from rich.console import Console

    from vec import cli

    monkeypatch.setattr(cli, "console", Console(emoji=False, highlight=False, width=10_000))  # nothing wraps or crops
    return cli


def _api(routes, seen):
    """Loopback API: a GET is answered by path from `routes` (path -> answers, the last one repeated), a request of
    any other method is refused; every request goes into `seen` as (method, path)."""

    class Api(BaseHTTPRequestHandler):
        def log_message(self, *args):
            pass

        def _answer(self, code, body):
            raw = json.dumps(body).encode()
            self.send_response(code)
            self.send_header("Content-Length", str(len(raw)))
            self.end_headers()
            self.wfile.write(raw)

        def do_GET(self):
            seen.append(("GET", self.path))
            answers = routes.get(self.path) or [(404, {"detail": "not found"})]
            self._answer(*(answers.pop(0) if len(answers) > 1 else answers[0]))

        def do_POST(self):
            seen.append(("POST", self.path))
            self._answer(405, {"detail": "refused"})

    srv = ThreadingHTTPServer(("127.0.0.1", 0), Api)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    return srv


TRAIN_ROWS = [{"key": "T2/E9.5.h5ad", "task": "T2", "role": "train", "stage": "E9.5", "modality": "MERFISH",
               "bytes": 5}]
VAL_ROWS = [{"key": "T2/E8.5.h5ad", "task": "T2", "role": "validation", "stage": "E8.5", "modality": "MERFISH",
             "bytes": 7, "sha256": "ab12"}]
UNLOCKED = {"unlocked": True, "available": True, "lock_scope": "team"}


def test_validation_answers_are_listed_and_linked_from_their_own_endpoints_once_unlocked(monkeypatch):
    import vec.data as vd

    seen = []
    srv = _api({"/challenge/data/manifest": [(200, {"files": TRAIN_ROWS})],
                "/challenge/data/validation/state": [(200, UNLOCKED)],
                "/challenge/data/validation/manifest": [(200, {"files": VAL_ROWS})],
                "/challenge/data/validation/link?key=T2/E8.5.h5ad": [(200, {"url": "https://signed/val"})],
                "/challenge/data/link?key=T2/E9.5.h5ad": [(200, {"url": "https://signed/train"})]}, seen)
    try:
        monkeypatch.setattr(vd, "API_BASE", f"http://127.0.0.1:{srv.server_port}")
        files, state = vd.released_files(TOKEN)
        assert [(f.key, f.source, f.bytes, f.extra) for f in files] == [
            ("T2/E9.5.h5ad", "train", 5, {}), ("T2/E8.5.h5ad", "validation", 7, {"sha256": "ab12"})]
        assert state == UNLOCKED
        assert vd.signed_link(files[1].key, TOKEN, files[1].source) == "https://signed/val"
        assert vd.signed_link(files[0].key, TOKEN, files[0].source) == "https://signed/train"
    finally:
        srv.shutdown()
    assert [m for m, _ in seen] == ["GET"] * 5 and not any("unlock" in p for _, p in seen)


@pytest.mark.parametrize("state", [{"unlocked": False, "available": True}, {"unlocked": True, "available": False}, {}])
def test_a_locked_or_unavailable_validation_release_is_never_listed(monkeypatch, state):
    import vec.data as vd

    seen = []
    srv = _api({"/challenge/data/manifest": [(200, {"files": TRAIN_ROWS})],
                "/challenge/data/validation/state": [(200, state)]}, seen)
    try:
        monkeypatch.setattr(vd, "API_BASE", f"http://127.0.0.1:{srv.server_port}")
        files, _ = vd.released_files(TOKEN)
    finally:
        srv.shutdown()
    assert [f.source for f in files] == ["train"]
    assert seen == [("GET", "/challenge/data/manifest"), ("GET", "/challenge/data/validation/state")]


def test_a_validation_key_that_is_also_a_training_key_is_refused(monkeypatch):
    import vec.data as vd

    srv = _api({"/challenge/data/manifest": [(200, {"files": TRAIN_ROWS})],
                "/challenge/data/validation/state": [(200, UNLOCKED)],
                "/challenge/data/validation/manifest": [(200, {"files": [dict(VAL_ROWS[0], key="T2/E9.5.h5ad")]})]},
               [])
    try:
        monkeypatch.setattr(vd, "API_BASE", f"http://127.0.0.1:{srv.server_port}")
        with pytest.raises(RuntimeError, match="both as training and validation"):
            vd.released_files(TOKEN)
    finally:
        srv.shutdown()


def test_no_code_can_unlock_the_validation_release():
    """The unlock permanently ends P2 submissions for the whole team: only a person on the site's Data page makes it."""
    from pathlib import Path

    import vec

    for path in Path(vec.__file__).parent.rglob("*.py"):
        src = path.read_text()
        assert "validation/unlock" not in src and "confirm_end_p2" not in src, path


def test_fetch_one_asks_the_validation_link_endpoint_and_keeps_a_finished_unsized_file(tmp_path, monkeypatch):
    import vec.data as vd

    asked, ran = [], []
    monkeypatch.setattr(vd, "signed_link", lambda key, token, source="train": asked.append((key, source)) or "https://x")

    def curl_then_mv(cmd, input, text, check):  # what the shell line does once curl succeeds
        ran.append(input)
        (tmp_path / "T2").mkdir(exist_ok=True)
        (tmp_path / "T2" / "E8.5.h5ad").write_bytes(b"1234567")

    monkeypatch.setattr(vd.subprocess, "run", curl_then_mv)
    f = DataFile(key="T2/E8.5.h5ad", task="T2", role="validation", stage="E8.5", modality="MERFISH", bytes=0,
                 source=vd.VALIDATION)
    assert vd.fetch_one(f, TOKEN, str(tmp_path), None) == "downloaded"
    assert asked == [("T2/E8.5.h5ad", "validation")] and ran == ['url = "https://x"\n']
    assert vd.fetch_one(f, TOKEN, str(tmp_path), None) == "present"  # no size listed: the finished file is kept
    sized = DataFile(key="T2/E8.5.h5ad", task="T2", role="validation", stage="E8.5", modality="MERFISH", bytes=9,
                     source=vd.VALIDATION)
    with pytest.raises(RuntimeError, match="got 7 bytes, manifest says 9"):  # a listed size still decides
        vd.fetch_one(sized, TOKEN, str(tmp_path), None)
    assert len(asked) == 2


def test_server_text_is_never_markup_in_the_data_commands(tmp_path, monkeypatch):
    from typer.testing import CliRunner

    import vec.data as vd

    files = [{"key": "T1/[red]E8.5[/].h5ad", "task": "[bold]T1", "role": "train", "stage": "E8.5", "modality": "[/]rna",
              "bytes": 10}]
    val = [{"key": "T1/[bold]E10.5.h5ad", "task": "T1", "role": "[/]val", "stage": "E10.5", "modality": "rna",
            "note": "[red]x"}]
    srv = _api({"/challenge/data/manifest": [(401, {"detail": "expired [/]"}), (200, {"files": files})],
                "/challenge/data/validation/state": [(200, {**UNLOCKED, "lock_scope": "[red]team"})],
                "/challenge/data/validation/manifest": [(200, {"files": val})]}, [])
    try:
        cli = _wide_cli(monkeypatch)
        monkeypatch.setattr(vd, "API_BASE", f"http://127.0.0.1:{srv.server_port}")
        monkeypatch.setattr(vd, "fetch_one", lambda f, token, dest, host: "present")  # nothing is downloaded
        monkeypatch.setenv("VEC_TOKEN", TOKEN)
        runner, unused = CliRunner(), ["--token-file", str(tmp_path / "unused")]
        res = runner.invoke(cli.app, ["data", "manifest", *unused])  # an AuthError quoting the server's detail
        assert res.exit_code == 2 and "auth error HTTP 401" in res.output and "expired [/]" in res.output
        res = runner.invoke(cli.app, ["data", "manifest", *unused])
        assert res.exit_code == 0, res.output
        assert "T1/[red]E8.5[/].h5ad" in res.output and "[bold]T1" in res.output and "[/]rna" in res.output
        assert "T1/[bold]E10.5.h5ad" in res.output and "[/]val" in res.output and '{"note": "[red]x"}' in res.output
        assert '"lock_scope": "[red]team"' in res.output and "1 validation files" in res.output
        res = runner.invoke(cli.app, ["data", "fetch", "--host", "", *unused])
        assert res.exit_code == 0, res.output
        assert "[1/2] present    T1/[red]E8.5[/].h5ad (0 MB, train)" in res.output
        assert "[2/2] present    T1/[bold]E10.5.h5ad (size not listed, validation)" in res.output
        assert TOKEN not in res.output
    finally:
        srv.shutdown()


def test_file_names_and_labels_are_never_markup_in_inventory_and_qc(tmp_path, monkeypatch):
    import anndata
    from typer.testing import CliRunner

    import vec.data as vd
    from vec import spatial

    cli = _wide_cli(monkeypatch)
    f = tmp_path / "raw" / "T1" / "[red]E8.5.h5ad"  # named by the server's key
    f.parent.mkdir(parents=True)
    f.write_bytes(b"")
    monkeypatch.setattr(vd, "inventory", lambda path: {"n_obs": 1, "n_vars": 2, "nonzero_fraction": 0.5,
                                                       "expm1_rowsum_median": 1.0, "layers": ["[/]counts"],
                                                       "obsm_keys": ["[bold]xy"]})
    res = CliRunner().invoke(cli.app, ["data", "inventory", "--root", str(tmp_path / "raw"), "--out",
                                       str(tmp_path / "inventory.json")])
    assert res.exit_code == 0, res.output
    assert "T1/[red]E8.5.h5ad" in res.output and "[/]counts" in res.output and "[bold]xy" in res.output
    row = {**dict.fromkeys(["non_integer_count_rows", "z_raw_values", "z_levels", "z_dense_levels", "z_step", "z_gaps",
                            "spearman_row_vs_z", "z_runs", "dup_redundant", "dup_by_gap", "dup_any_xy_redundant",
                            "identical_rows_min_library"], ""),
           "cells": 1, "obs_names_duplicated": 0, "identical_rows": 0, "section_labels": "[/]sec:ok",
           "centred_on": "[red]origin"}
    monkeypatch.setattr(anndata, "read_h5ad", lambda path: None)
    monkeypatch.setattr(spatial, "qc_stage", lambda a, **kw: {"cells": 1})
    monkeypatch.setattr(spatial, "summary_row", lambda file, qc: row)
    res = CliRunner().invoke(cli.app, ["data", "qc", str(f), "--out", str(tmp_path / "qc.json")])
    assert res.exit_code == 0, res.output
    assert "T1/[red]E8.5.h5ad" in res.output and "[/]sec:ok" in res.output and "[red]origin" in res.output
