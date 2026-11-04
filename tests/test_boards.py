import hashlib
import json

import pytest

from vec.boards import IntakeError, load_boards, sanitise_index, update_boards, update_boards_from_file


def _spec(key, split, anchors=True):
    s = {"key": key, "task": key.split(":")[0], "split": split, "genes_file": key.replace(":", "__") + ".genes.txt",
         "needs_coords": False, "min_cells": 1000, "max_cells": 5000, "truth_cells": 1500, "ref_cells": 1700,
         "floor_model": "copy_last", "secret_extra": 1}
    if anchors:
        s["anchors"] = {"mmd_u": {"floor": 0.1, "ceiling": 0.01}}
    return s


def test_sanitise_strips_everything_but_the_spec_from_non_validation_boards():
    index = {"T1:val": _spec("T1:val", "val"), "T1:test": _spec("T1:test", "test"),
             "T3:gata4": _spec("T3:gata4", "test"), "T2:heart:new": _spec("T2:heart:new", "val_extrap")}
    clean, reduced = sanitise_index(index)
    assert clean["T1:val"] == index["T1:val"]
    assert sorted(reduced) == ["T1:test", "T2:heart:new", "T3:gata4"]
    for k in reduced:
        assert "anchors" not in clean[k] and "secret_extra" not in clean[k]
        assert clean[k]["truth_cells"] == 1500 and clean[k]["max_cells"] == 5000


def test_update_boards_writes_no_test_anchors_and_quarantines_unreleased_composition(tmp_path):
    genes = b"Gata4\nNkx2-5\n"
    index = {"T1:val": _spec("T1:val", "val"), "T1:test": _spec("T1:test", "test")}
    for s in index.values():
        s["genes_sha256"] = hashlib.sha256(genes).hexdigest()
    comp = {"stages": {"E8.5": {}, "E9.5": {}, "E12.5": {"x": 1}}}
    files = {"index.json": json.dumps(index).encode(), "t1_composition.json": json.dumps(comp).encode(),
             "T1__val.genes.txt": genes, "T1__test.genes.txt": genes}
    summary = update_boards(tmp_path, fetch=lambda url: files[url.rsplit("/", 1)[1]])
    written = json.loads((tmp_path / "index.json").read_text())
    assert "anchors" not in written["T1:test"] and written["T1:val"]["anchors"]
    assert summary["reduced"] == ["T1:test"]
    assert not (tmp_path / "t1_composition.json").exists() and summary["t1_composition_unreleased_stages"] == 1
    boards = load_boards(tmp_path)
    assert boards["T1:test"].anchors == {} and boards["T1:test"].truth_cells == 1500


def test_update_boards_rejects_a_bad_gene_hash(tmp_path):
    index = {"T1:val": _spec("T1:val", "val")}
    index["T1:val"]["genes_sha256"] = "0" * 64
    files = {"index.json": json.dumps(index).encode(), "T1__val.genes.txt": b"A\n"}
    with pytest.raises(ValueError, match="sha256"):
        update_boards(tmp_path, fetch=lambda url: files[url.rsplit("/", 1)[1]])


def _phase(split):
    body = {"phase": "p3" if split == "test" else "p2", "split": split, "accepts_submissions": True}
    return json.dumps(body).encode()


@pytest.mark.parametrize("phase, state", [(_phase("val"), "closed"), (_phase("test"), "open"), (None, "unreadable"),
                                          (b"<html>", "unreadable"), (b"[]", "unreadable")])
def test_e10_5_composition_waits_for_the_final_phase(tmp_path, phase, state):
    """E10.5 is the T1 validation answer, released with the final phase (moved to 2026-10-27): its composition is
    written only when the public phase endpoint reports the test split; an unreachable or odd answer keeps it out."""
    genes = b"Gata4\n"
    index = {"T1:val": _spec("T1:val", "val")}
    index["T1:val"]["genes_sha256"] = hashlib.sha256(genes).hexdigest()
    comp = {"stages": {"E8.5": {}, "E9.5": {}, "E10.5": {}}}
    files = {"index.json": json.dumps(index).encode(), "t1_composition.json": json.dumps(comp).encode(),
             "T1__val.genes.txt": genes}
    if phase is not None:
        files["phase"] = phase

    def fetch(url):
        name = url.rsplit("/", 1)[1]
        if name not in files:
            raise OSError("unreachable")
        return files[name]

    summary = update_boards(tmp_path, fetch=fetch)
    written = state == "open"
    assert (tmp_path / "t1_composition.json").exists() is written
    assert summary["t1_composition_written"] is written and summary["final_phase"] == state
    assert summary["t1_composition_unreleased_stages"] == (0 if written else 1)


def test_released_only_composition_is_written_without_asking_the_phase(tmp_path):
    genes = b"Gata4\n"
    index = {"T1:val": _spec("T1:val", "val")}
    index["T1:val"]["genes_sha256"] = hashlib.sha256(genes).hexdigest()
    comp = {"stages": {"E8.5": {}, "E9.5": {}}}
    files = {"index.json": json.dumps(index).encode(), "t1_composition.json": json.dumps(comp).encode(),
             "T1__val.genes.txt": genes}
    asked = []

    def fetch(url):
        asked.append(url)
        return files[url.rsplit("/", 1)[1]]

    summary = update_boards(tmp_path, fetch=fetch)
    assert summary["t1_composition_written"] and summary["final_phase"] is None
    assert not any(u.endswith("/challenge/phase") for u in asked)


# ---------------------------------------------------------------------------------------------------------------------
# Offline intake of a board index the user saved with the browser (2026-10-28). The sentinel stands for a test-board
# anchor: it must never reach the terminal or any file in boards/, on success or on any refusal.

SENTINEL = 0.918273645
TAG = "918273645"
GENES = b"Gata4\nNkx2-5\n"


def _final_phase_index(genes=GENES, test_sha=True):
    """A final-phase index: a validation board with its anchors and a test board whose anchors and unknown field
    carry the sentinel."""
    val, test = _spec("T1:val", "val"), _spec("T1:test", "test")
    test["anchors"] = {"mmd_u": {"floor": SENTINEL, "ceiling": SENTINEL / 10}}
    test["secret_extra"] = f"hidden-{TAG}"
    val["genes_sha256"] = hashlib.sha256(GENES).hexdigest()
    if test_sha:
        test["genes_sha256"] = hashlib.sha256(genes).hexdigest()
    return {"T1:val": val, "T1:test": test}


def _save(folder, index, genes=GENES, comp=None, gene_files=("T1__val.genes.txt", "T1__test.genes.txt")):
    folder.mkdir(parents=True, exist_ok=True)
    raw = index if isinstance(index, bytes) else json.dumps(index).encode()
    (folder / "index.json").write_bytes(raw)
    for name in gene_files:
        (folder / name).write_bytes(GENES if name.startswith("T1__val") else genes)
    if comp is not None:
        (folder / "t1_composition.json").write_bytes(comp if isinstance(comp, bytes) else json.dumps(comp).encode())
    return folder / "index.json"


def _tree(d):
    return {p.name: p.read_bytes() for p in sorted(d.iterdir())} if d.exists() else {}


def _intake_cli(src, boards_dir, *extra):
    from typer.testing import CliRunner

    from vec import cli

    return CliRunner().invoke(cli.app, ["boards", "update", "--from-file", str(src), "--boards-dir", str(boards_dir),
                                        *extra])


def test_file_intake_writes_what_the_online_intake_writes(tmp_path):
    index = _final_phase_index()
    comp = {"stages": {"E8.5": {}, "E9.5": {}, "E10.5": {}}}
    files = {"index.json": json.dumps(index).encode(), "t1_composition.json": json.dumps(comp).encode(),
             "T1__val.genes.txt": GENES, "T1__test.genes.txt": GENES, "phase": _phase("test")}
    online, offline = tmp_path / "online", tmp_path / "offline"
    online.mkdir()
    offline.mkdir()
    a = update_boards(online, fetch=lambda url: files[url.rsplit("/", 1)[1]])
    src = _save(tmp_path / "saved", index, comp=comp)
    b = update_boards_from_file(src, offline, fetch=lambda url: files[url.rsplit("/", 1)[1]])
    assert _tree(online) == _tree(offline) and set(_tree(online)) == {"index.json", "t1_composition.json",
                                                                      "T1__val.genes.txt", "T1__test.genes.txt"}
    assert {k: b[k] for k in a} == a and a["reduced"] == ["T1:test"] and a["final_phase"] == "open"
    assert b["t1_composition_supplied"] and b["genes_reused"] == [] and b["deleted"] == [] and src.exists()
    written = json.loads((offline / "index.json").read_text())
    assert written["T1:val"]["anchors"] == index["T1:val"]["anchors"]  # validation anchors stay (calibration)
    assert TAG not in (offline / "index.json").read_text()


def test_file_intake_keeps_the_sentinel_out_of_the_terminal_and_boards(tmp_path):
    src = _save(tmp_path / "saved", _final_phase_index(), comp={"stages": {"E8.5": {}, "E12.5": {"cm": 0.5}}})
    boards_dir = tmp_path / "boards"
    boards_dir.mkdir()
    res = _intake_cli(src, boards_dir)
    assert res.exit_code == 0, res.output
    assert TAG not in res.output and all(TAG.encode() not in data for data in _tree(boards_dir).values())
    out = " ".join(res.output.split())
    assert "offline intake" in out and "submission spec (1): T1:test" in out and "still hold the test anchors" in out
    assert "T1:test: fields not kept: anchors, secret_extra" in out  # names only, for decision point D1
    assert "not written: it describes 1 unreleased" in out and not (boards_dir / "t1_composition.json").exists()
    assert load_boards(boards_dir)["T1:test"].anchors == {} and src.exists()


def _bad_inputs():
    index = _final_phase_index()
    wrapped = {"version": 2, "boards": index}
    path_name = _final_phase_index()
    path_name["T1:test"]["genes_file"] = f"../{TAG}.txt"
    no_coords = _final_phase_index()
    del no_coords["T1:test"]["needs_coords"]
    odd_sha = _final_phase_index()
    odd_sha["T1:test"]["genes_sha256"] = SENTINEL
    not_genes, nested, fraction, split_list = (_final_phase_index() for _ in range(4))
    not_genes["T1:test"]["genes_file"] = "t1_composition.json"
    nested["T1:test"]["floor_model"] = {"name": "copy_last", "mmd_u": {"floor": SENTINEL}}
    fraction["T1:test"]["max_cells"] = SENTINEL
    split_list["T1:test"]["split"] = ["test", SENTINEL]
    return {
        "a list": (json.dumps([index]).encode(), None, "not a JSON object of boards"),
        "a wrapped index": (json.dumps(wrapped).encode(), None, "board 'version'"),
        "a spec that is a number": (json.dumps({"T1:test": SENTINEL}).encode(), None, "not a JSON object"),
        "a saved web page": (f"<html><body>Just a moment... {TAG}</body></html>".encode(), None, "not JSON"),
        "bytes that are not text": (b"\xff\xfe" + TAG.encode(), None, "not JSON"),
        "a genes_file path": (json.dumps(path_name).encode(), None, "not a plain gene-list name"),
        "a missing field": (json.dumps(no_coords).encode(), None, "no needs_coords"),
        "a sha256 that is not text": (json.dumps(odd_sha).encode(), None, "genes_sha256 is not a string"),
        "a composition page": (json.dumps(index).encode(), f"<html>{TAG}</html>".encode(), "not JSON"),
        "composition stages that are a number": (json.dumps(index).encode(), json.dumps({"stages": SENTINEL}).encode(),
                                                 "no object of stages"),
        # found by the 2026-10-28 review: each of these used to reach boards/ with the sentinel
        "a composition that is the index": (json.dumps(index).encode(), json.dumps(index).encode(),
                                            "no object of stages"),
        "a composition without its stages wrapper": (json.dumps(index).encode(),
                                                     json.dumps({"E8.5": {}, "E12.5": {"cm": SENTINEL}}).encode(),
                                                     "no object of stages"),
        "a composition with an extra field": (json.dumps(index).encode(),
                                              json.dumps({"stages": {"E8.5": {}}, "held_out": {"E12.5": SENTINEL}})
                                              .encode(), "unexpected top-level field(s) held_out"),
        "a genes_file that is not a gene list": (json.dumps(not_genes).encode(), None, "not a plain gene-list name"),
        "a nested floor model": (json.dumps(nested).encode(), None, "unexpected value type for floor_model (dict)"),
        "a fractional cell limit": (json.dumps(fraction).encode(), None, "max_cells (float)"),
        "a split that is a list": (json.dumps(split_list).encode(), None, "split (list)"),
    }


@pytest.mark.parametrize("case", sorted(_bad_inputs()))
def test_file_intake_refuses_without_echoing_the_file(tmp_path, case):
    raw, comp, said = _bad_inputs()[case]
    src = _save(tmp_path / "saved", raw, comp=comp)
    boards_dir = tmp_path / "boards"
    boards_dir.mkdir()
    (boards_dir / "index.json").write_text("{}\n")  # the old specs stay as they are
    res = _intake_cli(src, boards_dir)
    out = " ".join(res.output.split())
    assert res.exit_code == 1 and said in out and "Nothing was written" in out and "kept for a re-run" in out
    assert TAG not in res.output and _tree(boards_dir) == {"index.json": b"{}\n"} and src.exists()


def test_an_unexpected_error_in_the_file_intake_is_reported_by_type_only(tmp_path, monkeypatch):
    def boom(index):
        raise RuntimeError(f"bug while reading {index}")  # its message would quote the saved file

    monkeypatch.setattr("vec.boards.sanitise_index", boom)
    boards_dir = tmp_path / "boards"
    boards_dir.mkdir()
    res = _intake_cli(_save(tmp_path / "saved", _final_phase_index()), boards_dir)
    assert res.exit_code == 1 and "RuntimeError" in res.output and "withheld" in res.output
    assert TAG not in res.output and _tree(boards_dir) == {}


def test_the_real_entry_point_never_echoes_the_saved_file(tmp_path):
    """The installed `vec` entry point (Typer's own exception printing, no CliRunner), on a refusal and a success."""
    import subprocess
    import sys

    cases = {"page": (f"<html>{TAG}</html>".encode(), 1), "wrapped": (json.dumps({"b": _final_phase_index()}).encode(),
                                                                       1), "good": (_final_phase_index(), 0)}
    for name, (raw, code) in cases.items():
        src, boards_dir = _save(tmp_path / name / "saved", raw), tmp_path / name / "boards"
        boards_dir.mkdir()
        run = subprocess.run([sys.executable, "-c", "from vec.cli import app; app()", "boards", "update",
                              "--from-file", str(src), "--boards-dir", str(boards_dir)], capture_output=True, text=True,
                             timeout=120)
        assert run.returncode == code, (name, run.stdout, run.stderr)
        assert TAG not in run.stdout + run.stderr and "Traceback" not in run.stderr
        assert all(TAG.encode() not in data for data in _tree(boards_dir).values())


def test_tracebacks_never_show_local_variables():
    from vec import cli

    assert cli.app.pretty_exceptions_show_locals is False  # a local could hold the saved index or the team token


def test_file_intake_refuses_a_file_inside_the_repository_before_reading_it(tmp_path, monkeypatch):
    repo = tmp_path / "repo"
    monkeypatch.setattr("vec.boards.REPO_ROOT", repo)
    inside = _save(repo / "downloads", b"not json at all")
    link = tmp_path / "outside" / "index.json"
    link.parent.mkdir()
    link.symlink_to(inside)
    for src in (inside, link):
        with pytest.raises(IntakeError, match="inside a git work tree.*it was not read"):
            update_boards_from_file(src, tmp_path / "boards")


def test_file_intake_takes_gene_lists_from_boards_by_published_sha256_only(tmp_path):
    boards_dir = tmp_path / "boards"
    boards_dir.mkdir()
    (boards_dir / "T1__val.genes.txt").write_bytes(GENES)
    src = _save(tmp_path / "saved", _final_phase_index(), gene_files=())  # no gene list saved
    s = update_boards_from_file(src, boards_dir, fetch=lambda url: b"")
    assert s["genes_reused"] == ["T1__val.genes.txt", "T1__test.genes.txt = T1__val.genes.txt"]
    assert (boards_dir / "T1__test.genes.txt").read_bytes() == GENES

    other = b"Ctnnb1\nAxin2\n"
    for case, index, gene_files, said in (
            ("no file has the hash", _final_phase_index(genes=other), (), "no file in boards/ has its published"),
            ("the saved list is wrong", _final_phase_index(genes=other), ("T1__test.genes.txt",), "does not match"),
            ("no hash and not saved", _final_phase_index(test_sha=False), (), "no sha256 published")):
        fresh = tmp_path / case
        fresh.mkdir()
        (fresh / "T1__val.genes.txt").write_bytes(GENES)
        src = _save(tmp_path / f"saved-{case}", index, gene_files=gene_files)
        with pytest.raises(IntakeError, match=said):
            update_boards_from_file(src, fresh, fetch=lambda url: b"")
        assert _tree(fresh) == {"T1__val.genes.txt": GENES}


def test_file_intake_leaves_the_composition_alone_unless_one_is_saved(tmp_path):
    boards_dir = tmp_path / "boards"
    boards_dir.mkdir()
    (boards_dir / "t1_composition.json").write_text('{"stages": {"E8.5": {}}}\n')
    asked = []
    s = update_boards_from_file(_save(tmp_path / "saved", _final_phase_index()), boards_dir,
                                fetch=lambda url: asked.append(url) or b"")
    assert not s["t1_composition_supplied"] and not s["t1_composition_written"] and asked == []
    assert (boards_dir / "t1_composition.json").read_text() == '{"stages": {"E8.5": {}}}\n'


def test_delete_source_removes_the_saved_index_and_composition_after_success_only(tmp_path):
    boards_dir = tmp_path / "boards"
    boards_dir.mkdir()
    comp = {"stages": {"E8.5": {}, "E9.5": {}}}
    src = _save(tmp_path / "bad", _final_phase_index(), comp=comp, gene_files=("T1__val.genes.txt",))
    other = _final_phase_index(genes=b"Ctnnb1\n")
    src.write_text(json.dumps(other))
    res = _intake_cli(src, boards_dir, "--delete-source")
    assert res.exit_code == 1 and src.exists() and (src.parent / "t1_composition.json").exists()

    src = _save(tmp_path / "good", _final_phase_index(), comp=comp)
    res = _intake_cli(src, boards_dir, "--delete-source")
    assert res.exit_code == 0, res.output
    assert not src.exists() and not (src.parent / "t1_composition.json").exists()
    assert (src.parent / "T1__test.genes.txt").exists() and "deleted the saved file(s)" in res.output
    assert (boards_dir / "t1_composition.json").exists() and TAG not in res.output


def test_delete_source_needs_from_file():
    from typer.testing import CliRunner

    from vec import cli

    res = CliRunner().invoke(cli.app, ["boards", "update", "--delete-source"])
    assert res.exit_code == 2 and "goes with --from-file" in res.output


def test_a_saved_gene_list_that_is_the_index_is_refused(tmp_path):
    """Found by the 2026-10-28 review: with no published sha256, whatever sat under the gene-list name was copied."""
    for how in ("copy", "symlink"):
        folder = tmp_path / how
        src = _save(folder, _final_phase_index(test_sha=False), gene_files=("T1__val.genes.txt",))
        if how == "copy":
            (folder / "T1__test.genes.txt").write_bytes(src.read_bytes())
        else:
            (folder / "T1__test.genes.txt").symlink_to(src)
        boards_dir = tmp_path / f"boards-{how}"
        boards_dir.mkdir()
        res = _intake_cli(src, boards_dir)
        out = " ".join(res.output.split())
        assert res.exit_code == 1 and ("not a plain gene list" in out or "the saved index itself" in out), out
        assert TAG not in res.output and _tree(boards_dir) == {}


def test_a_composition_linked_to_the_index_is_refused(tmp_path):
    src = _save(tmp_path / "saved", _final_phase_index())
    (src.parent / "t1_composition.json").symlink_to(src)
    boards_dir = tmp_path / "boards"
    boards_dir.mkdir()
    res = _intake_cli(src, boards_dir, "--delete-source")
    assert res.exit_code == 1 and "the saved index itself" in " ".join(res.output.split())
    assert TAG not in res.output and _tree(boards_dir) == {} and src.exists()


def test_the_online_refresh_refuses_a_composition_that_is_not_one(tmp_path):
    """The same review finding on the online route: `comp.get("stages", {})` used to wave such a file through."""
    index = _final_phase_index()
    for comp in (index, {"E8.5": {}, "E12.5": {"cm": SENTINEL}}, {"stages": {"E8.5": {}}, "x": [SENTINEL]}):
        files = {"index.json": json.dumps(index).encode(), "t1_composition.json": json.dumps(comp).encode(),
                 "T1__val.genes.txt": GENES, "T1__test.genes.txt": GENES}
        boards_dir = tmp_path / f"b{len(list(tmp_path.iterdir()))}"
        boards_dir.mkdir()
        with pytest.raises(IntakeError):
            update_boards(boards_dir, fetch=lambda url, f=files: f[url.rsplit("/", 1)[1]])
        assert _tree(boards_dir) == {}


def test_two_boards_naming_one_gene_list_must_agree(tmp_path):
    index = _final_phase_index(genes=b"Ctnnb1\nAxin2\n")
    index["T1:test"]["genes_file"] = "T1__val.genes.txt"  # same name, different published sha256
    boards_dir = tmp_path / "boards"
    boards_dir.mkdir()
    (boards_dir / "old.genes.txt").write_bytes(b"Ctnnb1\nAxin2\n")
    src = _save(tmp_path / "saved", index, gene_files=("T1__val.genes.txt",))
    with pytest.raises(IntakeError, match="another board names T1__val.genes.txt"):
        update_boards_from_file(src, boards_dir, fetch=lambda url: b"")
    assert _tree(boards_dir) == {"old.genes.txt": b"Ctnnb1\nAxin2\n"}


def test_a_file_inside_any_git_work_tree_is_refused_before_reading(tmp_path):
    """Found by the 2026-10-28 review: on this machine the repository's parent folder is a git work tree too."""
    tree = tmp_path / "home"
    (tree / ".git").mkdir(parents=True)
    src = _save(tree / "vec-panels", b"not json at all")
    with pytest.raises(IntakeError, match="inside a git work tree"):
        update_boards_from_file(src, tmp_path / "boards")


def test_an_unreadable_phase_keeps_the_saved_files_for_the_rerun(tmp_path, monkeypatch):
    monkeypatch.setattr("vec.boards.final_phase_state", lambda fetch: "unreadable")
    src = _save(tmp_path / "saved", _final_phase_index(), comp={"stages": {"E8.5": {}, "E9.5": {}, "E10.5": {}}})
    boards_dir = tmp_path / "boards"
    boards_dir.mkdir()
    res = _intake_cli(src, boards_dir, "--delete-source")
    out = " ".join(res.output.split())
    assert res.exit_code == 0 and "kept for that re-run" in out and "--from-file" in out and "deleted" not in out
    assert src.exists() and (src.parent / "t1_composition.json").exists()
    assert not (boards_dir / "t1_composition.json").exists() and (boards_dir / "index.json").exists()


def test_a_failed_delete_is_reported_with_its_own_exit_code(tmp_path):
    src = _save(tmp_path / "saved", _final_phase_index())
    boards_dir = tmp_path / "boards"
    boards_dir.mkdir()
    src.parent.chmod(0o500)  # the folder cannot be written, so its files cannot be deleted
    try:
        res = _intake_cli(src, boards_dir, "--delete-source")
    finally:
        src.parent.chmod(0o700)
    out = " ".join(res.output.split())
    assert res.exit_code == 3 and "boards/ was updated" in out and "could not be deleted" in out, out
    assert src.exists() and (boards_dir / "index.json").exists() and TAG not in res.output


def test_a_failed_write_replaces_nothing(tmp_path):
    boards_dir = tmp_path / "boards"
    boards_dir.mkdir()
    (boards_dir / "index.json").write_text("{}\n")
    (boards_dir / ".T1__test.genes.txt.intake-tmp").mkdir()  # the second temporary write fails
    res = _intake_cli(_save(tmp_path / "saved", _final_phase_index()), boards_dir)
    assert res.exit_code == 1 and "the intake stopped" in res.output and TAG not in res.output
    assert {p.name for p in boards_dir.iterdir()} == {"index.json", ".T1__test.genes.txt.intake-tmp"}
    assert (boards_dir / "index.json").read_text() == "{}\n"


def test_other_copies_left_in_the_folder_are_named(tmp_path):
    src = _save(tmp_path / "saved", _final_phase_index())
    (src.parent / "index (1).json").write_text(json.dumps(_final_phase_index()))
    boards_dir = tmp_path / "boards"
    boards_dir.mkdir()
    res = _intake_cli(src, boards_dir, "--delete-source")
    out = " ".join(res.output.split())
    assert res.exit_code == 0 and "may be copies of the saved index" in out and "index (1).json" in out
    assert not src.exists() and (src.parent / "index (1).json").exists() and TAG not in res.output


def test_online_403_reports_the_block_kind_from_headers_only(monkeypatch):
    """The 403 message names the kind of block (a bot challenge or a firewall rule) from response headers; the body,
    which on another day could be the board file itself, is never read."""
    from email.message import Message
    from urllib.error import HTTPError

    from typer.testing import CliRunner

    from vec import boards as vec_boards
    from vec import cli

    class Body:
        def read(self, *a):
            raise AssertionError("the 403 body must not be read")

    hdr = Message()
    hdr["server"], hdr["cf-mitigated"], hdr["content-type"] = "cloudflare", "challenge", "text/html"

    def blocked(**kw):
        raise HTTPError("https://example.invalid/panels/index.json", 403, "Forbidden", hdr, Body())

    monkeypatch.setattr(vec_boards, "update_boards", blocked)
    res = CliRunner().invoke(cli.app, ["boards", "update"])
    out = " ".join(res.output.split())
    assert res.exit_code == 1 and "nothing was written" in out
    assert "'cf-mitigated': 'challenge'" in out and "'server': 'cloudflare'" in out

    class Plain:  # a plain-text Cloudflare error page: only its numeric code is read and shown
        def read(self, n=-1):
            assert 0 < n <= 200
            return b"error code: 1010 and nothing else is shown"

    plain = Message()
    plain["server"], plain["content-type"] = "cloudflare", "text/plain; charset=UTF-8"

    def blocked_plain(**kw):
        raise HTTPError("https://example.invalid/panels/index.json", 403, "Forbidden", plain, Plain())

    monkeypatch.setattr(vec_boards, "update_boards", blocked_plain)
    out = " ".join(CliRunner().invoke(cli.app, ["boards", "update"]).output.split())
    assert "'cloudflare_error_code': '1010'" in out and "nothing else is shown" not in out


# ---------------------------------------------------------------------------------------------------------------------
# s1 build review (2026-11-01): the site's gene-list digest, null limits, required keys, hidden test-board counts.

def _site_digest(genes: bytes) -> str:
    """The digest the site publishes: the first 16 hex digits of sha256 over the symbols joined by newlines."""
    symbols = [g.strip() for g in genes.decode().splitlines() if g.strip()]
    return hashlib.sha256("\n".join(symbols).encode()).hexdigest()[:16]


def test_every_committed_gene_list_matches_its_published_digest():
    """The committed boards/index.json is the site's own file (2026-10-18): each list must verify, else the intake
    refuses every board."""
    from pathlib import Path

    from vec.boards import BOARDS_DIR, _genes_match

    index = json.loads((Path(BOARDS_DIR) / "index.json").read_text())
    for key, spec in index.items():
        data = (Path(BOARDS_DIR) / spec["genes_file"]).read_bytes()
        assert _genes_match(data, spec["genes_sha256"]), key
        assert not _genes_match(data + b"ExtraGene\n", spec["genes_sha256"]), key


@pytest.mark.parametrize("served", [GENES, GENES.rstrip(b"\n"), GENES.replace(b"\n", b"\r\n")])
def test_intakes_accept_the_site_digest_however_the_list_ends(tmp_path, served):
    index = {"T1:val": _spec("T1:val", "val"), "T1:test": _spec("T1:test", "test")}
    for s in index.values():
        s["genes_sha256"] = _site_digest(GENES)
    files = {"index.json": json.dumps(index).encode(), "T1__val.genes.txt": served, "T1__test.genes.txt": served,
             "t1_composition.json": json.dumps({"stages": {"E8.5": {}, "E9.5": {}}}).encode()}
    online = tmp_path / "online"
    online.mkdir()
    update_boards(online, fetch=lambda url: files[url.rsplit("/", 1)[1]])
    assert load_boards(online)["T1:test"].genes == ["Gata4", "Nkx2-5"]
    offline = tmp_path / "offline"
    offline.mkdir()
    src = _save(tmp_path / "saved", index, genes=served)
    update_boards_from_file(src, offline, fetch=lambda url: (_ for _ in ()).throw(OSError("offline")))
    assert load_boards(offline)["T1:test"].genes == ["Gata4", "Nkx2-5"]


def test_intake_reuses_a_boards_list_by_the_site_digest(tmp_path):
    boards_dir = tmp_path / "boards"
    boards_dir.mkdir()
    (boards_dir / "T1__val.genes.txt").write_bytes(GENES)
    index = {"T1:val": _spec("T1:val", "val"), "T1:test": _spec("T1:test", "test")}
    for s in index.values():
        s["genes_sha256"] = _site_digest(GENES)
    src = _save(tmp_path / "saved", index, gene_files=())
    summary = update_boards_from_file(src, boards_dir, fetch=lambda url: (_ for _ in ()).throw(OSError("offline")))
    assert sorted(summary["genes_reused"]) == ["T1__test.genes.txt = T1__val.genes.txt", "T1__val.genes.txt"]


def test_intakes_refuse_a_list_the_site_digest_does_not_name(tmp_path):
    index = {"T1:val": _spec("T1:val", "val")}
    index["T1:val"]["genes_sha256"] = _site_digest(b"Gata4\nTbx5\n")
    files = {"index.json": json.dumps(index).encode(), "T1__val.genes.txt": GENES,
             "t1_composition.json": json.dumps({"stages": {"E8.5": {}}}).encode()}
    with pytest.raises(IntakeError, match="sha256"):
        update_boards(tmp_path, fetch=lambda url: files[url.rsplit("/", 1)[1]])
    short = {"T1:val": _spec("T1:val", "val")}
    short["T1:val"]["genes_sha256"] = _site_digest(GENES)[:8]  # too short to name a list
    files["index.json"] = json.dumps(short).encode()
    with pytest.raises(IntakeError, match="sha256"):
        update_boards(tmp_path, fetch=lambda url: files[url.rsplit("/", 1)[1]])


def test_a_null_min_cells_loads_as_the_default_and_required_keys_load(tmp_path):
    spec = _spec("T1:test", "test", anchors=False)
    spec.update(min_cells=None, obs_required=[], obsm_required=["spatial_3D"])
    clean, _ = sanitise_index({"T1:test": spec})
    (tmp_path / "index.json").write_text(json.dumps(clean))
    (tmp_path / "T1__test.genes.txt").write_bytes(GENES)
    b = load_boards(tmp_path)["T1:test"]
    assert b.min_cells == 1000 and b.obs_required == [] and b.obsm_required == ["spatial_3D"]


# The published baselines (docs/CHALLENGE.md, "Metric panel"; the captain's guiding principle, 2026-11-04): copy_last on
# T1 and T2, wt_identity on T3. Until 2026-11-04 a spec without floor_model read as copy_last on every task (D1,
# ae6d6b3, lists copy_last for T3:bcatenin).

@pytest.mark.parametrize("key, published, floor", [
    ("T3:bcatenin", "omitted", "wt_identity"), ("T3:bcatenin", None, "wt_identity"),
    ("T1:test", "omitted", "copy_last"), ("T2:heart:test_extrap", None, "copy_last"),
    ("T2:embryo:test_interp", "omitted", "copy_last"), ("T3:gata4", "wt_identity", "wt_identity"),
    ("T3:bcatenin", "copy_last", "copy_last")])
def test_a_spec_without_a_floor_model_gets_its_tasks_published_floor(tmp_path, key, published, floor):
    spec = _spec(key, "test", anchors=False)
    if published == "omitted":
        del spec["floor_model"]
    else:
        spec["floor_model"] = published  # a published value is kept as published
    clean, _ = sanitise_index({key: spec})
    (tmp_path / "index.json").write_text(json.dumps(clean))
    (tmp_path / spec["genes_file"]).write_bytes(GENES)
    assert load_boards(tmp_path)[key].floor_model == floor


def test_every_task_has_a_published_floor():
    from vec.boards import GROUPS, TASK_FLOOR

    assert set(TASK_FLOOR) == set(GROUPS) and TASK_FLOOR["T3"] == "wt_identity"


@pytest.mark.parametrize("key, split, hidden", [
    ("T1:test", "test", True), ("T2:heart:test_extrap", None, True), ("T3:bcatenin", "test", True),
    ("T1:other", "test_x", True), ("T9:new", "weird", True), ("T1:val", "val", False),
    ("pseudo:T1/E9.5->E10.5", "pseudo", False), ("T1:synthetic", None, False)])
def test_counts_are_hidden_on_test_boards_only(key, split, hidden):
    from vec.boards import Board

    assert Board(key=key, task="T1", genes=["A"], needs_coords=False, split=split).counts_hidden is hidden
