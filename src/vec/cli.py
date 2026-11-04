"""`vec` command line: board listing, free format checks, local scoring and reference baselines."""

from __future__ import annotations

import csv
import json
import re
import shlex
import subprocess
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import typer
from rich.console import Console
from rich.markup import escape
from rich.table import Table
from typer.core import TyperGroup

from vec import baselines, spatial
from vec import data as vec_data
from vec import submit as vec_submit
from vec.boards import REPO_ROOT, get_board, load_boards
from vec.evaluate import Stage, score_board_seeds
from vec.io import check_submission, make_submission, read_stage, write_submission

app = typer.Typer(no_args_is_help=True, add_completion=False)
console = Console(emoji=False, highlight=False)
REGISTRY = REPO_ROOT / "results" / "registry.csv"


def _git_commit() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=REPO_ROOT, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def boards():
    """List the known boards and their constraints."""
    t = Table("board", "task", "genes", "coords", "cells", "floor")
    for key, b in load_boards().items():  # site text: table cells are markup, so every one is escaped
        cells = "hidden (test board)" if b.counts_hidden else f"{b.min_cells}-{b.max_cells}"
        t.add_row(*map(escape, [key, b.task, str(b.n_genes), str(b.needs_coords), cells, str(b.floor_model or "-")]))
    console.print(t)


boards_app = typer.Typer()
app.add_typer(boards_app, name="boards")


@boards_app.callback(invoke_without_command=True)
def boards_main(ctx: typer.Context):
    """List the known boards and their constraints (`vec boards update` refreshes them)."""
    if ctx.invoked_subcommand is None:
        boards()


@boards_app.command("update")
def boards_update(
    from_file: Path | None = typer.Option(
        None, "--from-file", help="offline intake of a board index saved from the site with the browser (Save "
        "link as, never opened), kept outside every git work tree; not cleared for use until the user decides "
        "(docs/DATA_POLICY.md); gene lists come from the same folder or from boards/ files with the published "
        "sha256"),
    delete_source: bool = typer.Option(False, "--delete-source", help="with --from-file: delete the saved index "
                                       "(and a saved t1_composition.json) after a successful intake"),
    boards_dir: Path | None = typer.Option(None, "--boards-dir", hidden=True),
):
    """Refresh boards/ from the public panel files, or from a saved board index (--from-file); test-board anchors
    are stripped in memory and never printed."""
    from urllib.error import HTTPError

    from vec import boards as vec_boards

    if delete_source and from_file is None:
        console.print("--delete-source goes with --from-file; nothing was done")
        raise typer.Exit(2)
    where = {} if boards_dir is None else {"boards_dir": boards_dir}
    try:
        if from_file is None:
            s = vec_boards.update_boards(**where)
        else:
            s = vec_boards.update_boards_from_file(from_file, delete_source=delete_source, **where)
    except vec_boards.IntakeError as e:  # its message names boards, fields and files only
        console.print(f"refused: {escape(str(e))}.")
        console.print("Nothing was written.")
        raise typer.Exit(1) from None
    except HTTPError as e:
        if from_file is not None or e.code != 403:
            raise
        # which kind of 403 it is (a bot challenge or browser-signature block a browser passes, or a firewall rule),
        # from the headers and, for a short plain-text Cloudflare error page only, its numeric code; nothing else of
        # the body is read or shown
        kind = {k: e.headers.get(k) for k in ("server", "cf-mitigated", "content-type")
                if e.headers is not None and e.headers.get(k)}
        if str(kind.get("content-type", "")).startswith("text/plain"):
            try:
                code = re.search(r"error code: (\d{3,5})", e.read(200).decode("ascii", "replace"))
            except Exception:  # noqa: BLE001 (diagnostic only)
                code = None
            if code:
                kind["cloudflare_error_code"] = code.group(1)
        console.print("the panel files answered HTTP 403; nothing was written. Not worked around: a copy saved in "
                      "the browser for --from-file would be a workaround too and needs the user's decision first "
                      f"(docs/DATA_POLICY.md). Response headers: {escape(str(kind or 'none'))}.")
        raise typer.Exit(1) from None
    except Exception as e:  # noqa: BLE001 (offline intake only: an error message could quote the saved file)
        if from_file is None:
            raise
        said = (escape(str(e)) if isinstance(e, OSError)
                else f"{type(e).__name__}; its message is withheld because it could quote the saved file")
        console.print(f"the intake stopped ({said}). Nothing in boards/ was replaced unless the error names a file "
                      "there. The saved index was kept and still holds the test anchors: do not open it to check.")
        raise typer.Exit(1) from None
    if from_file is not None:
        console.print(f"source: saved board index {escape(s['source'])} (offline intake)")
        if s["genes_reused"]:
            console.print("gene lists taken from boards/ by their published sha256: "
                          f"{escape(', '.join(s['genes_reused']))}")
    console.print(f"boards ({len(s['boards'])}): {escape(', '.join(s['boards']))}")  # keys from the site: never markup
    console.print(f"reduced to their submission spec ({len(s['reduced'])}): {escape(', '.join(s['reduced'])) or '-'}")
    for key, names in (s.get("dropped_fields") or {}).items():  # field names only (decision point D1), never values
        console.print(f"  {escape(key)}: fields not kept: {escape(', '.join(names))}")
    if from_file is not None and not s["t1_composition_supplied"]:
        console.print("T1 composition not supplied: boards/t1_composition.json left as it was")
    elif s["t1_composition_written"]:
        console.print("T1 composition written")
    else:
        why = {"closed": "; E10.5 counts once the final phase is open",
               "unreadable": "; the public phase endpoint could not be read, so E10.5 counts as unreleased: re-run"}
        console.print(f"T1 composition not written: it describes {escape(str(s['t1_composition_unreleased_stages']))} "
                      f"unreleased stage(s){why.get(s.get('final_phase'), '')}")
    if from_file is None:
        return
    saved = [s["source"]] + ([str(Path(s["source"]).parent / "t1_composition.json")]
                             if s["t1_composition_supplied"] else [])
    if s["kept_for_rerun"]:
        console.print("the saved file(s) were kept for that re-run and still hold the test anchors (do not open "
                      "them): uv run vec boards update --from-file "
                      f"{escape(shlex.quote(s['source']))} --delete-source")
    elif s["deleted"]:
        console.print(f"deleted the saved file(s), which held the test anchors: {escape(', '.join(s['deleted']))}")
    elif not s["undeletable"]:
        console.print("the saved file(s) still hold the test anchors: delete them now, e.g. rm -- "
                      f"{escape(' '.join(map(shlex.quote, saved)))}")
    if s["possible_copies"]:
        console.print("other files in that folder may be copies of the saved index (not opened, not deleted): "
                      f"{escape(', '.join(s['possible_copies']))}; delete them too")
    if s["undeletable"]:
        console.print(f"boards/ was updated, but these could not be deleted: {escape(', '.join(s['undeletable']))}. "
                      "They still hold the test anchors: delete them by hand.")
        raise typer.Exit(3)


@app.command()
def check(pred: Path, board: str = typer.Option(..., help="board key, e.g. T1:val")):
    """Free format check of a prediction file (no scoring)."""
    b = get_board(board)
    rep = check_submission(pred, b)
    for w in rep.warnings:  # the messages quote the site's gene panel and the file's contents: never markup
        console.print(f"[yellow]warning[/] {escape(w)}")
    for e in rep.errors:
        console.print(f"[red]error[/] {escape(e)}")
    count = ("max_cells" if rep.n_cells == b.max_cells else "count not shown") if b.counts_hidden else rep.n_cells
    console.print(f"{'[green]OK' if rep.ok else '[red]REJECTED'}[/] ({count} cells)")
    raise typer.Exit(0 if rep.ok else 1)


@app.command("eval")
def evaluate(
    pred: Path = typer.Option(..., help="prediction .h5ad"),
    target: Path = typer.Option(..., help="full held-out stage (e.g. released validation answers)"),
    ref: Path = typer.Option(..., help="DE reference: preceding stage (T1/T2) or matched wild type (T3)"),
    board: str = typer.Option(..., help="board key"),
    seeds: str = typer.Option("0,1,2", help="comma-separated subsampling seeds"),
    anchors: str = typer.Option("computed", help="computed | published"),
    gate: str = typer.Option("floor_normalised", help="T2 shape gate: none | multiplicative | floor_normalised"),
    frac: float = typer.Option(0.1, help="working-copy fraction (0.1 = server; 1.0 for small files)"),
    protocol: str = typer.Option("bake", help="bake (the server's draws, calibrated 2026-10-31) | emulator"),
    build: str = typer.Option("main", help="scorer build: main (compare candidates) | prefix (calibration only)"),
    method: str = typer.Option("", help="method name for the registry"),
    out: Path = typer.Option(None, help="write full JSON scores here"),
    register: bool = typer.Option(False, help="append a row to results/registry.csv"),
):
    """Score a prediction with the local board emulator (by default the calibrated bake protocol). Validation boards
    only: no test target exists locally, and a test board's counts are never written out."""
    from vec import harness, registry
    from vec.boards import VALIDATION_BOARDS

    if board not in VALIDATION_BOARDS:
        raise typer.BadParameter(f"{board} is not a validation board; vec eval scores validation boards only")
    if register and not out:
        raise typer.BadParameter("--register needs --out, so the row points at the full scores")
    b = get_board(board)
    stages = [Stage.from_anndata(read_stage(p, b.genes), b.genes, need_coords=b.needs_coords)
              for p in (pred, target, ref)]
    res = score_board_seeds(b, *stages, seeds=tuple(int(s) for s in seeds.split(",")), anchors=anchors, gate=gate,
                             frac=frac, protocol=protocol, build=build)
    console.print(f"{board}: score {res['score_mean']:.2f} +/- {res['score_sd']:.2f} over seeds {seeds}")
    first = res["per_seed"][0]
    for name, s in first["group_skill"].items():
        console.print(f"  {name}: {s:.3f} (seed {first['seed']})")
    if first["extras"].get("inverted_anchors"):
        console.print(f"[yellow]warning[/] ceiling no better than floor for {first['extras']['inverted_anchors']}: "
                      "scores on this board are not meaningful")
    if out:
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(res, indent=2, default=float))
    if register:
        registry.append([{
            "utc_time": datetime.now(UTC).isoformat(timespec="seconds"), "commit": harness._git("rev-parse", "HEAD"),
            "tree_clean": harness._git("status", "--porcelain", "--untracked-files=no") == "",
            "source": "vec eval", "run_id": "", "test_board": "", "pseudo_task": "", "scoring_board": board,
            "candidate": method or pred.stem,
            "params": json.dumps({k: {"path": str(v), "sha256": harness._sha256(v)}
                                  for k, v in (("pred", pred), ("target", target), ("ref", ref))}),
            "seeds": seeds, "score_mean": f"{res['score_mean']:.3f}", "score_sd": f"{res['score_sd']:.3f}",
            "protocol": protocol, "build": build, "anchors": anchors, "gate": gate, "frac": frac,
            "scores_json": str(out) if out else ""}], path=REGISTRY)


@app.command()
def baseline(
    name: str = typer.Argument(..., help="copy_last | pseudobulk_shift | wt_identity | shift_transfer | gene_ko"),
    board: str = typer.Option(..., help="board key (defines the gene panel)"),
    last: Path = typer.Option(None, help="last observed stage (copy_last, pseudobulk_shift)"),
    prev: Path = typer.Option(None, help="stage before last (pseudobulk_shift)"),
    wt: Path = typer.Option(None, help="matched wild type of the target (wt_identity, shift_transfer, gene_ko)"),
    ko_train: Path = typer.Option(None, help="observed knockout (shift_transfer, gene_ko)"),
    wt_train: Path = typer.Option(None, help="wild type matched to the observed knockout"),
    target_gene: str = typer.Option("", help="knocked-out gene (gene_ko)"),
    n_cells: int = typer.Option(0, help="keep this many random rows (0 = all), e.g. to meet the board's cell limits"),
    seed: int = typer.Option(0, help="seed of the --n-cells draw"),
    out: Path = typer.Option(..., help="output .h5ad"),
):
    """Write a reference-baseline prediction for a board."""
    b = get_board(board)

    def load(p: Path) -> Stage:
        return Stage.from_anndata(read_stage(p, b.genes), b.genes, need_coords=b.needs_coords)

    if name == "copy_last":
        pred = baselines.copy_last(load(last))
    elif name == "pseudobulk_shift":
        pred = baselines.pseudobulk_shift(load(prev), load(last))
    elif name == "wt_identity":
        pred = baselines.wt_identity(load(wt))
    elif name == "shift_transfer":
        pred = baselines.shift_transfer(load(ko_train), load(wt_train), load(wt))
    elif name == "gene_ko":
        pred = baselines.gene_ko(load(ko_train), load(wt_train), load(wt), b.genes, target_gene)
    else:
        raise typer.BadParameter(f"unknown baseline {name!r}")
    if n_cells:
        if not 0 < n_cells <= pred.n:
            raise typer.BadParameter(f"--n-cells must be between 1 and {pred.n} (the baseline's cells)")
        pred = pred.take(np.random.default_rng(seed).choice(pred.n, n_cells, replace=False))
    path = write_submission(make_submission(pred.X, b.genes, pred.coords), out)
    console.print(f"wrote {path} ({pred.n} cells)")


data_app = typer.Typer(no_args_is_help=True, help="Competition data: list and download released files.")


def _authed(call):
    try:
        return call()
    except vec_data.AuthError as e:  # carries the server's error text (token redacted): never markup
        console.print(f"[red]auth error[/] {escape(str(e))}")
        raise typer.Exit(2) from None
app.add_typer(data_app, name="data")


@data_app.command("manifest")
def data_manifest(token_file: Path = typer.Option(vec_data.DEFAULT_TOKEN_FILE, help="bearer token file (mode 600)")):
    """List the files your account can download: the training files, then the validation answers once unlocked."""
    files, state = _authed(lambda: vec_data.released_files(vec_data.load_token(token_file)))
    for source in ("train", vec_data.VALIDATION):
        part = [f for f in files if f.source == source]
        if source == vec_data.VALIDATION:
            _print_validation_state(state)
            if not part:
                continue
        t = Table("key", "task", "role", "stage", "modality", "MB", "other fields")
        for f in part:  # server text
            mb = f"{f.bytes / 1e6:.0f}" if f.bytes else "?"
            t.add_row(*map(escape, [f.key, f.task, f.role, f.stage, f.modality, mb,
                                    json.dumps(f.extra) if f.extra else "-"]))
        console.print(t)
        console.print(f"{len(part)} {source} files, {sum(f.bytes for f in part) / 1e9:.2f} GB")


def _print_validation_state(state: dict):
    shown = {k: v for k, v in state.items() if not isinstance(v, (dict, list))}  # server text
    if state.get("unlocked") and state.get("available"):
        console.print(f"validation release: unlocked {escape(json.dumps(shown))}")
    else:
        console.print(f"validation release: not listed {escape(json.dumps(shown))}; a person unlocks it on the "
                      "site's Data page, which permanently ends P2 submissions for the whole team (never this tool)")


@data_app.command("fetch")
def data_fetch(
    dest: str = typer.Option("~/v/virtual-embryo-challenge/data/raw", help="destination directory"),
    host: str = typer.Option("PrakashDGX_H2", help="ssh host to download onto; empty string = this machine"),
    match: str = typer.Option("", help="only keys containing this substring"),
    token_file: Path = typer.Option(vec_data.DEFAULT_TOKEN_FILE, help="bearer token file (mode 600)"),
):
    """Download released files, validation answers included once unlocked (one signed link at a time); already
    complete files are skipped."""
    token = _authed(lambda: vec_data.load_token(token_file))
    listed, state = _authed(lambda: vec_data.released_files(token))
    _print_validation_state(state)
    files = [f for f in listed if match in f.key]
    for i, f in enumerate(files, 1):
        status = _authed(lambda f=f: vec_data.fetch_one(f, token, dest, host or None))
        size = f"{f.bytes / 1e6:.0f} MB" if f.bytes else "size not listed"
        console.print(f"[{i}/{len(files)}] {status:10s} {escape(f.key)} ({size}, {f.source})")


@data_app.command("inventory")
def data_inventory(
    root: Path = typer.Option(Path("data/raw"), help="directory holding the downloaded h5ad files"),
    out: Path = typer.Option(Path("runs/inventory.json"), help="write the full JSON summary here"),
):
    """Summarise every h5ad under ROOT (run on the DGX, where the data lives)."""
    files = sorted(root.rglob("*.h5ad"))
    if not files:
        console.print(f"[red]no .h5ad files under {root}[/]")
        raise typer.Exit(1)
    rows = []
    t = Table("file", "cells", "genes", "nnz", "expm1 sum", "layers", "obsm")
    for f in files:  # file names (from the server's keys) and the names inside the files: never markup
        s = vec_data.inventory(f)
        rows.append(s)
        t.add_row(*map(escape, [str(f.relative_to(root)), str(s["n_obs"]), str(s["n_vars"]),
                                f"{s['nonzero_fraction']:.3f}", f"{s['expm1_rowsum_median']:.0f}",
                                ",".join(s["layers"]) or "-", ",".join(s["obsm_keys"]) or "-"]))
    console.print(t)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rows, indent=2))
    console.print(f"wrote {out}")


def _cell(v) -> str:
    return "-" if v in ("", None) else str(v)


@data_app.command("qc")
def data_qc(
    files: list[Path] = typer.Argument(..., help="h5ad files to check (run on the DGX, where the data lives)"),
    out: Path = typer.Option(Path("runs/qc.json"), help="write the full JSON report here"),
    csv_out: Path = typer.Option(None, "--csv", help="also write one summary row per file here"),
    xy_tol: float = typer.Option(0.01, help="duplicate records: largest |dx| and |dy| between copies (inf: any)"),
    min_library: float = typer.Option(20, help="duplicate records: never flag rows with fewer total counts"),
    max_level_gap: int = typer.Option(1, help="duplicate records: largest z-level distance between copies"),
):
    """Intake QC per file: z levels, step, gaps, sections, row order, duplicated records, coordinates (aggregates)."""
    import anndata as ad

    reports, rows = [], []
    t = Table("file", "cells", "dup names", "non-int rows", "z raw/levels", "step", "gaps", "sections", "rho/runs",
              "dup records", "identical rows", "centred")
    for f in files:
        qc = spatial.qc_stage(ad.read_h5ad(f), xy_tol=xy_tol, min_library=min_library, max_level_gap=max_level_gap)
        reports.append({"file": str(f), **qc})
        r = spatial.summary_row(str(f), qc)
        rows.append(r)
        has_z = r["z_levels"] != ""
        dense = f"dense levels {r['z_dense_levels']}" if r["z_dense_levels"] else ""
        lib = f" ({r['identical_rows_min_library']} >= min lib)" if r["identical_rows_min_library"] != "" else ""
        dups = (f"{r['dup_redundant']} ({r['dup_by_gap'] or 'none'}); any xy {r['dup_any_xy_redundant']}"
                if r["dup_redundant"] != "" else "-")
        t.add_row(*map(escape, [  # file names and the files' labels: never markup
            f"{f.parent.name}/{f.name}", _cell(r["cells"]), _cell(r["obs_names_duplicated"]),
            _cell(r["non_integer_count_rows"]), f"{r['z_raw_values']}/{r['z_levels']}" if has_z else "-",
            _cell(r["z_step"]), _cell(r["z_gaps"]), _cell("; ".join(filter(None, [r["section_labels"], dense]))),
            f"{r['spearman_row_vs_z']}/{r['z_runs']}" if has_z else "-", dups, f"{r['identical_rows']}{lib}",
            _cell(r["centred_on"])]))
    console.print(t)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(reports, indent=2, allow_nan=False))
    console.print(f"wrote {out}")
    if csv_out:
        csv_out.parent.mkdir(parents=True, exist_ok=True)
        with csv_out.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        console.print(f"wrote {csv_out}")


class _SubmitGroup(TyperGroup):
    """`vec submit FILE ...` runs `vec submit upload FILE ...`; `vec submit status` is the read-only view and
    `vec submit resolve` records a person's decision about an upload the ledger cannot settle."""

    def resolve_command(self, ctx, args):
        if args and args[0] not in self.commands:
            return "upload", self.commands["upload"], args
        return super().resolve_command(ctx, args)


submit_app = typer.Typer(cls=_SubmitGroup, no_args_is_help=True,
                         help="Submission client: free format checks by default, guarded official uploads, one ledger "
                              "row per attempt (results/submissions.csv). `vec submit FILE --board KEY` runs upload. "
                              "A file built on the DGX: `vec submit PATH --board KEY --from-dgx`, with PATH relative "
                              f"to {vec_submit.DGX_REPO} on {vec_submit.DGX_HOST}; it is pulled into a private "
                              "temporary copy for that call only, and the key never leaves this machine.")
app.add_typer(submit_app, name="submit")


def _print_plan(plan: vec_submit.Plan) -> None:
    if plan.official:
        console.print(f"[bold red]OFFICIAL upload[/] to {escape(plan.board)}: uses one of the board's "
                      f"{vec_submit.FINAL_PHASE_LIMIT} scored submissions for the whole phase; it cannot be withdrawn")
    else:
        console.print(f"Format check on {escape(plan.board)} (format_only=true): nothing is scored and no attempt is "
                      "spent")
    src = f"{vec_submit.DGX_HOST}:{plan.remote}, pulled for this call" if plan.remote else str(plan.file)
    console.print(f"  file      {escape(src)} (sent as {escape(plan.upload_name)}, "
                  f"{plan.file.stat().st_size / 1e6:.1f} MB)")
    console.print(f"  sha256    {plan.sha256}")
    console.print(f"  fields    {escape(' '.join(f'{k}={v}' for k, v in plan.fields.items()))}")
    if plan.official:
        console.print(f"  key       Idempotency-Key {escape(plan.idempotency_key or '-')}")
    console.print(f"  disclose  {escape(', '.join(plan.disclosure)) or '-'}")
    console.print(f"  commit    {plan.git_commit}")
    console.print(f"  build     {escape(plan.build_commit) or '-'}")
    if plan.remote:
        console.print(f"  origin    {escape(plan.origin)}")
    for w in plan.warnings:
        console.print(f"  [yellow]warning[/] {escape(w)}")
    for g in plan.guards:
        console.print(f"  {'[green]PASS[/]' if g.ok else '[red]FAIL[/]'} {g.name:<19} {escape(g.detail)}")


def _print_result(res: vec_submit.Result) -> None:
    p = res.plan
    detail = f" {escape(res.detail)}" if res.detail else ""
    status, sid = escape(res.server_status), escape(res.server_id or "?")  # server text: never markup
    if not p.official:
        if res.ok:
            console.print(f"[green]server format check passed[/]{detail}")
            src, dgx = (shlex.quote(p.remote), " --from-dgx") if p.remote else (str(p.file), "")
            gate = vec_submit.upload_gate_guard()
            if not gate.ok:  # the captain's upload gate (docs/PLAN.md section 8)
                console.print(f"[yellow]no official upload: {escape(gate.detail)}[/]")
            console.print(f"official upload (only with the captain's go-ahead, once his upload gate passes): vec "
                          f"submit {escape(src)} --board "
                          f"{escape(p.board)}{dgx} --official --confirm {p.sha256[:12]} --disclose "
                          f"{escape(','.join(p.disclosure)) or '<ids from EXTERNAL_DATA.md>'}")
        else:
            console.print(f"[red]server format check: {status}[/]{detail}")
    elif res.replayed:
        console.print(f"already sent: submission {sid} ({status}); nothing was uploaded")
    elif res.ok:
        console.print(f"[green]uploaded[/] submission {sid} ({status}){detail}")
    else:
        console.print(f"[red]upload {status}[/]{detail} (vec submit status shows the server's view)")


@submit_app.command("upload")
def submit_upload(
    file: Path = typer.Argument(..., help="prediction .h5ad (with --from-dgx, its path on the DGX)"),
    board: str = typer.Option(..., help="board key, e.g. T1:val"),
    from_dgx: bool = typer.Option(False, "--from-dgx", help=f"FILE is on {vec_submit.DGX_HOST}, relative to "
                                  f"{vec_submit.DGX_REPO} (or absolute inside it): pulled into a private temporary "
                                  "copy for this call only; the key stays here"),
    official: bool = typer.Option(False, "--official", help="OFFICIAL scored upload; also needs --confirm"),
    confirm: str = typer.Option("", help="first 12 hex characters of the file's sha256 (official uploads)"),
    disclose: str = typer.Option("", help="EXTERNAL_DATA.md ids, comma-separated (official uploads need one)"),
    note: str = typer.Option("", help="free-text note for the ledger"),
    model: str = typer.Option("", help="optional method name sent with the upload (shown with the team)"),
    architecture: str = typer.Option("", help="optional architecture name sent with the upload"),
):
    """Format check FILE for a board (default, free), or upload it officially with --official --confirm SHA12; add
    --from-dgx for a file built on the DGX (vec submit PATH --board KEY --from-dgx, also with --official)."""
    try:
        res = vec_submit.submit(file, board, official=official, confirm=confirm or None, disclose=disclose, note=note,
                                model=model or None, architecture=architecture or None, from_dgx=from_dgx,
                                on_plan=_print_plan, echo=lambda m: console.print(escape(m)))
    except vec_submit.SubmitRefused as e:
        console.print(f"[red]refused[/] {escape(str(e))}; nothing was sent")
        raise typer.Exit(1) from None
    except vec_data.AuthError as e:
        console.print(f"[red]auth error[/] {escape(str(e))}")
        raise typer.Exit(2) from None
    _print_result(res)
    raise typer.Exit(0 if res.ok else 1)


@submit_app.command("status")
def submit_status(board: str = typer.Option("", help="only this board key")):
    """Phase, quota and submissions from the server (read-only GETs) next to the ledger, per board."""
    rep = vec_submit.status(board or None)
    if rep.phase:
        ph = rep.phase
        console.print(escape(f"phase {ph.get('phase')} ({ph.get('label', '')}): split {ph.get('split')}, accepts "
                             f"submissions {ph.get('accepts_submissions')}"))
    for n in rep.notes:
        console.print(f"[yellow]note[/] {escape(n)}")
    daily = (rep.quota or {}).get("daily")
    if isinstance(daily, dict):
        console.print("daily: " + escape("; ".join(
            f"{k} scored {v.get('used')}/{v.get('limit')}, uploads {v.get('uploads')}/{v.get('upload_limit')}"
            if isinstance(v, dict) else f"{k} {v}" for k, v in daily.items())))
    t = Table("board", "server quota", "spent", "unlisted (ledger)", "ambiguous", "unresolved", "format ok",
              "ledger rows", "server entries (last 3)")
    for bs in rep.boards:  # table cells are markup too: everything from the server or the ledger is escaped
        b = bs.budget
        entries = ", ".join(f"{e.get('id')} {e.get('status')} {str(e.get('sha256') or '-')[:12]}"
                            for e in bs.entries[-3:])
        files = [" ".join(s[:12] + (f" x{n}" if n > 1 else "") for s, n in sorted(b.unlisted.items())),
                 " ".join(sorted(s[:12] for s in b.ambiguous)),
                 " ".join(sorted({r["sha256"][:12] for r in b.unresolved}))]
        t.add_row(escape(bs.board), f"{b.server_used}/{b.server_limit}" if b.server_used is not None else "unknown",
                  f"{b.spent}/{b.limit}", *(escape(c) or "-" for c in files), str(bs.format_ok), str(bs.ledger_rows),
                  escape(f"{len(bs.entries)}: {entries}") if bs.entries else "-")
    console.print(t)
    for bs in rep.boards:  # what each official upload in the ledger was sent at (commit) and built at (build)
        for r in bs.uploads:
            console.print(escape(f"official {bs.board}: {r.get('sha256', '')[:12]} {r.get('server_id') or '-'} "
                                 f"({r.get('server_status')}), commit {r.get('git_commit') or '-'}, build "
                                 f"{r.get('build_commit') or '-'}"))
    for bs in rep.boards:
        key, pending = escape(bs.board), vec_submit.resolve_commands(bs.budget.unresolved)
        for (_, sha), (kind, cmd) in sorted(pending.items()):
            console.print(f"[yellow]unresolved[/] {key}: a format check of {escape(sha[:12])} came back as "
                          f"{escape(kind)}; nothing is sent to any board until a person has checked the portal and "
                          f"run: {escape(cmd)}")
        for sha in sorted(bs.budget.ambiguous - {s for _, s in pending}):
            console.print(f"[yellow]ambiguous[/] {key}: {escape(sha[:12])} may be stored and counts until the server "
                          f"lists it; if the portal shows it is not: vec submit resolve {escape(sha[:12])} --board "
                          f"{escape(shlex.quote(bs.board))} --as not-stored --note TEXT")
    reasons = sorted({bs.budget.reason for bs in rep.boards if bs.budget.reason})
    for r in reasons:
        console.print(f"[yellow]server budget unknown[/] {escape(r)}: official uploads are refused")
    if reasons and isinstance(rep.quota, dict) and "official_per_board" in rep.quota:
        console.print(f"official_per_board as served: {escape(json.dumps(rep.quota['official_per_board'])[:600])}")


@submit_app.command("resolve")
def submit_resolve(
    sha: str = typer.Argument(..., help="first 12 (or more) hex characters of the file's sha256"),
    board: str = typer.Option(..., help="board key, e.g. T1:test"),
    decision: str = typer.Option(..., "--as", help="not-stored (the portal shows no submission of this file) or "
                                 "acknowledged (a format check created a submission: it keeps counting)"),
    note: str = typer.Option(..., help="what was checked (recorded in the ledger)"),
):
    """Record a person's decision about an upload the ledger cannot settle (local: nothing is sent, no token read)."""
    try:
        row = vec_submit.resolve(sha, board, decision, note)
    except vec_submit.SubmitRefused as e:
        console.print(f"[red]refused[/] {escape(str(e))}; nothing was recorded")
        raise typer.Exit(1) from None
    console.print(escape(f"recorded: {row['sha256'][:12]} on {board} resolved as {decision}; {row['detail']}"))


if __name__ == "__main__":
    app()
