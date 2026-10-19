"""`vec` command line: board listing, free format checks, local scoring and reference baselines."""

from __future__ import annotations

import csv
import json
import shlex
import subprocess
from datetime import UTC, datetime
from pathlib import Path

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
        t.add_row(*map(escape, [key, b.task, str(b.n_genes), str(b.needs_coords), f"{b.min_cells}-{b.max_cells}",
                                str(b.floor_model or "-")]))
    console.print(t)


boards_app = typer.Typer()
app.add_typer(boards_app, name="boards")


@boards_app.callback(invoke_without_command=True)
def boards_main(ctx: typer.Context):
    """List the known boards and their constraints (`vec boards update` refreshes them)."""
    if ctx.invoked_subcommand is None:
        boards()


@boards_app.command("update")
def boards_update():
    """Refresh boards/ from the public panel files; test-board anchors are stripped in memory and never printed."""
    from vec.boards import update_boards

    s = update_boards()  # board keys come from the site: never markup
    console.print(f"boards ({len(s['boards'])}): {escape(', '.join(s['boards']))}")
    console.print(f"reduced to their submission spec ({len(s['reduced'])}): {escape(', '.join(s['reduced'])) or '-'}")
    if s["t1_composition_written"]:
        console.print("T1 composition written")
    else:
        console.print(f"T1 composition not written: it describes {escape(str(s['t1_composition_unreleased_stages']))} "
                      "unreleased stage(s)")


@app.command()
def check(pred: Path, board: str = typer.Option(..., help="board key, e.g. T1:val")):
    """Free format check of a prediction file (no scoring)."""
    rep = check_submission(pred, get_board(board))
    for w in rep.warnings:  # the messages quote the site's gene panel and the file's contents: never markup
        console.print(f"[yellow]warning[/] {escape(w)}")
    for e in rep.errors:
        console.print(f"[red]error[/] {escape(e)}")
    console.print(f"{'[green]OK' if rep.ok else '[red]REJECTED'}[/] ({rep.n_cells} cells)")
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
    method: str = typer.Option("", help="method name for the registry"),
    out: Path = typer.Option(None, help="write full JSON scores here"),
    register: bool = typer.Option(False, help="append a row to results/registry.csv"),
):
    """Score a prediction with the local board emulator."""
    b = get_board(board)
    stages = [Stage.from_anndata(read_stage(p, b.genes), b.genes, need_coords=b.needs_coords)
              for p in (pred, target, ref)]
    res = score_board_seeds(b, *stages, seeds=tuple(int(s) for s in seeds.split(",")), anchors=anchors, gate=gate,
                             frac=frac)
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
        REGISTRY.parent.mkdir(parents=True, exist_ok=True)
        new = not REGISTRY.exists()
        with REGISTRY.open("a", newline="") as f:
            w = csv.writer(f)
            if new:
                w.writerow(["utc_time", "commit", "board", "method", "pred", "score_mean", "score_sd", "seeds",
                            "anchors", "gate", "frac", "scores_json"])
            w.writerow([datetime.now(UTC).isoformat(timespec="seconds"), _git_commit(), board,
                        method or pred.stem, str(pred), f"{res['score_mean']:.3f}", f"{res['score_sd']:.3f}",
                        seeds, anchors, gate, frac, str(out) if out else ""])


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
    """List the files your account can download."""
    files = _authed(lambda: vec_data.manifest(vec_data.load_token(token_file)))
    t = Table("key", "task", "role", "stage", "modality", "MB")
    for f in files:  # server text
        t.add_row(*map(escape, [f.key, f.task, f.role, f.stage, f.modality, f"{f.bytes / 1e6:.0f}"]))
    console.print(t)
    console.print(f"{len(files)} files, {sum(f.bytes for f in files) / 1e9:.2f} GB")


@data_app.command("fetch")
def data_fetch(
    dest: str = typer.Option("~/v/virtual-embryo-challenge/data/raw", help="destination directory"),
    host: str = typer.Option("PrakashDGX_H2", help="ssh host to download onto; empty string = this machine"),
    match: str = typer.Option("", help="only keys containing this substring"),
    token_file: Path = typer.Option(vec_data.DEFAULT_TOKEN_FILE, help="bearer token file (mode 600)"),
):
    """Download released files (one signed link at a time); already complete files are skipped."""
    token = _authed(lambda: vec_data.load_token(token_file))
    files = [f for f in _authed(lambda: vec_data.manifest(token)) if match in f.key]
    for i, f in enumerate(files, 1):
        status = _authed(lambda f=f: vec_data.fetch_one(f, token, dest, host or None))
        console.print(f"[{i}/{len(files)}] {status:10s} {escape(f.key)} ({f.bytes / 1e6:.0f} MB)")


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
                              "row per attempt (results/submissions.csv). `vec submit FILE --board KEY` runs upload.")
app.add_typer(submit_app, name="submit")


def _print_plan(plan: vec_submit.Plan) -> None:
    if plan.official:
        console.print(f"[bold red]OFFICIAL upload[/] to {escape(plan.board)}: uses one of the board's "
                      f"{vec_submit.FINAL_PHASE_LIMIT} scored submissions for the whole phase; it cannot be withdrawn")
    else:
        console.print(f"Format check on {escape(plan.board)} (format_only=true): nothing is scored and no attempt is "
                      "spent")
    console.print(f"  file      {escape(str(plan.file))} (sent as {escape(plan.upload_name)}, "
                  f"{plan.file.stat().st_size / 1e6:.1f} MB)")
    console.print(f"  sha256    {plan.sha256}")
    console.print(f"  fields    {escape(' '.join(f'{k}={v}' for k, v in plan.fields.items()))}")
    if plan.official:
        console.print(f"  key       Idempotency-Key {escape(plan.idempotency_key or '-')}")
    console.print(f"  disclose  {escape(', '.join(plan.disclosure)) or '-'}")
    console.print(f"  commit    {plan.git_commit}")
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
            console.print(f"official upload (only with the user's go-ahead): vec submit {escape(str(p.file))} --board "
                          f"{escape(p.board)} --official --confirm {p.sha256[:12]} --disclose "
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
    file: Path = typer.Argument(..., help="prediction .h5ad"),
    board: str = typer.Option(..., help="board key, e.g. T1:val"),
    official: bool = typer.Option(False, "--official", help="OFFICIAL scored upload; also needs --confirm"),
    confirm: str = typer.Option("", help="first 12 hex characters of the file's sha256 (official uploads)"),
    disclose: str = typer.Option("", help="EXTERNAL_DATA.md ids, comma-separated (official uploads need one)"),
    note: str = typer.Option("", help="free-text note for the ledger"),
    model: str = typer.Option("", help="optional method name sent with the upload (shown with the team)"),
    architecture: str = typer.Option("", help="optional architecture name sent with the upload"),
):
    """Format check FILE for a board (default, free), or upload it officially with --official --confirm SHA12."""
    try:
        res = vec_submit.submit(file, board, official=official, confirm=confirm or None, disclose=disclose, note=note,
                                model=model or None, architecture=architecture or None, on_plan=_print_plan,
                                echo=lambda m: console.print(escape(m)))
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
