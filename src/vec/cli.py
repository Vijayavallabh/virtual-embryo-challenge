"""`vec` command line: board listing, free format checks, local scoring and reference baselines."""

from __future__ import annotations

import csv
import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

from vec import baselines
from vec import data as vec_data
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


@app.command()
def boards():
    """List the known boards and their constraints."""
    t = Table("board", "task", "genes", "coords", "cells", "floor")
    for key, b in load_boards().items():
        t.add_row(key, b.task, str(b.n_genes), str(b.needs_coords), f"{b.min_cells}-{b.max_cells}", b.floor_model)
    console.print(t)


@app.command()
def check(pred: Path, board: str = typer.Option(..., help="board key, e.g. T1:val")):
    """Free format check of a prediction file (no scoring)."""
    rep = check_submission(pred, get_board(board))
    for w in rep.warnings:
        console.print(f"[yellow]warning[/] {w}")
    for e in rep.errors:
        console.print(f"[red]error[/] {e}")
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
    except vec_data.AuthError as e:
        console.print(f"[red]auth error[/] {e}")
        raise typer.Exit(2) from None
app.add_typer(data_app, name="data")


@data_app.command("manifest")
def data_manifest(token_file: Path = typer.Option(vec_data.DEFAULT_TOKEN_FILE, help="bearer token file (mode 600)")):
    """List the files your account can download."""
    files = _authed(lambda: vec_data.manifest(vec_data.load_token(token_file)))
    t = Table("key", "task", "role", "stage", "modality", "MB")
    for f in files:
        t.add_row(f.key, f.task, f.role, f.stage, f.modality, f"{f.bytes / 1e6:.0f}")
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
        console.print(f"[{i}/{len(files)}] {status:10s} {f.key} ({f.bytes / 1e6:.0f} MB)")


if __name__ == "__main__":
    app()
