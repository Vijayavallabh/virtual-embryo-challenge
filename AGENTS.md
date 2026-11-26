# AGENTS.md

Team Mithril's repository for the Virtual Embryo Challenge (NeurIPS 2026, Human Team track). `README.md` describes the
tasks and the `vec` toolkit, `METHOD_SUMMARY.md` the method, and `session-handoff.md` the current state.

## Guiding principle (the captain, 2026-11-04)

> Treat the virtual embryo not as an image to perfect, but as a process of **becoming**. Predict which cell states
> emerge, in what proportions, how their expression changes, and, where scored, how they organize in space. A
> population is not an average cell; an embryo is not expression without structure; a perturbation is not an
> unchanged wild type.
>
> Begin with the published copy-last-stage and unchanged-wild-type baselines across the boards you intend to enter.
> Improve the weakest weighted component before polishing what already works. Let developmental biology suggest
> testable hypotheses, not excuse complexity. Retain changes that generalize to unseen stages or knockouts—not merely
> those that flatter observed data or validation scores.
>
> Ask of every component: “What meaningful change beyond the baseline does this predict, and what held-out evidence
> shows it?”

## Setup and verification

```bash
uv sync      # local, CPU only (on the GPU machine: uv sync --group gpu --group jax)
./init.sh    # uv sync --locked, ruff, the full test suite
```

Run the repository's Python as `python -I -B`, and point `TMPDIR` at a disk with room before running tests.

## Rules

1. **Data stays out of git.** Competition data lives in `data/` (gitignored). Any public use of the unpublished data
   needs prior approval from Dr. Neil Chi and the data-generating team (`docs/DATA_NOTES.md`).
2. **External sources are registered first.** Every external dataset, model or piece of code gets a row in
   `EXTERNAL_DATA.md` before it enters a model, checked against the held-out windows in `docs/DATA_POLICY.md`.
3. **Scores choose, never compute.** Leaderboard returns may pick between built predictions, never solve for a
   property of a hidden target.
4. **Write it down first.** Register every try in `experiments/harness/TRIES.md` before it runs, date every decision
   in `docs/PLAN.md`, and record every built file with its sha256 in `experiments/builds/results/` before any score.
5. **Uploads are the captain's.** Agents prepare files and records; the team captain makes official uploads and
   final marks. The API token (`~/.config/vec/token`) is never committed or copied.
6. **Quarantine.** Do not read `docs/research/` or the paths `docs/DATA_POLICY.md` lists as quarantined.

## Done means

The tests pass (`./init.sh`), the records above are written, and the change is committed and pushed.
