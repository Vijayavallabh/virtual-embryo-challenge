# AGENTS.md

Team Mithril's repository for the Virtual Embryo Challenge (NeurIPS 2026, Human Team track). `README.md` covers the
toolkit and the best submission, `METHOD_SUMMARY.md` the method, and `session-handoff.md` the current state.

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

## Working here

- Set up with `uv sync` and verify with `./init.sh`. Run the repository's Python as `python -I -B`.
- Competition data stays in `data/` (gitignored). Public use of it needs the data owners' approval
  (`docs/DATA_NOTES.md`).
- Register an external source in `EXTERNAL_DATA.md` before it enters a model.
- Scores choose between built predictions; they never compute one.
- Write before you run: a try in `experiments/harness/TRIES.md`, a decision in `docs/PLAN.md`, and a built file's
  sha256 in `experiments/builds/results/`.
- Official uploads and final marks are the captain's. Never commit or copy the API token.
- Do not read `docs/research/` or the paths `docs/DATA_POLICY.md` quarantines.
- Done means `./init.sh` passes, the records are written, and the change is pushed.
