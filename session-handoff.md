# Session hand-off

State on 2026-11-26, the last upload day (uploads close at 11:59 UTC). The full previous hand-off is in git history:
`git show c10e079:session-handoff.md`.

- **Best submission and final picks:** the table in `README.md`, overall 205.62. Row R81 of
  `experiments/harness/TRIES.md` names the embryo pick. The captain marks one final submission per board on the
  account page by 2026-12-02, and winners are announced on 2026-12-11.
- **Records:**
  - uploads and headlines: `results/manual_uploads.csv`;
  - built files and their sha256: `experiments/builds/results/`;
  - decisions: `docs/PLAN.md`.
- **Open:** the final marks, which are the captain's. Row R80 (four embryo files) is ruled but unbuilt; nothing runs.
- **Kept after the cleanup of 2026-11-26:** the code (this repository and its history) and the four final files in
  `runs/manual_upload/`. The competition data, every other build output, the agent worktrees and every challenge
  folder on the DGX were deleted; `vec data fetch` downloads the data again.
- **Resume:** run `./init.sh`, then read `AGENTS.md` and the last entries of `docs/PLAN.md`.
