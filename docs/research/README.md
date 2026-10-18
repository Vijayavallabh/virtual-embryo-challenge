# Research notes (2026-10-18)

Background reports written during the initial research sweep. `docs/PLAN.md` is the decision document; these are the
evidence behind it. Paths starting with `SCRATCH/` refer to the session scratch area where probes and downloads were
kept; they are not part of the repo. Held-out-stage numbers found in the organisers' veckit 0.1.2 docstring were
redacted before these files were added (see `docs/DATA_POLICY.md`).

| Report | What it covers |
|---|---|
| `metrics.md` | Exact scorer behaviour, the skill map, where the points are, the structure gate, a per-task metric playbook |
| `landscape.md` | Live leaderboard (388 teams), organiser reference rows, Navigo/Spateo/dynamo, veckit issues and PRs, community resources |
| `community_digest.md` | Measured recipes from about 50 public participant repos (licences noted) |
| `slack_digest.md` | Every message and thread in the PantheonOS `#virtual-embryo-challenge` channel, with organiser rulings |
| `t1_methods.md` | Temporal forecasting methods, metric probes, candidate T1 pipelines, foundation-model legality |
| `t2_spatial.md` | Spatiotemporal OT and generation, coordinate facts, embryo and heart pipelines |
| `t3_methods.md` | Perturbation-prediction methods, why shift_transfer fails, T3 pipelines and validation folds |
| `t3_biology.md` | Mesp1-Cre beta-catenin, Gata4/Gata6 and Mab21l2 biology, and the rules question (now settled: no gene-level knockout findings) |
| `external_data.md` | Public dataset inventory with exact per-task removal rules |
| `architecture.md` | Shared-architecture options, cross-task data rules, portfolio and submission policy |
| `engineering.md` | uv environment on the DGX, compute and memory estimates, reproducibility and submission safety |
| `digest_2026-10-19.md` | Screened summary of the verified sweep (9 dimensions, 9 verifiers, critic, 5 gap-fills): the lead's source for PLAN v1.2 |
| `compliance_items_2026-10-19.md` | Window, size and judgement items from that screen, described without banned content |

## Read these notes with the data policy (2026-10-19)

Two meaning-level screens removed gene-level knockout content and held-out-stage values from these notes (marked
`[redacted 2026-10-19: A/B]`). Several recommendations that remain conflict with `docs/DATA_POLICY.md`; where they do,
the policy and `docs/PLAN.md` win:

- JAX embryo lists "1-24 and 50-74" or "after E13.5 by label": use embryos 1-22 and 54-74 only (23-24 and 50-53 come
  from boundary or in-window litters); embryos 1-12 are T1-only. Never filter atlases by curated stage labels alone.
- Label-only E9.5 sources, Theiler-stage filters and a pooled "E8.25" sample kept on its label: do not use them
  without the organisers' somite convention (follow-up question 7).
- The heart Visium dataset whose capture areas mix stages: skip it or remove spots by tissue position, and disclose.
- The data-generating lab's lineage dataset: never its E7.5/E7.75 samples, nor in any model that serves T2-embryo.
- Wnt-adjacent mutant datasets proposed as legal for T3: treat as possible phenocopies; do not use.
- In-window developmental literature used for E12.5 composition, emergent classes or niche placement: contingent on
  the organisers' answer; the plan uses the legal route meanwhile.
- Cross-task variants that feed released stages into another task's window: wait for the cross-task ruling (A15).
- The organisers' atlas "for orientation": prohibited at the target time points.
- Any advice to install or read veckit 0.1.2: never (pin git `cc3859a`).
- Raw research-workflow outputs and unscreened scratch copies: never open; use `digest_2026-10-19.md`.
