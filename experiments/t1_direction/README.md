# T1 direction readout

Pre-registered before the final-phase release (2026-10-19). `readout.py` with `config.yaml` runs once on release day
(`docs/P3_RUNBOOK.md` step 4.3) on the E9.5 -> E10.5 transition and prints one verdict, "build" or "do not build".

2026-10-28: not run yet. The E10.5 file (`data/raw/T1/E10.5_RNA.h5ad`) did not come on release day (2026-10-27); it
reached the DGX on 2026-10-28, after the user unlocked the validation answers (runbook step 2). The readout still runs
once, at runbook step 4.3, after the user agrees to start step 4.

"build" needs all three checks:
1. Feng's E9.5-to-E14.5 direction beats cross-fitted momentum by at least 0.05 in de_score.
2. It also does so in de_direction.
3. In cardiomyocytes, Feng's direction has the right sign for at least 0.70 of the genes that change.

The verdict decides two things (`docs/PLAN.md` section 4.1):
1. whether M2 (OT bracketing from our E10.5 to external E14.5) is built;
2. the direction source of M1, the main T1 pipeline: the Feng bracket on "build"; damped momentum on "do not build"
   (a fail drops Feng from M1).

The second use was in the plan but not in `config.yaml` or `readout.py`. It was written down here when the user signed
off the release-day rules on 2026-10-21 (`docs/RELEASE_RULES_SUMMARY.md`, answer 6: all three checks kept). The rules
were drafted by the AI assistant and approved by the user. This file changes no code or setting.
