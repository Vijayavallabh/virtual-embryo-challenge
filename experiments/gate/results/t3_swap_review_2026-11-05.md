# Adversarial review of the T3 second-upload swap (2026-11-05)

Workflow wf_6b41d823-d5a: three independent reviewers (statistics, biology from released data and committed results,
governance), each told to refute `experiments/gate/t3_swap.py`'s case for making T3's s2 the barred-branch s1 (today's
s2 cells and composition plus the Mab21l2 factor). Read-only; no target read and no score computed. Verdicts: modify
(P(mistake) 0.55), refute (0.60), modify (0.60). All three: keep today's s2, as the frozen rule says. The swap is not
made (`docs/PLAN.md` decision log 2026-11-05 (d)).

What they found, recomputed from committed files:
- **The factor that carried to Gata4 is mostly a stage offset.** `transfer:mab21l2_all` divides Mab21l2's E9.5
  knockout cells by the test's own WT E8.75 file. Against its own WT E9.5 (`transfer:mab21l2_vs95`), the same
  knockout's factor scores 52.04 / 46.24 on the Gata4 replicates, below plain resampling; `all` minus `vs95` is
  +7.36 / +15.91, the WT E9.5-against-E8.75 axis. The Gata4 knockout leans towards the E9.5 file (stage weights 0.22 /
  0.32), and in-sample stage mixing alone (`fac:stage`, 59.10 / 60.36) matches B (59.40 / 62.15).
- **So the two bets are coupled.** s1's Gata4 factor is also referenced to WT E8.75 and carries the same offset; K2a
  failed on a fold referenced to WT E9.5, where that offset cannot carry. The script's world KB (s1 fails, the swap's
  factor carries) splices two folds that rarely co-occur, and all of the priced value came from it (world A adds
  0.000). Coupled, the expected T3 gain is about -0.22 to +0.11.
- **The paired +9.90 (20 of 20) is two replicates times ten resamples of the same files:** +5.78 on rep1 (seed sd
  3.56) and +14.02 on rep2. 83% of it is response magnitude, because today's s2 stand-in sits at the severity floor in
  14 of 20 seeds and B in none; with neither at the floor the margin is +3.12. B's response identity is at chance.
- **The downside was floored.** The both-fail world puts s1 at K2a (50.08), only 0.66 under the WT copy; failing
  constructs on the committed folds score 46-48. With s1 four points lower, the expected T3 gain turns negative in
  four of six weightings, and in the most pessimistic one P(Overall >= our median) falls from 50% to 37%. "Never lowers
  P(Overall >= c)" held only under that floor.
- **Neither stand-in is the file to be uploaded.** B has no composition change, and the barred s1 (factor plus the
  central composition) was never scored; factor-by-composition interactions measured on Gata4 run from -1.0 to +5.6.
- **For a top-3 aim it is worth about nothing:** above our 90th percentile the largest change is +0.05 points of chance.
- **Governance.** Today's s2 is the captain's ruling of 2026-11-01 and rules.yaml note 6's pre-registered consequence of
  D11b's miss. The swap has no counting fold under the adoption rule: Mab21l2 -> Gata4 is the fold D11b read, and the
  other is the factor's own knockout. Three of the four T3 ruling switches would break, one by uploading the same file
  twice. A10's go/no-go (D11b >= 60) had failed on 2026-10-28, but the ledger still said open.
