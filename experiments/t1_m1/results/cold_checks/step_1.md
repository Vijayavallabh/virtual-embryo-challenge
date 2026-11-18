Record of a cold adversarial check (the devil-advocate procedure: a fresh agent is given the claim block only, never the
conversation, and returns a kill shot, three ranked attacks with checks, a cheaper alternative and a chance of a
mistake). Committed on 2026-11-18 so that the decision log's references can be checked; the text holds practice
figures and the claim's own numbers only. The Advocate's output is verbatim except that pronouns for the team
captain are replaced by bracketed words. What was done about each attack: `docs/PLAN.md` decision log, the entry
named below.

Decision-log entry: 2026-11-17 (b) to (d).

Cold adversarial check of the Task 1 step request (devil-advocate skill; agent ae07af21fa6c2ce63; 2026-11-17).
The claim block is in the session record; the Advocate's output, verbatim:

KILL SHOT
0.69 is supported only while the test's real share is at least 0.28 (measurement stops at 1.43 times real); the author's floor is 0.29. Practice fixes only pace times scale (0.254), so on the even course that floor is (scale - 0.254) / 2, and the scale was never measured over time: 0.83-0.96 comes from cell-type contrasts within E9.5, where sample-level effects (batch, dissection) cancel; d contrasts samples five days apart, where they do not. At a time scale of 0.70 the test holds 0.22, today's 0.40 already delivers 0.233, and 0.69 over-delivers 1.8 times on both uploads.

RANKED ATTACKS
1. Wrong-contrast scale (above); its 0.96 end is the noise-free limit, which a noisy d never reaches. | Holds if: a quarter of d is sample difference or noise. | Check: split d into its part common to all cell groups (cancelled in cell-type contrasts, so untested by 0.83-0.96) and the rest; a lower released-day share on the common part confirms it.
2. 0.35 is a middle course, not a floor. The ruler back-loads: a true first fifth reads 0.68 of nominal, so later fifths read above it and an even true pace puts days two-three under half the remainder. | Holds if: pace stops decaying after E10.5. | Check: read fifths two-three against four-five on the reviewer's synthetic counts, or the outside dataset's middle stages if any.
3. The shortfall story fails across truths: library 3 has the larger real share (0.326) yet flattens first; libraries 1-2 (0.254) still gain 0.4-0.7 per half-step at 1.4-1.8 times theirs. The over-delivery 0.69 leans on is one library's behaviour on one in-sample day. | Holds if: the flat point follows truth depth, not the ratio. | Check: thin library 3's truth to 0.3 depth; re-score the stored predictions.

CHEAPER ALTERNATIVE
Same build, smaller value: set 0.55 unless checks 1-2 clear 0.69. It delivers 0.32, inside the author's range; by library 3's curve it earns 0.7-0.8 points against the modelled 1.0-1.2 and stays inside the measured range down to a real share of 0.22 (time scale 0.70), where 0.69 needs 0.28 (0.81).

P(MISTAKE)
30% - I could not break the practice dose response or its placebo control, so the sign probably carries; the risk is an unmeasured time scale leaving 0.69 over-delivering on both uploads for no gain or a loss.
