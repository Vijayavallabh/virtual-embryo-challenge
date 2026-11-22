<!-- Written by a fresh score-blind analyst and corrected by a cold checker on 2026-11-22 (workflow embryo-so-ohi-question-and-build),
     from row E51's committed reading; every figure was checked against reading.json. The main session, recused on the embryo's
     shape forms and occupancy term, passes it to the captain unchanged; it fills the two sha256 blanks once the builds agree. -->

Embryo board: two new candidate files, SO and O_hi

O is the overlay file: the exact-position file P with the overlay shape form. SO keeps O's coordinates and takes its expression from the series file s2xs (P's positions with the series expression). O_hi keeps O's coordinates, rows, count, seed and size, and tops up every drawn cell's expression to the deeper bracket's median depth, with no per-gene factor. Both passed a practice read (nine practice tasks, 10 seeds) under pass rules fixed before the run. A pass only puts the question to you. The read never compares the two candidates directly and does not choose between them.

Practice figures in board points, candidate minus its control (seeds ahead, of 10). "Merged" means five practice stages at equal weight. "Test-like" means the four tasks on E8.0, E8.25 and E8.5.

| | Merged | Test-like |
|---|---|---|
| SO vs s2xs | +2.240 (10 ahead) | +3.147 (10 ahead) |
| O_hi vs O | +1.059 (10 ahead) | +1.838 (10 ahead) |

All four are diagnostics. The ledger's noise floor is 6.54 board points and every gain is below it, so none of them forecasts a test score. Reported, deciding nothing: SO vs O, merged -0.043 (5 ahead, 5 behind) and test-like +3.111 (10 ahead). O_hi vs s2xs, merged +3.341 and test-like +1.874 (10 ahead each).

At the test (E7.75, which is not a practice stage):
- SO needs E7.75 to share E8.0's detection profile and depth (judged 0.4 to 0.5). It also needs the series factor to hold at E7.75. An earlier check (E26) found that factor does not carry across the second profile jump.
- O_hi needs E7.75's depth to be near E8.0's. That depth is not known, and if it is shallow O_hi loses about 1 point against O. O_hi also needs E7.75 to have its own profile (judged 0.2 to 0.3).
- Both need O's neighbourhood gain to hold at the test pair, which has not been shown without the target. The live gate, scorer build and metric seed are not shown.

Price: a file that scores lower does not replace the board's best, so a miss costs one upload slot (at most two a day; the window closes 26 Nov 11:59 UTC). The embryo's order rule puts these after P, C and the sized copy-rows file. Among themselves the order is SO, then T2embryo_s2 (an earlier embryo file built on 13 Nov), then O_hi. You may use fewer, none, or another order.

[Added by the main session on 2026-11-22: both files were built after this text was written; for each, two builds agree byte for byte and the identity checks hold. The line that stood here said neither was built yet.]
Files: SO = T2embryo_s2xs_ov, sha256 77062082d97555961901c811d8a1aaf9b5b2403417c19a6988c816ff553f1165 (`runs/manual_upload/T2embryo_s2xs_ov.h5ad`); O_hi = T2embryo_s2x_ovhi, sha256 7e4d1eaafc40080d9d9c461d4b9396f0cee6353686538e79fb92ca31d1d2deeb (`runs/manual_upload/T2embryo_s2x_ovhi.h5ad`). Data used: released challenge data only; no outside data.

1. Do you name SO? Yes / No (default: No)
2. Do you name O_hi? A yes also lifts the embryo's depth line (the deeper expression), which the ledger had closed for files, for this one file only. Yes / No (default: No)

If you give no answer, neither file is uploaded.
