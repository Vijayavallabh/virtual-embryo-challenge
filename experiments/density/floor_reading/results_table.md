# Floor-case reading: G (next best gain) vs C (next count up), Monte Carlo results

Noise: per-seed paired SD fitted to the 20-seed shuffle50 review (A=0.172, B=1.115, rho=0.362); base 1-doubling paired SD 0.94 points. Paired SDs at the rule's levels: 0.25 vs 1 1.07, 0.5 vs 1 0.94, 0.75 vs 1 0.75, 0.5 vs 0.75 0.82.

D = true test score of C's choice (0.5 x max_T) minus G's (0.75 x max_T), points, in the disagreement case. Base geometry r_V = max_V/truth_V ~ U[4, 8], r_T ~ U[3, 4); no transfer tilt.

## Headline

- x mapping (score depends on n/truth), all x families: disagreement in 0.34% of releases; P(C better) 0.82; mean D +0.40 (sd 0.59). Smooth families only: 0.16%, 0.89, +0.26 (sd 0.22).
- abs mapping: P(C better) 0.93, mean D +0.35. f-mapping stress test: 0.38, -0.20 (G wins only through rugged or wiggly effects carried over by f).
- One-at-a-time sensitivities, x mapping: P(C better) 0.49 to 0.85, mean D +0.07 to +0.50. Break-even test tilt toward more cells: +0.45 (smooth) to +0.68 (all x) points per doubling.
- Expected gain of C over G per release (floor geometry assumed): +0.0014 points (all x), +0.0004 (smooth x).

## Base, per family

| family | P(floor) | P(disagree) | P(dis given floor) | n disagree | P(C better) | mean D | sd D | 5% | 95% | P(C better, one realised scoring) | P(true val g0.5 > g0.75) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| X0 more-is-better (offline) | 0.000237 | 1.7e-06 | 0.007 | 2 | 0.000 | -0.133 | 0.025 | -0.148 | -0.117 | 0.436 | 0.000 |
| X1 decreasing above 1x | 0.819 | 0.0037 | 0.0045 | 1113 | 1.000 | +0.359 | 0.135 | +0.160 | +0.598 | 0.667 | 1.000 |
| X1s decrease then plateau | 0.422 | 0.0024 | 0.0056 | 712 | 1.000 | +0.295 | 0.129 | +0.079 | +0.499 | 0.639 | 1.000 |
| X2q unimodal, log-quadratic | 0.473 | 0.0022 | 0.0047 | 665 | 0.680 | +0.072 | 0.195 | -0.267 | +0.359 | 0.534 | 1.000 |
| X2b unimodal bump to plateau | 0.369 | 0.0013 | 0.0036 | 399 | 0.820 | +0.285 | 0.348 | -0.224 | +0.911 | 0.625 | 1.000 |
| X3 flat plus penalty below 1x | 0.0187 | 0.00015 | 0.008 | 180 | 0.000 | -0.040 | 0.047 | -0.142 | -0.000 | 0.481 | 0.000 |
| X4 GP wiggle, flat trend | 0.171 | 0.0033 | 0.019 | 988 | 0.698 | +0.336 | 0.721 | -0.795 | +1.628 | 0.620 | 0.365 |
| X4d GP wiggle, decreasing trend | 0.606 | 0.015 | 0.024 | 4352 | 0.803 | +0.499 | 0.697 | -0.555 | +1.695 | 0.681 | 0.442 |
| F0 increasing in f | 8.08e-05 | 8.3e-07 | 0.01 | 1 | 0.000 | -0.112 | nan | -0.112 | -0.112 | 0.446 | 0.000 |
| F1 decreasing in f | 0.804 | 0.0035 | 0.0044 | 1050 | 1.000 | +0.347 | 0.140 | +0.150 | +0.608 | 0.662 | 1.000 |
| F1s decrease-plateau in f | 0.548 | 0.0025 | 0.0045 | 745 | 1.000 | +0.199 | 0.121 | +0.028 | +0.429 | 0.595 | 1.000 |
| F2q unimodal in f, quad | 0.355 | 0.0019 | 0.0054 | 580 | 1.000 | +0.403 | 0.158 | +0.151 | +0.681 | 0.685 | 1.000 |
| F2b unimodal bump in f | 0.288 | 0.0011 | 0.0037 | 320 | 0.991 | +0.298 | 0.176 | +0.007 | +0.569 | 0.639 | 0.991 |
| F3 rugged level effects | 0.25 | 0.023 | 0.093 | 6928 | 0.145 | -0.412 | 0.439 | -1.222 | +0.166 | 0.332 | 0.145 |
| F4 GP wiggle in f | 0.172 | 0.0035 | 0.021 | 1063 | 0.367 | -0.142 | 0.351 | -0.750 | +0.356 | 0.439 | 0.368 |
| A1 decreasing in n (abs) | 0.802 | 0.0035 | 0.0044 | 1064 | 1.000 | +0.357 | 0.148 | +0.145 | +0.616 | 0.666 | 1.000 |
| A2b unimodal bump in n (abs) | 0.289 | 0.0011 | 0.004 | 343 | 0.875 | +0.233 | 0.303 | -0.167 | +0.775 | 0.604 | 0.988 |
| A2q unimodal quad in n (abs) | 0.356 | 0.0018 | 0.0049 | 529 | 0.828 | +0.419 | 0.468 | -0.219 | +1.306 | 0.668 | 0.998 |
| A0 increasing in n (abs) | 8.92e-05 | 0 (of 1,200,000) | - | 0 | - | - | - | - | - | - | - |

## Pooled over families (posterior weight = prior x P(disagree | family))

| prior | P(disagree) per release | P(C better) | mean D | sd D | E[gain of C per release] | posterior weights |
|---|---|---|---|---|---|---|
| x families incl. GP wiggle, equal prior | 0.0034 | 0.820 | +0.395 | 0.589 | +0.0014 | X1 0.13; X1s 0.09; X2q 0.08; X2b 0.05; X3 0.01; X4 0.12; X4d 0.53 |
| x smooth families, equal prior | 0.0016 | 0.887 | +0.262 | 0.222 | +0.00043 | X1 0.38; X1s 0.24; X2q 0.23; X2b 0.14; X3 0.02 |
| x families, offline-heavy prior (X0 0.6, X3 0.1, fewer-is-better 0.075 each) | 0.00074 | 0.882 | +0.260 | 0.223 | +0.00019 | X1 0.38; X1s 0.24; X2q 0.23; X2b 0.14; X3 0.02 |
| fewer-is-better x families only | 0.0024 | 0.901 | +0.267 | 0.221 | +0.00064 | X1 0.39; X1s 0.25; X2q 0.23; X2b 0.14 |
| f-mapping stress families, equal prior | 0.0051 | 0.382 | -0.203 | 0.489 | -0.001 | F1 0.10; F1s 0.07; F2q 0.05; F2b 0.03; F3 0.65; F4 0.10 |
| abs-mapping stress families, equal prior | 0.0016 | 0.931 | +0.352 | 0.303 | +0.00057 | A1 0.55; A2b 0.18; A2q 0.27 |
| all families: mapping x 1/2, f 1/4, abs 1/4 | 0.0034 | 0.670 | +0.166 | 0.598 | +0.00057 | X1 0.07; X1s 0.04; X2q 0.04; X2b 0.02; X4 0.06; X4d 0.27; F1 0.04; F1s 0.03; F2q 0.02; F2b 0.01; F3 0.24; F4 0.04; A1 0.07; A2b 0.02; A2q 0.03 |

## One-at-a-time sensitivity (pooled, equal prior within each mapping)

| assumption changed | x: P(dis) | x: P(C better) | x: mean D | f: P(C better) | f: mean D | abs: P(C better) | abs: mean D |
|---|---|---|---|---|---|---|---|
| (base) | 0.00345 | 0.820 | +0.395 | 0.382 | -0.203 | 0.931 | +0.352 |
| noise sd1 0.45 (large heart truth half) | 0.00381 | 0.782 | +0.467 | 0.161 | -0.312 | 0.944 | +0.265 |
| noise sd1 1.35 (noisiest target) | 0.00278 | 0.845 | +0.418 | 0.503 | -0.093 | 0.918 | +0.488 |
| noise independent across levels | 0.00558 | 0.851 | +0.352 | 0.589 | -0.023 | 0.917 | +0.371 |
| noise strongly correlated (rho 0.6/doubling) | 0.00299 | 0.798 | +0.405 | 0.305 | -0.274 | 0.916 | +0.338 |
| noise t-tails (df 4) | 0.00473 | 0.811 | +0.383 | 0.424 | -0.148 | 0.939 | +0.375 |
| noise sd1 per world, log-uniform 0.4-1.35 | 0.00375 | 0.806 | +0.404 | 0.293 | -0.238 | 0.945 | +0.338 |
| r_V in [4, 5] | 0.00351 | 0.720 | +0.169 | - | - | - | - |
| r_V in [6, 8] | 0.00345 | 0.854 | +0.496 | - | - | - | - |
| r_V = 8 (worst case for C) | 0.0034 | 0.798 | +0.410 | - | - | - | - |
| r_V in [3, 4) (0.25x below truth on validation) | 0.00148 | 0.665 | +0.107 | - | - | - | - |
| r_T = 3 exactly | 0.00341 | 0.822 | +0.459 | - | - | - | - |
| r_T = 3.9 | 0.00351 | 0.796 | +0.331 | - | - | - | - |
| transfer tilt sd 0.3 pts/doubling | 0.00354 | 0.792 | +0.380 | 0.375 | -0.195 | 0.851 | +0.334 |
| transfer tilt sd 0.6 pts/doubling | 0.00349 | 0.757 | +0.395 | 0.403 | -0.187 | 0.816 | +0.407 |
| test favours more cells: tilt mean +0.3, sd 0.3 | 0.00338 | 0.658 | +0.236 | 0.254 | -0.382 | 0.736 | +0.191 |
| test favours more cells: tilt mean +0.6, sd 0.3 | 0.00341 | 0.492 | +0.067 | 0.118 | -0.561 | 0.459 | -0.016 |

## Family parameter ranges (base scenario)

| family variant | P(disagree) | n disagree | P(C better) | mean D | 5% | 95% |
|---|---|---|---|---|---|---|
| X1 shallow (b 0.2-0.6) | 0.0045 | 904 | 1.000 | +0.278 | +0.134 | +0.455 |
| X1 steep (b 1-3) | 0.0013 | 265 | 1.000 | +0.587 | +0.347 | +0.870 |
| X2q peak 1x-2x | 0.0027 | 536 | 0.976 | +0.186 | +0.025 | +0.420 |
| X2q peak 2x-4x | 0.0018 | 363 | 0.223 | -0.097 | -0.361 | +0.081 |
| X2q peak 2x-4x, shallow right side (DR 0.05-0.3) | 0.00081 | 163 | 0.227 | -0.086 | -0.344 | +0.060 |
| X2b peak 2x-4x | 0.0006 | 120 | 0.233 | -0.082 | -0.350 | +0.174 |
| X2b peak 1.5x-3x, narrow (w 0.25-0.6) | 0.00068 | 137 | 0.825 | +0.554 | -0.549 | +1.505 |
| X2b peak 1x-8x | 0.00097 | 194 | 0.840 | +0.324 | -0.198 | +0.991 |
| X4 GP, small wiggle (sig 0.05-0.3) | 0.0006 | 120 | 0.600 | +0.041 | -0.284 | +0.343 |
| X4 GP, large wiggle (sig 0.5-1.5) | 0.0061 | 1226 | 0.697 | +0.512 | -1.179 | +2.303 |
| X4 GP, short length (ell 0.1-0.3) | 0.011 | 2155 | 0.638 | +0.320 | -1.204 | +1.949 |
| X4 GP, long length (ell 1-2) | 0.001 | 201 | 0.861 | +0.172 | -0.061 | +0.417 |
| X4d GP on decreasing trend, small wiggle (sig 0.05-0.3) | 0.0081 | 1625 | 0.955 | +0.385 | +0.017 | +0.801 |
| F4 GP in f, small wiggle (sig 0.05-0.3) | 0.00064 | 129 | 0.512 | -0.010 | -0.347 | +0.328 |
| F3 rugged, small level effects (tau 0.1-0.4) | 0.031 | 6250 | 0.218 | -0.172 | -0.590 | +0.184 |
| F3 rugged, large level effects (tau 1-3) | 0.02 | 3908 | 0.085 | -0.847 | -2.450 | +0.098 |

## Shift scan (smooth abs-mapping shapes A1, A2b, A2q, equal prior)

s = test level positions minus validation level positions on the true curve, in doublings of count. x mapping: s = log2(r_T/r_V), -1.42 to 0; f mapping: s = 0; abs mapping: s = log2(max_T/max_V).

| s | P(disagree) | P(C better) | mean D |
|---|---|---|---|
| -1.4 | 0.0022 | 0.673 | +0.148 |
| -1.0 | 0.002 | 0.808 | +0.243 |
| -0.8 | 0.0021 | 0.854 | +0.248 |
| -0.5 | 0.002 | 0.941 | +0.295 |
| +0.0 | 0.0021 | 0.995 | +0.364 |
| +0.5 | 0.0021 | 0.996 | +0.410 |
| +1.0 | 0.0021 | 0.999 | +0.483 |

Extended to larger negative shifts (test levels far lower on the curve than validation levels; per-family P(C better)/mean D; A1 at -3.0 is a tie on its flat tail):

| s | P(disagree) | P(C better) | mean D | per family P(C better)/mean D |
|---|---|---|---|---|
| -1.7 | 0.0021 | 0.563 | +0.071 | A1 1.00/+0.41; A2b 0.06/-0.33; A2q 0.02/-0.35 |
| -2.0 | 0.0021 | 0.576 | +0.083 | A1 1.00/+0.52; A2b 0.00/-0.44; A2q 0.00/-0.55 |
| -2.5 | 0.0021 | 0.546 | -0.236 | A1 1.00/+0.13; A2b 0.00/-0.32; A2q 0.00/-0.87 |
| -3.0 | 0.0021 | 0.000 | -0.385 | A1 0.00/+0.00; A2b 0.00/-0.20; A2q 0.00/-1.15 |

## Diagnostics

- X0 more-is-better (offline): {"p_C_better_k_lt_1_29": 0.0, "p_C_better_k_ge_1_29": null, "share_k_ge_1_29": 0.0}
- X1 decreasing above 1x: {"p_C_better_k_lt_1_29": 1.0, "p_C_better_k_ge_1_29": 1.0, "share_k_ge_1_29": 0.017969451931716084, "k 0.00-0.50": {"n": 287, "p_C_better": 1.0, "mean_D": 0.3520389153576737}, "k 0.50-0.90": {"n": 453, "p_C_better": 1.0, "mean_D": 0.366715885105663}, "k 0.90-1.29": {"n": 353, "p_C_better": 1.0, "mean_D": 0.3535832347359628}, "k 1.29-1.50": {"n": 20, "p_C_better": 1.0, "mean_D": 0.39083706678481694}}
- X1s decrease then plateau: {"p_C_better_k_lt_1_29": 1.0, "p_C_better_k_ge_1_29": 1.0, "share_k_ge_1_29": 0.014044943820224719, "k 0.00-0.50": {"n": 221, "p_C_better": 1.0, "mean_D": 0.26551509537895335}, "k 0.50-0.90": {"n": 307, "p_C_better": 1.0, "mean_D": 0.2961492998091336}, "k 0.90-1.29": {"n": 174, "p_C_better": 1.0, "mean_D": 0.3274198249864254}}
- X2q unimodal, log-quadratic: {"p_C_better_k_lt_1_29": 0.6919504643962848, "p_C_better_k_ge_1_29": 0.2631578947368421, "share_k_ge_1_29": 0.02857142857142857, "peak_x_quartiles_in_disagreement": [1.3557909079171484, 1.7809291760354247, 2.300268952462134], "k 0.00-0.50": {"n": 135, "p_C_better": 0.9185185185185185, "mean_D": 0.23826396400064165}, "k 0.50-0.90": {"n": 261, "p_C_better": 0.7471264367816092, "mean_D": 0.09507998276743418}, "k 0.90-1.29": {"n": 250, "p_C_better": 0.512, "mean_D": -0.026557621935285467}}
- X2b unimodal bump to plateau: {"p_C_better_k_lt_1_29": 0.8350785340314136, "p_C_better_k_ge_1_29": 0.47058823529411764, "share_k_ge_1_29": 0.042606516290726815, "peak_x_quartiles_in_disagreement": [1.2488696810960978, 1.5078836969798526, 1.8974005750154859], "k 0.00-0.50": {"n": 67, "p_C_better": 0.9402985074626866, "mean_D": 0.3363899087169353}, "k 0.50-0.90": {"n": 163, "p_C_better": 0.8773006134969326, "mean_D": 0.29414522765340634}, "k 0.90-1.29": {"n": 152, "p_C_better": 0.743421052631579, "mean_D": 0.27263664230488943}}
- X3 flat plus penalty below 1x: {"p_C_better_k_lt_1_29": 0.0, "p_C_better_k_ge_1_29": 0.0, "share_k_ge_1_29": 0.027777777777777776, "k 0.00-0.50": {"n": 21, "p_C_better": 0.0, "mean_D": -0.012077241381438822}, "k 0.50-0.90": {"n": 80, "p_C_better": 0.0, "mean_D": -0.03309303371285601}, "k 0.90-1.29": {"n": 74, "p_C_better": 0.0, "mean_D": -0.053300465508674455}}
- X4 GP wiggle, flat trend: {"p_C_better_k_lt_1_29": 0.7034339229968782, "p_C_better_k_ge_1_29": 0.5185185185185185, "share_k_ge_1_29": 0.027327935222672066, "k 0.00-0.50": {"n": 231, "p_C_better": 0.4329004329004329, "mean_D": -0.09026839424581835}, "k 0.50-0.90": {"n": 385, "p_C_better": 0.7766233766233767, "mean_D": 0.40245057756123515}, "k 0.90-1.29": {"n": 345, "p_C_better": 0.8028985507246377, "mean_D": 0.5682744024743817}, "k 1.29-1.50": {"n": 27, "p_C_better": 0.5185185185185185, "mean_D": 0.055233864707966394}}
- X4d GP wiggle, decreasing trend: {"p_C_better_k_lt_1_29": 0.8050787679285211, "p_C_better_k_ge_1_29": 0.696969696969697, "share_k_ge_1_29": 0.02274816176470588, "k 0.00-0.50": {"n": 1112, "p_C_better": 0.585431654676259, "mean_D": 0.06291104830727894}, "k 0.50-0.90": {"n": 1683, "p_C_better": 0.8805704099821747, "mean_D": 0.612282424091164}, "k 0.90-1.29": {"n": 1458, "p_C_better": 0.8854595336076817, "mean_D": 0.7151099369964579}, "k 1.29-1.50": {"n": 99, "p_C_better": 0.696969696969697, "mean_D": 0.28083623035489674}}
- A1 decreasing in n (abs): {"p_C_better_delta_-1.5_-0.8": 1.0, "p_C_better_delta_-0.8_+0.0": 1.0, "p_C_better_delta_+0.0_+0.8": 1.0, "p_C_better_delta_+0.8_+1.5": 1.0}
- A2b unimodal bump in n (abs): {"p_C_better_delta_-1.5_-0.8": 0.5897435897435898, "p_C_better_delta_-0.8_+0.0": 0.9310344827586207, "p_C_better_delta_+0.0_+0.8": 1.0, "p_C_better_delta_+0.8_+1.5": 0.9444444444444444}
- A2q unimodal quad in n (abs): {"p_C_better_delta_-1.5_-0.8": 0.4251968503937008, "p_C_better_delta_-0.8_+0.0": 0.8741258741258742, "p_C_better_delta_+0.0_+0.8": 1.0, "p_C_better_delta_+0.8_+1.5": 1.0}
- x families incl. GP, pooled (equal prior), by k = log2(r_V/r_T) and by peak: {"k 0.00-0.50": {"share": 0.2487461478034926, "p_C_better": 0.7025385643143448, "mean_D": 0.12800167372291021}, "k 0.50-0.90": {"share": 0.39543174814188164, "p_C_better": 0.8798899755501224, "mean_D": 0.46288343827515965}, "k 0.90-1.29": {"share": 0.3324672185630552, "p_C_better": 0.8491457651762995, "mean_D": 0.5286150222187638}, "k 1.29-1.50": {"share": 0.023354885491570488, "p_C_better": 0.6520051746442433, "mean_D": 0.20488309773682706}, "peak 1x-1.5x (unimodal families)": {"share_of_unimodal": 0.3956766917293233, "p_C_better": 1.0, "mean_D": 0.33562363069884454}, "peak 1.5x-2x (unimodal families)": {"share_of_unimodal": 0.2894736842105264, "p_C_better": 0.9318181818181818, "mean_D": 0.18318614435464453}, "peak 2x-2.5x (unimodal families)": {"share_of_unimodal": 0.1588345864661654, "p_C_better": 0.4201183431952663, "mean_D": -0.020121965502649244}, "peak 2.5x-3x (unimodal families)": {"share_of_unimodal": 0.09680451127819546, "p_C_better": 0.0, "mean_D": -0.17579664986503704}, "peak 3x-4.01x (unimodal families)": {"share_of_unimodal": 0.05921052631578946, "p_C_better": 0.0, "mean_D": -0.2306904081341881}}

## What each reading picks across all floor cases (floor_composition.py, base, x families)

| family | P(floor) | floor: G picks 0.5/0.75/max | floor: C picks 0.5/0.75/max | P(dis given floor) | E[D over all floor cases] | break-even test tilt (pts/doubling toward more cells) |
|---|---|---|---|---|---|---|
| X0 more-is-better (offline) | 0.000258 | 0.074/0.049/0.877 | 0.074/0.049/0.877 | 0.0000 | +0.00000 | +nan |
| X1 decreasing above 1x | 0.818 | 0.732/0.035/0.233 | 0.736/0.031/0.233 | 0.0042 | +0.00156 | +0.64 |
| X1s decrease then plateau | 0.421 | 0.193/0.050/0.757 | 0.198/0.045/0.757 | 0.0057 | +0.00170 | +0.51 |
| X2q unimodal, log-quadratic | 0.472 | 0.907/0.026/0.067 | 0.912/0.021/0.067 | 0.0053 | +0.00034 | +0.11 |
| X2b unimodal bump to plateau | 0.369 | 0.451/0.038/0.511 | 0.455/0.034/0.511 | 0.0035 | +0.00105 | +0.51 |
| X3 flat plus penalty below 1x | 0.0189 | 0.126/0.070/0.803 | 0.133/0.064/0.803 | 0.0067 | -0.00025 | -0.06 |
| X4 GP wiggle, flat trend | 0.173 | 0.348/0.083/0.569 | 0.368/0.063/0.569 | 0.0202 | +0.00674 | +0.57 |
| X4d GP wiggle, decreasing trend | 0.606 | 0.551/0.082/0.367 | 0.575/0.058/0.367 | 0.0241 | +0.01224 | +0.87 |

## Rule cross-check (verify_rule.py)

Vectorised rule 7 and reading G compared with vec.t2_build.count_override and _passing_ranking on 4,000 random 10-seed tables: winner mismatches 0; floor cases 1,039; G-reading mismatches 0; disagreement cases 112.
