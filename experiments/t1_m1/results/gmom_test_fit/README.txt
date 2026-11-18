Two target-free listings of 2026-11-18 (decision log 2026-11-18 (f)), made for the cold check of the claim "raise
g_mom's damping to 1.0" and for a note the captain forwarded the same day. Each reads the three released Task 1
stages only; neither builds a file, reads a board, a target or a score, or compares a prediction with anything.

sizes.txt       g_mom's finished factors at the test pair's fit (start E9.5, reference E10.5, as the builder fits
                them) against the practice fit (start E8.5, reference E9.5), at dampings 0.5, 0.75 and 1.0.
                Script as run: gmom_fit_sizes.py.txt (run from the repository root on the local machine, 16 threads).
collection.txt  The share of each category's factors' squared size that goes with a metacell's share of cells from
                each library, with a permutation null, for both fits. Script as run: gmom_collection.py.txt (run in
                the DGX's first campaign checkout at 6214e28, 16 threads).

The scripts are kept as text beside their outputs: they are records of two listings, not maintained code.
