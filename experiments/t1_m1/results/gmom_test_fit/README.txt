Target-free listings of 2026-11-18 (decision log 2026-11-18 (f) to (h)), the first two made for the cold check of the claim "raise
g_mom's damping to 1.0" and for a note the captain forwarded the same day. Each reads the three released Task 1
stages only; neither builds a file, reads a board, a target or a score, or compares a prediction with anything.

sizes.txt       g_mom's finished factors at the test pair's fit (start E9.5, reference E10.5, as the builder fits
                them) against the practice fit (start E8.5, reference E9.5), at dampings 0.5, 0.75 and 1.0.
                Script as run: gmom_fit_sizes.py.txt (run from the repository root on the local machine, 16 threads).
collection.txt  The share of each category's factors' squared size that goes with a metacell's share of cells from
                each library, with a permutation null, for both fits. Script as run: gmom_collection.py.txt (run in
                the DGX's first campaign checkout at 6214e28, 16 threads).

The scripts are kept as text beside their outputs: they are records of two listings, not maintained code.

emulation.txt   (decision log 2026-11-18 (g)) The practice fit refitted with the smaller library of its reference
                thinned to two thirds and to one third: the library-linked share of the factors stays at its null.
                An independent review showed that this design could not have shown an effect (entry (h)): it
                says nothing about depth. Script as run: gmom_emulate.py.txt (the
                DGX's first campaign checkout at 6214e28, 16 threads).

single_kind.txt (decision log 2026-11-18 (g), which fixed its use before it ran, and (h)) The test pair's fit on the
                reference as built (ALL), on library 3 only (L3) and on libraries 1 and 2 only (L12), and the
                practice fit (BACK), compared by reading N21's measures. Script as run: gmom_single_kind.py.txt
                (the DGX's first campaign checkout at 6214e28, 16 threads).

within_metacell.txt (decision log 2026-11-18 (h), which fixed its use before it ran, and (i)) Inside the metacells of
                the fit as built that hold at least ten cells of each kind, library 3's cells against the cells of
                libraries 1 and 2, projected on the library slope; a joint regression on the two portion
                libraries' shares; the genes with the largest slope and the largest within-metacell contrast; the
                library statistic on a no-development fit of E10.5. Script as run: gmom_within.py.txt (the DGX's
                first campaign checkout at 6214e28, 16 threads).
