"""What each reading picks across all floor cases (base scenario, x families), and E[D] over floor cases.
Run from the repo root: uv run --no-sync python -I <this file>"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sim_floor as S  # noqa: E402

base = S.Scenario()
lines = ["family | P(floor) | floor: G picks 0.5/0.75/max | floor: C picks 0.5/0.75/max | P(dis given floor) | "
         "E[D over all floor cases] | break-even test tilt (pts/doubling toward more cells)"]
for i, fam in enumerate(S.X_FAMILIES):
    n = 1_200_000 if fam.shape in ("inc", "knee") else 150_000
    acc = S.run(fam, base, n, 909 + i)
    nf = acc["n_floor"]
    if not nf:
        continue
    g = acc["n_floor_G"][1:] / nf
    c = acc["n_floor_C"][1:] / nf
    fd = acc["floor_D"].mean()
    be = acc["D"].mean() / 0.585 if acc["n_dis"] else float("nan")
    lines.append(f"{fam.name} | {nf / n:.3g} | {g[0]:.3f}/{g[1]:.3f}/{g[2]:.3f} | {c[0]:.3f}/{c[1]:.3f}/{c[2]:.3f} | "
                 f"{acc['n_dis'] / nf:.4f} | {fd:+.5f} | {be:+.2f}")
out = "\n".join(lines)
(Path(__file__).resolve().parent / "floor_composition.txt").write_text(out + "\n")
print(out)
