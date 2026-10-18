"""Start the v2 raw directory from the first run: copy every row whose cells the block-centre fix does not change
(random, jitter, refstage), drop block/head/refblock/refhead so that run_density.py --patterns block head refblock
refhead re-scores them. head and refhead do not depend on the centre; re-scoring them checks that nothing else moved.

Usage (DGX, repo root): uv run python experiments/density/seed_v2.py --old runs/density/raw --new runs/density/raw_v2
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

RERUN = {"block", "head", "refblock", "refhead"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", default="runs/density/raw")
    ap.add_argument("--new", default="runs/density/raw_v2")
    args = ap.parse_args()
    new = Path(args.new)
    new.mkdir(parents=True, exist_ok=True)
    for p in sorted(Path(args.old).glob("*.jsonl")):
        out = new / p.name
        if out.exists():
            raise SystemExit(f"{out} exists; seed into an empty directory")
        rows = [line for line in p.read_text().splitlines() if line.strip()]
        keep = [line for line in rows if json.loads(line)["pattern"] not in RERUN]
        out.write_text("".join(line + "\n" for line in keep))
        print(f"{p.name}: kept {len(keep)} of {len(rows)} rows")


if __name__ == "__main__":
    main()
