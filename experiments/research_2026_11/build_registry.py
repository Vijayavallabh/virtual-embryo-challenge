"""Build the global citation registry (deep-research P3) from the audited notes in final/ only.

Usage: python -I -B build_registry.py <final_dir> <out_md>
Reads final/task-*.md, parses each '## Sources' section, merges sources by URL (keeping the higher authority),
assigns global numbers in task order, maps each task's local [n] to the global number, and prints quality stats.
Lines marked [removed: X] are skipped and counted.
"""
import re
import sys
from pathlib import Path
from urllib.parse import urlparse

final_dir, out_md = Path(sys.argv[1]), Path(sys.argv[2])
LINE = re.compile(r"^\[(\d+)\]\s*(.+)$")
reg, order, local_map, removed, bad = {}, [], {}, 0, []
for f in sorted(final_dir.glob("task-*.md")):
    tid = f.stem.split("-", 1)[1]
    text = f.read_text()
    m = re.search(r"^## Sources\s*$(.*?)^## ", text, re.S | re.M)
    if not m:
        bad.append(f"{f.name}: no Sources section")
        continue
    for raw in m.group(1).splitlines():
        raw = raw.strip()
        if not raw:
            continue
        if "[removed" in raw:
            removed += 1
            continue
        lm = LINE.match(raw)
        if not lm:
            continue
        n, rest = lm.group(1), lm.group(2)
        parts = [p.strip() for p in rest.split("|")]
        url = next((p for p in parts if p.startswith("http")), None)
        if not url:
            bad.append(f"{f.name} [{n}]: no URL")
            continue
        fields = {}
        for p in parts:
            if ":" in p and not p.startswith("http"):
                k, v = p.split(":", 1)
                fields[k.strip().lower()] = v.strip()
        auth = re.match(r"(\d+)", fields.get("authority", "0"))
        auth = int(auth.group(1)) if auth else 0
        key = url.rstrip("/")
        entry = {"title": parts[0], "url": url, "type": fields.get("source-type", "other"),
                 "access": fields.get("accessibility", "public"), "asof": fields.get("as of", "undated"),
                 "auth": auth, "tasks": {tid}}
        if key in reg:
            reg[key]["tasks"].add(tid)
            if auth > reg[key]["auth"]:
                reg[key].update({k: entry[k] for k in ("title", "type", "access", "asof", "auth")})
        else:
            reg[key] = entry
            order.append(key)
        local_map[(tid, n)] = key
gnum = {k: i + 1 for i, k in enumerate(order)}
types = {}
for k in order:
    types[reg[k]["type"]] = types.get(reg[k]["type"], 0) + 1
domains = {urlparse(reg[k]["url"]).netloc.lower().removeprefix("www.") for k in order}
dom_count = {}
for k in order:
    d = urlparse(reg[k]["url"]).netloc.lower().removeprefix("www.")
    dom_count[d] = dom_count.get(d, 0) + 1
lines = ["# Citation registry (built from the audited notes in final/ only)", "", "## Approved", ""]
for k in order:
    e = reg[k]
    lines.append(f"[{gnum[k]}] {e['title']} | {e['url']} | Source-Type: {e['type']} | Accessibility: {e['access']} | "
                 f"As Of: {e['asof']} | Auth: {e['auth']} | From: {','.join(sorted(e['tasks']))}")
lines += ["", "## Local-to-global map", ""]
for (tid, n), k in sorted(local_map.items(), key=lambda x: (x[0][0], int(x[0][1]))):
    lines.append(f"task-{tid} [{n}] -> [{gnum[k]}]")
total = len(order)
off = types.get("official", 0) + types.get("academic", 0)
lines += ["", "## Stats", "",
          f"Approved: {total}; removed source lines skipped: {removed}; unique domains: {len(domains)}",
          "Source types: " + ", ".join(f"{t} {c}" for t, c in sorted(types.items())),
          f"Official+academic share: {100 * off / max(total, 1):.0f}%",
          f"Max single-domain share: {100 * max(dom_count.values(), default=0) / max(total, 1):.0f}% "
          f"({max(dom_count, key=dom_count.get) if dom_count else '-'})",
          "Parse problems: " + ("; ".join(bad) if bad else "none")]
out_md.write_text("\n".join(lines) + "\n")
print("\n".join(lines[-6:]))
