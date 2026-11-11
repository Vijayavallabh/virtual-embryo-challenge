#!/usr/bin/env python
"""Print one window of a web page's visible text from its HTML on stdin, for a screener (decision 2026-11-11 (c)).

Stdlib only; reads stdin, writes stdout, never a file. The agent pipes one plain curl request into it:

  curl -sS -w '\n[HTTP %{http_code}]\n' URL | python -I -B experiments/intake/page_text.py OFFSET LENGTH

Scripts, styles, templates, comments and markup are dropped, so no inline application payload or other page's data is
printed. It prints the HTTP status curl appended, the total length of the visible text, and the characters
[offset, offset + length), so a long page is read in windows small enough that no tool output is saved to a file.
"""

from __future__ import annotations

import html
import re
import sys

STATUS = re.compile(r"\[HTTP (\d{3})\]\s*$")


def visible_text(page: str) -> str:
    t = re.sub(r"(?is)<(script|style|noscript|template)\b[^>]*>.*?</\1\s*>", " ", page)
    t = re.sub(r"(?is)<(script|style|noscript|template)\b[^>]*>.*\Z", " ", t)  # an unclosed one runs to the end
    t = re.sub(r"(?s)<!--.*?-->", " ", t)
    t = re.sub(r"(?i)<br\s*/?>|</(p|div|li|h[1-6]|tr|section|article|header|footer)\s*>", "\n", t)
    t = re.sub(r"(?s)<[^>]+>", " ", t)
    t = html.unescape(t)
    t = re.sub(r"[ \t\r\f\v]+", " ", t)
    return re.sub(r"\n\s*\n+", "\n", t).strip()


def main() -> None:
    if len(sys.argv) != 3:
        sys.exit("usage: page_text.py OFFSET LENGTH < page.html")
    offset, length = int(sys.argv[1]), int(sys.argv[2])
    raw = sys.stdin.read()
    m = STATUS.search(raw)
    status = m.group(1) if m else "unknown"
    out = visible_text(raw[: m.start()] if m else raw)
    print(f"[HTTP {status}; visible text length {len(out)}; showing {offset}-{min(offset + length, len(out))}]")
    print(out[offset: offset + length])


if __name__ == "__main__":
    main()
