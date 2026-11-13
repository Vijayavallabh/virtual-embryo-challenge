#!/bin/bash
set -e

echo "=== Harness Initialization ==="

echo "=== uv sync --locked ==="
uv sync --locked

echo "=== uv run ruff check src tests ==="
uv run ruff check src tests

echo "=== uv run pytest -q ==="
uv run pytest -q

echo "=== Verification Complete ==="
echo ""
echo "Next steps:"
echo "0. Read session-handoff.md: it says what waits. In the final phase the uploads follow docs/P3_RUNBOOK.md"
echo "   (section 5), each on the captain's go-ahead, and only files in the current reference list"
echo "1. Read feature_list.json to see current feature state"
echo "2. Pick ONE unfinished feature to work on"
echo "3. Before choosing or changing work on an entry: the guiding principle in AGENTS.md, then the headroom"
echo "   tables and rules in experiments/harness/TRIES.md (register the try before its run)"
echo "4. Implement only that feature"
echo "5. Re-run verification before claiming done"
