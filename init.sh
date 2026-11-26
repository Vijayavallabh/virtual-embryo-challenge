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
echo "1. Read session-handoff.md (current state) and AGENTS.md (rules)"
echo "2. Register a try in experiments/harness/TRIES.md before running it; date decisions in docs/PLAN.md"
echo "3. Re-run ./init.sh before claiming done"
