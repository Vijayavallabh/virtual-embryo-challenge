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
echo "1. Read feature_list.json to see current feature state"
echo "2. Pick ONE unfinished feature to work on"
echo "3. Implement only that feature"
echo "4. Re-run verification before claiming done"
