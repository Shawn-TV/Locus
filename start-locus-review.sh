#!/usr/bin/env bash
set -e

cd "$(dirname "$0")"

echo "Starting Locus Review Editor..."

PYTHON_CMD="${PYTHON:-python3}"

if ! command -v "$PYTHON_CMD" >/dev/null 2>&1; then
  echo "Python 3 was not found. Please install Python 3.10 or newer, then run this file again."
  exit 1
fi

if [ ! -x ".venv/bin/python" ]; then
  echo "Creating local Python environment..."
  "$PYTHON_CMD" -m venv .venv
fi

if [ ! -f ".venv/.locus_review_deps_ok" ]; then
  echo "Installing required packages. This may take a few minutes the first time..."
  ".venv/bin/python" -m pip install --upgrade pip
  ".venv/bin/python" -m pip install -r requirements-review.txt
  touch ".venv/.locus_review_deps_ok"
fi

".venv/bin/python" tools/launch_review_editor.py
