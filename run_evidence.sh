#!/usr/bin/env bash
set -euo pipefail
here="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
base="${PAPER3_ROOT:-$HOME/git/repository}"
py="${PAPER3_PYTHON:-}"
if [[ -z "$py" ]]; then
  known="$base/paper1_public_evidence/paper3_mechanism_source_103eceef0bec/.venv/bin/python"
  if [[ -x "$known" ]]; then py="$known"; else py=python; fi
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
opts=()
if [[ "${PAPER3_OFFLINE:-0}" != 1 ]]; then opts+=(--allow-downloads); fi
exec "$py" -u "$here/reproduce.py" evidence --root "$base" "${opts[@]}" "$@"
