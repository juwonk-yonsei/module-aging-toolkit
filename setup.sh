#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# setup.sh -- prepare the module-aging-toolkit for use.
#
# v0.3.0 checks the checksums of the author-licensed inputs (Supplementary
# Table S20) and does not retrain module clocks. The v0.2.0 retraining path
# is withdrawn. Public GEO/GTEx raw data are not downloaded here.
#
# Usage:  MAT_EXPRESSION=... MAT_CLOCKS=... MAT_TABLE5B=... bash setup.sh --verify-only
# ---------------------------------------------------------------------------
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

if [ "${1:-}" = "--verify-only" ]; then
  python scripts/setup_v030.py --verify-only
  exit $?
fi

TAGE_DIR="${TAGE_DIR:-$ROOT/third_party/tAge}"
ZENODO_RECORD="18763485"

say() { printf '\n\033[1;34m==>\033[0m %s\n' "$*"; }

# 1) Python deps -------------------------------------------------------------
say "Installing Python dependencies (requirements.txt)"
python -m pip install -r requirements.txt
# make the module_aging package importable
python -m pip install -e . || echo "  (editable install skipped; using src/ via PYTHONPATH)"

# 2) tAge dependency (MGB Open Access License 1.0) ---------------------------
if [ ! -d "$TAGE_DIR/.git" ]; then
  say "Cloning Gladyshev-Lab/tAge into third_party/tAge"
  git clone --depth 1 https://github.com/Gladyshev-Lab/tAge.git "$TAGE_DIR"
else
  say "tAge already present at $TAGE_DIR (skipping clone)"
fi
echo "  NOTE: tAge is licensed under the MGB Open Access License 1.0"
echo "        (non-commercial / academic). See THIRD_PARTY_NOTICES.md."

# 3) Module-clock training data from Zenodo ----------------------------------
mkdir -p data/rodent data/supp data/clock_models
say "Downloading module-clock training data (Zenodo record $ZENODO_RECORD)"
if command -v zenodo_get >/dev/null 2>&1; then
  ( cd data/rodent && zenodo_get "$ZENODO_RECORD" )
else
  echo "  'zenodo_get' not found. Install it (pip install zenodo_get) OR download"
  echo "  record $ZENODO_RECORD manually from:"
  echo "     https://zenodo.org/records/$ZENODO_RECORD"
  echo "  and place these files as follows:"
  echo "     data/rodent/Expression_data_relative_rodents_Scaled.csv"
  echo "     data/rodent/Data_annotation_relative_rodents.xlsx"
  echo "     data/supp/module_membership_rodent.csv"
  echo "     data/clock_models/EN_*.pkl   (base tAge clocks)"
fi

# 4) Checksums of the author-licensed inputs. No retraining. -----------------
say "Checking S20 checksums (no download, no clock extraction)"
python scripts/setup_v030.py --verify-only

# 5) Sanity check ------------------------------------------------------------
say "Verifying shipped results against the v0.3.0 claim table"
python verify/verify_claims.py --r2 || true

say "Done. Next:"
echo "  * Reuse the method:  python examples/reuse_quickstart.py"
echo "  * Reproduce figures: see docs/REPRODUCE.md (needs GEO/GTEx data)"
