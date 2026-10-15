#!/bin/zsh
# hybrid4 curvature-simplex basis: production update chain (2026-10-14)
set -e
cd /Users/dkipping/Storage1/Work/Documents/Transit_Work/CODES/LDhybrid
PY=.venv/bin/python
echo "== bands $(date)"; $PY -W ignore scripts/p3_update_law.py hybrid4 phoenix mps1 mps2 --workers 6
echo "== spectra $(date)"; $PY -W ignore scripts/p3_spectra_update.py hybrid4 --workers 6
echo "== corrections $(date)"; $PY -W ignore scripts/p3_corrections.py
echo "== assemble $(date)"; $PY -W ignore scripts/p3_assemble.py phoenix mps1 mps2
echo "== cds $(date)"; $PY -W ignore scripts/p5_cds.py
echo "== g3 $(date)"; $PY -W ignore studies/g3/g3_checks.py
echo "== stats $(date)"; $PY -W ignore studies/rev/production_stats.py
echo "== CHAIN DONE $(date)"
