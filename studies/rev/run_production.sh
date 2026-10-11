#!/bin/zsh
# Full v1 production with resample="exact" (revision), then profiles, corrections, assembly, CDS.
set -e
cd /Users/dkipping/Storage1/Work/Documents/Transit_Work/CODES/LDhybrid
PY=".venv/bin/python -W ignore"
for g in phoenix mps1 mps2; do
  eval $PY scripts/p3_tables.py $g --workers 5
done
eval $PY scripts/p3_extend.py mps2 --workers 5        # profiles only (bands complete)
eval $PY scripts/p3_corrections.py
eval $PY scripts/p3_assemble.py phoenix mps1 mps2
eval $PY scripts/p5_cds.py
echo "ALL DONE"
