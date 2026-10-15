#!/bin/zsh
# hybrid4 curvature-simplex basis: reruns of the studies that feed the paper (after the production chain)
set -e
cd /Users/dkipping/Storage1/Work/Documents/Transit_Work/CODES/LDhybrid
PY=.venv/bin/python
echo "== phoenix refit with face snapping $(date)"
mv data/products/bands/phoenix/.hybrid4_updated data/products/bands/phoenix/.hybrid4_updated.presnap
$PY -W ignore scripts/p3_update_law.py hybrid4 phoenix --workers 6
$PY -W ignore scripts/p3_assemble.py phoenix
$PY -W ignore scripts/p5_cds.py
echo "== standalone $(date)"; $PY -W ignore scripts/make_standalone.py
echo "== kipping $(date)"; $PY -W ignore studies/rev3/h4_kipping_check.py > studies/rev3/h4_kipping_check.log 2>&1
echo "== bound $(date)"; $PY -W ignore studies/rev4/rv_bound.py > studies/rev4/rv_bound.log 2>&1
echo "== validate $(date)"; $PY -W ignore studies/rev4/rv_validate.py > studies/rev4/rv_validate.log 2>&1
echo "== p2 weighting $(date)"; $PY -W ignore studies/rev4/rv_p2.py weighting > studies/rev4/rv_p2_weighting.log 2>&1
echo "== p2 limb $(date)"; $PY -W ignore studies/rev4/rv_p2.py limb > studies/rev4/rv_p2_limb.log 2>&1
echo "== bias $(date)"; $PY -W ignore studies/rev4/rv_bias.py > studies/rev4/rv_bias.log 2>&1
echo "== fisher $(date)"; $PY -W ignore studies/rev4/rv_fisher.py > studies/rev4/rv_fisher.log 2>&1
echo "== lc section $(date)"; $PY -W ignore studies/lc/lc_section.py > studies/rev4/rv_lc_section.log 2>&1
echo "== REV4 DONE $(date)"
