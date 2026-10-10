#!/bin/sh
# P2 studies run one after another (each uses all cores)
cd "$(dirname "$0")/../.."
for s in p2d_limb p2b_constraints p2c_hybrid2_pole p2f_fisher; do
  echo "=== $s $(date -u +%FT%TZ)"
  .venv/bin/python -u studies/p2/$s.py > studies/p2/$s.log 2>&1 || echo "FAILED $s"
done
echo "=== queue done $(date -u +%FT%TZ)"
