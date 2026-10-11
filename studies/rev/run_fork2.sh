#!/bin/zsh
# Fork-2 revision analyses, exact quadrature, in sequence (5 workers each).
cd ${0:a:h}
PY=../../.venv/bin/python
for s in "p2_exact.py weighting" "p2_exact.py pole" "p2_exact.py limb" "p2s_study.py" "t3_study.py" \
         "d13_photon_energy.py" "d12b_limb_excluded.py" "d12c_mu_min.py" "e1_bias.py" \
         "d5b_interp_mu.py" "d2_d6_full.py"; do
  echo "=== $s $(date +%T)"
  eval $PY $s > /dev/null 2> run_fork2_err.log || { echo "FAILED $s"; cat run_fork2_err.log; }
done
echo "=== all done $(date +%T)"
