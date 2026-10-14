# hybrid4 vs the Kipping (2013) conditions

* poles: (0.028601207249486824, 0.38582417584735224)
* T_j(0) = [1. 1. 1. 1.], T_j(1) = [0. 0. 0. 0.]; max dT_j/dxi over [0, 1] = 0 (must be <= 0)
* simplex vertices: 5 sets; violate (I) 0, violate (II) 0; min I = -2.22e-16, min dI/dxi = -0
* edge midpoints: 25 sets; violate (I) 0, violate (II) 0; min I = -2.22e-16, min dI/dxi = -0
* random interior (flat Dirichlet): 100000 sets; violate (I) 0, violate (II) 0; min I = 5.17e-07, min dI/dxi = 0.00073
* random points on the face sum w = 1: 100000 sets; violate (I) 0, violate (II) 0; min I = -4.44e-16, min dI/dxi = 0.0011
* random points on the faces w_j = 0: 100000 sets; violate (I) 0, violate (II) 0; min I = 0.00211, min dI/dxi = 3.13e-05
* q-map: min w = 0, max sum w = 1.000000000000000
* q-map images of q in [0,1]^4 (incl. corners): 50016 sets; violate (I) 0, violate (II) 0; min I = -2.22e-16, min dI/dxi = -0
* phoenix table: 134560 rows; min w = -0, max sum w = 1.000001, rows with sum w = 1 (zero limb intensity): 6
* phoenix tabulated hybrid4 (r and mu methods): 134560 sets; violate (I) 1, violate (II) 0; min I = -1e-06, min dI/dxi = 0.0156
* mps1 table: 1391868 rows; min w = 0.0085, max sum w = 0.922830, rows with sum w = 1 (zero limb intensity): 0
* mps1 tabulated hybrid4 (r and mu methods): 1391868 sets; violate (I) 0, violate (II) 0; min I = 0.0772, min dI/dxi = 0.0188
* mps2 table: 1391868 rows; min w = 0.00901, max sum w = 0.914385, rows with sum w = 1 (zero limb intensity): 0
* mps2 tabulated hybrid4 (r and mu methods): 1391868 sets; violate (I) 0, violate (II) 0; min I = 0.0856, min dI/dxi = 0.0211
* mps2 + ssd: outside simplex 0 of 645012
* mps2 + ssd corrected sets: 645012 sets; violate (I) 0, violate (II) 0; min I = 0.0915, min dI/dxi = 0.033
* mps2 + B100G: outside simplex 0 of 645012
* mps2 + B100G corrected sets: 645012 sets; violate (I) 0, violate (II) 0; min I = 0.092, min dI/dxi = 0.0318
* mps2 + B200G: outside simplex 0 of 645012
* mps2 + B200G corrected sets: 645012 sets; violate (I) 0, violate (II) 0; min I = 0.0954, min dI/dxi = 0.0292
* mps2 + B300G: outside simplex 0 of 645012
* mps2 + B300G corrected sets: 645012 sets; violate (I) 0, violate (II) 0; min I = 0.0977, min dI/dxi = 0.0246
* phoenix interpolated (random points): 3000 sets; violate (I) 0, violate (II) 0; min I = 0.0569, min dI/dxi = 0.0634
* mps2 interpolated (random points): 3000 sets; violate (I) 0, violate (II) 0; min I = 0.241, min dI/dxi = 0.064
* standalone files (random points): 4000 sets; violate (I) 0, violate (II) 0; min I = 0.0426, min dI/dxi = 0.108
