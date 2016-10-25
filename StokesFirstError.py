""" Compute error norms for Stokes first problem """
from math import sqrt
import numpy as np
from numpy import linalg as la
from scipy.special import erf

MU = 1e-4

print("N\t||error||_1\t\t||error||_2\t\t||error||_inf")
for i in [100, 200, 400, 1000]:
    r = np.genfromtxt('./Results/StokesFirstProblem_1DX_Nx' + str(i) +
                     '_Ny1_100.out', delimiter='\t', dtype=float, names=True)
    analytic = 0.1 * erf(r['x'] / (2 * sqrt(MU)))
    eps = r['v'] - analytic
    one_norm = la.norm(eps, 1) / la.norm(analytic, 1)
    two_norm = la.norm(eps) / la.norm(analytic)
    inf_norm = la.norm(eps, np.inf) / la.norm(analytic, np.inf)
    print(i, "\t", one_norm, "\t", two_norm, "\t", inf_norm)
