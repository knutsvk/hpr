""" Compute error norms for Convergence studies """
from math import sqrt
import numpy as np
from numpy import linalg as la
from scipy.special import erf

print("N\t||error||_1\t\t||error||_2\t\t||error||_inf")
for i in [100]:
    exact = np.genfromtxt('./Results/ConvergenceStudies_2D_Nx' + str(i) + '_Ny'
                          + str(i) + '_Exact.out', delimiter='\t', dtype=float, names=True)
    compd = np.genfromtxt('./Results/ConvergenceStudies_2D_Nx' + str(i) + '_Ny'
                          + str(i) + '_100.out', delimiter='\t', dtype=float, names=True)
    eps = exact['rho'] - compd['rho']
    one_norm = la.norm(eps, 1) / la.norm(exact, 1)
    two_norm = la.norm(eps) / la.norm(exact)
    inf_norm = la.norm(eps, np.inf) / la.norm(exact, np.inf)
    print(i, "\t", one_norm, "\t", two_norm, "\t", inf_norm)
