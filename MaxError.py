import numpy as np
from scipy import special
from numpy import linalg as la

for i in [1, 2, 4, 8]:
    r = np.genfromtxt('./StokesFirstProblem_1DX_Nx' + str(10 * i) + '_Ny' +
            str(i) + '_100.out', delimiter='\t', dtype=float, names=True)
    eps = abs(r['v'] - (0.1 * special.erf(5.0 * r['x'])) )
    m1 = eps[eps.argmax()]
    m2 = np.mean(eps)
    m3 = la.norm(eps)
    print i, "\t", 1.0/i, "\t", m1, "\t", m2, "\t", m3
