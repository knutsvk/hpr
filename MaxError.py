import numpy as np
from scipy import special

for i in [1, 2, 4, 8]:
    r = np.genfromtxt('./StokesFirstProblem_1DX_Nx' + str(10 * i) + '_Ny' +
            str(i) + '_100.out', delimiter='\t', dtype=float, names=True)
    eps = abs(r['v'] - 0.1 * special.erf(5.0 * r['x']))
    m = eps[eps.argmax()]
    m = np.mean(eps)
    print i, "\t", m
