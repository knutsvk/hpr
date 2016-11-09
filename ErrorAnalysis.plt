set xlabel "no of cells"
set ylabel "error"
set logscale xy

error_file = "./Results/ConvergenceErrors.out"
n = 4               # 2: 1norm, 3: 2norm, 4: infnorm
FIT_LIMIT = 1e-6

f(x) = a*x**b
fit f(x) error_file u 1:n via a,b

plot error_file u 1:n, f(x)
pause(-1)
