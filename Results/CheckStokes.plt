set terminal qt size 1728, 972
set multiplot layout 2,3
set xrange[-0.5:0.5]
set yrange[-0.12:0.12]

plot \
    "full_mu2_N100.out" u 1:4 w l, \
    "full_mu2_N200.out" u 1:4 w l, \
    "full_mu2_N400.out" u 1:4 w l, \
    "full_mu2_N1000.out" u 1:4 w l, \
    0.1*erf(x/(2*sqrt(1e-2)))

plot \
    "full_mu3_N100.out" u 1:4 w l, \
    "full_mu3_N200.out" u 1:4 w l, \
    "full_mu3_N400.out" u 1:4 w l, \
    "full_mu3_N1000.out" u 1:4 w l, \
    0.1*erf(x/(2*sqrt(1e-3)))

plot \
    "full_mu4_N100.out" u 1:4 w l, \
    "full_mu4_N200.out" u 1:4 w l, \
    "full_mu4_N400.out" u 1:4 w l, \
    "full_mu4_N1000.out" u 1:4 w l, \
    0.1*erf(x/(2*sqrt(1e-4)))

plot \
    "simp_mu2_N100.out" u 1:4 w l, \
    "simp_mu2_N200.out" u 1:4 w l, \
    "simp_mu2_N400.out" u 1:4 w l, \
    "simp_mu2_N1000.out" u 1:4 w l, \
    0.1*erf(x/(2*sqrt(1e-2)))

plot \
    "simp_mu3_N100.out" u 1:4 w l, \
    "simp_mu3_N200.out" u 1:4 w l, \
    "simp_mu3_N400.out" u 1:4 w l, \
    "simp_mu3_N1000.out" u 1:4 w l, \
    0.1*erf(x/(2*sqrt(1e-3)))

plot \
    "simp_mu4_N100.out" u 1:4 w l, \
    "simp_mu4_N200.out" u 1:4 w l, \
    "simp_mu4_N400.out" u 1:4 w l, \
    "simp_mu4_N1000.out" u 1:4 w l, \
    0.1*erf(x/(2*sqrt(1e-4)))

pause(-1)
