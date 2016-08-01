set xrange [0.4:20]
set yrange [0.0008:0.06]
set xlabel "no of cells"
set ylabel "max error"
set logscale xy

FIT_LIMIT = 1e-6
f(x) = a*x**b
fit f(x) "ErrorStrang.out" via a,b

g(x) = x*c**d
fit g(x) "ErrorNoStrang.out" via c,d


plot "ErrorStrang.out" u 1:2, f(x), "ErrorNoStrang.out" u 1:2, g(x)
pause(-1)
