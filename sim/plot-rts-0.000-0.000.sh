# Gnuplot: plot RTS behaviour for a given region of parameter space

reset
r = 0.2
unset key
set cbrange [-r:r]
set palette defined (-r "#8888FF", 0 "#FFFFFF", r "#FF8888")
set pm3d map #corners2color c1
set xrange [-0.024:0.021]
set xtics 0.01
set yrange [0.5:2.]
lambda1 = 0.3
delta = -0.005
set ylabel "Selective threshold P"
set xlabel "Relative foreign selfishness (λ₂ - λ₁)"

splot "mtdna-mito-0.000-0.000-zoom-overall-0.300.txt" u ($1-lambda1-delta/2.):2:(2*($3-0.5))
