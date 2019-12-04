# Gnuplot: plot RTS behaviour for a given region of parameter space

reset
r = 0.4
unset key
set cbrange [-r:r]
set palette defined (-r "#8888FF", 0 "#FFFFFF", r "#FF8888")
set pm3d map #corners2color c1
set xrange [-0.08:0.05]
set yrange [0.:2]
lambda1 = 0.3
delta = -0.01
set ylabel "Selectivity"
set xlabel "Relative foreign selfishness"

splot "mtdna-mito-zoom-overall-0.300.txt" u ($1-lambda1-delta/2.):2:(2*($3-0.5))
