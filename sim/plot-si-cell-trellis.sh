# Gnuplot: plot trellis of cell-level RTS behaviour for different parameterisations
# Suggested terminal: png 1200 800

reset
set multiplot
set size 0.33,0.5
r = 0.2
set cbrange [-r:r]
set palette defined (-r "#AAAAFF", 0 "white", r "#FFAAAA")
set pm3d map corners2color c1
set xlabel "Relative foreign\nselfishness"
set ylabel "Cell-level selectivity"
unset key
set xtics 0.1
set origin 0,0

set origin 0,0.5

set label 1 at graph 0.05, graph 0.9 "λ1 = 0.2" front
splot "mtdna-cell-1-10-0.000-0.2.txt" u ($1-0.2):2:($3 == -1 ? 1/0 : $4-$3)
set origin 0.33,0.5
set label 1 at graph 0.05, graph 0.9 "λ1 = 0.5" front
splot "mtdna-cell-1-10-0.000-0.5.txt" u ($1-0.5):2:($3 == -1 ? 1/0 : $4-$3)
set origin 0.66,0.5
set label 1 at graph 0.05, graph 0.9 "λ1 = 0.8" front
splot "mtdna-cell-1-10-0.000-0.8.txt" u ($1-0.8):2:($3 == -1 ? 1/0 : $4-$3)

set origin 0,0.
set label 1 at graph 0.05, graph 0.9 "λ1 = 0.2" front
splot "mtdna-cell-0-10-0.000-0.2.txt" u ($1-0.2):2:($3 == -1 ? 1/0 : $4-$3)
set origin 0.33,0.
set label 1 at graph 0.05, graph 0.9 "λ1 = 0.5" front
splot "mtdna-cell-0-10-0.000-0.5.txt" u ($1-0.5):2:($3 == -1 ? 1/0 : $4-$3)
set origin 0.66,0.
set label 1 at graph 0.05, graph 0.9 "λ1 = 0.8" front
splot "mtdna-cell-0-10-0.000-0.8.txt" u ($1-0.8):2:($3 == -1 ? 1/0 : $4-$3)
