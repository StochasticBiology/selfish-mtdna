# Gnuplot: plot trellis of RTS behaviour with different lambda1 values
# suggested terminal: png 1600 1000

reset
set multiplot
set size 0.33

unset key

set cbrange [-1:1]
set palette defined (-1 "#AAAAFF", 0 "#FFFFFF", 1 "#FFAAAA")
set pm3d map
lambda1 = 0.4
set ylabel "Selective threshold P"
set xlabel "Relative foreign selfishness (λ₂ - λ₁)"

offset = 0

set object rectangle from graph 0,0 to graph 1,1 behind fillcolor rgbcolor "#FFFFDD"
set yrange [0:2.8]
set origin 0,0
lambda1 = 0.1
set label 1 at graph 0.05, graph 0.9 "λ₁ = 0.1" front
set xrange [-lambda1-offset:0.9-lambda1-offset]
splot "mtdna-mito-0.000-0.010-overall-0.100.txt" u ($1-lambda1):2:(2*($3-0.5))

set origin 0.33,0
lambda1 = 0.2
set label 1 at graph 0.05, graph 0.9 "λ₁ = 0.2" front
set xrange [-lambda1-offset:0.9-lambda1-offset]
splot "mtdna-mito-0.000-0.010-overall-0.200.txt" u ($1-lambda1):2:(2*($3-0.5))

set origin 0.66,0
lambda1 = 0.3
#set label 2 at graph 0.92,graph 0.85 "*" front
set label 1 at graph 0.05, graph 0.9 "λ₁ = 0.3" front
set xrange [-lambda1-offset:0.9-lambda1-offset]
splot "mtdna-mito-0.000-0.010-overall-0.300.txt" u ($1-lambda1):2:(2*($3-0.5))


set origin 0,0.33
lambda1 = 0.4
set label 1 at graph 0.05, graph 0.9 "λ₁ = 0.4" front
set xrange [-lambda1-offset:0.9-lambda1-offset]
splot "mtdna-mito-0.000-0.010-overall-0.400.txt" u ($1-lambda1):2:(2*($3-0.5))

set origin 0.33,0.33
lambda1 = 0.5
set label 1 at graph 0.05, graph 0.9 "λ₁ = 0.5" front
set xrange [-lambda1-offset:0.9-lambda1-offset]
splot "mtdna-mito-0.000-0.010-overall-0.500.txt" u ($1-lambda1):2:(2*($3-0.5))

set origin 0.66,0.33
lambda1 = 0.6
set label 1 at graph 0.05, graph 0.9 "λ₁ = 0.6" front
set xrange [-lambda1-offset:0.9-lambda1-offset]
splot "mtdna-mito-0.000-0.010-overall-0.600.txt" u ($1-lambda1):2:(2*($3-0.5))



set origin 0,0.66
lambda1 = 0.7
set label 1 at graph 0.05, graph 0.9 "λ₁ = 0.7" front
set xrange [-lambda1-offset:0.9-lambda1-offset]
splot "mtdna-mito-0.000-0.010-overall-0.700.txt" u ($1-lambda1):2:(2*($3-0.5))

set origin 0.33,0.66
lambda1 = 0.8
set label 1 at graph 0.05, graph 0.9 "λ₁ = 0.8" front
set xrange [-lambda1-offset:0.9-lambda1-offset]
splot "mtdna-mito-0.000-0.010-overall-0.800.txt" u ($1-lambda1):2:(2*($3-0.5))

set origin 0.66,0.66
lambda1 = 0.9
set label 1 at graph 0.05, graph 0.9 "λ₁ = 0.9" front
set xrange [-lambda1-offset:0.9-lambda1-offset]
splot "mtdna-mito-0.000-0.010-overall-0.900.txt" u ($1-lambda1):2:(2*($3-0.5))
