# Gnuplot: plot time series of RTS behaviour

reset
set multiplot
set size 0.5,0.33
set origin 0,0.66
set yrange [0:1]
set ylabel "Type 2 proportion"
set xlabel "Time"
set key top left
plot "mtdna-mito-trace-0.300-0.200.txt" u ($1 == 0 ? $2 : 1/0):3 w l lw 4 lc rgbcolor "#FF0000" title "P = 0", "" u ($1 == 0.5 ? $2 : 1/0):3 w l lw 4 lc rgbcolor "#FF8888" title "P = 0.5", "" u ($1 == 1 ? $2 : 1/0):3 w l lw 4 lc rgbcolor "#AAAAAA" title "P = 1", "" u ($1 == 1.5 ? $2 : 1/0):3 w l lw 4 lc rgbcolor "#8888FF" title "P = 1.5", "" u ($1 == 2 ? $2 : 1/0):3 w l lw 4 lc rgbcolor "#0000FF" title "P = 2"

unset key
set origin 0.5,0.66
plot "mtdna-mito-trace-0.300-0.400.txt" u ($1 == 0 ? $2 : 1/0):3 w l lw 4 lc rgbcolor "#FF0000" title "P = 0", "" u ($1 == 0.5 ? $2 : 1/0):3 w l lw 4 lc rgbcolor "#FF8888" title "P = 0.5", "" u ($1 == 1 ? $2 : 1/0):3 w l lw 4 lc rgbcolor "#AAAAAA" title "P = 1", "" u ($1 == 1.5 ? $2 : 1/0):3 w l lw 4 lc rgbcolor "#8888FF" title "P = 1.5", "" u ($1 == 2 ? $2 : 1/0):3 w l lw 4 lc rgbcolor "#0000FF" title "P = 2"

set yrange [0:1000]
set ylabel "Total n"
set origin 0,0.33
plot "mtdna-mito-trace-0.300-0.200.txt" u ($1 == 0 ? $2 : 1/0):5 w l lw 4 lc rgbcolor "#FF0000" title "P = 0", "" u ($1 == 0.5 ? $2 : 1/0):5 w l lw 4 lc rgbcolor "#FF8888" title "P = 0.5", "" u ($1 == 1 ? $2 : 1/0):5 w l lw 4 lc rgbcolor "#AAAAAA" title "P = 1", "" u ($1 == 1.5 ? $2 : 1/0):5 w l lw 4 lc rgbcolor "#8888FF" title "P = 1.5", "" u ($1 == 2 ? $2 : 1/0):5 w l lw 4 lc rgbcolor "#0000FF" title "P = 2"

set origin 0.5,0.33
plot "mtdna-mito-trace-0.300-0.400.txt" u ($1 == 0 ? $2 : 1/0):5 w l lw 4 lc rgbcolor "#FF0000" title "P = 0", "" u ($1 == 0.5 ? $2 : 1/0):5 w l lw 4 lc rgbcolor "#FF8888" title "P = 0.5", "" u ($1 == 1 ? $2 : 1/0):5 w l lw 4 lc rgbcolor "#AAAAAA" title "P = 1", "" u ($1 == 1.5 ? $2 : 1/0):5 w l lw 4 lc rgbcolor "#8888FF" title "P = 1.5", "" u ($1 == 2 ? $2 : 1/0):5 w l lw 4 lc rgbcolor "#0000FF" title "P = 2"

set logscale 
set yrange [*:*]
set ylabel "V'(h)"
set origin 0,0
plot "mtdna-mito-trace-0.300-0.200.txt" u ($1 == 0 ? $2 : 1/0):(($4**2)/($3*(1-$3))) w l lw 4 lc rgbcolor "#FF0000" title "P = 0", "" u ($1 == 0.5 ? $2 : 1/0):(($4**2)/($3*(1-$3))) w l lw 4 lc rgbcolor "#FF8888" title "P = 0.5", "" u ($1 == 1 ? $2 : 1/0):(($4**2)/($3*(1-$3))) w l lw 4 lc rgbcolor "#AAAAAA" title "P = 1", "" u ($1 == 1.5 ? $2 : 1/0):(($4**2)/($3*(1-$3))) w l lw 4 lc rgbcolor "#8888FF" title "P = 1.5", "" u ($1 == 2 ? $2 : 1/0):(($4**2)/($3*(1-$3))) w l lw 4 lc rgbcolor "#0000FF" title "P = 2", x/1000.

set origin 0.5,0
plot "mtdna-mito-trace-0.300-0.400.txt" u ($1 == 0 ? $2 : 1/0):(($4**2)/($3*(1-$3))) w l lw 4 lc rgbcolor "#FF0000" title "P = 0", "" u ($1 == 0.5 ? $2 : 1/0):(($4**2)/($3*(1-$3))) w l lw 4 lc rgbcolor "#AA4444" title "P = 0.5", "" u ($1 == 1 ? $2 : 1/0):(($4**2)/($3*(1-$3))) w l lw 4 lc rgbcolor "#888888" title "P = 1", "" u ($1 == 1.5 ? $2 : 1/0):(($4**2)/($3*(1-$3))) w l lw 4 lc rgbcolor "#4444FF" title "P = 1.5", "" u ($1 == 2 ? $2 : 1/0):(($4**2)/($3*(1-$3))) w l lw 4 lc rgbcolor "#0000FF" title "P = 2", x/1000.
