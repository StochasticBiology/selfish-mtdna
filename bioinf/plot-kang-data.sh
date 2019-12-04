# Gnuplot: plot classification of reverting vs non-reverting sequences from analysis of Kang et al. data

reset
set xlabel "CSBII difference"
set ylabel "TAS difference"
set xtics 1
set ytics 1
unset key
set xrange [-2.4:1.4]
set yrange [-2.4:1.4]
set palette defined (0 "#EEEEEE" , 0.5 "#FFFFFF", 1 "#FFDDDD")
set view map
set pm3d
unset surface
set isosamples 40
set cbrange [0:1]
unset colorbox

# functions to jitter datapoints according to their label
fx(x) = -0.05*((x/5)-3.)
fy(x) = -0.065*(x-5*floor(x/5)-2.5)

# parameters for (absence of) shift, and pointsize
myshift = 0
myps = 2

# parameters from logistic regression (see accompanying R code)
intercept = -2.463
l_300_diff = 1.661
l_16345_diff = -1.701

# first plot is the inferred regression function; subsequent plots are the data. 6 plots: 3 each for reverting and non-reverting types ($5 == 1 and $5 == 0 respectively); one for point, one for outline, one for label
splot exp(l_300_diff*x + l_16345_diff*y + intercept)/(1. + exp(l_300_diff*x + l_16345_diff*y + intercept)) lw 0, "kang-competition.txt" u ($5 == 0 ? $2+fx($1) : 1/0):($4 + fy($1)):(2):("") w labels point pt 7 ps myps lc rgbcolor "#AAAAAA", "" u ($5 == 0 ? $2+fx($1) : 1/0):($4 + fy($1)):(2):("") w labels point pt 6 ps myps lc rgbcolor "#888888",  "" u ($5 == 1 ? $2+fx($1) : 1/0):($4 + fy($1)):(2):("") w labels point pt 11 ps myps lc rgbcolor "#FFCCCC", "" u ($5 == 1 ? $2+fx($1) : 1/0):($4 + fy($1)):(2):("") w labels point pt 10 ps myps lc rgbcolor "#FFAAAA", "" u ($5 == 0 ? $2+fx($1)-myshift : 1/0):($4 + fy($1)-myshift):(2):($1) w labels font "Arial, 10" tc rgbcolor "#FFFFFF", "" u ($5 == 1 ? $2+fx($1)-myshift : 1/0):($4 + fy($1)-myshift):(2):($1) w labels font "Arial, 10" tc rgbcolor "#FF0000"
