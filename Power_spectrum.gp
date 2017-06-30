#!/usr/bin/gnuplot
fname   = 'RF_wht_window.csv'
set palette model RGB defined (0.0 1.0 1.0 1.0, 1/256. 1.0 1.0 1.0, 1/255. 0.0 0.0 0.51, 0.34 0.0 0.81 1.0, 0.61 0.87 1.0 0.12, 0.84 1.0 0.2 0.0, 1.0 0.51 0.0 0.0) positive

set output "RF_wht_window.eps"
set terminal postscript eps enhanced solid color "Helvetica"

set grid
set multiplot

set xlabel '[MHz]'
set xrange [0:200*10]
set ylabel '[dBm]'
set yrange [-150:30]
unset key

set label 1 'Thermal noise {/Symbol m}=0V, {/Symbol s}=0.05mV' at 195*10,-135 right front
set label 2 'Hanning window, FFT size=512, loops=8000' at 195*10,-145 right front 
set label 3 'peak {/Symbol m}=0.3185V, {/Symbol s}=1.36455e-06V' at 60*10,27 center front
plot fname u ($1*10):2:3 w p pt 5 lc palette
unset xlabel
unset ylabel
unset label 1
unset label 2
unset label 3
show label
unset grid

set origin 0.23,0.52
set size 0.2,0.4

set bmargin 0; set tmargin 0; set lmargin 0; set rmargin 0
clear
set xrange[129.6*10:129.77*10]
set yrange[20.0633:20.0636]
set xtics(129.6*10,129.77*10)
set ytics(20.0633,20.0636)
unset key
unset xlabel
unset ylabel
plot fname index 0 u ($1*10):2:3 w p pt 5 lc palette

unset multiplot
unset output
set term x11
