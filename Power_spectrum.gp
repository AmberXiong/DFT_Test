#!/usr/bin/gnuplot
fname='RF_power_spectrum.csv'
set palette model RGB defined (0.0 1.0 1.0 1.0, 1/256. 1.0 1.0 1.0, 1/255. 0.0 0.0 0.51, 0.34 0.0 0.81 1.0, 0.61 0.87 1.0 0.12, 0.84 1.0 0.2 0.0, 1.0 0.51 0.0 0.0) positive

set output "Power_spectrum.eps"
set terminal postscript eps enhanced solid color "Helvetica"

set xlabel '[Hz]'
set xrange [0:1500000]
set ylabel '[dBm]'
#set yrange [-200:50]

#set key bottom right

#set label 1 "Response to 100mV tail-pulse" at first 0.52,1.0 front
plot 'RF_power_spectrum.csv' u 1:2:3 w p pt 5 lc palette

unset output
#set term x11
