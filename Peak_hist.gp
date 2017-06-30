#!/usr/bin/gnuplot
fname1 = 'Average_amp_dist_y_wht_w.csv'
fname2 = 'Average_amp_dist_y_qnt_w.csv'
fname3 = 'Average_amp_dist_y_phs_w.csv'
fname4 = 'Average_amp_dist_y_all_w.csv'

set output "wht_w_dist.eps"
set terminal postscript eps enhanced solid color "Helvetica"

set bmargin 0
set tmargin 0
set lmargin 1
set rmargin 1
set ylabel 'Counts'

set multiplot layout 4,1 title "Peak amplitude distribution\n" font ",12"

unset xlabel
set format x ""
unset ylabel
set yrange[0:16]
set key top right

set label 1 'average number:160, repeat times:50' at first 0.318546,13 right front
set label 2 '{/Symbol m}=0.318546V, {/Symbol s}=9.042e-08V' at first 0.318546,11 right front
set style fill solid border -1

plot fname1 u 1:2 with boxes title 'with thermal noise' linecolor rgb "#ff82ab"
#set boxwidth 0.7
unset label 1
unset label 2

unset xlabel
set format x ""
set yrange[0:22]
set key top right
set label 3 '{/Symbol m}=0.318539V, {/Symbol s}=2.9893e-06V' at first 0.318539,19 right front

plot fname2 u 1:2 with boxes title 'with quantization noise' linecolor rgb "#ffc125"
unset label 3

unset xlabel
set format x ""
set yrange[0:16]
set label 4 '{/Symbol m}=0.318546V, {/Symbol s}=1.0847e-08V' at first 0.31854603,13 center front
set key top right

plot fname3 u 1:2 with boxes title 'with phase noise' linecolor rgb "#b4eeb4"
unset label 4

set xlabel 'V'
set yrange[0:14]
set label 5 '{/Symbol m}=0.318537V, {/Symbol s}=3.54137e-06V' at first 0.318537,11 right front
set key top right

plot fname4 u 1:2 with boxes title 'with all noise' linecolor rgb "#9f79ee"
unset label 5

unset multiplot
unset output
