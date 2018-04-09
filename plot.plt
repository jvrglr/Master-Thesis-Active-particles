set key left
set ylabel "y/R"
set xlabel "x/R"
#set log xy
#set format y "%2.0t{/Symbol \264}10^{%L}"
#set xrange [1:2]
#set yrange [0:1]
set term postscript eps enhanced color 16
set output "Picture.eps"
#set arrow from 3.1,0 to 3.1,10 lt 2 nohead
plot 'self_U.dat' u 1:2:($3/100):($4/100) w vectors notitle
