set key left
set ylabel "y/R"
set xlabel "x/R"
#set log xy
#set format y "%2.0t{/Symbol \264}10^{%L}"
set term postscript eps enhanced color 16
set output "Picture.eps"
#set arrow from 3.1,0 to 3.1,10 lt 2 nohead 
set parametric
fx(t)=3.1
fy(t)=t
set style line 1 lt 10 lw 1 pt 3
plot 'data.dat' pt 7 ps 0.1 lt rgb "#8A2BE2" notitle

