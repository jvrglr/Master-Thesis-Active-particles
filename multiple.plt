set term png size 1600,1200
do for [t=1:500]{
set key left
set ylabel "y/R"
set xlabel "x/R"
set xrange [0:10]
set yrange [0:10]
set output "image.".t.".png"
plot "force.".t.".dat" u 1:2 pt 7 ps 1 lt rgb "#8A2BE2" notitle,\
"force.".t.".dat" u 1:2:($5/10.0):($6/10.0) w vectors head filled title "Self propulsion",\
"force.".t.".dat" u 1:2:($3/10.0):($4/10.0) w vectors head filled title "Repulsive force",\
}
