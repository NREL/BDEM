set terminal pngcairo dashed size 800,500
set output "compare.png"
set termoption font "Helvetica,18"
set key spacing 1.2
set lmargin 9
set rmargin 5
set bmargin 3.5
set xlabel "Time (sec)"
set ylabel "Position (m)"
plot 'p1_refsoln' u 1:($2*0.01) w lp lw 1 ps 2 pt 5 title "P1 Reference",\
     'p2_refsoln' u 1:($2*0.01) w lp lw 1 ps 2 pt 5 title "P2 Reference" axis x1y2,\
     'particle_positions.dat' u 1:2 w l lw 2 title "P1 computed",\
     'particle_positions.dat' u 1:4 w l lw 2 title "P2 computed" axis x1y2
