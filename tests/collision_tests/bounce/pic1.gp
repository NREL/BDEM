set terminal pngcairo dashed size 800,500
set output "compare.png"
set termoption font "Helvetica,18"
set key spacing 1.2
set lmargin 9
set rmargin 5
set bmargin 3.5
set xlabel "Time (sec)"
set ylabel "Position (m)"
plot 'ref_ysoln' u 1:($2*0.01) w lp lw 1 ps 2 pt 5 title "Reference pos soln",\
     'ref_velsoln' u 1:($2*0.01) w lp lw 1 ps 2 pt 5 title "Reference vel soln" axis x1y2,\
     'particle_data.dat' u 1:2 w l lw 2 title "Computed pos soln ",\
     'particle_data.dat' u 1:3 w l lw 2 title "Computed vel soln" axis x1y2
