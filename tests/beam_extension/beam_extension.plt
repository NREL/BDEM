#/bin/gnuplot -p
set terminal epslatex color font 10 dashed
 
set style line 1 pt 6 lw 3 lt 1 lc rgb '#e41a1c' # red
set style line 2 pt 6 lw 3 lt 1 lc rgb '#377eb8' # blue
set style line 3 pt 6 lw 3 lt 1 lc rgb '#4daf4a' # green
set style line 4 pt 6 lw 3 lt 1 lc rgb '#984ea3'   # purple
set style line 5 pt 6 lw 3 lt 1 lc rgb '#ff7f00'   # orange
set style line 6 pt 6 lw 3 lt 1 lc rgb 'black'
set style line 7 pt 6 lw 3 lt 1 lc rgb '#377eb8' # blue 

set style line 11 pt 6 lw 3 lt 2 lc rgb '#e41a1c' # red
set style line 12 pt 6 lw 3 lt 2 lc rgb '#377eb8' # blue
set style line 13 pt 6 lw 3 lt 2 lc rgb '#4daf4a' # green
set style line 14 pt 6 lw 3 lt 2 lc rgb '#984ea3'   # purple
set style line 15 pt 6 lw 3 lt 2 lc rgb '#ff7f00'   # orange
set style line 16 pt 6 lw 3 lt 2 lc rgb 'black'
set style line 17 pt 6 lw 3 lt 2 lc rgb '#377eb8' # blue 

set style line 21 pt 6 lw 3 lt 3 lc rgb '#e41a1c' # red
set style line 22 pt 6 lw 3 lt 3 lc rgb '#377eb8' # blue
set style line 23 pt 6 lw 3 lt 3 lc rgb '#4daf4a' # green
set style line 24 pt 6 lw 3 lt 3 lc rgb '#984ea3'   # purple
set style line 25 pt 6 lw 3 lt 3 lc rgb '#ff7f00'   # orange
set style line 26 pt 6 lw 3 lt 3 lc rgb 'black'
set style line 27 pt 6 lw 3 lt 3 lc rgb '#377eb8' # blue 

# set style line 3 lt 1 lw 1.5 pt 186  ps 1.5*1.25  lc rgb '#66c2a5'   # aqua
# set style line 6 lt 1 lw 1.5 pt 16   ps 1.5*1.5   lc rgb 'black'     # black
# set style line 7 lt 1 lw 1.5 pt 17   ps 1.5*1.5   lc rgb '#377eb8'     # blue
# set style line 8 lt 1 lw 1.5 pt 19   ps 1.5*1.5   lc rgb '#e41a1c'     # red
#
set pointsize 1.5
 
set size square
set output 'beam_extension.tex'
  
set xrange[0:1]
set yrange[0:1]
   
set xlabel '$x/L$'
set ylabel '$\delta_n / \delta_{n,\text{max}}$' offset 8,12.0 rotate by 0


#set logscale y
#set logscale y2

#set format y "10^{%L}"
#set format x "10^{%L}"

#set format y2 "10^{%L}"

set xtics auto
set ytics auto
    
set key top center
 
f1='particle_data.dat'
f2='ref_soln_displacement.dat'

s=2

p f1 u ($1):($2) with points ls 1 title 'BDEM',\
  f2 u ($1):($2) w l ls 12 title 'ref.';
