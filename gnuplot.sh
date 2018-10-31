set xrange [1e-6:1e+6]
set yrange [1e+4:1e+8]
set terminal png
set out "Trho.png"
set logscale xy
set title "T v rho"
plot "profile.dat" u 2:3 w l 
unset xrange 
set title "T v radius"
set out "T.png"
plot "profile.dat" u 1:3 w l
set yrange [1e-6:1e+6]
set title "rho v radius"
plot "profile.dat" u 1:2 w l



