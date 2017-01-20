
set  fit errorvariables
set encoding iso_8859_1
set key right bottom right outside



set terminal jpeg
set term jpeg size 700, 600
set output "ToandL.jpg"
set xlabel 'Number of sites'
set ylabel 'Auto Correlation Time'
 




f(x)= x**a 
fit f(x) 'OUTPUT' using 1:2  via a

set label sprintf("z=%.3g \261 %.3g", a,a_err) at screen .20,.80 tc rgb "black"


plot  [0:80] [0:1000] 'OUTPUT' using 1:2:3 lt rgb "violet" lw 6 title 'data'  with yerrorbars,  f(x)  lt rgb "black" lw 4 title 'L^z'




        


 exit
