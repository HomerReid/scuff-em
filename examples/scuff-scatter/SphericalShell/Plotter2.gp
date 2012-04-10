Analytical(x) = 36*x/(7.0*x*x+ 22.0*x+ 7.0);

set xlabel 'Epsilon / Epsilon_0'
set ylabel 'E_z(0,0,0) / E_0'
set title 'Z-component of E-field at origin vs. shell permittivity'
set key at 90, 0.7
plot  Analytical(x) t 'Analytical result: 36x/(7x^2 + 22x + 7)' w l lw 2.0 \
    ,'EzVsEps.dat' every 2 u 1:2 t 'scuff-scatter calculation' w p pt 7 ps 1.5
     

call 'png' 'EzVsEps.png'
