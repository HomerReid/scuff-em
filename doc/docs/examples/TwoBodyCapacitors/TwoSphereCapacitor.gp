COARSE='TwoSphereCapacitor_501.out'
FINE='TwoSphereCapacitor_2604.out'

set logscale xy

set key at 0.5, 100 spacing 5

set xlabel 'Surface-surface separation $d$'
set ylabel 'Capacitance $C/\epsilon_0$'

C11Short(d) = 0.981755 + 0.17057*d  - 0.006183*d*d - 0.25*(1.0+d/6.0 - d*d/180.0)*log(d);
C12Short(d) = 0.288608 + 0.096712*d - 0.003570*d*d - 0.25*(1.0+d/6.0 - d*d/180.0)*log(d);

C11Long(d)  = 4.0*pi*abs((1.0 - 1/d + 1.0/(d*d-1.0) - 1.0/(d*d*d-2.0*d)));
C12Long(d)  = (4.0*pi/d)*(1.0 + 1/(d*d-2.0) + 1/(d*d*d*d - 4.0*d*d + 3.0));

ShortC11(d) = d > 1.0 ? 1/0 : 4.0*pi*C11Short(d);
ShortC12(d) = d > 1.0 ? 1/0 : 4.0*pi*C12Short(d);
LongC11(d)  = d < 3.0 ? 1/0 : C11Long(d)+C12Long(d);
LongC12(d)  = d < 3.0 ? 1/0 : C12Long(d);

set yrange [0.02:]
set logscale xy
set xlabel 'Surface-surface separation $d$'
set ylabel 'Capacitance coefficient $C_{ij}/\epsilon_0$'
plot COARSE       u 1:(abs($2))  t '$C_{11}$ (SCUFF, $N=501$)'   w p pt 7 ps 1 lc RED \
    ,FINE         u 1:(abs($2))  t '$C_{11}$ (SCUFF, $N=2604$)'  w p pt 6 ps 2 lc RED \
    ,ShortC11(x)                 t '$C_{11}$ (small-d series)'   w l lw 2      lc RED \
    ,LongC11(x+2.0)              t '$C_{11}$ (large-d series)'   w l lt 0 lw 5 lc RED \
    ,COARSE       u 1:(abs($3))  t '$C_{12}$ (SCUFF, $N=501$)'   w p pt 7 ps 1 lc BLUE \
    ,FINE         u 1:(abs($3))  t '$C_{12}$ (SCUFF, $N=2604$)'  w p pt 6 ps 2 lc BLUE \
    ,ShortC12(x)                 t '$C_{12}$ (small-d series)'   w l lw 2      lc BLUE \
    ,LongC12(x+2.0)              t '$C_{12}$ (large-d series)'   w l lt 0 lw 5 lc BLUE

set key at 0.5,4 spacing 4
call 'latex' 'TwoSphereCapacitorData'
