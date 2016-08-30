COARSE='TwoPlateCapacitor_280.out'
FINE='TwoPlateCapacitor_1160.out'

set logscale xy

set key at 0.5, 100 spacing 5

set xlabel 'Plate-plate separation $d$'
set ylabel 'Capacitance $C/\epsilon_0$'

plot COARSE u 1:(abs($3)) t 'SCUFF, $N=280$'   w p pt 7 ps 1 \
    ,FINE u 1:(abs($3))   t 'SCUFF, $N=1160$'  w p pt 6 ps 2 \
    ,COARSE u 1:(abs($2)) t 'SCUFF, $N=280$'   w p pt 7 ps 1 \
    ,FINE u 1:(abs($2))   t 'SCUFF, $N=1160$'  w p pt 6 ps 2 \
    ,1/x                  t '$\frac{A^2}{d}$'  w l lw 2

call 'latex' 'TwoPlateCapacitorData'
