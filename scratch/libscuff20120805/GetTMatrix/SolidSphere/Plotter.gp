set title 'T-matrix elements for a sphere of epsilon=10'
set xlabel 'Omega (units of 2pic/R)'
set ylabel 'magnitude of T-matrix element'
set xrange [0.1:5]
set key 1,1e-5

plot  'CONST_EPS_10.SphereTMatrix' u 1:(D2($4,$5)) t 'exact' w l \
     ,'E10_654.E11.E11'            u 1:(D2($2,$3)) t 'scuff' w p pt 7 ps 1.5 \
     ,'CONST_EPS_10.SphereTMatrix' u 1:(D2($6,$7)) t 'exact' w l \
     ,'E10_654.M20.M20'            u 1:(D2($2,$3)) t 'scuff' w p pt 7 ps 1.5 \
