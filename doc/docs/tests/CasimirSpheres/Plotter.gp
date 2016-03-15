###################################################
# theoretical predictions
###################################################
CP7=143.0/(16.0*pi);
CP9=7947.0/(160.0*pi);
CP10=2065.0/(32.0*pi);
EPEC(x)=CP7/(x**7) + CP9/(x**9) + CP10/(x**10);
FPEC(x)=7.0*CP7/(x**8) + 9.0*CP9/(x**10) + 10.0*CP10/(x**11);

Eps=10;
CE7=(23.0/4.0) * ( ((Eps-1.0)/(Eps+1.0))**2 ) / pi;
EE10(x)=CE7/(x**7);
FE10(x)=7.0*CE7/(x**8);

###################################################
# data files ######################################
###################################################
PECDATA='PECSpheres_327.out'
E10DATA='E10Spheres_327.out'

set format '%g'

set xlabel  'Center-center separation d'

set ylabel 'Force (\hbar c/R)'
set xrange [3:]
set xtics (3,5,10,30,50,100)

set key at 50,0.01 spacing 5

set terminal x11 1
set title   'Casimir force between R=1 micron spheres'
plot PECDATA u 1:(abs($2)) t 'PEC (SCUFF)'  w p pt 7 ps 1.5 \
    ,FPEC(x) lw 2          t 'PEC (theory)' \
    ,E10DATA u 1:(abs($2)) t 'Epsilon=10 (SCUFF)'  w p pt 7 ps 1.5 \
    ,FE10(x) lw 2          t 'Epsilon=10 (theory)'
