#!/bin/bash 

##################################################
##################################################
##################################################
W0P5N4L10='./Run5/StriplineCapacitorFull.W0.5.N4.CapVsT'
W0P5N2L10='./Run5/StriplineCapacitorFull.W0.5.N2.CapVsT'
W1P0N4L10='./Run6/StriplineCapacitorFull.W1.0.N4.CapVsT'
W1P0N2L10='./Run6/StriplineCapacitorFull.W1.0.N2.CapVsT'
W1P0N4L5='./Run3/StriplineCapacitor.W1.0.L5.N4.CapVsT'
W1P0N2L5='./Run3/StriplineCapacitor.W1.0.L5.N2.CapVsT'
W1P0N2L20='./Run3/StriplineCapacitor.W1.0.L20.N2.CapVsT'
W1P0N4L20='./Run3/StriplineCapacitor.W1.0.L20.N4.CapVsT'
W0P5N4L5B20='./Run4/StriplineCapacitor.W0.5.L5.B20.N4.CapVsT'
W0P5N2L5B20='./Run4/StriplineCapacitor.W0.5.L5.B20.N2.CapVsT'
W1P0N4L5B20='./Run4/StriplineCapacitor.W1.0.L5.B20.N4.CapVsT'
W1P0N2L5B20='./Run4/StriplineCapacitor.W1.0.L5.B20.N2.CapVsT'

##################################################
# capacitance per unit length from analog devices memo
##################################################
Eps0=8.85e-3 # pf per mm
InchesPerMM=(1.0/25.4)
EpsR=4.0
W=1.0;
T=0.0;
CPUL(H,W) = InchesPerMM*0.67*(EpsR+1.41) / log( 5.98*H / (0.8*W + T) )
CPULV(H,W) = InchesPerMM*0.67*(1.0+1.41) / log( 5.98*H / (0.8*W + T) )

##################################################
# scuff-static data files ########################
##################################################

MDATA='StriplineCapacitor_1868.CapVsT'
FDATA='StriplineCapacitor.CapVsT'

MC(C11,C12,C22)=1.0 / ( 1.0/C11 + 1.0/C22 - 2.0/C12 )

set xlabel 'PCB dielectric thickness (mm)'
set ylabel 'Capacitance per unit length (pf/mm)'
#plot CPUL(x,1.0)  w l                            t 'W=1, L=20 theory'   \
#    ,W1P0N4L5     u 1:(FF*MC($2,$3,$4)*Eps0/5)   t 'W=1, L=5  SCUFF N4' w p pt 7 ps 1 \
#    ,W1P0N4L5B20  u 1:(FF*MC($2,$3,$4)*Eps0/5)   t 'W=1, L=5  SCUFF N4 B 20' w p pt 4 ps 1 \
#    ,W1P0N4L10    u 1:(FF*MC($2,$3,$4)*Eps0/10)  t 'W=1, L=5  SCUFF N4' w p pt 6 ps 2 \
#    ,W1P0N4L20    u 1:(FF*MC($2,$3,$4)*Eps0/20)  t 'W=1, L=20 SCUFF N4' w p pt 4 ps 2 \
#    ,CPUL(x,0.5)  w l                            t 'W=0.5, L=20 theory'   \
#    ,W0P5N4L5B20  u 1:(FF*MC($2,$3,$4)*Eps0/5)   t 'W=0.5, L=5  SCUFF N4 B 20' w p pt 4 ps 1 \
#    ,W0P5N4L10    u 1:(FF*MC($2,$3,$4)*Eps0/10)  t 'W=0.5, L=5  SCUFF N4' w p pt 6 ps 2 \

plot CPUL(x,1.0)  w l                            t 'W=1, L=10 theory'   \
    ,W1P0N2L10    u 1:(FF*MC($2,$3,$4)*Eps0/10)  t 'W=1, L=10 SCUFF N2' w p pt 7 ps 1 \
    ,W1P0N4L10    u 1:(FF*MC($2,$3,$4)*Eps0/10)  t 'W=1, L=510  SCUFF N4' w p pt 6 ps 2 \
    ,CPUL(x,0.5)  w l                            t 'W=0.5, L=20 theory'   \
    ,W0P5N2L10    u 1:(FF*MC($2,$3,$4)*Eps0/10)  t 'W=0.5, L=10  SCUFF N2' w p pt 7 ps 1 \
    ,W0P5N4L10    u 1:(FF*MC($2,$3,$4)*Eps0/10)  t 'W=0.5, L=10  SCUFF N4' w p pt 6 ps 2
