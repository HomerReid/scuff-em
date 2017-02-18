#!/bin/bash 

##################################################
# capacitance per unit length from analog devices memo
##################################################
Eps0=8.85e-3 # pf per mm
InchesPerMM=(1.0/25.4)
EpsR=4.0
W=1.0;
T=0.0;
CPUL(H) = InchesPerMM*0.67*(EpsR+1.41) / log( 5.98*H / (0.8*W + T) )

##################################################
# scuff-static data files ########################
##################################################
MDATA='StriplineCapacitor_1868.CapVsT'
FDATA='StriplineCapacitor.CapVsT'
L=10.0; # length of trace in mm

MC(C11,C12,C22)=1.0 / ( 1.0/C11 + 1.0/C22 - 2.0/C12 )

set xlabel 'PCB dielectric thickness (mm)'
set ylabel 'Capacitance per unit length (pf/mm)'
plot CPUL(x) w l t 'Theory'                                                      \
    ,MDATA   u 1:(FF2*MC($2,$3,$4)*Eps0/L) t 'SCUFF (medium mesh)' w p pt 7 ps 1 \
    ,FDATA   u 1:(FF2*MC($2,$3,$4)*Eps0/L) t 'SCUFF (fine   mesh)' w p pt 6 ps 2
