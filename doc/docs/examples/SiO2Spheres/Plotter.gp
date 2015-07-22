REF='SiO2Spheres_Kruger.SIFlux'
DATA0='SiO2Sphere_501.SIFlux'
DATA1='SiO2Spheres_501.SIFlux'
DATA2='SiO2Spheres_1479.SIFlux'

ifeq(x,y,z)=(x==y ? z : 1/0)

set xlabel 'Omega (3e14 rad/sec)'

set logscale y
set format y "%.0e"
set format x "%.1f"

##################################################
##################################################
##################################################
unset logscale x
set yrange [1e-6:]
set xrange [0.01:1]
set terminal x11 1
set title 'Power radiation'
set ylabel 'Radiated Power flux (dimensionless)'
set key at 0.85,1e-3
plot REF   u 2:(abs($3)) t 'Krueger'      w l lw 2        \
    ,DATA0 u 2:4 t 'SCUFF, N=501' w p pt 7 ps 1.25
call 'png' 'SiO2Sphere_PowerRadiation.png'

##################################################
##################################################
##################################################
unset logscale x
unset yrange
set autoscale y
set xrange [0.1:1]
set terminal x11 2
set title 'Power transfer'
set ylabel 'Transferred Power flux (dimensionless)'
set key at 0.85,1e-8
plot REF   u 2:4                t 'Krueger'      w l lw 2        \
    ,DATA1 u (ifeq($3,12,$2)):4 t 'SCUFF, N=501' w p pt 7 ps 1.25\
    ,DATA2 u (ifeq($3,12,$2)):4 t 'SCUFF, N=1479' w p pt 6 ps 1.25
call 'png' 'SiO2Spheres_PowerTransfer.png'

##################################################
##################################################
##################################################
set terminal x11 3
set title 'Force flux (1-->2)'
set key at 0.85,1e-6
set ylabel 'Force density (nanonewtons/watts)'
plot REF   u 2:(abs($5))                t 'Krueger'      w l lw 2         \
    ,DATA1 u (ifeq($3,12,$2)):(abs($5)) t 'SCUFF, N=501' w p pt 7 ps 1.00 \
    ,DATA2 u (ifeq($3,12,$2)):(abs($5)) t 'SCUFF, N=1479' w p pt 6 ps 1.25 
call 'png' 'SiO2Spheres_F12.png'

##################################################
##################################################
##################################################
set terminal x11 4
set key at 0.85,1e-6
set title 'Force flux 2-->2'
set ylabel 'Force density (nanonewtons/watts)'
plot REF   u 2:(+$6)                    t 'Krueger (positive)' w l lw 2   \
    ,REF   u 2:(-$6)                    t 'Krueger (negative)' w l lw 2   \
    ,DATA1 u (ifeq($3,22,$2)):(abs($5)) t 'SCUFF, N=501' w p pt 7 ps 1.00 \
    ,DATA2 u (ifeq($3,22,$2)):(abs($5)) t 'SCUFF, N=1479' w p pt 6 ps 1.25
call 'png' 'SiO2Spheres_F22.png'
