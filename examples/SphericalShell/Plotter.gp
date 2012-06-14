Eps=10.0;

Denom(Eps) = 7.0*Eps*Eps + 22.0*Eps + 7.0;

A=(7.0/2.0)*(Eps-1.0)*(1.0+2.0*Eps) / Denom(Eps);
B=-12.0*(1.0+2.0*Eps) / Denom(Eps);
C=-(3.0/2.0)*(Eps-1.0) / Denom(Eps);
D=-36.0*Eps / Denom(Eps);

EzOut(z) = 1 + 2*A/(z*z*z);
EzMid(z) = -B + 2*C/(z*z*z);
EzIn(z)  = -D;

Ez(z) = z<0.5 ? EzIn(z) : z<1.0 ? EzMid(z) : EzOut(z);

set xlabel 'z'
set ylabel 'E_z / E_0'
set title 'Z-component of E-field along the Z axis'
set key at 1.0, 2.0 
plot  Ez(x) t 'Analytical result' w l lw 1.5 , \
     'LineOfPoints.total' every 2 u 3:8 t 'scuff-scatter calculation' w p pt 7 ps 1.5

call 'png' 'EzVsZ.png'

