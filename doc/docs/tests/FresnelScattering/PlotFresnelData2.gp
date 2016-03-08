D2R=pi/180.0;

sind(x) = sin(x*D2R);
cosd(x) = cos(x*D2R);

Radical(Eps,Theta)   = sqrt( Eps - sind(Theta)*sind(Theta) );

DenomTE(Eps,Theta) = cosd(Theta) + Radical(Eps,Theta);
tTE(Eps,Theta)     = 2.0*cosd(Theta) / DenomTE(Eps,Theta);
rTE(Eps,Theta)     = (cosd(Theta)-Radical(Eps,Theta)) / DenomTE(Eps,Theta);

DenomTM(Eps,Theta) = Eps*cosd(Theta) + Radical(Eps,Theta);
tTM(Eps,Theta)     = 2.0*sqrt(Eps)*cosd(Theta) / DenomTM(Eps,Theta);
rTM(Eps,Theta)     = (Eps*cosd(Theta) - Radical(Eps,Theta)) / DenomTM(Eps,Theta);

FILE='E10HalfSpace_40.transmission'
Epsilon=10

set key spacing 5
set key at 120,1
set xlabel '$\theta$'
set ylabel 'Transmission/reflection coefficient magnitude'
plot FILE u 2:7          t '$t^{\text{\tiny{TE}}}$ ({\sc scuff})' w p pt 7 ps 1.5 \
    ,abs(tTE(Epsilon,x)) t '$t^{\text{\tiny{TE}}}$ (theory)' w l lw 2 \
    ,FILE u 2:9          t '$t^{\text{\tiny{TM}}}$ ({\sc scuff})' w p pt 7 ps 1.5 \
    ,abs(tTM(Epsilon,x)) t '$t^{\text{\tiny{TM}}}$ (theory)' w l lw 2 \
    ,FILE u 2:11         t '$r^{\text{\tiny{TE}}}$ ({\sc scuff})' w p pt 7 ps 1.5 \
    ,abs(rTE(Epsilon,x)) t '$r^{\text{\tiny{TE}}}$ (theory)' w l lw 2 \
    ,FILE u 2:13         t '$r^{\text{\tiny{TM}}}$ ({\sc scuff})' w p pt 7 ps 1.5 \
    ,abs(rTM(Epsilon,x)) t '$r^{\text{\tiny{TM}}}$ (theory)' w l lw 2 

call 'latex' 'Fresnel'
#!evince Fresnel.pdf
