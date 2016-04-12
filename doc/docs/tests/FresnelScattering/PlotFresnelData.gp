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

plot FILE u 2:7          t 'tTE (SCUFF)'  w p pt 7 ps 1.5 \
    ,abs(tTE(Epsilon,x)) t 'tTE (theory)' w l lw 2        \
    ,FILE u 2:9          t 'tTM (SCUFF)'  w p pt 7 ps 1.5 \
    ,abs(tTM(Epsilon,x)) t 'tTM (theory)' w l lw 2        \
    ,FILE u 2:11         t 'rTE (SCUFF)'  w p pt 7 ps 1.5 \
    ,abs(rTE(Epsilon,x)) t 'rTE (theory)' w l lw 2        \
    ,FILE u 2:13         t 'rTM (SCUFF)'  w p pt 7 ps 1.5 \
    ,abs(rTM(Epsilon,x)) t 'rTM (theory)' w l lw 2
