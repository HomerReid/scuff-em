Eps=10.0;
Mu=1.0;
n = sqrt(Eps*Mu);
g =sqrt(Eps/Mu); # gamma
i = {0,1};
r(x) = (g*g-1)*sin(n*x) / ( (1+g*g)*sin(n*x) + 2*i*g*cos(n*x))
t(x) = 2*i*g*exp(-i*x) / ( (1+g*g)*sin(n*x) + 2*i*g*cos(n*x))

plot 'E10Slab_104.transmission.1' u 1:3 w p pt 7 ps 1.5  t 'rPerp' \
   ,                           '' u 1:5 w p pt 7 ps 1.5  t 'rPar'  \
   ,                    abs(r(x))       w l t 'r, theory'          \
   , 'E10Slab_104.transmission.1' u 1:7 w p pt 7 ps 1.5  t 'tPerp' \
   ,                           '' u 1:9 w p pt 7 ps 1.5  t 'tPar' \
   ,                    abs(t(x))       w l t 'r, theory' 
