i={0,1};
nn={1,0.1};

Fa3(x,k) = 6 / ( ((-i*k*x)**4)*x );
Fb3(x,k) = exp(i*k*x) * (   -6/((x*k)**4)                    \
                          +6*i/((x*k)**3)                    \
                          +  3/((x*k)**2)    		     \
                          -  i/((x*k)**1) ) / x;
F3(x,k) = Fa3(x,k) + Fb3(x,k);

Fa2(x,k) = -2*i/ ( ((k*x)**3) * x );

Fb2(x,k) = exp(i*k*x) * (  2*i/((x*k)**3) 	             \
                            +2/((x*k)**2)                    \
                            -i/((x*k)**1) ) / x;

F2(x,k) = Fa2(x,k) + Fb2(x,k);

set terminal x11 1 
set title 'P=2'
plot [0.1:100] abs(Fa2(10,nn*x)), abs(Fb2(10,nn*x)), abs(F2(10,nn*x)), \
               abs(Fa2(1,nn*x)), abs(Fb2(1,nn*x)), abs(F2(1,nn*x))

set terminal x11 2 
set title 'P=3'
plot [0.1:100] abs(Fa3(10,nn*x)), abs(Fb3(10,nn*x)), abs(F3(10,nn*x)), \
               abs(Fa3(1,nn*x)), abs(Fb3(1,nn*x)), abs(F3(1,nn*x))
