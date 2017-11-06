##################################################
# L=1 T-matrix elements for homogeneous sphere
# Eps, Mu = relative permittivity, permeability
# a0 = k0*R
# note \overline{f(x)} == f(x) + x*df/dx
##################################################
I={0.0,1.0}

# spherical bessel function j_1(x) and \overline{j_1}(x)
j1(x)    = ( abs(x) < 1.0e-2 ? x*(1.0-x*x/10.0)/3.0 : (sin(x)-x*cos(x))/(x*x) )
j1Bar(x) = ( abs(x) < 1.0e-2 ? (2.0/3.0)*x*(1.0-x*x/5.0) : sin(x)-j1(x) )

# spherical Hankel function h_1(x)
h1(x)    = ( abs(x) < 1.0e-2 ? (-I/(x*x) - I/2.0) : exp(I*x)*(-I-x)/(x*x) )
h1Bar(x) = ( abs(x) < 1.0e-2 ? (I-0.5*I*x*x+2*x*x*x/3.0)/(x*x) : exp(I*x)*(I+x-I*x*x)/(x*x) )

T1M(Eps,Mu,a0,a1) = -1.0*( Mu*j1(a1)*j1Bar(a0) - j1(a0)*j1Bar(a1) )         \
                       / ( Mu*j1(a1)*h1Bar(a0) - h1(a0)*j1Bar(a1) );

T1N(Eps,Mu,a0,a1) = -1.0*( Eps*j1(a1)*j1Bar(a0) - j1(a0)*j1Bar(a1) )         \
                       / ( Eps*j1(a1)*h1Bar(a0) - h1(a0)*j1Bar(a1) );

TM(Eps,Mu,a0) = T1M( Eps, Mu, a0, sqrt(Eps*Mu)*a0 );
TN(Eps,Mu,a0) = T1N( Eps, Mu, a0, sqrt(Eps*Mu)*a0 );

##################################################
##################################################
##################################################
E10COARSE='Run1/E10Sphere_327.TMatrix'
E10FINE='Run1/E10Sphere_1362.TMatrix'
E10M5COARSE='Run1/E10M5Sphere_327.TMatrix'
E10M5FINE='Run1/E10M5Sphere_1362.TMatrix'

Alpha1=0  # I will plot the diagonal T-matrix element for the
Alpha2=1  # spherical waves with indices (L,M,P)=(1,-1,0)
          #                          and (L,M,P)=(1,-1,1)

set logscale xy
set terminal x11 1 
set yrange [:2.0]

set xlabel '$\omega R/c$'
set ylabel '$|T|$'

##################################################
##################################################
##################################################
set key at 4,0.001 spacing 6
set title '$T_1^M$ for $\epsilon_r=10$ sphere'
plot E10COARSE   u (ifeq($2,Alpha1,ifeq($6,Alpha1,$1))):(D2($10,$11)) \
                 t '$T_1^M$ {\sc scuff}, $N$=327'       \
                 w p pt 7 ps 1.5                                      \
    ,E10FINE     u (ifeq($2,Alpha1,ifeq($6,Alpha1,$1))):(D2($10,$11)) \
                 t '$T_1^M$ {\sc scuff}, $N$=1362'      \
                 w p pt 6 ps 1.5                                      \
   ,abs(TM(10,1,x)) t '$T_1^M$ (theory)' w l lw 2  

#call 'latex' 'E10_T1M'

set terminal x11 2
set title '$T_1^N$ for $\epsilon_r=10$ sphere'

set key at 4,0.01 spacing 6
plot E10COARSE   u (ifeq($2,Alpha2,ifeq($6,Alpha2,$1))):(D2($10,$11))  \
                  t '$T_1^N$ {\sc scuff}, $N$=327'      \
                  w p pt 7 ps 1.5                                      \
    ,E10FINE     u (ifeq($2,Alpha2,ifeq($6,Alpha2,$1))):(D2($10,$11))  \
                  t '$T_1^N$ {\sc scuff}, $N$=1362'      \
                  w p pt 6 ps 2                                        \
    ,abs(TN(10,1,x)) t '$T_1^N$ (theory)' w l lw 2  
#call 'latex' 'E10_T1N'

##################################################
##################################################
##################################################
set key at 4,0.01 spacing 6
set title '$T_1^M$ for $\{\epsilon_r,\mu_r\}=\{10,5\}$ sphere'
plot E10M5COARSE   u (ifeq($2,Alpha1,ifeq($6,Alpha1,$1))):(D2($10,$11)) \
                   t '$T_1^M$ {\sc scuff}, $N$=327'       \
                   w p pt 7 ps 1.5                                      \
    ,E10M5FINE     u (ifeq($2,Alpha1,ifeq($6,Alpha1,$1))):(D2($10,$11)) \
                   t '$T_1^M$ {\sc scuff}, $N$=1362'      \
                   w p pt 6 ps 1.5                                      \
   ,abs(TM(10,5,x)) t '$T_1^M$ (theory)' w l lw 2  

call 'latex' 'E10M5_T1M'

set terminal x11 2
set title '$T_1^N$ for $\{\epsilon_r,\mu_r\}=\{10,5\}$ sphere'

set key at 4,0.01 spacing 6
plot E10M5COARSE   u (ifeq($2,Alpha2,ifeq($6,Alpha2,$1))):(D2($10,$11))  \
                    t '$T_1^N$ {\sc scuff}, $N$=327'      \
                    w p pt 7 ps 1.5                                      \
    ,E10M5FINE     u (ifeq($2,Alpha2,ifeq($6,Alpha2,$1))):(D2($10,$11))  \
                    t '$T_1^N$ {\sc scuff}, $N$=1362'      \
                    w p pt 6 ps 2                                        \
    ,abs(TN(10,5,x)) t '$T_1^N$ (theory)' w l lw 2  
call 'latex' 'E10M5_T1N'
