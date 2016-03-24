#######################################################################
# KrugerFormulas.jl -- a julia code to evaluate the T-matrix trace
#                   -- formulas of Krueger et al. for 
#                   -- (a) the radiation of a single SiO2 sphere and 
#                   -- (b) the radiative heat transfer and non-equilibrium 
#                   --     Casimir force between two SiO2 spheres
#
#                   -- Ref: Krueger et al, Phys. Rev. B 86 115423 (2012)
#                   -- (Equation numbers below refer to this paper)
#
# Julia code by Homer Reid, 6/2015
#######################################################################

##################################################
# \hbar * \omega_0 in units of joules (\omega_0 = 3e14 rad/sec)
##################################################
global const HBAROMEGA0 = 3.16371518e-20

##################################################
# \hbar * \omega_0^2 in units of watts
##################################################
global const HBAROMEGA02 = 9.491145534e-06

##################################################
# boltzmann constant k / (\hbar\omega_0)
# so that BOLTZMANNK * HBAROMEGA0 * (T in kelvin) = thermal energy in joules
##################################################
global const BOLTZMANNK = 4.36763e-4

##################################################
# inverse speed of light in units of nanoNewtons/watt
##################################################
global const ONE_OVER_C = (10.0/3.0)

###################################################
# SiO2 dielectric function 
# [(u is dimensionless angular frequency, u = omega / (3e14 rad/sec)]
###################################################
function EpsSIO2(u)

  w=3.0e14*u;
  A1              = 8.2736e+13;
  w01             = 8.54484e+13;
  G1              = 8.46448e+12;
  A2              = 1.58004e+14;
  w02             = 2.029e+14;
  G2              = 1.06449e+13;
  A3              = 3.39786e+13;
  w03             = 1.51198e+14;
  G3              = 8.33205e+12;
  EpsInf          = 2.03843;

  EpsInf + A1*A1/(w01*w01 - w*w - im*w*G1) + A2*A2/(w02*w02 - w*w - im*w*G2) + A3*A3/(w03*w03 - w*w - im*w*G3);
end


###################################################
function EpsGold(u)
    w=3.0e14*u;
    wp = 1.37e16; 
    gamma = 5.32e13;
    1 - wp^2 / (w * (w + im*gamma));
end
###################################################

###################################################
# Get the L=1 T-matrix elements for an SiO2 sphere.
# Returns (T_1^M, T_1^N).
# u is angular frequency in units of 3e14 rad/sec
# R is sphere radius in microns
###################################################
function GetTMatrixElements(u, Material, R=1.0)

  kR = u*R;

  kR3 = kR*kR*kR;
  kR5 = kR3*kR*kR;
  kR6 = kR5*kR;

  Eps = Material=="Gold" ? EpsGold(u) : EpsSiO2(u);
  Mu  = 1.0;
 
  Factor1 = 2.0*(Eps-1.0)/ ( 3.0*(Eps+2.0 ));
  Factor2 = (2.0 - 3.0*Eps + Eps*Eps*(1.0+Mu) ) / (5.0*(2.0+Eps)*(2.0+Eps));
  Factor3 = -4.0*(Eps-1.0)*(Eps-1.0) / (9.0 * (2.0+Eps)*(2.0+Eps));

  TN = im*Factor1*kR3 + 2.0*im*Factor2*kR5 + Factor3*kR6;

  Mu  = Eps;
  Eps = 1.0;

  Factor1 = 2.0*(Eps-1.0)/ ( 3.0*(Eps+2.0) );
  Factor2 = (2.0 - 3.0*Eps + Eps*Eps*(1.0+Mu) ) / (5.0*(2.0+Eps)*(2.0+Eps));
  Factor3 = -4.0*(Eps-1.0)*(Eps-1.0) / (9.0 * (2.0+Eps)*(2.0+Eps));

  TM = im*Factor1*kR3 + 2.0*im*Factor2*kR5 + Factor3*kR6;

  (TM, TN)

end 

###################################################
# Get the 'flux' quantities for radiated power,
# power transfer, and force
#
# The 'flux' is the quantity multiplied by Delta Theta(\omega,T) and 
# integrated over \omega to yield the full power or force (with
# no extra factor of 2 or pi or anything).
# 
# DeltaTheta = Theta(\omega, T1) - Theta(\omega, TEnv)
# 
# Theta = \hbar\omega / (e^{\hbar\omega/kT} - 1)
#
# inputs: 
#  Omega = angular frequency in units of 3e14 rad/sec
#      d = center-center separation (microns)
#      R = sphere radius (microns)
#
# returns ( Phi_{PRad}, Phi_{PTransfer}, Phi_F^{1->2}, Phi_F^{2->2} )
###################################################
function GetFluxes(Omega, d, Material="Gold", R=1)

  ##################################################
  ##################################################
  ##################################################
  kd  = Omega*d;
  kd2 = kd*kd;
  kd3 = kd2*kd;
  kd4 = kd3*kd;
  kd5 = kd4*kd;
  kd6 = kd5*kd;
  kd7 = kd6*kd;
  
  ##################################################
  ##################################################
  ##################################################
  (TM,TN) = GetTMatrixElements(Omega, Material);
  T=[TM; TN];

  ##################################################
  ##################################################
  ##################################################
  ReTpT2    = zeros(2);
  ReTpT2[1] = real(TM) + (abs(TM))^2;
  ReTpT2[2] = real(TN) + (abs(TN))^2;

  ##################################################
  # single-sphere radiated power flux, Kruger equation (131)
  ##################################################
  PRad = -2.0*3.0*( ReTpT2[1] + ReTpT2[2] ) / pi; 
  
  ##################################################
  # two-sphere power transfer, Kruger equation (137)
  ##################################################
  Ua=9.0/(2.0*kd2) + 9.0/(2.0*kd4);
  Ub=27.0/(2.0*kd6);
  PTransfer = 0.0;
  for p=1:2, pPrime=1:2
    U = ( (p==pPrime) ? Ua+Ub : Ua );
    PTransfer += 2.0 * ReTpT2[p] * ReTpT2[pPrime] * U / pi;
  end
  
  ##################################################
  # force on body 2 due to body 1, Krueger eq. (149)
  ##################################################
  F12=0.0; 
  for p=1:2, pP=1:2

    Delta = (p==pP) ? 1.0 : 0.0;

    TpBar   = T[3-p];
    TpTpBar = T[p]*conj(TpBar);

    U1 = (9.0/kd2) * ( real(T[pP]) + Delta*real(TpTpBar) );
    U2 = imag(T[pP])*( 9.0/kd3 + Delta*81.0/kd7 );
    U3 = (imag(T[pP]) - 0.5*Delta*imag(TpTpBar)) * 18.0/kd5;

    F12 += ReTpT2[p]*(U1 + U2 + U3);
  end
  F12 *= -1.0*ONE_OVER_C/pi;

  ##################################################
  # force on body 2 due to body 2, Krueger eq. (150)
  ##################################################
  F22=0.0; 
  ExpFac = exp(2.0*im*kd);
  for p=1:2

    Tp   = T[p];
    TpBar= T[3-p];

    U1 = (Tp-TpBar)*(9.0/kd2 + im*27.0/kd3);
    U2 = -1.0*(Tp-TpBar/2.0)*72.0/kd4;
    U3 = -1.0*(Tp-TpBar/8.0)*im*144.0/kd5;
    U4 =  Tp*(162.0/kd6 + im*81.0/kd7);

    F22 += ReTpT2[p]*real((U1 + U2 + U3 + U4)*ExpFac);

  end
  F22 *= 1.0*ONE_OVER_C/pi;

  ##################################################
  ##################################################
  ##################################################
  (PRad, PTransfer, F12, F22)
end

###################################################
###################################################
###################################################
function Theta(Omega, T)

  if (T==0)
    return 0.0;
  end
 
  if (Omega==0.0)
    return BOLTZMANNK*T;
  end

  return Omega / ( exp( Omega/(BOLTZMANNK*T) ) - 1.0 );

end

###################################################
###################################################
###################################################
function Integrand(x, v, d, T1, T2, TEnv, f)
  
  if ( (x==0.0) || (x==1.0) )
    v[:]=zeros(4);
    return;
  end

  Omega = x / (1.0-x);
  J     = 1.0 / ( (1.0-x)*(1.0-x) );
  DeltaTheta1 = Theta(Omega,T1); #// - Theta(Omega,TEnv);
  DeltaTheta2 = Theta(Omega,T2); #// - Theta(Omega,TEnv);

  (PR, PT, F12, F22)=GetFluxes(Omega, d);

  HBAROMEGA02 = 9.491145534e-06;
  v[:] = HBAROMEGA02*[ PR*DeltaTheta1
                       PT*(DeltaTheta1-DeltaTheta2)
                       F12*DeltaTheta1
                       F22*DeltaTheta2
                     ];

  if (f!=0)
    @printf(f,"%e %e %e %e %e %e \n",d,Omega,
                v[1],v[2],v[3],v[4]);
  end
   
  v*=J;

end

#################################################################
#################################################################
#################################################################
function WriteFilePreamble(f)
     @printf(f,"# 1: center--center separation\n");
     @printf(f,"# 2: omega \n");
     @printf(f,"# 3: radiated power integrand\n");
     @printf(f,"# 4: power transfer integrand\n");
     @printf(f,"# 5: force (1->2) integrand\n");
     @printf(f,"# 6: force (2->2) integrand\n");
end

###################################################
# if ByOmega==true then FixedParm is interpreted as 
# the center-center separation and the frequency is 
# scanned
# if ByOmega==false then FixedParm is interpreted as 
# Omega and the center-center separation is scanned
###################################################
function PlotFlux(FixedParm, Material="Gold", ByOmega=true, FileName=0)

  if (FileName==0)
   FileName=@sprintf("%sSpheres_Kruger.SIIntegrand",Material);
  end

  ScanVector=ByOmega ? logspace(-2,2,200) : linspace(5,100,96);

  f=open(FileName,"a");
  WriteFilePreamble(f);
  for ScanParm in ScanVector
    Omega = ByOmega ? ScanParm  : FixedParm;
    d     = ByOmega ? FixedParm : ScanParm;
    (PR, PT, F12, F22)=GetFluxes(Omega, d, Material)
    @printf(f,"%e %e %e %e %e %e \n",d,Omega,PR,PT,F12,F22);
  end
  close(f);
end 

###################################################
# Integrate over frequency for a given separation and
# temperatures
###################################################
using Cubature;
function GetPF(d, T1, T2, TEnv, WriteSIIntegrandFile=false)

 ##################################################
 # open .SIIntegrand file and write preamble if necessary
 ##################################################
 f=0;
 if (WriteSIIntegrandFile)
   FileName="SiO2Spheres_Kruger.SIIntegrand";
   WritePreamble = !isfile(FileName);
   f=open(FileName,"a");
   if (WritePreamble)
     @printf(f,"# (T1, T2, TEnv) = (%e, %e, %e)\n",T1,T2,TEnv);
     WriteFilePreamble(f);
   end
 end
 
 ##################################################
 # do the integration
 ##################################################
 (PFT,err) = pquadrature(4, (x,v)->Integrand(x,v,d,T1,T2,TEnv,f), 0.0, 1.0,
                         reltol=1.0e-4)
 
 WriteSIIntegrandFile ? close(f) : 0;

 ##################################################
 # write results to NEQPFT file
 ##################################################
 FileName="SiO2Spheres_Kruger.NEQPFT";
 WritePreamble = !isfile(FileName);
  
 f=open(FileName,"a");
 if (WritePreamble)
   @printf(f,"# (T1, T2, TEnv) = (%e, %e, %e)\n",T1,T2,TEnv);
   @printf(f,"# 1: center--center separation\n");
   @printf(f,"# 2: radiated power (watts)\n");
   @printf(f,"# 3: power transfer (1->2) (watts)\n");
   @printf(f,"# 4: force (1->2) (nanoNewtons)\n");
   @printf(f,"# 5: force (2->2) (nanoNewtons)\n");
 end
 @printf(f,"%e %e %e %e %e\n",d,PFT[1],PFT[2],PFT[3],PFT[4]);

 close(f);

end
