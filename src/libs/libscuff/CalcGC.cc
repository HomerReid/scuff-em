/*
 * CalcGC.cc     -- routines for calculating the cartesian components
 *               -- of the homogeneous dyadic electromagnetic green's 
 *               -- functions, and their derivatives, at real or 
 *               -- imaginary frequency. 
 *
 * homer reid   -- 2/2008
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <libhrutil.h>

#define II cdouble(0,1)

/***************************************************************/
/* three-dimensional kronecker delta and levi-civita symbol    */
/***************************************************************/
static double Delta[3][3]=
 {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

static double LeviCivita[3][3][3]=
{ { { 0.0,  0.0,  0.0 }, {  0.0, 0.0, +1.0 }, {  0.0, -1.0, +0.0 }  },
  { { 0.0,  0.0, -1.0 }, {  0.0, 0.0,  0.0 }, { +1.0, +0.0, +0.0 }  },
  { { 0.0, +1.0,  0.0 }, { -1.0, 0.0,  0.0 }, {  0.0,  0.0,  0.0 }  }
};

/***************************************************************/
/* calculate tensor (dyadic) green's functions and their       */
/* derivatives.                                                */
/*                                                             */
/* inputs:                                                     */
/*                                                             */
/*  R[0..2]    ... cartesian components of evaluation point    */
/*  Frequency  ... either \omega or \xi                        */
/*  RealFreq   ... 1 for real frequency, 0 for imag frequency  */
/*  Eps        ... relative permittivity of medium             */
/*  Mu         ... relative permeability of medium             */
/*                                                             */
/* outputs:                                                    */
/*                                                             */
/*  GMuNu[Mu][Nu]         ... G_{\mu\nu}                       */
/*  CMuNu[Mu][Nu]         ... C_{\mu\nu}                       */
/*  GMuNuRho[Mu][Nu][Rho] ... d/dR_\rho G_{\mu\nu}             */
/*  CMuNuRho[Mu][Nu][Rho] ... d/dR_\rho C_{\mu\nu}             */
/*                                                             */
/*                                                             */
/* note:                                                       */
/*                                                             */
/*  G_\mu\nu(R) = f1(r) * [   f2(r)*Delta_\mu\nu               */
/*                          + f3*R_\mu*R\nu / r^2              */
/*                        ]                                    */
/*                                                             */
/*  C_\mu\nu(R) = \LeviCivita_{\mu\nu\sigma}                   */
/*                 * R_\sigma * f4(r) * f5(r)                  */
/***************************************************************/
void CalcGC(double R[3], double Frequency, int RealFreq,
            double EpsR, double MuR,
            cdouble GMuNu[3][3], cdouble CMuNu[3][3],
            cdouble GMuNuRho[3][3][3], cdouble CMuNuRho[3][3][3])
{ 
  int Mu, Nu, Rho, Sigma;
  double k;
  double r, r2;
  cdouble ik, ikr, ikr2, ikr3, ExpFac;
  cdouble f1, f2, f3, f4, f5, f1p, f2p, f3p, f4p, f5p;

  r2=R[0]*R[0] + R[1]*R[1] + R[2]*R[2];

  if (r2==0.0)
   return;

  r=sqrt(r2);

  k=sqrt(EpsR*MuR)*Frequency;
  if(RealFreq)
   ik=II*k;
  else
   ik=-1.0*k;
  ExpFac=exp(ik*r);
  ikr=ik*r;
  ikr2=ikr*ikr;
  ikr3=ikr2*ikr;

  /* scalar factors */
  f1=MuR*ExpFac / (4.0*M_PI*ikr2*r);
  f1p=MuR*(ikr-3.0) * ExpFac / (4.0*M_PI*ikr2*r*r);

  f2=1.0 - ikr + ikr2;
  f2p=-1.0*ik + 2.0*ik*ik*r;

  f3=-3.0 + 3.0*ikr - ikr2;
  f3p=3.0*ik - 2.0*ik*ik*r;

  f4=(-1.0+ikr)/(ik);
  f4p=1.0;

  f5=MuR*ExpFac / (4.0*M_PI*r*r*r);
  f5p=MuR*(ikr-3.0)*ExpFac/(4.0*M_PI*r*r*r*r);

  /* computation of G_{\mu\nu} */
  if (GMuNu)
   { for(Mu=0; Mu<3; Mu++)
      for(Nu=0; Nu<3; Nu++)
       GMuNu[Mu][Nu]=f1*(f2*Delta[Mu][Nu] + f3*R[Mu]*R[Nu]/r2);
   };

  /* computation of C_{\mu\nu} */
  if (CMuNu)
   { for(Mu=0; Mu<3; Mu++)
      for(Nu=0; Nu<3; Nu++)
       for( CMuNu[Mu][Nu]=0.0, Rho=0; Rho<3; Rho++ )
        CMuNu[Mu][Nu]+=LeviCivita[Mu][Nu][Rho]*R[Rho]*f4*f5;
   };
   
  /* computation of G_{\mu\nu\rho} = d/dR_\rho G_\mu\nu  */
  if (GMuNuRho)
   { 
     for(Mu=0; Mu<3; Mu++)
      for(Nu=0; Nu<3; Nu++)
       for(Rho=0; Rho<3; Rho++)
        GMuNuRho[Mu][Nu][Rho]
         = (  (f1p*f2 + f1*f2p)*Delta[Mu][Nu] 
             +(f1p*f3 + f1*f3p - 2.0*f1*f3/r)*R[Mu]*R[Nu]/r2
           ) * R[Rho] / r
           + f1*f3*(Delta[Mu][Rho]*R[Nu] + R[Mu]*Delta[Nu][Rho])/r2;
   };

  /* computation of C_{\mu\nu\rho} = d/dR_\rho C_\mu\nu  */
  if (CMuNuRho)
   { 
     for(Mu=0; Mu<3; Mu++)
      for(Nu=0; Nu<3; Nu++)
       for(Rho=0; Rho<3; Rho++)
        for(CMuNuRho[Mu][Nu][Rho]=0.0, Sigma=0; Sigma<3; Sigma++)
         CMuNuRho[Mu][Nu][Rho]
          += LeviCivita[Mu][Nu][Sigma] *(   Delta[Rho][Sigma]*f4*f5
                                          + R[Sigma]*R[Rho]*(f4p*f5+f4*f5p)/r
                                        );
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void CalcGC(double R1[3], double R2[3], double Frequency, int RealFreq,
            double EpsR, double MuR, cdouble GMuNu[3][3], cdouble CMuNu[3][3],
            cdouble GMuNuRho[3][3][3], cdouble CMuNuRho[3][3][3])
{ 
  double R[3];

  R[0]=R1[0]-R2[0];
  R[1]=R1[1]-R2[1];
  R[2]=R1[2]-R2[2];

  CalcGC(R,Frequency,RealFreq,EpsR,MuR,GMuNu,CMuNu,GMuNuRho,CMuNuRho);
}
