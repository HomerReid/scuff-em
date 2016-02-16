/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This file is part of SCUFF-EM.
 *
 * SCUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * SCUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

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

namespace scuff {

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
/*  Omega      ... angular frequency                           */
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
void CalcGC(double R[3], cdouble Omega,
            cdouble EpsR, cdouble MuR,
            cdouble GMuNu[3][3], cdouble CMuNu[3][3],
            cdouble GMuNuRho[3][3][3], cdouble CMuNuRho[3][3][3])
{ 
  cdouble k=csqrt2(EpsR*MuR)*Omega;
  cdouble ik=II*k;

  double r2=R[0]*R[0] + R[1]*R[1] + R[2]*R[2];

  if (r2==0.0)
   { 
     cdouble GDiag  = II*real(k)/(6.0*M_PI);
     cdouble dCDiag = -real(k*k)/(12.0*M_PI);

     for(int Mu=0; Mu<3; Mu++)
      for(int Nu=0; Nu<3; Nu++)
       { GMuNu[Mu][Nu] = CMuNu[Mu][Nu] = 0.0;
         for(int Rho=0; Rho<3; Rho++)
          GMuNuRho[Mu][Nu][Rho] = CMuNuRho[Mu][Nu][Rho] = 0.0;
       };

     GMuNu[0][0] = GMuNu[1][1] = GMuNu[2][2] = GDiag;
     CMuNuRho[0][1][2] = CMuNuRho[1][2][0] = CMuNuRho[2][0][1] = dCDiag;
     CMuNuRho[0][2][1] = CMuNuRho[1][0][2] = CMuNuRho[2][1][0] = -dCDiag;
     return;
   };

  double r=sqrt(r2);

  cdouble ExpFac=exp(ik*r);
  cdouble ikr=ik*r;
  cdouble ikr2=ikr*ikr;
  cdouble ikr3=ikr2*ikr;

  /* scalar factors */
  cdouble f1=ExpFac / (4.0*M_PI*ikr2*r);
  cdouble f1p=(ikr-3.0) * ExpFac / (4.0*M_PI*ikr2*r*r);

  cdouble f2=1.0 - ikr + ikr2;
  cdouble f2p=-1.0*ik + 2.0*ik*ik*r;

  cdouble f3=-3.0 + 3.0*ikr - ikr2;
  cdouble f3p=3.0*ik - 2.0*ik*ik*r;

  cdouble f4=(-1.0+ikr)/(ik);
  cdouble f4p=1.0;

  cdouble f5=ExpFac / (4.0*M_PI*r*r*r);
  cdouble f5p=(ikr-3.0)*ExpFac/(4.0*M_PI*r*r*r*r);

  /* computation of G_{\mu\nu} */
  if (GMuNu)
   { for(int Mu=0; Mu<3; Mu++)
      for(int Nu=0; Nu<3; Nu++)
       GMuNu[Mu][Nu]=f1*(f2*Delta[Mu][Nu] + f3*R[Mu]*R[Nu]/r2);
   };

  /* computation of C_{\mu\nu} */
  if (CMuNu)
   { int Rho;
     for(int Mu=0; Mu<3; Mu++)
      for(int Nu=0; Nu<3; Nu++)
       for( CMuNu[Mu][Nu]=0.0, Rho=0; Rho<3; Rho++ )
        CMuNu[Mu][Nu]+=LeviCivita[Mu][Nu][Rho]*R[Rho]*f4*f5;
   };
   
  /* computation of G_{\mu\nu\rho} = d/dR_\rho G_\mu\nu  */
  if (GMuNuRho)
   { 
     for(int Mu=0; Mu<3; Mu++)
      for(int Nu=0; Nu<3; Nu++)
       for(int Rho=0; Rho<3; Rho++)
        GMuNuRho[Mu][Nu][Rho]
         = (  (f1p*f2 + f1*f2p)*Delta[Mu][Nu] 
             +(f1p*f3 + f1*f3p - 2.0*f1*f3/r)*R[Mu]*R[Nu]/r2
           ) * R[Rho] / r
           + f1*f3*(Delta[Mu][Rho]*R[Nu] + R[Mu]*Delta[Nu][Rho])/r2;
   };

  /* computation of C_{\mu\nu\rho} = d/dR_\rho C_\mu\nu  */
  if (CMuNuRho)
   { 
     int Sigma;
     for(int Mu=0; Mu<3; Mu++)
      for(int Nu=0; Nu<3; Nu++)
       for(int Rho=0; Rho<3; Rho++)
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
void CalcGC(double R1[3], double R2[3], cdouble Omega, cdouble EpsR, cdouble MuR, 
            cdouble GMuNu[3][3], cdouble CMuNu[3][3],
            cdouble GMuNuRho[3][3][3], cdouble CMuNuRho[3][3][3])
{ 
  double R[3];

  R[0]=R1[0]-R2[0];
  R[1]=R1[1]-R2[1];
  R[2]=R1[2]-R2[2];

  CalcGC(R,Omega,EpsR,MuR,GMuNu,CMuNu,GMuNuRho,CMuNuRho);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
#if 0
  cdouble GijNew[3][3], CijNew[3][3];
  cdouble k2=k*k, ik=II*k;
  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    { GijNew[Mu][Nu] = (Mu==Nu) ? G0 : 0.0;
      GijNew[Mu][Nu] += ddGBar[3*Mu + Nu]/k2;
    };

  for(int Mu=0; Mu<3; Mu++)
   { int MP1=(Mu+1)%3, MP2=(Mu+2)%3;
     CijNew[Mu][Mu]=0.0;
     CijNew[Mu][MP1]=dGBar[MP2]/(ik);
     CijNew[MP1][Mu]=-CijNew[Mu][MP1];
   };
#endif
/***************************************************************/
cdouble GetG(double R[3], cdouble k, cdouble *dG, cdouble *ddG)
{
  double r2 = R[0]*R[0] + R[1]*R[1] + R[2]*R[2];

  if (r2==0.0)
   {
     if (dG) 
      memset(dG,  0, 3*sizeof(cdouble));
     if (ddG) 
      { memset(ddG, 0, 9*sizeof(cdouble));
        cdouble ddGii = -II*(k*k*k/(12.0*M_PI));
        ddG[3*0+0] = ddG[3*1+1] = ddG[3*2+2] = ddGii;
      };
     return II*k/(4.0*M_PI);
   };

  double r     = sqrt(r2);
  cdouble ikr  = II*k*r;

  cdouble Phi = exp(ikr)/(4.0*M_PI*r);

  if (dG || ddG)
   {  
     cdouble Psi = Phi * (ikr - 1.0) / r2;
     if (dG)
      { dG[0] = R[0] * Psi;
        dG[1] = R[1] * Psi;
        dG[2] = R[2] * Psi;
      };

     if (ddG)
      { double r4     = r2*r2;
        cdouble ikr2  = ikr*ikr;
        cdouble Zeta  = Phi * (ikr2 - 3.0*ikr + 3.0) / r4;
        ddG[3*0+0] = Psi + R[0]*R[0]*Zeta;
        ddG[3*1+1] = Psi + R[1]*R[1]*Zeta;
        ddG[3*2+2] = Psi + R[2]*R[2]*Zeta;
        ddG[3*0+1] = ddG[3*1+0] = R[0]*R[1]*Zeta;
        ddG[3*0+2] = ddG[3*2+0] = R[0]*R[2]*Zeta;
        ddG[3*1+2] = ddG[3*2+1] = R[1]*R[2]*Zeta;
      };
   };

  return Phi;
}

/***************************************************************/
/* Compute the G and C dyadic Green's functions, retaining only*/
/* far-field contributions.                                    */
/***************************************************************/
void CalcGCFarField(double R[3], cdouble k,
                    cdouble GMuNu[3][3], cdouble CMuNu[3][3])
{ 
  double r2=R[0]*R[0] + R[1]*R[1] + R[2]*R[2];

  if (r2==0.0)
   return;

  double r=sqrt(r2);

  cdouble ExpFac=exp(II*k*r)/(4.0*M_PI*r);

  GMuNu[0][0] =               ExpFac * ( 1.0 - R[0]*R[0]/r2 );
  GMuNu[0][1] = GMuNu[1][0] = ExpFac * (     - R[0]*R[1]/r2 );
  GMuNu[0][2] = GMuNu[2][0] = ExpFac * (     - R[0]*R[2]/r2 );
  GMuNu[1][1] =               ExpFac * ( 1.0 - R[1]*R[1]/r2 );
  GMuNu[1][2] = GMuNu[2][1] = ExpFac * (     - R[1]*R[2]/r2 );
  GMuNu[2][2] =               ExpFac * ( 1.0 - R[2]*R[2]/r2 );

  CMuNu[0][0]=CMuNu[1][1]=CMuNu[2][2]=0.0;

  CMuNu[0][1] = ExpFac * R[2]/r;
  CMuNu[1][2] = ExpFac * R[0]/r;
  CMuNu[2][0] = ExpFac * R[1]/r;

  CMuNu[1][0] = -CMuNu[0][1];
  CMuNu[2][1] = -CMuNu[1][2];
  CMuNu[0][2] = -CMuNu[2][0];

}

} // namespace scuff
