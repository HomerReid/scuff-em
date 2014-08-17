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
#include <stdlib.h>
#include <math.h>

#include <libhrutil.h>
#include "/home/homer/work/scuff-em/src/libs/libMDInterp/libMDInterp.h"
#include <libscuff.h>
#include <libscuffInternals.h>

#include "GBarAccelerator.h"

#define II cdouble(0,1)

namespace scuff{

/***************************************************************/
/* prototypes for routines in GBarVDEwald                      */
/***************************************************************/
void GBarVDEwald(double *R, cdouble k, double *kBloch,
                 double *LBV[2], int LDim,
                 double E, bool ExcludeInnerCells,
                 cdouble *GBarVD);

/***************************************************************/
/* entry point for GBarVD that has the proper prototype for    */
/* passage to the Interp2D() initialization routine.           */
/***************************************************************/
void GBarVDPhi2D(double x, double Rho, void *UserData, double *PhiVD)
{
  GBarAccelerator *GBA = (GBarAccelerator *)UserData;
 
  double R[3];
  R[0]=x;
  R[1]=Rho;
  R[2]=0.0;

  cdouble GBarVD[8];
  GBarVDEwald(R, GBA->k, GBA->kBloch, GBA->LBV, GBA->LDim,
              -1.0, GBA->ExcludeInnerCells, GBarVD);
 
  PhiVD[0] = real(GBarVD[0]); // real(G)
  PhiVD[1] = real(GBarVD[1]); // real(dGdX)
  PhiVD[2] = real(GBarVD[2]); // real(dGdY) = real(dGdRho)
  PhiVD[3] = real(GBarVD[4]); // real(dG2dXdY) = real(dG2dXdRho)

  PhiVD[4] = imag(GBarVD[0]); // imag(G)
  PhiVD[5] = imag(GBarVD[1]); // imag(dGdX)
  PhiVD[6] = imag(GBarVD[2]); // imag(dGdRho)
  PhiVD[7] = imag(GBarVD[4]); // imag(dG2dXdRho)

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetOptimalGridSpacing1D(GBarAccelerator *GBA, double x, double Rho,
                             double RelTol, double OptimalDelta[2])
{
  double Delta[2];

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double L0 = GBA->LBV[0][0];
  Delta[0] = L0 / 10.0;
  Delta[1] = (GBA->RhoMax - GBA->RhoMin)/10.0;
  if (Delta[1]<=0.0)
   Delta[1] = Delta[0];

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  Interp2D *I2D=new Interp2D(x,   x+Delta[0],   2,
                             Rho, Rho+Delta[1], 2,
                             2, GBarVDPhi2D, (void *)GBA);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double Phi[8];
  cdouble GExact, GInterp;
  double RelError, ScaleFactor;

  // estimate relative error in X direction
  GBarVDPhi2D( x+0.5*Delta[0], Rho, (void *)GBA, Phi);
  GExact=cdouble(Phi[0], Phi[4]);
  I2D->Evaluate(x+0.5*Delta[0], Rho, Phi);
  GInterp=cdouble(Phi[0], Phi[1]);
  RelError = abs(GInterp-GExact) / abs(GExact);
  OptimalDelta[0] = Delta[0] * pow( RelTol/RelError, 0.25 );

  // estimate relative error in Rho direction
  GBarVDPhi2D( x, Rho+0.5*Delta[1], (void *)GBA, Phi);
  GExact=cdouble(Phi[0], Phi[4]);
  I2D->Evaluate(x, Rho+0.5*Delta[1], Phi);
  GInterp=cdouble(Phi[0], Phi[1]);
  RelError = abs(GInterp-GExact) / abs(GExact);
  OptimalDelta[1] = Delta[1] * pow( RelTol/RelError, 0.25 );

  delete I2D;

  Log("Optimal spacing at (%e,%e)=(%e,%e)",
       x,Rho,OptimalDelta[0],OptimalDelta[1]);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
GBarAccelerator *CreateGBarAccelerator(int LDim, double *LBV[2],
                                       double RhoMin, double RhoMax,
                                       cdouble k, double *kBloch,
                                       double RelTol, bool ExcludeInnerCells)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  GBarAccelerator *GBA = (GBarAccelerator *)mallocEC( sizeof(GBarAccelerator) );

  GBA->k                 = k;
  GBA->kBloch            = kBloch;
  GBA->ExcludeInnerCells = ExcludeInnerCells;
  GBA->RhoMin            = RhoMin;
  GBA->RhoMax            = RhoMax;

  GBA->LDim = LDim;
  GBA->LBV1[0]=LBV[0][0];
  GBA->LBV1[1]=LBV[0][1];
  if (LDim==2)
   { GBA->LBV2[0]=LBV[1][0];
     GBA->LBV2[1]=LBV[1][1];
   };
  GBA->LBV[0] = GBA->LBV1;
  GBA->LBV[1] = GBA->LBV2;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (LDim==1)
   {
      // estimate the optimal grid spacings at a few points
      // in the domain of interest and take the smallest 
      double L0 = LBV[0][0], MinDelta[2], Delta[2];
      GetOptimalGridSpacing1D(GBA, 0.0*L0, RhoMin, RelTol, MinDelta);

      GetOptimalGridSpacing1D(GBA, 0.5*L0, RhoMin, RelTol, Delta);
      MinDelta[0]=fmin(MinDelta[0], Delta[0]);
      MinDelta[1]=fmin(MinDelta[1], Delta[1]);

      GetOptimalGridSpacing1D(GBA, 1.0*L0, RhoMin, RelTol, Delta);
      MinDelta[0]=fmin(MinDelta[0], Delta[0]);
      MinDelta[1]=fmin(MinDelta[1], Delta[1]);

      GetOptimalGridSpacing1D(GBA, 0.0*L0, RhoMax, RelTol, MinDelta);
      MinDelta[0]=fmin(MinDelta[0], Delta[0]);
      MinDelta[1]=fmin(MinDelta[1], Delta[1]);

      GetOptimalGridSpacing1D(GBA, 0.5*L0, RhoMax, RelTol, Delta);
      MinDelta[0]=fmin(MinDelta[0], Delta[0]);
      MinDelta[1]=fmin(MinDelta[1], Delta[1]);

      GetOptimalGridSpacing1D(GBA, 1.0*L0, RhoMax, RelTol, Delta);
      MinDelta[0]=fmin(MinDelta[0], Delta[0]);
      MinDelta[1]=fmin(MinDelta[1], Delta[1]);

      int nx   = ceil( L0 / MinDelta[0] );
      int nRho = ceil( (RhoMax-RhoMin) / MinDelta[1] );

      if (nx<2) nx=2;
      if (nRho<2) nRho=2;

      printf("Creating %ix%i interpolation table...\n",nx,nRho);
      GBA->I3D=0;
      GBA->I2D=new Interp2D(0.0,    L0,     nx,
                            RhoMin, RhoMax, nRho,
                            2, GBarVDPhi2D, (void *)GBA);
   }
  else // LDim==2
   {
      GBA->I2D=0;
#if 0
      GBA->I3D=new Interp3D(0.0, LBV[0][0], nx,
                            RhoMin, RhoMax, nRho,
                            2, GBarVDPhi2D, (void *)GBA);
#endif
   };

  return GBA;
   
}

/***************************************************************/
/* Similar to AddGFull, but does not accumulate, and also      */
/* computes unmixed second partials.                           */
/* T[0] = GFull                                                */
/* T[1] = d_{x} GFull                                          */
/* T[2] = d_{y} GFull                                          */
/* T[3] = d_{z} GFull                                          */
/* T[4] = d_{xx} GFull                                         */
/* T[5] = d_{xy} GFull                                         */
/* T[6] = d_{xz} GFull                                         */
/* T[7] = d_{yy} GFull                                         */
/* T[8] = d_{yz} GFull                                         */
/* T[9] = d_{zz} GFull                                         */
/***************************************************************/
void GetGFullTerm(double R[3], cdouble k, double kBloch[2],
                  double Lx, double Ly, cdouble T[10])
{
  memset(T, 0, 10*sizeof(cdouble));
  double RmL[3];
  double r2, r;
  cdouble PhaseFactor, IKR, Phi, Psi, Zeta, Upsilon;
   
  PhaseFactor=exp( II*(Lx*kBloch[0] + Ly*kBloch[1]) );

  RmL[0]=R[0] - Lx;
  RmL[1]=R[1] - Ly;
  RmL[2]=R[2];

  r2=RmL[0]*RmL[0] + RmL[1]*RmL[1] + RmL[2]*RmL[2];
  r=sqrt(r2);
  if ( r < 1.0e-8 )
   return;
  IKR=II*k*r;
  Phi=exp(IKR)/(4.0*M_PI*r);
  T[0] = PhaseFactor * Phi;

  Psi=(IKR-1.0)*Phi/r2;
  Zeta=(3.0 + IKR*(-3.0 + IKR))*Phi/(r2*r2);
  Upsilon=(-15.0 + IKR*(15.0 + IKR*(-6.0 + IKR)))*Phi/(r2*r2*r2);

  T[1] += PhaseFactor * RmL[0] * Psi;
  T[2] += PhaseFactor * RmL[1] * Psi;
  T[3] += PhaseFactor * RmL[2] * Psi;
  T[4] += PhaseFactor * (Psi + RmL[0] * RmL[0] * Zeta);
  T[5] += PhaseFactor * (      RmL[0] * RmL[1] * Zeta);
  T[6] += PhaseFactor * (      RmL[0] * RmL[2] * Zeta);
  T[7] += PhaseFactor * (Psi + RmL[1] * RmL[1] * Zeta);
  T[8] += PhaseFactor * (      RmL[1] * RmL[2] * Zeta);
  T[9] += PhaseFactor * (Psi + RmL[2] * RmL[2] * Zeta);
 
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
cdouble GetGBar_1D(double R[3], GBarAccelerator *GBA,
                   cdouble *dGBar, cdouble *ddGBar)
{
  cdouble k              = GBA->k;
  double *kBloch         = GBA->kBloch;
  bool ExcludeInnerCells = GBA->ExcludeInnerCells;

  /*--------------------------------------------------------------*/
  /* if we have no interpolation grid, or we have one but the     */
  /* transverse coordinate lies outside it, then we just do the   */
  /* calculation directly                                         */
  /*--------------------------------------------------------------*/
  double Rho2 = R[1]*R[1] + R[2]*R[2];
  double Rho = sqrt(Rho2);
  if ( GBA->I2D==0 || Rho<GBA->RhoMin || Rho>GBA->RhoMax )
   { cdouble G[8];
     GBarVDEwald(R, k, kBloch, GBA->LBV, 1, -1.0, ExcludeInnerCells, G);
     if (dGBar) 
      { dGBar[0]=G[1];
        dGBar[1]=G[2];
        dGBar[2]=G[3];
      };
     if (ddGBar)
      memset(ddGBar, 0, 9*sizeof(cdouble));
     return G[0];
   };
 
  /*--------------------------------------------------------------*/
  /* get xBar = unit-cell representative of x coordinate          */
  /*--------------------------------------------------------------*/
  double L0 = GBA->LBV[0][0];
  int m = (int)(floor( R[0] / L0 )); 
  double xBar = R[0] - m*L0;

  /*--------------------------------------------------------------*/
  /* get GBar(xBar, Rho) and as many derivatives as necessary     */
  /*--------------------------------------------------------------*/
  cdouble G, dGdx=0, dGdRho=0.0, d2Gdx2=0.0, d2GdxdRho=0.0, d2GdRho2=0.0;
  if (ddGBar)
   { 
     double Phi[12];
     GBA->I2D->EvaluatePlusPlus(xBar, Rho, Phi);
     G         = cdouble(Phi[0], Phi[ 6]);
     dGdx      = cdouble(Phi[1], Phi[ 7]);
     dGdRho    = cdouble(Phi[2], Phi[ 8]);
     d2Gdx2    = cdouble(Phi[3], Phi[ 9]);
     d2GdxdRho = cdouble(Phi[4], Phi[10]);
     d2GdRho2  = cdouble(Phi[5], Phi[11]);
   }
  else if (dGBar)
   { 
     double Phi[8];
     GBA->I2D->EvaluatePlus(xBar, Rho, Phi);
     G         = cdouble(Phi[0], Phi[4]);
     dGdx      = cdouble(Phi[1], Phi[5]);
     dGdRho    = cdouble(Phi[2], Phi[6]);
   }
  else
   GBA->I2D->Evaluate(xBar, Rho, (double *)&G);

  /*--------------------------------------------------------------*/
  /* correct for the fact that we may have needed to              */
  /* translate the evaluation point into the unit cell            */
  /*--------------------------------------------------------------*/
  if (m!=0)
   { 
     cdouble PhaseFactor = exp( II*((double) m)*kBloch[0]*L0 );

     if (ExcludeInnerCells)
      { 
        cdouble TP[10], TM[10];
        if (m==1) 
         { 
           GetGFullTerm(R, k, kBloch,  1.0*L0, 0.0, TP);   // T_{+1}
           GetGFullTerm(R, k, kBloch, -2.0*L0, 0.0, TM);   // T_{+2}
         }
        else if (m==-1) 
         { 
           GetGFullTerm(R, k, kBloch, -1.0*L0, 0.0, TP);   // T_{+1}
           GetGFullTerm(R, k, kBloch,  2.0*L0, 0.0, TM);   // T_{+2}
         }
        else 
         ErrExit("evaluation point too far from origin in GetGBar");

        G         += TP[0]-TM[0];
        dGdx      += TP[1]-TM[1];
        dGdRho    += TP[2]-TM[2];
        d2Gdx2    += TP[4]-TM[4];
        d2GdxdRho += TP[5]-TM[5];
        d2GdRho2  += TP[7]-TM[7];

      }; // if (ExcludeInnerCells)

     G         *= PhaseFactor;
     dGdx      *= PhaseFactor;
     dGdRho    *= PhaseFactor;
     d2Gdx2    *= PhaseFactor;
     d2GdxdRho *= PhaseFactor;
     d2GdRho2  *= PhaseFactor;
        
   }; // if (m!=0)

  /*--------------------------------------------------------------*/
  /*- convert derivatives from the (x,Rho) system to the (x,y,z) system  */
  /*--------------------------------------------------------------*/
  if (Rho==0.0)
   { dGdRho=d2GdxdRho=0.0;
     Rho=1.0; // avoid division-by-zero
   };

  if (dGBar)
   { dGBar[0]  = dGdx;
     dGBar[1]  = (R[1]/Rho) * dGdRho;  // dG/dy = (y/Rho) dG/dRho
     dGBar[2]  = (R[2]/Rho) * dGdRho;  // dG/dz = (z/Rho) dG/dRho
   };

  if (ddGBar) 
   { 
     ddGBar[ 3*0 + 0 ] = d2Gdx2;

     if (Rho2==0.0)
      { ddGBar[ 3*1 + 1 ] = ddGBar[3*2 + 2] = d2GdRho2;
        ddGBar[ 3*1 + 2 ] = ddGBar[3*2 + 1] = 0.0;
      }
     else 
      { double Rho3 = Rho*Rho2;
        ddGBar[ 3*1 + 1 ] = (R[2]*R[2]*dGdRho + R[1]*R[1]*Rho*d2GdRho2)/Rho3;
        ddGBar[ 3*2 + 2 ] = (R[1]*R[1]*dGdRho + R[2]*R[2]*Rho*d2GdRho2)/Rho3;
        ddGBar[ 3*1 + 2 ] = ddGBar[ 3*2 + 1 ] = R[1]*R[2]*(-dGdRho + Rho*d2GdRho2)/Rho3;
      };

     ddGBar[ 3*0 + 1 ] = ddGBar[ 3*1 + 0 ] = (R[1]/Rho) * d2GdxdRho;
     ddGBar[ 3*0 + 2 ] = ddGBar[ 3*2 + 0 ] = (R[2]/Rho) * d2GdxdRho;
   };

  return G;
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
cdouble GetGBar_2D(double R[3], GBarAccelerator *GBA, 
                   cdouble *dGBar, cdouble *ddGBar)
{
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
cdouble GetGBar(double R[3], GBarAccelerator *GBA,
                cdouble *dGBar, cdouble *ddGBar)
             
{
  if (GBA->LDim==1)
   return GetGBar_1D(R, GBA, dGBar, ddGBar);
  else
   return GetGBar_2D(R, GBA, dGBar, ddGBar);
}

} // namespace scuff
