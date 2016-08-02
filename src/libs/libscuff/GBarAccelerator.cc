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
 * GBarAccelerator.h -- 
 */

#include <stdlib.h>
#include <math.h>

#include <libhrutil.h>
#include <libMDInterp.h>
#include <libhmat.h>
#include <libscuff.h>

#include "GBarAccelerator.h"

#define II cdouble(0,1)

namespace scuff{

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
/* this is an entry point for GBarVD that has the proper       */
/* prototype for passage to the Interp3D() initialization      */
/* routine; the structure GBarData is defined in libscuffInternals.h.*/
/***************************************************************/
void GBarVDPhi3D(double X1, double X2, double X3, void *UserData, double *PhiVD)
{
  GBarAccelerator *GBA = (GBarAccelerator *)UserData;

  double R[3];
  R[0]=X1;
  R[1]=X2;
  R[2]=X3;

  cdouble GBarVD[8];
  GBarVDEwald(R, GBA->k, GBA->kBloch, GBA->LBV, GBA->LDim,
              -1.0, GBA->ExcludeInnerCells, GBarVD);
 
  PhiVD[ 0] = real(GBarVD[0]);
  PhiVD[ 1] = real(GBarVD[1]);
  PhiVD[ 2] = real(GBarVD[2]);
  PhiVD[ 3] = real(GBarVD[3]);
  PhiVD[ 4] = real(GBarVD[4]);
  PhiVD[ 5] = real(GBarVD[5]);
  PhiVD[ 6] = real(GBarVD[6]);
  PhiVD[ 7] = real(GBarVD[7]);
  PhiVD[ 8] = imag(GBarVD[0]);
  PhiVD[ 9] = imag(GBarVD[1]);
  PhiVD[10] = imag(GBarVD[2]);
  PhiVD[11] = imag(GBarVD[3]);
  PhiVD[12] = imag(GBarVD[4]);
  PhiVD[13] = imag(GBarVD[5]);
  PhiVD[14] = imag(GBarVD[6]);
  PhiVD[15] = imag(GBarVD[7]);

} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetOptimalGridSpacing1D(GBarAccelerator *GBA, double x, double Rho,
                             double RelTol, double OptimalDelta[2],
                             double *MinDelta)
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
  double xMin = x, xMax=x+Delta[0];
  if (EqualFloat(x,L0))
   { xMin=x-Delta[0]; xMax=x; }
  double RhoMin = Rho, RhoMax=Rho+Delta[1];
  if (EqualFloat(Rho,GBA->RhoMax))
   { RhoMin = Rho-Delta[1]; RhoMax=Rho;
   };

  Interp2D *I2D=new Interp2D(xMin,     xMax, 2,
                             RhoMin, RhoMax, 2,
                             2, GBarVDPhi2D, (void *)GBA, LMDI_LOGLEVEL_NONE);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double Phi[8];
  cdouble GExact, GInterp;
  cdouble dGExact[3], dGInterp[3];
  double RelError;

  // estimate relative error in Mu direction
  //int Worst[2];
  for(int Mu=0; Mu<2; Mu++)
   { 
     double R[2];
     R[0] = xMin;
     R[1] = RhoMin;
     R[Mu] += 0.5*Delta[Mu];
     GBarVDPhi2D( R[0], R[1], (void *)GBA, Phi);
     GExact     = cdouble(Phi[0], Phi[4]);
     dGExact[0] = cdouble(Phi[1], Phi[5]);
     dGExact[1] = cdouble(Phi[2], Phi[6]);
     I2D->EvaluatePlus(R[0], R[1], Phi);
     GInterp     = cdouble(Phi[0], Phi[4]);
     dGInterp[0] = cdouble(Phi[1], Phi[5]);
     dGInterp[1] = cdouble(Phi[2], Phi[6]);
     RelError = abs(GInterp-GExact) / abs(GExact);
     OptimalDelta[Mu] = Delta[Mu] * pow( RelTol/RelError, 0.25 );
     //Worst[Mu]=0;
     for(int Nu=0; Nu<2; Nu++)
      { if ( abs(dGExact[Nu]) < 1.0e-6*abs(GExact) ) continue;
        RelError = abs(dGInterp[Nu]-dGExact[Nu]) / abs(dGExact[Nu]);
        double OptDeltaNu = Delta[Mu] * pow( RelTol/RelError, 0.33 );
        if (OptDeltaNu < OptimalDelta[Mu] )
         { OptimalDelta[Mu] = OptDeltaNu;
           //Worst[Mu]=1+Nu;
         };
      };
   };

  delete I2D;

#if 0
  if (RWGGeometry::LogLevel >= SCUFF_VERBOSE2)
   Log("Optimal spacing at (%e,%e)=(%e,%e) (worst: %s, %s)",
        x,Rho,OptimalDelta[0],OptimalDelta[1],
        Worst[0]==0 ? "value" : Worst[0]==1 ? "xDeriv" : "RhoDeriv",
        Worst[1]==0 ? "value" : Worst[1]==1 ? "xDeriv" : "RhoDeriv");
#endif

  if (MinDelta)
   { MinDelta[0] = fmin(MinDelta[0], OptimalDelta[0]);
     MinDelta[1] = fmin(MinDelta[1], OptimalDelta[1]);
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetOptimalGridSpacing2D(GBarAccelerator *GBA,
                             double x, double y, double z,
                             double RelTol, double OptimalDelta[3],
                             double *MinDelta)
{
  double Delta[3];

  if (GBA->LBV[0][1]!=0.0 || GBA->LBV[1][0]!=0.0)
   ErrExit("%s:%i: non-square lattice not yet supported",__FILE__,__LINE__);
  double L0x=GBA->LBV[0][0], L0y=GBA->LBV[1][1];

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  Delta[0] = L0x / 10.0;
  Delta[1] = L0y / 10.0;
  Delta[2] = (GBA->RhoMax - GBA->RhoMin)/10.0;
  if (Delta[2]<=0.0)
   Delta[2] = Delta[0];

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  Interp3D *I3D=new Interp3D(x, x+Delta[0],   2,
                             y, y+Delta[1],   2,
                             z, z+Delta[2],   2,
                             2,   GBarVDPhi3D, (void *)GBA,
                             LMDI_LOGLEVEL_NONE);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double R[3], Phi[16];
  cdouble GExact, GInterp;
  cdouble dGExact[3], dGInterp[3];
  double RelError;
  
  for(int Mu=0; Mu<3; Mu++)
   { 
    // estimate relative error in Mu direction
    R[0] = x;
    R[1] = y;
    R[2] = z;
    R[Mu] += 0.5*Delta[Mu];
    GBarVDPhi3D( R[0], R[1], R[2], (void *)GBA, Phi);
    GExact     = cdouble(Phi[0], Phi[8]);
    dGExact[0] = cdouble(Phi[1], Phi[9]);
    dGExact[1] = cdouble(Phi[2], Phi[10]);
    dGExact[2] = cdouble(Phi[3], Phi[11]);
    I3D->EvaluatePlus(R[0], R[1], R[2], Phi);
    GInterp     = cdouble(Phi[0], Phi[8]);
    dGInterp[0] = cdouble(Phi[1], Phi[9]);
    dGInterp[1] = cdouble(Phi[2], Phi[10]);
    dGInterp[2] = cdouble(Phi[3], Phi[11]);
    RelError = abs(GInterp-GExact) / abs(GExact);
    OptimalDelta[Mu] = Delta[Mu] * pow( RelTol/RelError, 0.25 );
    for(int Nu=0; Nu<3; Nu++)
     { if ( abs(dGExact[Nu]) < 1.0e-6*abs(GExact) ) continue;
       RelError = abs(dGInterp[Nu]-dGExact[Nu]) / abs(dGExact[Nu]);
       double OptDeltaNu = Delta[Mu] * pow( RelTol/RelError, 0.33 );
       OptimalDelta[Mu] = fmin(OptimalDelta[Mu], OptDeltaNu);
     };
   };

  delete I3D;

#if 0
  Log("Optimal spacing at (%e,%e,%e)=(%e,%e,%e)",
       x,y,z,OptimalDelta[0],OptimalDelta[1],OptimalDelta[2]);
#endif

  if (MinDelta)
   { MinDelta[0] = fmin(MinDelta[0], OptimalDelta[0]);
     MinDelta[1] = fmin(MinDelta[1], OptimalDelta[1]);
     MinDelta[2] = fmin(MinDelta[2], OptimalDelta[2]);
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
GBarAccelerator *CreateGBarAccelerator(HMatrix *LBasis,
                                       double RhoMin, double RhoMax,
                                       cdouble k, double *kBloch,
                                       double RelTol, bool ExcludeInnerCells,
                                       int LMDILogLevel)
{
  CheckLattice(LBasis);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  GBarAccelerator *GBA = (GBarAccelerator *)mallocEC( sizeof(GBarAccelerator) );

  GBA->k                  = k;
  GBA->kBloch             = kBloch;
  GBA->ExcludeInnerCells  = ExcludeInnerCells;
  GBA->RhoMin             = RhoMin;
  GBA->RhoMax             = RhoMax;
  GBA->I2D                = 0;
  GBA->I3D                = 0;

  int LDim = GBA->LDim = LBasis->NC;
  for(int nd=0; nd<LDim; nd++)
   for(int j=0; j<3; j++)
    GBA->LBV[nd][j] = LBasis->GetEntryD(j,nd);

  /***************************************************************/
  /* If the wavenumber has an imaginary component, LMax is the   */
  /* distance from the origin at which GBar has decayed to 1e-6  */
  /* of its value at the origin. If LMax lies within the Wigner- */
  /* Seitz cell then we can truncate the range of our            */
  /* interpolation table accordingly.                            */
  /***************************************************************/
  double LMax=1.0e87;
  double Kappa=imag(k);
  if (Kappa!=0.0)
   { 
     double GoK = 14.0/Kappa; // note 14.0 \approx log(1.0e6)
     LMax = sqrt( GoK*GoK + 2.0*RhoMin*GoK );
   };
  GBA->LMax=LMax;

  GBA->ForceFullEwald = false;
  if (RhoMin > RhoMax) 
   GBA->ForceFullEwald=true;
  else 
   { char *str=getenv("SCUFF_EWALD_FULL");
     if ( str && str[0]=='1' )
     { Log("Forcing full Ewald summation.");
       GBA->ForceFullEwald=true;
     };
   };
  
  if (GBA->ForceFullEwald)
   return GBA;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (LDim==1)
   {
      // sample the optimal grid spacing at a few points
      // in the domain of interest and take the smallest
      double L0 = GBA->LBV[0][0], MinDelta[2], Delta[2];
      if ( L0 > LMax )
       { L0=LMax;
         Log("  Cutting off interpolation table at L=LMax=%e.",LMax);
       };
      GetOptimalGridSpacing1D(GBA,  0.0*L0, RhoMin, RelTol, MinDelta, 0);
      GetOptimalGridSpacing1D(GBA, -0.5*L0, RhoMin, RelTol, Delta, MinDelta);
      GetOptimalGridSpacing1D(GBA,  0.0*L0, RhoMax, RelTol, Delta, MinDelta);
      GetOptimalGridSpacing1D(GBA, -0.5*L0, RhoMax, RelTol, Delta, MinDelta);

      double DeltaX = MinDelta[0];
      int nx = ceil( L0 / DeltaX );
      if (nx<2) nx=2;
      DeltaX = L0 / ((double)(nx-1));
      double *XPoints = new double[nx];
      for(int n=0; n<nx; n++)
       XPoints[n] = ((double)n)*DeltaX - 0.5*L0;

      #define MAXRHOPOINTS 1000
      double RhoPoints[MAXRHOPOINTS];
      double MinDeltaRho = (RhoMax-RhoMin) / MAXRHOPOINTS;
      RhoPoints[0]=RhoMin;
      int nRho=0;
      bool Done=false;
      while(!Done)
       { 
         double Rho = RhoPoints[nRho];
         GetOptimalGridSpacing1D(GBA,  0.0*L0, Rho, RelTol, MinDelta, 0);
         GetOptimalGridSpacing1D(GBA, -0.5*L0, Rho, RelTol, Delta, MinDelta);
         double DeltaRho = fmax(MinDeltaRho,MinDelta[1]);
         if (Rho + DeltaRho >= RhoMax)
          { DeltaRho = RhoMax - Rho;
            Done=true;
          };
         RhoPoints[++nRho] = Rho + DeltaRho;
       };
      nRho++; 

      Log("  Initializing %ix%i interpolation table with ",nx,nRho);
      Log("  X points at [%g:%g:%g] ",0.0,DeltaX,L0);
      Log("  Rho points at (");
      for(int n=0; n<(nRho-1); n++)
       LogC("%.2e, ",RhoPoints[n]);
      LogC("%.2e)",RhoPoints[nRho-1]);

      GBA->I3D=0;
      GBA->I2D=new Interp2D(XPoints, nx, RhoPoints, nRho,
                            2, GBarVDPhi2D, (void *)GBA, LMDILogLevel);

      delete[] XPoints;

   }
  else // LDim==2
   {
      double Lx = GBA->LBV[0][0], Ly=GBA->LBV[1][1];
      if ( LMax<Lx ) 
       { Lx=LMax;
         Log("  Cutting off interpolation table at Lx=LMax=%e.",LMax);
       };
      if ( LMax<Ly ) 
       { Ly=LMax;
         Log("  Cutting off interpolation table at Ly=LMax=%e.",LMax);
       };

      // estimate the optimal grid spacings at a few points
      // in the domain of interest and take the smallest
      double MinDelta[3], Delta[3];
      GetOptimalGridSpacing2D(GBA,  0.0*Lx,  0.0*Ly, RhoMin, RelTol, MinDelta, 0);
      GetOptimalGridSpacing2D(GBA,  0.0*Lx, -0.5*Ly, RhoMin, RelTol, Delta, MinDelta);
      GetOptimalGridSpacing2D(GBA, -0.5*Lx,  0.0*Ly, RhoMin, RelTol, Delta, MinDelta);
      GetOptimalGridSpacing2D(GBA, -0.5*Lx, -0.5*Ly, RhoMin, RelTol, Delta, MinDelta);
      if (RhoMax!=RhoMin)
       { GetOptimalGridSpacing2D(GBA,  0.0*Lx,  0.0*Ly, RhoMax, RelTol, Delta, MinDelta);
         GetOptimalGridSpacing2D(GBA,  0.0*Lx, -0.5*Ly, RhoMax, RelTol, Delta, MinDelta);
         GetOptimalGridSpacing2D(GBA, -0.5*Lx,  0.0*Ly, RhoMax, RelTol, Delta, MinDelta);
         GetOptimalGridSpacing2D(GBA, -0.5*Lx, -0.5*Ly, RhoMax, RelTol, Delta, MinDelta);
       };

      int nx   = ceil( Lx / MinDelta[0] );
      if (nx<2) nx=2;

      int ny   = ceil( Ly / MinDelta[1] );
      if (ny<2) ny=2;

      int nRho; 
      if (RhoMax<=RhoMin)
       { RhoMax = RhoMin + MinDelta[0];
         nRho=2;
       }
      else
       { nRho = ceil( (RhoMax-RhoMin) / MinDelta[1] );
         if (nRho<2) nRho=2;
       };

      GBA->I2D=0;
      GBA->I3D=new Interp3D(-0.5*Lx, 0.5*Lx, nx, -0.5*Ly, 0.5*Ly, ny,
                            RhoMin, RhoMax, nRho,
                            2, GBarVDPhi3D, (void *)GBA, LMDILogLevel);
   };

  return GBA;
   
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void DestroyGBarAccelerator(GBarAccelerator *GBA)
{ 
  if (GBA->I2D) delete GBA->I2D;
  if (GBA->I3D) delete GBA->I3D;
  free(GBA);
}

/***************************************************************/
/* Similar to AddGFull, but computes unmixed second partials.  */
/* T[0] += GFull                                               */
/* T[1] += d_{x} GFull                                         */
/* T[2] += d_{y} GFull                                         */
/* T[3] += d_{z} GFull                                         */
/* T[4] += d_{xx} GFull                                        */
/* T[5] += d_{xy} GFull                                        */
/* T[6] += d_{xz} GFull                                        */
/* T[7] += d_{yy} GFull                                        */
/* T[8] += d_{yz} GFull                                        */
/* T[9] += d_{zz} GFull                                        */
/***************************************************************/
void AddGFullTerm(double R[3], cdouble k, double kBloch[2],
                  double Lx, double Ly, cdouble T[10])
{
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
  T[0] += PhaseFactor * Phi;

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
cdouble GetGBarFullEwald(double R[3], GBarAccelerator *GBA,
                         cdouble *dGBar, cdouble *ddGBar)
{
  cdouble G[8];

  cdouble k              = GBA->k;
  double *kBloch         = GBA->kBloch;
  int LDim               = GBA->LDim;
  bool ExcludeInnerCells = GBA->ExcludeInnerCells;

  GBarVDEwald(R, k, kBloch, GBA->LBV, LDim, -1.0, ExcludeInnerCells, G);

  if (dGBar) 
   { dGBar[0]=G[1];
     dGBar[1]=G[2];
     dGBar[2]=G[3];
   };

  if (ddGBar)
   { 
     ddGBar[3*0 + 1] = ddGBar[3*1 + 0] = G[4];
     ddGBar[3*0 + 2] = ddGBar[3*2 + 0] = G[5];
     ddGBar[3*1 + 2] = ddGBar[3*2 + 1] = G[6];
    
     // finite-differencing to get unmixed second partials
     for(int Mu=0; Mu<3; Mu++)
      { 
         double Delta = (R[Mu]==0.0) ? 1.0e-4 : 1.0e-4*fabs(R[Mu]);
         double RR[3];
         RR[0]=R[0]; RR[1]=R[1]; RR[2]=R[2]; 

         RR[Mu] += Delta;
         cdouble Gp[8];
         GBarVDEwald(RR, k, kBloch, GBA->LBV, LDim, -1.0, ExcludeInnerCells, Gp);

         RR[Mu] -= 2.0*Delta;
         cdouble Gm[8];
         GBarVDEwald(RR, k, kBloch, GBA->LBV, LDim, -1.0, ExcludeInnerCells, Gm);

         ddGBar[3*Mu + Mu] = (Gp[0] + Gm[0] - 2.0*G[0]) / (Delta*Delta);
      };

   };

  return G[0];
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
bool InM101(int x) { return (-1<=x) && (x<=1); }

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
  /* transverse coordinate lies outside it, or if the caller      */
  /* requested that we do the calculation directly, then we just  */
  /* do the calculation directly via Ewald summation, skipping    */
  /* the interpolation step                                       */
  /*--------------------------------------------------------------*/
  double Rho2 = R[1]*R[1] + R[2]*R[2];
  double Rho = sqrt(Rho2);
  if (GBA->I2D==0 || Rho<GBA->RhoMin || Rho>GBA->RhoMax)
   return GetGBarFullEwald(R, GBA, dGBar, ddGBar);

  /*--------------------------------------------------------------*/
  /* get xBar, which is the periodic image of the x-coordinate    */
  /* lying within the Wigner-Seitz cell of the 1D lattice         */
  /*--------------------------------------------------------------*/
  double L0 = GBA->LBV[0][0];
  int m = (int)(lround( R[0] / L0 ));
  double xBar = R[0] - m*L0;

  if ( fabs(xBar) > GBA->LMax )
   { if (dGBar) memset(dGBar, 0, 3*sizeof(cdouble));
     if (ddGBar) memset(ddGBar, 0, 9*sizeof(cdouble));
     return 0.0;
   };

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
        memset(TP, 0, 10*sizeof(cdouble));
        memset(TM, 0, 10*sizeof(cdouble));
        double RBar[3];
        RBar[0] = xBar;
        RBar[1] = Rho;
        RBar[2] = 0.0;
        for(int n=-1; n<=1; n++)
         { if ( !InM101(n+m) )
            AddGFullTerm(RBar, k, kBloch, n*L0, 0.0, TP);
           if ( !InM101(n-m) )
            AddGFullTerm(RBar, k, kBloch, (n-m)*L0, 0.0, TM);
         };

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
  cdouble k              = GBA->k;
  double *kBloch         = GBA->kBloch;
  bool ExcludeInnerCells = GBA->ExcludeInnerCells;

  /*--------------------------------------------------------------*/
  /* if we have no interpolation grid, or we have one but the     */
  /* transverse coordinate lies outside it, or if the caller      */
  /* requested that we do the calculation directly, then we just  */
  /* do the calculation directly via Ewald summation, skipping    */
  /* the interpolation step                                       */
  /*--------------------------------------------------------------*/
  double Rho = fabs(R[2]);
  if (GBA->I3D==0 || Rho<GBA->RhoMin || Rho>GBA->RhoMax)
   return GetGBarFullEwald(R, GBA, dGBar, ddGBar);
 
  /*--------------------------------------------------------------*/
  /* get (xBar, yBar) = representative of (x,y) in the Wigner-    */
  /* Seitz cell                                                   */
  /*--------------------------------------------------------------*/
  if (GBA->LBV[0][1]!=0.0 || GBA->LBV[1][0]!=0.0)
   ErrExit("%s:%i: non-square lattice not yet supported",__FILE__,__LINE__);
  
  double L0x = GBA->LBV[0][0];
  int mx = (int)(lround( R[0] / L0x ));
  double xBar = R[0] - mx*L0x;

  double L0y = GBA->LBV[1][1];
  int my = (int)(lround( R[1] / L0y ));
  double yBar = R[1] - my*L0y;

  if ( fabs(xBar)>GBA->LMax || fabs(yBar)>GBA->LMax )
   { if (dGBar) memset(dGBar, 0, 3*sizeof(cdouble));
     if (ddGBar) memset(ddGBar, 0, 9*sizeof(cdouble));
     return 0.0;
   };

  /*--------------------------------------------------------------*/
  /* get GBar(RBar, Rho) and as many derivatives as necessary     */
  /*--------------------------------------------------------------*/
  cdouble GBar;
  if (ddGBar)
   { 
     double Phi[20];
     GBA->I3D->EvaluatePlusPlus(xBar, yBar, Rho, Phi);
     GBar            = cdouble(Phi[0], Phi[10]);
     dGBar[0]        = cdouble(Phi[1], Phi[11]);
     dGBar[1]        = cdouble(Phi[2], Phi[12]);
     dGBar[2]        = cdouble(Phi[3], Phi[13]);
     ddGBar[3*0 + 0] = cdouble(Phi[4], Phi[14]);
     ddGBar[3*0 + 1] = cdouble(Phi[5], Phi[15]);
     ddGBar[3*0 + 2] = cdouble(Phi[6], Phi[16]);
     ddGBar[3*1 + 1] = cdouble(Phi[7], Phi[17]);
     ddGBar[3*1 + 2] = cdouble(Phi[8], Phi[18]);
     ddGBar[3*2 + 2] = cdouble(Phi[9], Phi[19]);
   }
  else if (dGBar)
   { 
     double Phi[16];
     GBA->I3D->EvaluatePlus(xBar, yBar, Rho, Phi);
     GBar      = cdouble(Phi[0], Phi[8]);
     dGBar[0]  = cdouble(Phi[1], Phi[9]);
     dGBar[1]  = cdouble(Phi[2], Phi[10]);
     dGBar[2]  = cdouble(Phi[3], Phi[11]);
   }
  else
   GBA->I3D->Evaluate(xBar, yBar, Rho, (double *)&GBar);

  /*--------------------------------------------------------------*/
  /* correct for the fact that we may have needed to              */
  /* translate the evaluation point into the unit cell            */
  /*--------------------------------------------------------------*/
  if (mx!=0 || my!=0)
   { 
     cdouble ExpArg = kBloch[0]*( mx*(GBA->LBV[0][0]) + my*(GBA->LBV[1][0]) )
                     +kBloch[1]*( mx*(GBA->LBV[0][1]) + my*(GBA->LBV[1][1]) );
     cdouble PhaseFactor = exp( II*ExpArg );

     if (ExcludeInnerCells)
      { cdouble TP[10], TM[10];
        memset(TP, 0, 10*sizeof(cdouble));
        memset(TM, 0, 10*sizeof(cdouble));
        double RBar[3];
        RBar[0] = xBar;
        RBar[1] = yBar;
        RBar[2] = Rho;
        for(int nx=-1; nx<=1; nx++)
         for(int ny=-1; ny<=1; ny++)
         { if ( ! (InM101(nx+mx) && InM101(ny+my) ) )
            AddGFullTerm(RBar, k, kBloch, nx*L0x, ny*L0y, TP);
           if ( ! (InM101(nx-mx) && InM101(ny-my) ) )
            AddGFullTerm(RBar, k, kBloch, (nx-mx)*L0x, (ny-my)*L0y, TM);
         };
        GBar += TP[0]-TM[0];
        if (dGBar)
         for(int Mu=0; Mu<3; Mu++) dGBar[Mu] += TP[Mu+1] - TM[Mu+1];
        if (ddGBar)
         { ddGBar[3*0 + 0] += TP[4] - TM[4];
           ddGBar[3*0 + 1] += TP[5] - TM[5];
           ddGBar[3*0 + 2] += TP[6] - TM[6];
           ddGBar[3*1 + 1] += TP[7] - TM[7];
           ddGBar[3*1 + 2] += TP[8] - TM[8];
           ddGBar[3*2 + 2] += TP[9] - TM[9];
         };

      }; // if (ExcludeInnerCells)

     GBar      *= PhaseFactor;
     if (dGBar) for(int Mu=0; Mu<3; Mu++) dGBar[Mu]*=PhaseFactor;
     if (ddGBar) for(int Mu=0; Mu<9; Mu++) ddGBar[Mu]*=PhaseFactor;
        
   }; // if (m!=0)

  /*--------------------------------------------------------------*/
  /*- convert derivatives from the (x,y,Rho) system to the (x,y,z) system  */
  /*--------------------------------------------------------------*/
  double Sign = R[2] > 0.0 ? 1.0 : -1.0;
  if (dGBar)
   dGBar[2] *= Sign;
  if (ddGBar)
   { 
     ddGBar[3*0 + 2] *= Sign;
     ddGBar[3*1 + 2] *= Sign;

     ddGBar[3*1 + 0] = ddGBar[3*0 + 1];
     ddGBar[3*2 + 0] = ddGBar[3*0 + 2];
     ddGBar[3*2 + 1] = ddGBar[3*1 + 2];
   };

  return GBar;
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
cdouble GetGBar(double R[3], GBarAccelerator *GBA,
                cdouble *dGBar, cdouble *ddGBar, 
                bool ForceFullEwald)
{
  /***************************************************************/
  /* there are two ways to force the code to bypass the          */
  /* interpolation table in favor of doing full Ewald summation: */
  /*  (a) set the ForceFullEwald parameter to this function      */
  /*      to true                                                */
  /*  (b) set the ForceFullEwald field in the GBA structure      */
  /*      to true                                                */
  /* The reason we need both mechanisms is that (b) does not     */
  /* allow for thread-specific selection of full Ewald summation */
  /* (because the GBA structure is typically shared among many   */
  /*  threads).                                                  */
  /***************************************************************/
  if (ForceFullEwald || GBA->ForceFullEwald)
   return GetGBarFullEwald(R, GBA, dGBar, ddGBar);
  if (GBA->LDim==1)
   return GetGBar_1D(R, GBA, dGBar, ddGBar);
  else
   return GetGBar_2D(R, GBA, dGBar, ddGBar);
}

/***************************************************************/
/* Create a GBar accelerator suitable for computing GBar at    */
/* points R in the three-dimensional box with corners RMin and */
/* RMax.                                                       */
/***************************************************************/
GBarAccelerator *RWGGeometry::CreateRegionGBA(int nr, cdouble Omega, double *kBloch,
                                              double RMin[3], double RMax[3],
                                              bool ExcludeInnerCells)
{
  double Buffer[3];
  HMatrix LBasis1D(3,1,LHM_REAL,LHM_NORMAL,Buffer);
  HMatrix *GBA_LBasis;

  /***************************************************************/
  /* the lattice with respect to which the region is periodic    */
  /* is not necessarily the same as the full lattice of the      */
  /* geometry; it may be a proper sublattice.                    */
  /***************************************************************/
  int NumExtendedDimensions=0;
  if (RegionIsExtended[0][nr]) NumExtendedDimensions++;
  if (LDim>1 && RegionIsExtended[1][nr]) NumExtendedDimensions++;

  if (NumExtendedDimensions == 0 ) 
   return 0; // region is not extended
  else if (NumExtendedDimensions == LDim )
   GBA_LBasis=LBasis; // region has full periodicity of geometry
  else // region periodicity lattice is sublattice of full geometry lattice
   { 
     int ExtendedDim = RegionIsExtended[0][nr] ? 0 : 1;

     GBA_LBasis = &LBasis1D;
     for(int j=0; j<3; j++)
      GBA_LBasis->SetEntry(j,0,LBasis->GetEntryD(j,ExtendedDim));
   };
   
  /***************************************************************/
  /* do not bother to create an interpolation table if we are in */
  /* the regime in which GBar will be exponentially tiny at all  */
  /* evaluation points                                           */
  /***************************************************************/
  cdouble k = Omega * RegionMPs[nr]->GetRefractiveIndex(Omega);
  if ( ExcludeInnerCells && imag(k)>0.0 )
   { 
     double minL
      =2.0*sqrt(  pow( GBA_LBasis->GetEntryD(0,0), 2.0 )
                 +pow( GBA_LBasis->GetEntryD(1,0), 2.0 )
               );

     if (GBA_LBasis->NC >1)
      minL = fmin(minL, 2.0*sqrt(  pow( GBA_LBasis->GetEntryD(0,1), 2.0 )
                                  +pow( GBA_LBasis->GetEntryD(1,1), 2.0 )
                                )
                 );

     if ( imag(k)*minL > 14.0 )
      return 0;
   };

  /***************************************************************/
  /* determine the minimum and maximum values of the transverse  */
  /* coordinate Rho that our interpolation table will need to    */
  /* support. This requires a little care if RMin < 0, because   */
  /* we have to distinguish between cases in which the interval  */
  /* [RMin, RMax] does or does not straddle the origin.          */
  /* If RMin[0] > RMax[0], then we create a dummy 'accelerator'  */
  /* in which no interpolation table is created and GBar is      */
  /* always evaluated by full Ewald summation.                   */
  /***************************************************************/
  double RhoMin, RhoMax;
  if (RMin[0] > RMax[0])
   { RhoMin=1.0;
     RhoMax=0.0;
   }
  else
   { double AbsRMin[3], AbsRMax[3];
     for(int Mu=0; Mu<3; Mu++)
      { if ( (RMin[Mu]>=0.0) )
         { AbsRMin[Mu] = fabs(RMin[Mu]);
           AbsRMax[Mu] = fabs(RMax[Mu]);
         }
        else if ( (RMin[Mu]<0.0) && (RMax[Mu]>0.0) )
         { AbsRMin[Mu] = 0.0;
           AbsRMax[Mu] = fmax( fabs(RMax[Mu]), fabs(RMin[Mu]) );
         }
        else // ( (RMin[Mu]<0.0) && (RMax[Mu]<0.0) )
         { AbsRMin[Mu] = fabs(RMax[Mu]);
           AbsRMax[Mu] = fabs(RMin[Mu]);
         };
      };
    
     if (LDim==1)
      { RhoMin = sqrt(AbsRMin[1]*AbsRMin[1] + AbsRMin[2]*AbsRMin[2]);
        RhoMax = sqrt(AbsRMax[1]*AbsRMax[1] + AbsRMax[2]*AbsRMax[2]);
      }
     else
      { RhoMin = AbsRMin[2];
        RhoMax = AbsRMax[2];
      };
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double RelTol = 1.0e-3;
  char *str=getenv("SCUFF_INTERPOLATION_TOLERANCE");
  if ( str )
   { sscanf(str,"%le",&RelTol);
     Log("Setting interpolation tolerance to %e.",RelTol);
   };

  int LMDILogLevel = LMDI_LOGLEVEL_TERSE;
  if (LogLevel>=SCUFF_VERBOSE2)
   LMDILogLevel = LMDI_LOGLEVEL_VERBOSE;
  return CreateGBarAccelerator(GBA_LBasis, RhoMin, RhoMax,
                               k, kBloch, RelTol, ExcludeInnerCells,
                               LMDILogLevel);
}

/***************************************************************/
/* Create a GBarAccelerator object appropriate for computing   */
/* the periodic contribution to the scattered fields at a list */
/* of evaluation points.                                       */
/***************************************************************/
GBarAccelerator *RWGGeometry::CreateRegionGBA(int nr, cdouble Omega, double *kBloch, HMatrix *XMatrix)
{
  if (XMatrix->NR < 8)
   { double RMin[3]={1.0, 1.0, 1.0};
     double RMax[3]={0.0, 0.0, 0.0};
     return CreateRegionGBA(nr, Omega, kBloch, RMin, RMax, false);
   };

  if (LogLevel>=SCUFF_VERBOSE2)
   Log("Creating GBar accelerator for %i points in region %s...",XMatrix->NR,RegionLabels[nr]);

  /*--------------------------------------------------------------*/
  /* get bounding box enclosing all eval points in this region    */
  /*--------------------------------------------------------------*/
  int NumPointsInRegion=0;
  double PointXMax[3]={-1.0e89, -1.0e89, -1.0e89};
  double PointXMin[3]={+1.0e89, +1.0e89, +1.0e89}; 
  for(int np=0; np<XMatrix->NR; np++)
   { 
     double X[3];
     X[0] = XMatrix->GetEntryD(np,0);
     X[1] = XMatrix->GetEntryD(np,1);
     X[2] = XMatrix->GetEntryD(np,2);

     if ( PointInRegion(nr, X ) )
      { 
        NumPointsInRegion++;
        PointXMax[0] = fmax(PointXMax[0], X[0] );
        PointXMax[1] = fmax(PointXMax[1], X[1] );
        PointXMax[2] = fmax(PointXMax[2], X[2] );
        PointXMin[0] = fmin(PointXMin[0], X[0] );
        PointXMin[1] = fmin(PointXMin[1], X[1] );
        PointXMin[2] = fmin(PointXMin[2], X[2] );
      };
   };

  if (NumPointsInRegion==0) 
   return 0;

  /*--------------------------------------------------------------*/
  /* get the maximum and minimum values of X-Y where              */
  /* X runs over all vertices on all surfaces bounding the region */
  /* and Y runs over all evaluation points in the region          */
  /*--------------------------------------------------------------*/
  double RMax[3]={-1.0e89, -1.0e89, -1.0e89};
  double RMin[3]={+1.0e89, +1.0e89, +1.0e89};
  for(int ns=0; ns<NumSurfaces; ns++)
   { 
     if (    Surfaces[ns]->RegionIndices[0] == nr
          || Surfaces[ns]->RegionIndices[1] == nr
        )
      { 
        double *SurfaceXMax = Surfaces[ns]->RMax;
        double *SurfaceXMin = Surfaces[ns]->RMin;
        RMax[0] = fmax(RMax[0], SurfaceXMax[0] - PointXMin[0] );
        RMax[1] = fmax(RMax[1], SurfaceXMax[1] - PointXMin[1] );
        RMax[2] = fmax(RMax[2], SurfaceXMax[2] - PointXMin[2] );
        RMin[0] = fmin(RMin[0], SurfaceXMin[0] - PointXMax[0] );
        RMin[1] = fmin(RMin[1], SurfaceXMin[1] - PointXMax[1] );
        RMin[2] = fmin(RMin[2], SurfaceXMin[2] - PointXMax[2] );
      };
   };
  if ( RMin[0]>RMax[0] || RMin[1]>RMax[1] || RMin[2]>RMax[2] )
   ErrExit("%s:%i: internal error (%e,%e,%e) (%e,%e,%e)",
           __FILE__,__LINE__,RMin[0],RMin[1],RMin[2],RMax[0],RMax[1],RMax[2]);

  return CreateRegionGBA(nr, Omega, kBloch, RMin, RMax, false);

}

/***************************************************************/
/* Create a GBarAccelerator object appropriate for computing   */
/* the periodic contribution to the BEM matrix block for two   */
/* RWG surfaces.                                               */
/***************************************************************/
GBarAccelerator *RWGGeometry::CreateRegionGBA(int nr, cdouble Omega, double *kBloch, int ns1, int ns2)
{
  if (LogLevel>=SCUFF_VERBOSE2)
   Log("Creating GBar accelerator for region %s, surfaces %s,%s",
        RegionLabels[nr],Surfaces[ns1]->Label,Surfaces[ns2]->Label);

  /***************************************************************/
  /* get the maximum and minimum values of the cartesian         */
  /* coordinates of R=x1-x2 as x1, x2 range over surfaces 1,2    */
  /***************************************************************/
  double RMax[3], RMin[3];
  VecSub(Surfaces[ns1]->RMax, Surfaces[ns2]->RMin, RMax);
  VecSub(Surfaces[ns1]->RMin, Surfaces[ns2]->RMax, RMin);

  return CreateRegionGBA(nr, Omega, kBloch, RMin, RMax, true);

}

} // namespace scuff
