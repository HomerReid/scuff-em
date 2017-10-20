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
 * FullWave_SurfaceCurrents.cc -- compute the full-wave green's function
 *                             -- above a multi-layer dielectric substrate
 *                             -- using the method of induced surface currents
 *
 * homer reid                  -- 3/2017-9/2017
 *
 */
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>

#include "libhrutil.h"
#include "libMDInterp.h"
#include "libSGJC.h"
#include "libSpherical.h"
#include "libSubstrate.h"

#define II cdouble(0.0,1.0)
#ifndef ZVAC
#define ZVAC 376.73031346177
#endif

#define SQRT2 1.41421356237309504880

/***************************************************************/
/* compute the 3x3 submatrices that enter the quadrants of the */
/* Fourier transform of the (6x6) homogeneous DGF              */
/***************************************************************/
void GetGC0Twiddle(cdouble k2, double q2D[2], cdouble qz, double Sign,
                   cdouble GT[3][3], cdouble CT[3][3])
{
  double qx  = q2D[0];
  double qy  = q2D[1];

  GT[0][0] = 1.0 - qx*qx/k2;
  GT[1][1] = 1.0 - qy*qy/k2;
  GT[2][2] = 1.0 - qz*qz/k2;

  GT[0][1] = GT[1][0] =      -1.0*qx*qy/k2;

  GT[1][2] = GT[2][1] = -1.0*Sign*qy*qz/k2;

  CT[0][0] = CT[1][1] = CT[2][2] = 0.0;

  if (qz==0.0)
   CT[0][1]=CT[1][0]=CT[0][2]=CT[2][0]=CT[1][2]=CT[2][1]=0.0;
  else
   { CT[0][1] = -1.0*Sign;         CT[1][0]=-1.0*CT[0][1];
     CT[0][2] = +1.0*Sign*qy/qz;   CT[2][0]=-1.0*CT[0][2];
     CT[1][2] = -1.0*Sign*qx/qz;   CT[2][1]=-1.0*CT[1][2];
   };
}

/***************************************************************/
/* Get the (Fourier-space) 6x6 homogeneous DGF for the         */
/* appropriate region of the multilayer substrate.             */
/***************************************************************/
void LayeredSubstrate::GetScriptG0Twiddle(cdouble Omega, double q2D[2],
                                          double zDest, double zSource,
                                          cdouble ScriptG0Twiddle[6][6],
                                          int nr, double Sign, bool Accumulate)
{
  if (!Accumulate)
   memset( (cdouble *) ScriptG0Twiddle, 0, 6*sizeof(cdouble) );

  // autodetect region if not specified 
  if (nr==-1)
   { nr = GetRegionIndex(zDest);
     if (nr!=GetRegionIndex(zSource) ) 
      Warn("%s:%i: internal inconsistency",__FILE__,__LINE__);
   };

  // autodetect sign if not specified 
  if (Sign==0.0) 
   Sign = ( (zDest > zSource) ? 1.0 : -1.0 );

  cdouble EpsRel=EpsLayer[nr], MuRel=MuLayer[nr];
  cdouble k2=EpsRel*MuRel*Omega*Omega;
  cdouble qz = sqrt(k2 - q2D[0]*q2D[0] - q2D[1]*q2D[1]);

  cdouble GT[3][3], CT[3][3];
  cdouble ExpFac = exp(II*qz*fabs(zDest-zSource));
  GetGC0Twiddle(k2, q2D, qz, Sign, GT, CT);
  cdouble Factor = (qz==0.0 ? 0.0 : -0.5*Omega*ExpFac/qz);
  cdouble EEPrefac = Factor*MuRel*ZVAC;
  cdouble EMPrefac = +0.5*ExpFac;
  cdouble MEPrefac = -0.5*ExpFac;
  cdouble MMPrefac = Factor*EpsRel/ZVAC;
  for(int i=0; i<3; i++)
   for(int j=0; j<3; j++)
    { ScriptG0Twiddle[0+i][0+j] += EEPrefac * GT[i][j];
      ScriptG0Twiddle[0+i][3+j] += EMPrefac * CT[i][j];
      ScriptG0Twiddle[3+i][0+j] += MEPrefac * CT[i][j];
      ScriptG0Twiddle[3+i][3+j] += MMPrefac * GT[i][j];
    };

  if ( nr<NumInterfaces || isinf(zGP) ) return;

  /***************************************************************/
  /* add image contribution if we are in the bottommost layer    */
  /* and a ground plane is present                               */
  /***************************************************************/
  if (Sign>0.0) GetGC0Twiddle(k2, q2D, qz, -1.0, GT, CT);
  ExpFac = exp(II*qz*fabs(zDest + zSource - 2.0*zGP));
  Factor = (qz==0.0 ? 0.0 : -0.5*Omega*ExpFac/qz);
  EEPrefac = Factor*MuRel*ZVAC;
  EMPrefac = +0.5*ExpFac;
  MEPrefac = -0.5*ExpFac;
  MMPrefac = Factor*EpsRel/ZVAC;
  const static double ImageSign[3] = {1.0, 1.0, -1.0};
  for(int i=0; i<3; i++)
   for(int j=0; j<3; j++)
    { ScriptG0Twiddle[0+i][0+j] += ImageSign[j] * EEPrefac * GT[i][j];
      ScriptG0Twiddle[0+i][3+j] -= ImageSign[j] * EMPrefac * CT[i][j];
      ScriptG0Twiddle[3+i][0+j] += ImageSign[j] * MEPrefac * CT[i][j];
      ScriptG0Twiddle[3+i][3+j] -= ImageSign[j] * MMPrefac * GT[i][j];
    };

}

/**********************************************************************/
/* assemble the matrix that operates on the vector of external-field  */
/* Fourier coefficients in all regions to yield the vector of         */
/* surface-current Fourier coefficients on all material interfaces    */
/*                                                                    */
/* note: UpdateCachedEpsMu() should be called before this routine! We */
/* don't call it here to allow the possibility of passing Omega=0     */
/* to evaluate W in the DC limit (but using Epsilon, Mu values at     */
/* the frequency in question)                                         */
/**********************************************************************/
void LayeredSubstrate::ComputeW(cdouble Omega, double q2D[2], HMatrix *W)
{
  W->Zero();
  for(int a=0; a<NumInterfaces; a++)
   {
     int RowOffset = 4*a;
     double za     = zInterface[a];

     for(int b=a-1; b<=a+1; b++)
      { 
        if (b<0 || b>=NumInterfaces) continue;

        int ColOffset = 4*b;
        double zb     = zInterface[b];

        // contributions of surface currents on interface z_b
        // to tangential-field matching equations at interface z_a
        cdouble ScriptG0Twiddle[6][6];
        double Sign=-1.0;
        if ( b==(a-1) )
         GetScriptG0Twiddle(Omega, q2D, za, zb, ScriptG0Twiddle, a);
        else if ( b==(a+1) )
         GetScriptG0Twiddle(Omega, q2D, za, zb, ScriptG0Twiddle, b);
        else // (b==a)
         { GetScriptG0Twiddle(Omega, q2D, za, za, ScriptG0Twiddle, a,   +1.0);
           GetScriptG0Twiddle(Omega, q2D, za, za, ScriptG0Twiddle, a+1, -1.0, true);
           Sign=1.0;
         };

        for(int EH=0; EH<2; EH++)
         for(int KN=0; KN<2; KN++)
          for(int i=0; i<2; i++)
           for(int j=0; j<2; j++)
            W->AddEntry(RowOffset+2*EH+i, ColOffset+2*KN+j, Sign*ScriptG0Twiddle[3*EH+i][3*KN+j]);
      };
   };
  W->LUFactorize();

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void LayeredSubstrate::GetScriptGTwiddle(cdouble Omega, double q2D[2],
                                         double zDest, double zSource,
                                         HMatrix *WMatrix, HMatrix *GTwiddle)
{
  UpdateCachedEpsMu(Omega);

  cdouble OmegaArg = (StaticLimit ? 0.0 : Omega);
  cdouble OmegaFac = ((StaticLimit && Omega!=0.0) ? 1.0/Omega : 1.0);

  /**********************************************************************/
  /**********************************************************************/
  /**********************************************************************/
  if (ForceFreeSpace)
   { 
     cdouble ScriptG0Twiddle[6][6];
     GetScriptG0Twiddle(OmegaArg, q2D, zDest, zSource, ScriptG0Twiddle, 0);
     for(int Mu=0; Mu<6; Mu++)
      for(int Nu=0; Nu<6; Nu++)
       GTwiddle->SetEntry(Mu,Nu,OmegaFac*ScriptG0Twiddle[Mu][Nu]);
     return;
   };
 
  /**********************************************************************/
  /* assemble W matrix ***********************************************/
  /**********************************************************************/
  ComputeW(OmegaArg, q2D,  WMatrix);
  WMatrix->LUInvert();
  
  /**********************************************************************/
  /* assemble RHS vector for each (source point, polarization, orientation) */
  /**********************************************************************/
  int DestRegion  = GetRegionIndex(zDest);
  int pMin        = (DestRegion==0               ? 1 : 0 );
  int pMax        = (DestRegion==NumInterfaces   ? 0 : 1 );
  cdouble ScriptG0TDest[2][6][6];
  for(int p=pMin, nrDest=DestRegion-1+p; p<=pMax; p++, nrDest++)
   GetScriptG0Twiddle(OmegaArg, q2D, zDest, zInterface[nrDest], ScriptG0TDest[p], DestRegion);

  int SourceRegion  = GetRegionIndex(zSource);
  int qMin          = (SourceRegion==0             ? 1 : 0 );
  int qMax          = (SourceRegion==NumInterfaces ? 0 : 1 );
  cdouble ScriptG0TSource[2][6][6];
  for(int q=qMin, nrSource=SourceRegion-1+q; q<=qMax; q++, nrSource++)
   GetScriptG0Twiddle(OmegaArg, q2D, zInterface[nrSource], zSource, ScriptG0TSource[q], SourceRegion);

  /**********************************************************************/
  /**********************************************************************/
  /**********************************************************************/
  GTwiddle->Zero();
  for(int p=pMin, nrDest=DestRegion-1+p; p<=pMax; p++, nrDest++)
   for(int q=qMin, nrSource=SourceRegion-1+q; q<=qMax; q++, nrSource++)
    { 
      double Sign = (p==q) ? 1.0 : -1.0;
      Sign*=-1.0; // EXPLAIN ME 
      int RowOffset = 4*nrDest, ColOffset = 4*nrSource;
      for(int i=0; i<6; i++)
       for(int j=0; j<6; j++)
        for(int k=0; k<4; k++)
         for(int l=0; l<4; l++)
          { int kk = (k>=2 ? k+1 : k);
            int ll = (l>=2 ? l+1 : l);
            GTwiddle->AddEntry(i,j, Sign*OmegaFac
                                        *ScriptG0TDest[p][i][kk]
                                        *WMatrix->GetEntry(RowOffset+k, ColOffset+l)
                                        *ScriptG0TSource[q][ll][j]
                              );
         };
    };
}

void LayeredSubstrate::GetScriptGTwiddle(cdouble Omega, double qx, double qy,
                                         double zDest, double zSource,
                                         HMatrix *WMatrix, HMatrix *GTwiddle)
{ double q2D[2];
  q2D[0]=qx;
  q2D[1]=qy;
  GetScriptGTwiddle(Omega, q2D, zDest, zSource, WMatrix, GTwiddle);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void LayeredSubstrate::Getg0112(cdouble Omega, double qMag,
                                double zDest, double zSource,
                                HMatrix *WMatrix, HMatrix *g0112[4])
{
  // qx=qMag, qy=0
  HMatrix GT10(6,6,LHM_COMPLEX);
  GetScriptGTwiddle(Omega, qMag,       0.0,        zDest, zSource, WMatrix, &GT10);

  // qx=qMag/sqrt[2], qy=qMag/sqrt[2]
  HMatrix GT11(6,6,LHM_COMPLEX);
  GetScriptGTwiddle(Omega, qMag/SQRT2, qMag/SQRT2, zDest, zSource, WMatrix, &GT11);

  // qx=0, qy=qMag
  HMatrix GT01(6,6,LHM_COMPLEX);
  GetScriptGTwiddle(Omega, 0.0,        qMag,       zDest, zSource, WMatrix, &GT01);

  HMatrix *g0 = g0112[0], *g1x=g0112[1], *g1y=g0112[2], *g2=g0112[3];
  g0->Zero();
  g1x->Zero();
  g1y->Zero();
  g2->Zero();
  
  for(int p=0; p<2; p++)
   for(int q=0; q<2; q++)
    { 
      int ROfs = 3*p;
      int COfs = 3*q;

      g0->SetEntry(ROfs+0,COfs+0,GT01.GetEntry(ROfs+0,COfs+0));
      g0->SetEntry(ROfs+0,COfs+1,GT01.GetEntry(ROfs+0,COfs+1));
      g0->SetEntry(ROfs+1,COfs+0,GT01.GetEntry(ROfs+1,COfs+0));
      g0->SetEntry(ROfs+1,COfs+1,GT10.GetEntry(ROfs+1,COfs+1));

      g0->SetEntry(ROfs+0,COfs+2,GT10.GetEntry(ROfs+0,COfs+2) + GT01.GetEntry(ROfs+0, COfs+2) - SQRT2*GT11.GetEntry(ROfs+0, COfs+2));
      g0->SetEntry(ROfs+2,COfs+0,GT10.GetEntry(ROfs+2,COfs+0) + GT01.GetEntry(ROfs+2, COfs+0) - SQRT2*GT11.GetEntry(ROfs+2, COfs+0));
      g0->SetEntry(ROfs+1,COfs+2,GT10.GetEntry(ROfs+1,COfs+2) + GT01.GetEntry(ROfs+1, COfs+2) - SQRT2*GT11.GetEntry(ROfs+1, COfs+2));
      g0->SetEntry(ROfs+2,COfs+1,GT10.GetEntry(ROfs+2,COfs+1) + GT01.GetEntry(ROfs+2, COfs+1) - SQRT2*GT11.GetEntry(ROfs+2, COfs+1));

      g0->SetEntry(ROfs+2,COfs+2,GT01.GetEntry(ROfs+2,COfs+2));

      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      g1x->SetEntry(ROfs+0,COfs+2, GT10.GetEntry(ROfs+0, COfs+2) - g0->GetEntry(ROfs+0, COfs+2) );
      g1x->SetEntry(ROfs+2,COfs+0, GT10.GetEntry(ROfs+2, COfs+0) - g0->GetEntry(ROfs+2, COfs+0) );
      g1x->SetEntry(ROfs+1,COfs+2, GT10.GetEntry(ROfs+1, COfs+2) - g0->GetEntry(ROfs+1, COfs+2) );
      g1x->SetEntry(ROfs+2,COfs+1, GT10.GetEntry(ROfs+2, COfs+1) - g0->GetEntry(ROfs+2, COfs+1) );

      g1y->SetEntry(ROfs+0,COfs+2, GT01.GetEntry(ROfs+0, COfs+2) - g0->GetEntry(ROfs+0, COfs+2) );
      g1y->SetEntry(ROfs+2,COfs+0, GT01.GetEntry(ROfs+2, COfs+0) - g0->GetEntry(ROfs+2, COfs+0) );
      g1y->SetEntry(ROfs+1,COfs+2, GT01.GetEntry(ROfs+1, COfs+2) - g0->GetEntry(ROfs+1, COfs+2) );
      g1y->SetEntry(ROfs+2,COfs+1, GT01.GetEntry(ROfs+2, COfs+1) - g0->GetEntry(ROfs+2, COfs+1) );

      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      g2->SetEntry(ROfs+0,COfs+0,       GT10.GetEntry(ROfs+0, COfs+0) - g0->GetEntry(ROfs+0, COfs+0));
      g2->SetEntry(ROfs+0,COfs+1, 2.0*( GT11.GetEntry(ROfs+0, COfs+1) - g0->GetEntry(ROfs+0, COfs+1)));
      g2->SetEntry(ROfs+1,COfs+0, 2.0*( GT11.GetEntry(ROfs+1, COfs+0) - g0->GetEntry(ROfs+1, COfs+0)));
      g2->SetEntry(ROfs+1,COfs+1,       GT01.GetEntry(ROfs+1, COfs+1) - g0->GetEntry(ROfs+1, COfs+1));

    };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct IntegrandData 
 {
   LayeredSubstrate *Substrate;
   HMatrix *XMatrix;
   HMatrix *WMatrix;
   cdouble Omega;
   double rkLayer2;
   double q0;
   FILE *LogFile;
   bool Propagating;
   int nCalls;

 } IntegrandData;

int SubstrateDGFIntegrand_SC(unsigned ndim, const double *uVector,
                             void *UserData, unsigned fdim, double *fval)
{
  (void )ndim; //unused

  IntegrandData *Data         = (IntegrandData *)UserData;
  Data->nCalls++;

  LayeredSubstrate *Substrate = Data->Substrate;
  HMatrix *XMatrix            = Data->XMatrix;
  HMatrix *WMatrix            = Data->WMatrix;
  cdouble Omega               = Data->Omega;
  double rkLayer2             = Data->rkLayer2;
  FILE *LogFile               = Data->LogFile;
  double q0                   = Data->q0;
  bool Propagating            = Data->Propagating;

  int EntryOnly               = Substrate->EntryOnly;
  bool EEOnly                 = Substrate->EEOnly;
  bool XYOnly                 = Substrate->XYOnly;
  int NumInterfaces           = Substrate->NumInterfaces;
  double *zInterface          = Substrate->zInterface;
  cdouble *EpsLayer           = Substrate->EpsLayer;

  cdouble *Integrand = (cdouble *)fval;
  memset(Integrand, 0, fdim*sizeof(double));

  double Jacobian = 1.0, u=uVector[0], qMag;
#define UMIN 1.0e-3
  if (Propagating)
   { 
     if ( u<=UMIN )
      u=UMIN;
     if ( u>=(1.0-UMIN) )
      u=1.0-UMIN;
     double qz = u*q0;
     Jacobian  = q0*qz/(2.0*M_PI);
     qMag      = sqrt(rkLayer2 - qz*qz);
   }
  else // evanescent
   { 
     if (u>1.0-UMIN)
      return 0;    // integrand vanishes at qMag=iqz=infinity
     if (u<=UMIN)
      u=UMIN;
     double Denom = 1.0 / (1.0-u);
     qMag         = q0 + u*Denom;
     Jacobian     = qMag*Denom*Denom/(2.0*M_PI);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int nx=0, ni=0; nx<XMatrix->NR; nx++)
   { 
     double Rho[2];
     Rho[0]          = XMatrix->GetEntryD(nx,0) - XMatrix->GetEntryD(nx,3);
     Rho[1]          = XMatrix->GetEntryD(nx,1) - XMatrix->GetEntryD(nx,4);
     double RhoMag   = sqrt( Rho[0]*Rho[0] + Rho[1]*Rho[1] );
     double CosTheta = (RhoMag==0.0) ? 0.0 : Rho[0] / RhoMag;
     double SinTheta = (RhoMag==0.0) ? 0.0 : Rho[1] / RhoMag;
     double zDest    = XMatrix->GetEntryD(nx,2);
     double zSource  = XMatrix->GetEntryD(nx,5);

     // fetch scalar quantities
     HMatrix g0(6, 6, LHM_COMPLEX);
     HMatrix g1x(6, 6, LHM_COMPLEX);
     HMatrix g1y(6, 6, LHM_COMPLEX);
     HMatrix g2(6, 6, LHM_COMPLEX);
     HMatrix *g0112[4];
     g0112[0]=&g0;
     g0112[1]=&g1x;
     g0112[2]=&g1y;
     g0112[3]=&g2;
     Substrate->Getg0112(Omega, qMag, zDest, zSource, WMatrix, g0112);

     // fetch bessel functions
     cdouble J[3], J1oqRho;
     double qRho = qMag*RhoMag;
     if (qRho<1.0e-4) // series expansions for small arguments
      { double qRho2=qRho*qRho;
        J[0]    = 1.0 - qRho2/4.0;
        J1oqRho = 0.5*(1.0-qRho2/8.0);
        J[1]    = J1oqRho*qRho;
        J[2]    = 0.125*qRho2*(1.0-qRho2/12.0);
      }
     else if (qRho>1.0e3) // asymptotic forms for large argument
      { double JPreFac = sqrt(2.0/(M_PI*qRho));
        J[0] = JPreFac * cos(qRho - 0.25*M_PI);
        J[1] = JPreFac * cos(qRho - 0.75*M_PI);
        J1oqRho = J[1]/qRho;
        J[2] = JPreFac * cos(qRho - 1.25*M_PI);
      }
     else
      { double Workspace[12];
        AmosBessel('J', qRho, 0.0, 3, false, J, Workspace);
        J1oqRho = J[1]/qRho;
      };
     cdouble J0, J1[2], J2[2][2];
     J0       = J[0];
     J1[0]    = II*CosTheta*J[1];
     J1[1]    = II*SinTheta*J[1];
     J2[0][0] = -1.0*CosTheta*CosTheta*J[2] + J1oqRho;
     J2[0][1] = -1.0*CosTheta*SinTheta*J[2];
     J2[1][0] = -1.0*SinTheta*CosTheta*J[2];
     J2[1][1] = -1.0*SinTheta*SinTheta*J[2] + J1oqRho;

     // check for zDest, zSource both on a layer boundary
     // note that this calculation assumes Rho=(x,0))
     cdouble Gxx0=0.0, Gyy0=0.0;
     bool OnInterface=false;
     if ( EqualFloat(zDest, zSource) )
      for(int nz=0; nz<NumInterfaces && !OnInterface; nz++)
       if ( fabs(zDest-zInterface[nz]) < 1.0e-6 )
        { cdouble EpsA = EpsLayer[nz], EpsB=EpsLayer[nz+1];
          cdouble Factor = (-0.5*II*ZVAC/Omega)*(EpsA-EpsB)/(EpsA+EpsB);
          Gxx0 = Factor*qMag*(J0 - J1oqRho);
          Gyy0 = Factor*qMag*J1oqRho;
          OnInterface=true;
        };
     
     // assemble GTwiddle
     cdouble GTwiddle[6][6];
     memset( (cdouble *)GTwiddle, 0, 36*sizeof(cdouble));
     int pqMax = EEOnly ? 1 : 2;
     int ijMax = XYOnly ? 2 : 3;
     for(int p=0; p<pqMax; p++)
      for(int q=0; q<pqMax; q++)
       for(int i=0; i<ijMax; i++)
        for(int j=0; j<ijMax; j++)
         { 
           int Mu=3*p+i, Nu=3*q+j;
           if (EntryOnly!=-1 && (EntryOnly!=(6*Mu+Nu)) )
            continue;
            
           cdouble g0MuNu  = g0. GetEntry(Mu,Nu);
           cdouble g1xMuNu = g1x.GetEntry(Mu,Nu);
           cdouble g1yMuNu = g1y.GetEntry(Mu,Nu);
           cdouble g2MuNu  = g2. GetEntry(Mu,Nu);

           GTwiddle[Mu][Nu]   = g0MuNu * J0;

           if (i<2 && j<2)
            GTwiddle[Mu][Nu] += g2MuNu * J2[i][j];
           else if ( (i<2 && j==2) || (i==2 && j<2) )
            GTwiddle[Mu][Nu] += (g1xMuNu*J1[0] + g1yMuNu*J1[1]);
         };
     if (OnInterface)
      { GTwiddle[0][0] -= Gxx0;
        GTwiddle[1][1] -= Gyy0;
      };

     // stamp into integrand vector
     for(int Mu=0; Mu<6; Mu++)
      for(int Nu=0; Nu<6; Nu++)
       Integrand[ni++] = Jacobian * GTwiddle[Mu][Nu];

     if (LogFile)
      { fprintf(LogFile,"%e ",qMag);
        fprintf(LogFile,"%e %e %e %e ",Rho[0],Rho[1],zDest,zSource);
        fprintf(LogFile,"%e %e ",real(J[0]),imag(J[0]));
        fprintf(LogFile,"%e %e ",real(J1oqRho),imag(J1oqRho));
        fprintf(LogFile,"%e %e ",real(J[2]),imag(J[2]));
        fprintVec(LogFile,(cdouble *)GTwiddle, 36);
        fprintf(LogFile,"%e %e %e %e ",real(Gxx0), imag(Gxx0), real(Gyy0), imag(Gyy0));
        fprintf(LogFile,"\n");
        fflush(LogFile);
      };

   }; // for(int nx=0, ni=0; nx<XMatrix->NR; nx++)
  return 0;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void LayeredSubstrate::GetSubstrateDGF_SurfaceCurrent(cdouble Omega,
                                                      HMatrix *XMatrix,
                                                      HMatrix *GMatrix)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int NX = XMatrix->NR;
  int IDim = 36*NX;

  HMatrix *WMatrix  = new HMatrix(4*NumInterfaces, 4*NumInterfaces, LHM_COMPLEX);
  HMatrix *kLayer   = new HMatrix(NumLayers, 2, LHM_REAL);
  cdouble *Integral = new cdouble[IDim];
  cdouble *Error    = new cdouble[IDim];

  UpdateCachedEpsMu(Omega);
 
  // compute wavevectors in all layers and sort by real part
  for(int n=0; n<NumLayers; n++)
   { cdouble k=sqrt(EpsLayer[n]*MuLayer[n])*Omega;
     kLayer->SetEntry(n, 0, real(k));
     kLayer->SetEntry(n, 1, imag(k));
   };
  kLayer->Sort(0);

  double Lower=0.0, Upper=1.0;
  IntegrandData MyData, *Data=&MyData;
  Data->Substrate  = this;
  Data->XMatrix    = XMatrix;
  Data->WMatrix    = WMatrix;
  Data->Omega      = Omega;
  Data->LogFile    = 0;

  GMatrix->Zero();

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int n=0; n<NumLayers; n++)
   {
     if (LayerOnly!=-1 && LayerOnly!=n)
      continue;

     double knm1       = (n==0) ? 0.0 : kLayer->GetEntryD(n-1,0);
     double kn         = kLayer->GetEntryD(n,0);
     Data->rkLayer2    = kn*kn;
     Data->q0          = sqrt(kn*kn - knm1*knm1);
     Data->nCalls      = 0;
     Data->Propagating = true;
     Data->LogFile     = vfopen("/tmp/Propagating%i%s.log","w",n+1,StaticLimit ? "_SL" : "");
     pcubature(2*IDim, SubstrateDGFIntegrand_SC, (void *)Data, 1,
               &Lower, &Upper, qMaxEval, qAbsTol, qRelTol, 
               ERROR_PAIRED, (double *)Integral, (double *)Error);
     if (Data->LogFile) fclose(Data->LogFile);
     Data->LogFile=0;
     Log("Integral %i [%g,%g]: %i calls, GEExx,GEEzz,GMMxx,GMMzz=%s,%s",n+1,knm1,kn,Data->nCalls,
          CD2S(Integral[0]),CD2S(Integral[14]), CD2S(Integral[21]),CD2S(Integral[35]));

     for(int nx=0; nx<NX; nx++)
      for(int ng=0; ng<36; ng++)
       GMatrix->AddEntry(nx, ng, Integral[36*nx + ng]);
   };

  if (LayerOnly==-1 || LayerOnly==NumLayers)
   { Data->nCalls      = 0;
     Data->Propagating = false;
     Data->LogFile     = vfopen("/tmp/Evanescent%s.log","w", StaticLimit ? "_SL" : "");
     Data->q0          = kLayer->GetEntryD(NumLayers-1,0);
     pcubature(2*IDim, SubstrateDGFIntegrand_SC, (void *)Data, 1,
               &Lower, &Upper, qMaxEval, qAbsTol, qRelTol,
               ERROR_PAIRED, (double *)Integral, (double *)Error);
     Log("Integral %i [%g,inf]: %i calls, GEExx,GEEzz,GMMxx,GMMzz=%s,%s",NumLayers,
          Data->q0, Data->nCalls,CD2S(Integral[0]),CD2S(Integral[14]), CD2S(Integral[21]),CD2S(Integral[35]));
     if (Data->LogFile) fclose(Data->LogFile);
     Data->LogFile=0;

     for(int nx=0; nx<NX; nx++)
      for(int ng=0; ng<36; ng++)
       GMatrix->AddEntry(nx, ng, Integral[36*nx + ng]);
   };
  
  delete[] Error;
  delete[] Integral;
  delete kLayer;
  delete WMatrix;
  
}
