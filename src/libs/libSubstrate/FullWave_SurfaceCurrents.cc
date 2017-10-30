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

const char *TimeNames[]={"G0", "BESSEL", "W", "SOLVE", "STAMP"};

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
   memset( (cdouble *) ScriptG0Twiddle, 0, 36*sizeof(cdouble) );

  // autodetect region if not specified 
  if (nr==-1)
   { nr = GetRegionIndex(zDest);
     if (nr!=GetRegionIndex(zSource) ) 
      Warn("%s:%i: internal inconsistency",__FILE__,__LINE__);
   };

  // autodetect sign if not specified 
  if (Sign==0.0) 
   Sign = ( (zDest >= zSource) ? 1.0 : -1.0 );

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
  // only need to refetch if Sign was +1 before
  if (Sign>0.0) GetGC0Twiddle(k2, q2D, qz, -1.0, GT, CT);
  ExpFac = exp(II*qz*fabs(zDest + zSource - 2.0*zGP));
  Factor = (qz==0.0 ? 0.0 : -0.5*Omega*ExpFac/qz);
  EEPrefac = Factor*MuRel*ZVAC;
  EMPrefac = +0.5*ExpFac;
  MEPrefac = -0.5*ExpFac;
  MMPrefac = Factor*EpsRel/ZVAC;
  const static double ImageSign[3] = {-1.0, -1.0, +1.0};
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
/**********************************************************************/
void LayeredSubstrate::ComputeW(cdouble Omega, double q2D[2], HMatrix *W)
{
  UpdateCachedEpsMu(Omega);

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
  if (LogLevel >= LIBSUBSTRATE_VERBOSE) 
   Log("LU factorizing...");
  W->LUFactorize();
}

/***************************************************************/
/* Get the (Fourier-space) 6x6 *inhomogeneous* DGF, i.e. the   */
/* fields due to point sources in the presence of the substrate*/
/***************************************************************/
void LayeredSubstrate::GetScriptGTwiddle(cdouble Omega, double q2D[2],
                                         double zDest, double zSource,
                                         HMatrix *RTwiddle,
                                         HMatrix *WMatrix,
                                         HMatrix *STwiddle,
                                         HMatrix *GTwiddle)
{
  UpdateCachedEpsMu(Omega);

  /**********************************************************************/
  /**********************************************************************/
  /**********************************************************************/
  if (ForceFreeSpace)
   { 
     cdouble ScriptG0Twiddle[6][6];
     GetScriptG0Twiddle(Omega, q2D, zDest, zSource, ScriptG0Twiddle);
     for(int Mu=0; Mu<6; Mu++)
      for(int Nu=0; Nu<6; Nu++)
       GTwiddle->SetEntry(Mu,Nu,ScriptG0Twiddle[Mu][Nu]);
     return;
   };
  
  /**********************************************************************/
  /* assemble RHS vector for each (source point, polarization, orientation) */
  /**********************************************************************/
  double TT=Secs();
  int DestRegion  = GetRegionIndex(zDest);
  int pMin        = (DestRegion==0               ? 1 : 0 );
  int pMax        = (DestRegion==NumInterfaces   ? 0 : 1 );
  cdouble ScriptG0TDest[2][6][6];
  for(int p=pMin, nrDest=DestRegion-1+p; p<=pMax; p++, nrDest++)
   GetScriptG0Twiddle(Omega, q2D, zDest, zInterface[nrDest], ScriptG0TDest[p], DestRegion);

  int SourceRegion  = GetRegionIndex(zSource);
  int qMin          = (SourceRegion==0             ? 1 : 0 );
  int qMax          = (SourceRegion==NumInterfaces ? 0 : 1 );
  cdouble ScriptG0TSource[2][6][6];
  for(int q=qMin, nrSource=SourceRegion-1+q; q<=qMax; q++, nrSource++)
   GetScriptG0Twiddle(Omega, q2D, zInterface[nrSource], zSource, ScriptG0TSource[q], SourceRegion);

  Times[G0TIME] += Secs() - TT;
 
  /**********************************************************************/
  /* assemble W matrix **************************************************/
  /**********************************************************************/
  TT=Secs();
  ComputeW(Omega, q2D,  WMatrix);
  Times[WTIME]+=Secs()-TT;

  /**********************************************************************/
  /**********************************************************************/
  /**********************************************************************/
  TT=Secs();
  RTwiddle->Zero();
  STwiddle->Zero();
  double Sign[2]={-1.0, 1.0};
  for(int p=pMin, nrDest=DestRegion-1+p; p<=pMax; p++, nrDest++)
   for(int i=0; i<6; i++)
    { RTwiddle->SetEntry(i, 4*nrDest+0, Sign[p]*ScriptG0TDest[p][i][0]);
      RTwiddle->SetEntry(i, 4*nrDest+1, Sign[p]*ScriptG0TDest[p][i][1]);
      RTwiddle->SetEntry(i, 4*nrDest+2, Sign[p]*ScriptG0TDest[p][i][3]);
      RTwiddle->SetEntry(i, 4*nrDest+3, Sign[p]*ScriptG0TDest[p][i][4]);
    };
  for(int q=qMin, nrSource=SourceRegion-1+q; q<=qMax; q++, nrSource++)
   for(int j=0; j<6; j++)
    { STwiddle->SetEntry(4*nrSource+0, j, -1.0*Sign[q]*ScriptG0TSource[q][0][j]);
      STwiddle->SetEntry(4*nrSource+1, j, -1.0*Sign[q]*ScriptG0TSource[q][1][j]);
      STwiddle->SetEntry(4*nrSource+2, j, -1.0*Sign[q]*ScriptG0TSource[q][3][j]);
      STwiddle->SetEntry(4*nrSource+3, j, -1.0*Sign[q]*ScriptG0TSource[q][4][j]);
    };
  WMatrix->LUSolve(STwiddle);
  RTwiddle->Multiply(STwiddle, GTwiddle);
  Times[SOLVETIME]+=Secs()-TT;
}

void LayeredSubstrate::GetScriptGTwiddle(cdouble Omega, double qx, double qy,
                                         double zDest, double zSource,
                                         HMatrix *RTwiddle, HMatrix *WMatrix,
                                         HMatrix *STwiddle, HMatrix *GTwiddle)
{ double q2D[2];
  q2D[0]=qx;
  q2D[1]=qy;
  GetScriptGTwiddle(Omega, q2D, zDest, zSource, RTwiddle, WMatrix, STwiddle, GTwiddle);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void LayeredSubstrate::Getg0112(cdouble Omega, double qMag,
                                double zDest, double zSource,
                                HMatrix *RTwiddle, HMatrix *WMatrix, 
                                HMatrix *STwiddle, HMatrix *g0112[4])
{
  // qx=qMag, qy=0
  cdouble GTBuf[3*36];
  HMatrix GT10(6,6,LHM_COMPLEX, GTBuf+0);
  GetScriptGTwiddle(Omega, qMag,       0.0,        zDest, zSource, RTwiddle, WMatrix, STwiddle, &GT10);

  // qx=qMag/sqrt[2], qy=qMag/sqrt[2]
  HMatrix GT11(6,6,LHM_COMPLEX, GTBuf+36);
  GetScriptGTwiddle(Omega, qMag/SQRT2, qMag/SQRT2, zDest, zSource, RTwiddle, WMatrix, STwiddle, &GT11);

  // qx=0, qy=qMag
  HMatrix GT01(6,6,LHM_COMPLEX, GTBuf+72);
  GetScriptGTwiddle(Omega, 0.0,        qMag,       zDest, zSource, RTwiddle, WMatrix, STwiddle, &GT01);

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
/* fetch Bessel-function factors                               */
/* J[0,1,2,3] = J0, J1, J1/qRho, J_2                           */
/* dJ = dJ/dRho                                                */
/***************************************************************/
void GetJFactors(double q, double Rho, cdouble *J, cdouble *dJ)
{
  double qRho = q*Rho;
  cdouble J0, J1oqRho, J1, J2;
  cdouble dJ0, dJ1oqRho, dJ1, dJ2;
  if (qRho<1.0e-4) // series expansions for small arguments
   { double qRho2=qRho*qRho;
     J0    = 1.0 - qRho2/4.0;
     J1oqRho = 0.5*(1.0-qRho2/8.0);
     J1    = J1oqRho*qRho;
     J2    = 0.125*qRho2*(1.0-qRho2/12.0);
     if (dJ)
      { dJ0      = -0.5*q*qRho;
        dJ1oqRho = -0.125*q*qRho;
        dJ1      = q*J1oqRho + qRho*dJ1oqRho;
        dJ2      = q*qRho/4.0 - qRho2*qRho*Rho/24.0;
      };
   }
  else if (qRho>1.0e3) // asymptotic forms for large argument
   { double JPreFac = sqrt(2.0/(M_PI*qRho));
     J0   = JPreFac * cos(qRho - 0.25*M_PI);
     J1   = JPreFac * cos(qRho - 0.75*M_PI);
     J1oqRho = J1/qRho;
     J2   = JPreFac * cos(qRho - 1.25*M_PI);
     if (dJ)
      { dJ0      = -0.5*J0/Rho - q*JPreFac*sin(qRho-0.25*M_PI);
        dJ1      = -0.5*J1/Rho - q*JPreFac*sin(qRho-0.75*M_PI);
        dJ1oqRho =    dJ1/qRho - J1oqRho/Rho;
        dJ2      = -0.5*J2/Rho - q*JPreFac*sin(qRho-1.25*M_PI);
      };
   }
  else
   { double Workspace[12];
     AmosBessel('J', qRho, 0.0, 3, false, J, Workspace);
     J0=J[0];
     J1=J[1];
     J1oqRho = J[1]/qRho;
     J2=J[2];
     if (dJ)
      { dJ0 = -q*J1;
        dJ1 = 0.5*q*(J0 - J2);
        dJ1oqRho = dJ1/qRho - J1oqRho/Rho;
        dJ2 = q*(J1 - 2.0*J2/qRho);
      };
   };

  J[0] = J0;
  J[1] = J1;
  J[2] = J1oqRho;
  J[3] = J2;

  if (dJ==0) return;

  dJ[0] = dJ0;
  dJ[1] = dJ1;
  dJ[2] = dJ1oqRho;
  dJ[3] = dJ2;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct qFunctionSCData
 {
   HMatrix *XMatrix;
   HMatrix *RTwiddle;
   HMatrix *WMatrix;
   HMatrix *STwiddle;
   FILE *byqFile;
 } qFunctionSCData;

void qFunctionSC(LayeredSubstrate *Substrate,
                 double q2D[2], cdouble Omega,
                 void *UserData, cdouble *Integrand)
{
  qFunctionSCData *Data = (qFunctionSCData *)UserData;

  HMatrix *XMatrix    = Data->XMatrix;
  HMatrix *RTwiddle   = Data->RTwiddle;
  HMatrix *WMatrix    = Data->WMatrix;
  HMatrix *STwiddle   = Data->STwiddle;
  FILE *byqFile       = Data->byqFile;

  int EntryOnly       = Substrate->EntryOnly;
  bool EEOnly         = Substrate->EEOnly;
  bool XYOnly         = Substrate->XYOnly;
  int NumInterfaces   = Substrate->NumInterfaces;
  double *zInterface  = Substrate->zInterface;
  cdouble *EpsLayer   = Substrate->EpsLayer;

  double qMag = q2D[0];

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
     cdouble gBuf[4*36];
     HMatrix g0( 6, 6, LHM_COMPLEX, gBuf +  0  );
     HMatrix g1x(6, 6, LHM_COMPLEX, gBuf + 36  );
     HMatrix g1y(6, 6, LHM_COMPLEX, gBuf + 72  );
     HMatrix g2( 6, 6, LHM_COMPLEX, gBuf + 108 );
     HMatrix *g0112[4];
     g0112[0]=&g0;
     g0112[1]=&g1x;
     g0112[2]=&g1y;
     g0112[3]=&g2;
     Substrate->Getg0112(Omega, qMag, zDest, zSource, RTwiddle, WMatrix, STwiddle, g0112);

     // fetch bessel functions
     double TT=Secs();
bool NeedRhoDerivatives=false;
     cdouble J[4], dJ[4];
     GetJFactors(qMag, RhoMag, J, NeedRhoDerivatives ? dJ : 0);
     Substrate->Times[BESSELTIME] += (Secs()-TT);
  
cdouble J0, J1oqRho, J1[2], J2[2][2]; 
  J0       = J[0];
  J1[0]    = II*CosTheta*J[1];
  J1[1]    = II*SinTheta*J[1];
  J1oqRho  = J[2];
  J2[0][0] = -1.0*CosTheta*CosTheta*J[3] + J1oqRho;
  J2[0][1] = -1.0*CosTheta*SinTheta*J[3];
  J2[1][0] = -1.0*SinTheta*CosTheta*J[3];
  J2[1][1] = -1.0*SinTheta*SinTheta*J[3] + J1oqRho;
#if 0
dJ0, dJ1oqRho, dJ1[2], dJ2[2][2];
  dJData->J0       = dJ[0];
  dJData->J1oqRho  = dJ1oqRho;
  dJData->J1[0]    = II*CosTheta*dJ[1];
  dJData->J1[1]    = II*SinTheta*dJ[1];
  dJData->J2[0][0] = -1.0*CosTheta*CosTheta*dJ[2] + dJ1oqRho;
  dJData->J2[0][1] = -1.0*CosTheta*SinTheta*dJ[2];
  dJData->J2[1][0] = -1.0*SinTheta*CosTheta*dJ[2];
  dJData->J2[1][1] = -1.0*SinTheta*SinTheta*dJ[2] + dJ1oqRho;
#endif

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
     int pqMax = (EEOnly ? 1 : 2);
     int ijMax = (XYOnly ? 2 : 3);
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
       Integrand[ni++] = GTwiddle[Mu][Nu];

     Substrate->Times[STAMPTIME]+=Secs() - TT;

     if (byqFile)
      { fprintf(byqFile,"%e ",qMag);
        fprintf(byqFile,"%e %e %e %e ",Rho[0],Rho[1],zDest,zSource);
 //      fprintf(byqFile,"%e %e ",real(J[0]),imag(J[0]));
 //       fprintf(byqFile,"%e %e ",real(J1oqRho),imag(J1oqRho));
 //       fprintf(byqFile,"%e %e ",real(J[2]),imag(J[2]));
        fprintVec(byqFile,(cdouble *)GTwiddle, 36);
        fprintf(byqFile,"%e %e %e %e ",real(Gxx0), imag(Gxx0), real(Gyy0), imag(Gyy0));
        fprintf(byqFile,"\n");
        fflush(byqFile);
      };

   }; // for(int nx=0, ni=0; nx<XMatrix->NR; nx++)

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
  int FDim = 36*NX;

  UpdateCachedEpsMu(Omega);

  HMatrix *RTwiddle  = new HMatrix(6,               4*NumInterfaces, LHM_COMPLEX);
  HMatrix *WMatrix   = new HMatrix(4*NumInterfaces, 4*NumInterfaces, LHM_COMPLEX);
  HMatrix *STwiddle  = new HMatrix(4*NumInterfaces, 6,               LHM_COMPLEX);

  qFunctionSCData MyData, *Data=&MyData;
  Data->XMatrix    = XMatrix;
  Data->RTwiddle   = RTwiddle;
  Data->WMatrix    = WMatrix;
  Data->STwiddle   = STwiddle;
  Data->byqFile    = 0;

  qIntegrate(Omega, qFunctionSC, (void *)Data, GMatrix->ZM, FDim);
  
  delete RTwiddle;
  delete WMatrix;
  delete STwiddle;
}
