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

#define EE_XX 0
#define EE_XY 1
#define EE_XZ 2
#define ME_XX 3
#define ME_XY 4
#define ME_XZ 5
#define EE_YX 6
#define EE_YY 7
#define EE_YZ 8
#define ME_YX 9
#define ME_YY 10
#define ME_YZ 11
#define EE_ZX 12
#define EE_ZY 13
#define EE_ZZ 14
#define ME_ZX 15
#define ME_ZY 16
#define ME_ZZ 17
#define EM_XX 18
#define EM_XY 19
#define EM_XZ 20
#define MM_XX 21
#define MM_XY 22
#define MM_XZ 23
#define EM_YX 24
#define EM_YY 25
#define EM_YZ 26
#define MM_YX 27
#define MM_YY 28
#define MM_YZ 29
#define EM_ZX 30
#define EM_ZY 31
#define EM_ZZ 32
#define MM_ZX 33
#define MM_ZY 34
#define MM_ZZ 35

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
  GT[0][2] = GT[2][0] = -1.0*Sign*qx*qz/k2;
  GT[1][2] = GT[2][1] = -1.0*Sign*qy*qz/k2;

  CT[0][0] = CT[1][1] = CT[2][2] = 0.0;

  if (qz==0.0)
   CT[0][1]=CT[1][0]=CT[0][2]=CT[2][0]=CT[1][2]=CT[2][1]=0.0;
  else
   { CT[0][1] = -1.0*Sign;    CT[1][0]=-1.0*CT[0][1];
     CT[0][2] = +1.0*qy/qz;   CT[2][0]=-1.0*CT[0][2];
     CT[1][2] = -1.0*qx/qz;   CT[2][1]=-1.0*CT[1][2];
   };
}

/***************************************************************/
/* Get the (Fourier-space) 6x6 homogeneous DGF for the         */
/* medium occupying substrate layer #nl.                       */
/***************************************************************/
void LayeredSubstrate::GetScriptG0Twiddle(cdouble Omega, double q2D[2],
                                          double zDest, double zSource,
                                          cdouble ScriptG0Twiddle[6][6],
                                          int nl, double Sign, bool Accumulate,
                                          bool dzDest, bool dzSource)
{
  if (!Accumulate)
   memset( (cdouble *) ScriptG0Twiddle, 0, 36*sizeof(cdouble) );

  // autodetect region if not specified 
  if (nl==-1)
   { nl = GetLayerIndex(zDest);
     if (nl!=GetLayerIndex(zSource) ) 
      Warn("%s:%i: internal inconsistency",__FILE__,__LINE__);
   };

  // autodetect sign if not specified 
  if (Sign==0.0) 
   Sign = ( (zDest >= zSource) ? 1.0 : -1.0 );

  cdouble EpsRel=EpsLayer[nl], MuRel=MuLayer[nl];
  cdouble k2=EpsRel*MuRel*Omega*Omega;
  cdouble qz = sqrt(k2 - q2D[0]*q2D[0] - q2D[1]*q2D[1]);

  cdouble GT[3][3], CT[3][3];
  cdouble ExpFac = exp(II*qz*fabs(zDest-zSource));
  if (dzDest)
   ExpFac *= +1.0*Sign*II*qz;
  if (dzSource)
   ExpFac *= -1.0*Sign*II*qz;
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
  if ( nl<NumInterfaces || isinf(zGP) ) return;

  /***************************************************************/
  /* add image contribution if we are in the bottommost layer    */
  /* and a ground plane is present                               */
  /***************************************************************/
  // only need to refetch if Sign was +1 before
  if (Sign>0.0) GetGC0Twiddle(k2, q2D, qz, -1.0, GT, CT);
  ExpFac = exp(II*qz*fabs(zDest + zSource - 2.0*zGP));
  if (dzDest)
   ExpFac *= +1.0*II*qz;
  if (dzSource)
   ExpFac *= +1.0*II*qz;
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
/* (Fourier components of) the fields due to point sources in  */
/* the presence of the substrate.                              */
/***************************************************************/
void LayeredSubstrate::GetScriptGTwiddle(cdouble Omega, double q2D[2],
                                         double zDest, double zSource,
                                         HMatrix *RTwiddle,
                                         HMatrix *WMatrix,
                                         HMatrix *STwiddle,
                                         HMatrix *GTwiddle,
                                         bool dzDest, bool dzSource)
{
  UpdateCachedEpsMu(Omega);

  /**********************************************************************/
  /**********************************************************************/
  /**********************************************************************/
  if (ForceFreeSpace)
   { 
     cdouble ScriptG0Twiddle[6][6];
     GetScriptG0Twiddle(Omega, q2D, zDest, zSource, ScriptG0Twiddle, dzDest, dzSource);
     for(int Mu=0; Mu<6; Mu++)
      for(int Nu=0; Nu<6; Nu++)
       GTwiddle->SetEntry(Mu,Nu,ScriptG0Twiddle[Mu][Nu]);
     return;
   };
  
  /**********************************************************************/
  /* assemble RHS vector for each (source point, polarization, orientation).*/
  /* loops over p and q run over the interfaces above and below the dest*/
  /* or source point.                                                   */
  /**********************************************************************/
  double TT=Secs();
  int DestLayer   = GetLayerIndex(zDest);
  int pMin        = (DestLayer==0               ? 1 : 0 );
  int pMax        = (DestLayer==NumInterfaces   ? 0 : 1 );
  cdouble ScriptG0TDest[2][6][6];
  for(int p=pMin, nlDest=DestLayer-1+p; p<=pMax; p++, nlDest++)
   GetScriptG0Twiddle(Omega, q2D, zDest, zInterface[nlDest], ScriptG0TDest[p], DestLayer, 0.0, false, dzDest, false);

  int SourceLayer  = GetLayerIndex(zSource);
  int qMin         = (SourceLayer==0             ? 1 : 0 );
  int qMax         = (SourceLayer==NumInterfaces ? 0 : 1 );
  cdouble ScriptG0TSource[2][6][6];
  for(int q=qMin, nlSource=SourceLayer-1+q; q<=qMax; q++, nlSource++)
   GetScriptG0Twiddle(Omega, q2D, zInterface[nlSource], zSource, ScriptG0TSource[q], SourceLayer, 0.0, false, false, dzSource);

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
  for(int p=pMin, nlDest=DestLayer-1+p; p<=pMax; p++, nlDest++)
   for(int i=0; i<6; i++)
    { RTwiddle->SetEntry(i, 4*nlDest+0, Sign[p]*ScriptG0TDest[p][i][0]);
      RTwiddle->SetEntry(i, 4*nlDest+1, Sign[p]*ScriptG0TDest[p][i][1]);
      RTwiddle->SetEntry(i, 4*nlDest+2, Sign[p]*ScriptG0TDest[p][i][3]);
      RTwiddle->SetEntry(i, 4*nlDest+3, Sign[p]*ScriptG0TDest[p][i][4]);
    };
  for(int q=qMin, nlSource=SourceLayer-1+q; q<=qMax; q++, nlSource++)
   for(int j=0; j<6; j++)
    { STwiddle->SetEntry(4*nlSource+0, j, -1.0*Sign[q]*ScriptG0TSource[q][0][j]);
      STwiddle->SetEntry(4*nlSource+1, j, -1.0*Sign[q]*ScriptG0TSource[q][1][j]);
      STwiddle->SetEntry(4*nlSource+2, j, -1.0*Sign[q]*ScriptG0TSource[q][3][j]);
      STwiddle->SetEntry(4*nlSource+3, j, -1.0*Sign[q]*ScriptG0TSource[q][4][j]);
    };
  WMatrix->LUSolve(STwiddle);
  RTwiddle->Multiply(STwiddle, GTwiddle);
  Times[SOLVETIME]+=Secs()-TT;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct qFunctionData
 {
   HMatrix *XMatrix;
   HMatrix *RTwiddle;
   HMatrix *WMatrix;
   HMatrix *STwiddle;
   HMatrix *GTwiddle;
   bool dRho, dzDest, dzSource;
   FILE *byqFile;
   bool *RetainIScalar;
 } qFunctionData;

// integrand for "full," i.e. slow brute-force,
// computation of substrate DGF 
void qFunctionFullSC(LayeredSubstrate *Substrate,
                     double q2D[2], cdouble Omega,
                     void *UserData, cdouble *Integrand)
{
  qFunctionData *Data = (qFunctionData *)UserData;

  HMatrix *XMatrix    = Data->XMatrix;
  HMatrix *RTwiddle   = Data->RTwiddle;
  HMatrix *WMatrix    = Data->WMatrix;
  HMatrix *STwiddle   = Data->STwiddle;
  HMatrix *GTwiddle   = Data->GTwiddle;
  FILE *byqFile       = Data->byqFile;

  int EntryOnly       = Substrate->EntryOnly;
  bool EEOnly         = Substrate->EEOnly;
  bool XYOnly         = Substrate->XYOnly;
  int NumInterfaces   = Substrate->NumInterfaces;
  double *zInterface  = Substrate->zInterface;
  //cdouble *EpsLayer   = Substrate->EpsLayer;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double LastzDest=HUGE_VAL, LastzSource=HUGE_VAL;
  cdouble Gxx0=0.0, Gyy0=0.0;
  bool OnInterface=false;
  for(int nx=0, ni=0; nx<XMatrix->NR; nx++)
   { 
     double Rho[2];
     Rho[0]          = XMatrix->GetEntryD(nx,0) - XMatrix->GetEntryD(nx,3);
     Rho[1]          = XMatrix->GetEntryD(nx,1) - XMatrix->GetEntryD(nx,4);
     double zDest    = XMatrix->GetEntryD(nx,2);
     double zSource  = XMatrix->GetEntryD(nx,5);
     cdouble ExpFac  = exp(II*(q2D[0]*Rho[0] + q2D[1]*Rho[1]));

     /*--------------------------------------------------------------*/
     /* re-fetch G_0^\twiddle only if zDest or zSource has changed   */
     /*--------------------------------------------------------------*/
     if ( zDest!=LastzDest || zSource!=LastzSource )
      { LastzDest=zDest;
        LastzSource=zSource;
        Substrate->GetScriptGTwiddle(Omega, q2D, zDest, zSource,
                                     RTwiddle, WMatrix, STwiddle, GTwiddle);

        // check for zDest, zSource both on a layer boundary
        // FIXME this is temporarily disabled
        OnInterface=false;
        Gxx0=Gyy0=0.0;
        if ( EqualFloat(zDest, zSource) )
         for(int nz=0; nz<NumInterfaces && !OnInterface; nz++)
          if ( fabs(zDest-zInterface[nz]) < 1.0e-6 )
           { //cdouble EpsA = EpsLayer[nz], EpsB=EpsLayer[nz+1];
             //cdouble Factor = (-0.5*II*ZVAC/Omega)*(EpsA-EpsB)/(EpsA+EpsB);
             Gxx0 = 0.0; // Factor*qMag*(J0 - J1oqRho);
             Gyy0 = 0.0; // Factor*qMag*J1oqRho;
             OnInterface=true;
           };
        if (OnInterface)
         { GTwiddle->AddEntry(0,0,-1.0*Gxx0);
           GTwiddle->AddEntry(0,0,-1.0*Gyy0);
         };
        if (EEOnly)
         for(int Mu=0; Mu<6; Mu++)
          for(int Nu = (Mu<2 ? 2 : 0); Nu<6; Nu++)
           GTwiddle->SetEntry(Mu,Nu,0.0);
        if (XYOnly)
         for(int Mu=0; Mu<6; Mu++)
          { GTwiddle->SetEntry(Mu,2,0.0);
            GTwiddle->SetEntry(Mu,5,0.0);
            GTwiddle->SetEntry(2,Mu,0.0);
            GTwiddle->SetEntry(5,Mu,0.0);
          };
        if (EntryOnly!=-1)
         { cdouble GEntry = GTwiddle->ZM[EntryOnly];
           GTwiddle->Zero();
           GTwiddle->ZM[EntryOnly]=GEntry;
         };
      };

     // stamp into integrand vector
     for(int Mu=0; Mu<6; Mu++)
      for(int Nu=0; Nu<6; Nu++)
       Integrand[ni++] = GTwiddle->GetEntry(Mu, Nu) * ExpFac;

     if (byqFile)
      { fprintf(byqFile,"%e %e ",q2D[0],q2D[1]);
        fprintf(byqFile,"%e %e %e %e ",Rho[0],Rho[1],zDest,zSource);
        fprintVec(byqFile,(cdouble *)GTwiddle, 36);
        fprintf(byqFile,"%e %e %e %e ",real(Gxx0), imag(Gxx0), real(Gyy0), imag(Gyy0));
        fprintf(byqFile,"\n");
        fflush(byqFile);
      };

   }; // for(int nx=0, ni=0; nx<XMatrix->NR; nx++)

}

/***************************************************************/
/* Get the substrate Green's function by evaluating the full   */
/* 2D Fourier integral with no fancy accelerations. This is    */
/* too slow for use in practical calculations but offers a     */
/* helpful sanity check for debugging, etc.                    */
/***************************************************************/
void LayeredSubstrate::GetSubstrateDGF_FullSurfaceCurrent(cdouble Omega,
                                                          HMatrix *XMatrix,
                                                          HMatrix *GMatrix)
{
  UpdateCachedEpsMu(Omega);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  qFunctionData MyData, *Data=&MyData;
  Data->XMatrix    = XMatrix;
  Data->RTwiddle   = new HMatrix(6, 4*NumInterfaces, LHM_COMPLEX);
  Data->WMatrix    = new HMatrix(4*NumInterfaces, 4*NumInterfaces, LHM_COMPLEX);
  Data->STwiddle   = new HMatrix(4*NumInterfaces, 6,               LHM_COMPLEX);
  Data->GTwiddle   = new HMatrix(6,               6,               LHM_COMPLEX);
  Data->byqFile    = WritebyqFiles ? fopen("/tmp/q2D.log","w") : 0;

  int FDim = 36*XMatrix->NR;
  bool ThetaSymmetric=false;
  qIntegrate(Omega, qFunctionFullSC, (void *)Data, GMatrix->ZM, FDim, ThetaSymmetric);
  
  if (Data->byqFile) fclose(Data->byqFile);

  delete Data->RTwiddle;
  delete Data->WMatrix;
  delete Data->STwiddle;
  delete Data->GTwiddle;
}

/***************************************************************/
/* fetch Bessel-function factors                               */
/*  JdJFactors[0][i] = J[i]                                    */
/*  JdJFactors[1][i] = d/dRho J_i                              */
/*                                                             */
/* J[_J0] = J_0(q*Rho)                                         */
/* J[_J1] = I*J_1(q*Rho)                                       */
/* J[_J2] = -1*J_2(q*Rho)                                      */
/* J[_JJ] = J_1(q*Rho) / (q*Rho)                               */
/*                                                             */
/***************************************************************/
#define _J0 0
#define _J1 1
#define _J2 2
#define _JJ 3
void GetJdJFactors(double q, double Rho, cdouble JdJFactors[2][4], bool NeedRhoDerivatives)
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
     if (NeedRhoDerivatives)
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
     if (NeedRhoDerivatives)
      { dJ0      = -0.5*J0/Rho - q*JPreFac*sin(qRho-0.25*M_PI);
        dJ1      = -0.5*J1/Rho - q*JPreFac*sin(qRho-0.75*M_PI);
        dJ1oqRho =    dJ1/qRho - J1oqRho/Rho;
        dJ2      = -0.5*J2/Rho - q*JPreFac*sin(qRho-1.25*M_PI);
      };
   }
  else
   { double Workspace[12];
     cdouble JFactors[3];
     AmosBessel('J', qRho, 0.0, 3, false, JFactors, Workspace);
     J0=JFactors[0];
     J1=JFactors[1];
     J1oqRho = JFactors[1]/qRho;
     J2=JFactors[2];
     if (NeedRhoDerivatives)
      { dJ0 = -q*J1;
        dJ1 = 0.5*q*(J0 - J2);
        dJ1oqRho = dJ1/qRho - J1oqRho/Rho;
        dJ2 = q*(J1 - 2.0*J2/qRho);
      };
   };

  JdJFactors[0][_J0] = J0;
  JdJFactors[0][_J1] = II*J1;
  JdJFactors[0][_J2] = -1.0*J2;
  JdJFactors[0][_JJ] = J1oqRho;

  if (!NeedRhoDerivatives) return;

  JdJFactors[1][_J0] = dJ0;
  JdJFactors[1][_J1] = II*dJ1;
  JdJFactors[1][_J2] = -1.0*dJ2;
  JdJFactors[1][_JJ] = dJ1oqRho;

}

/***************************************************************/
/* Extract the q_\theta-indepedent scalar quantities from the  */
/* full 36-component Fourier-space Green's function.           */
/* In what follows, 'Reduced' refers to the case in which both */
/* source and destination pnts lie in the same substrate layer,*/
/* in which case some of the gScalars may be computed from the */
/* others, reducing the number of gScalars and IScalars we     */
/* need to compute by about 33%. However, at the moment we     */
/* don't bother to take advantage of this.                     */
/***************************************************************/
#define _EE0P  0
#define _EE0Z  1
#define _EE1   2
#define _EE2   3
#define _EM0   4
#define _EM1A  5
#define _EM1B  6
#define _EM2   7
#define _MM0P  8
#define _MM0Z  9
#define _MM1  10
#define _MM2  11
#define NUMGSCALARS_REDUCED 12

#define _EE1X  12
#define _ME0   13
#define _ME1A  14
#define _ME1B  15
#define _ME2   16
#define _MM1X  17

#define NUMGSCALARS_FULL 18
#define NUMGSCALARS NUMGSCALARS_FULL

// gdgFactors[0][0] = g scalar factors
// gdgFactors[0][1] = dg/dzDest   (if dzDest=true)
// gdgFactors[1][0] = dg/dzSource (if dzSource=true)
// gdgFactors[1][1] = d2g/dzDestdzSource (if dzDest=dzSource=true)
void GetgScalars(LayeredSubstrate *Substrate, cdouble Omega,
                 double qMag, double zDest, double zSource,
                 HMatrix *RTwiddle, HMatrix *WMatrix,
                 HMatrix *STwiddle, cdouble gdgFactors[2][2][NUMGSCALARS],
                 bool dzDest, bool dzSource)
{
  cdouble gBuf[4*36];
  HMatrix   GTwiddle(6,6,LHM_COMPLEX,gBuf+0*36);
  HMatrix     dGTdzD(6,6,LHM_COMPLEX,gBuf+1*36);
  HMatrix     dGTdzS(6,6,LHM_COMPLEX,gBuf+2*36);
  HMatrix d2GTdzDdzS(6,6,LHM_COMPLEX,gBuf+3*36);

  double q2D[2];
  q2D[0]=qMag;
  q2D[1]=0.0;

  for(int ndzDest=0; ndzDest <= (dzDest ? 1 : 0); ndzDest++)
   for(int ndzSource=0; ndzSource <= (dzSource ? 1 : 0); ndzSource++)
    {   
     Substrate->GetScriptGTwiddle(Omega, q2D, zDest, zSource,
                                  RTwiddle, WMatrix, STwiddle, 
                                  &GTwiddle, 
                                  (ndzDest==1), (ndzSource==1));

     cdouble *gScalars = gdgFactors[ndzDest][ndzSource];

     gScalars[_EE0P] = GTwiddle.GetEntry(0+1,0+1);
     gScalars[_EE0Z] = GTwiddle.GetEntry(0+2,0+2);
     gScalars[_EE1]  = GTwiddle.GetEntry(0+0,0+2);
     gScalars[_EE2]  = GTwiddle.GetEntry(0+0,0+0) - gScalars[_EE0P];

     gScalars[_EM0]  = GTwiddle.GetEntry(0+0,3+1);
     gScalars[_EM1A] = GTwiddle.GetEntry(0+1,3+2);
     gScalars[_EM1B] = GTwiddle.GetEntry(0+2,3+1);
     gScalars[_EM2]  = -1.0*(GTwiddle.GetEntry(0+1,3+0) + gScalars[_EM0]);

     gScalars[_MM0P] = GTwiddle.GetEntry(3+1,3+1);
     gScalars[_MM0Z] = GTwiddle.GetEntry(3+2,3+2);
     gScalars[_MM1]  = GTwiddle.GetEntry(3+0,3+2);
     gScalars[_MM2]  = GTwiddle.GetEntry(3+0,3+0) - gScalars[_MM0P];

     // the following are only independent of the preceding
     // if zDest, zSource lie in different substrate layers
     // so in principle their calculation could be omitted
     // for the same-layer case
     gScalars[_EE1X] = GTwiddle.GetEntry(0+2,0+0);
     gScalars[_ME0 ] = GTwiddle.GetEntry(3+1,0+0);
     gScalars[_ME1A] = GTwiddle.GetEntry(3+2,0+1);
     gScalars[_ME1B] = GTwiddle.GetEntry(3+1,0+2);
     gScalars[_ME2 ] = -1.0*(GTwiddle.GetEntry(3+0,0+1) + gScalars[_ME0]);
     gScalars[_MM1X] = GTwiddle.GetEntry(3+2,3+0);
   };

}

/***************************************************************/
/* 1D Fourier integrand for computing 'I scalar' integrals     */
/***************************************************************/
#define _EE0P_J0 0
#define _EE0Z_J0 1
#define _EE1_J1  2
#define _EE2_J2  3
#define _EE2_JJ  4
#define _EM0_J0  5
#define _EM1A_J1 6
#define _EM1B_J1 7
#define _EM2_J2  8
#define _EM2_JJ  9
#define _MM0P_J0 10
#define _MM0Z_J0 11
#define _MM1_J1  12
#define _MM2_J2  13
#define _MM2_JJ  14
#define NUMISCALARS_REDUCED 15

#define _EE1X_J1 15
#define _ME0_J0  16
#define _ME1A_J1 17
#define _ME1B_J1 18
#define _ME2_J2  19
#define _ME2_JJ  20
#define _MM1X_J1 21
#define NUMISCALARS_FULL 22
#define NUMISCALARS NUMISCALARS_FULL

typedef struct IScalarPair
 { int WhichG, WhichJ;
 } IScalarPair;

const IScalarPair IPairs[NUMISCALARS]={
 {_EE0P, _J0}, // EE0P_J0
 {_EE0Z, _J0}, // EE0Z_J0
 {_EE1,  _J1}, // EE1_J1
 {_EE2,  _J2}, // EE2_J2
 {_EE2,  _JJ}, // EE2_JJ
 {_EM0,  _J0}, // EM0_J0
 {_EM1A, _J1}, // EM1A_J1
 {_EM1B, _J1}, // EM1B_J1
 {_EM2,  _J2}, // EM2_J2
 {_EM2,  _JJ}, // EM2_JJ
 {_MM0P, _J0}, // MM0P_J0
 {_MM0Z, _J0}, // MM0Z_J0
 {_MM1,  _J1}, // MM1_J1
 {_MM2,  _J2}, // MM2_J2
 {_MM2,  _JJ}, // MM2_JJ
 {_EE1X, _J1}, // EE1X_J1
 {_ME0,  _J0}, // ME0_J0
 {_ME1A, _J1}, // ME1A_J1
 {_ME1B, _J1}, // ME1B_J1
 {_ME2,  _J2}, // ME2_J2
 {_ME2,  _JJ}, // ME2_JJ
 {_MM1X, _J1}  // MM1X_J1
};

// return the index, within a vector of IScalars, of
// the nVDth value/derivative of IScalar #nI at evaluation point #nx
// nVD = 0 --> I
// nVD = 1 --> dI/dRho
// nVD = 2 --> dI/dz
// nVD = 3 --> dI/dzS          or d^2I/dRhodzDS
// nVD = 4 --> d^2I/dRho dzD
// nVD = 5 --> d^2I/dRho dzS
// nVD = 6 --> d^2I/dzD dzS
// nVD = 7 --> d^3I/dRho dzD dzS
#define ISCALAR_INDEX(nx,nI,nVD) ( nx*NIVD + nI*NVD + nVD )

void qFunction_IScalars(LayeredSubstrate *Substrate,
                        double q2D[2], cdouble Omega,
                        void *UserData, cdouble *Integrand)
{
  qFunctionData *Data = (qFunctionData *)UserData;

  HMatrix *XMatrix     = Data->XMatrix;
  HMatrix *RTwiddle    = Data->RTwiddle;
  HMatrix *WMatrix     = Data->WMatrix;
  HMatrix *STwiddle    = Data->STwiddle;
  FILE *byqFile        = Data->byqFile;

  bool dRho            = Data->dRho;
  bool dzDest          = Data->dzDest;
  bool dzSource        = Data->dzSource;

  bool *RetainIScalar  = Data->RetainIScalar;

  //bool EEOnly         = Substrate->EEOnly;
  //bool XYOnly         = Substrate->XYOnly;
  int NumInterfaces   = Substrate->NumInterfaces;
  double *zInterface  = Substrate->zInterface;
  cdouble *EpsLayer   = Substrate->EpsLayer;
  cdouble *MuLayer    = Substrate->MuLayer;

  double qMag = q2D[0];
  
  int NVD=1;            // number of values and derivatives per IScalar
  if (dRho)     NVD*=2;
  if (dzDest)   NVD*=2;
  if (dzSource) NVD*=2;
  int NIVD = NVD*NUMISCALARS; // total number of IScalar values and derivatives per eval point

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double LastzDest=HUGE_VAL, LastzSource=HUGE_VAL;
  cdouble gdgFactors[2][2][NUMGSCALARS];
  cdouble GEEStatic, GMMStatic;
  bool OnInterface=false;
  for(int nx=0; nx<XMatrix->NR; nx++)
   { 
     double Rho[2];
     Rho[0]          = XMatrix->GetEntryD(nx,0) - XMatrix->GetEntryD(nx,3);
     Rho[1]          = XMatrix->GetEntryD(nx,1) - XMatrix->GetEntryD(nx,4);
     double zDest    = XMatrix->GetEntryD(nx,2);
     double zSource  = XMatrix->GetEntryD(nx,5);
     double RhoMag   = sqrt(Rho[0]*Rho[0] + Rho[1]*Rho[1]);

     /*--------------------------------------------------------------*/
     /* fetch Bessel-function factors -------------------------------*/
     /*--------------------------------------------------------------*/
     cdouble JdJFactors[2][4];
     GetJdJFactors(qMag, RhoMag, JdJFactors, dRho);

     /*--------------------------------------------------------------*/
     /* re-fetch g scalars only if zDest or zSource has changed      */
     /*--------------------------------------------------------------*/
     if ( zDest!=LastzDest || zSource!=LastzSource )
      { 
        LastzDest=zDest;
        LastzSource=zSource;
        GetgScalars(Substrate, Omega, qMag, zDest, zSource,
                    RTwiddle, WMatrix, STwiddle, gdgFactors,
                    dzDest, dzSource);

        // check for zDest, zSource both on a layer boundary
        OnInterface=false;
        GEEStatic = GMMStatic = 0.0;
        if ( EqualFloat(zDest, zSource) )
         for(int nz=0; nz<NumInterfaces && !OnInterface; nz++)
          if ( fabs(zDest-zInterface[nz]) < 1.0e-6 )
           { OnInterface=true;
             cdouble EpsA = EpsLayer[nz], EpsB=EpsLayer[nz+1];
             cdouble MuA  = MuLayer[nz],  MuB=MuLayer[nz+1];
             GEEStatic = qMag*(-0.5*II*ZVAC/Omega)*(EpsA-EpsB)/(EpsA+EpsB);
             GMMStatic = qMag*(-0.5*II*Omega/ZVAC)*(MuA-MuB)/(MuA+MuB);
           };
      };

     /*--------------------------------------------------------------*/
     /* assemble integrand ------------------------------------------*/
     /*--------------------------------------------------------------*/
     for(int nVD=0; nVD<NVD; nVD++)
      { 
        int d1Order = (nVD>>0) % 2;
        int d2Order = (nVD>>1) % 2;
        int d3Order = (nVD>>2) % 2;

        cdouble *J = JdJFactors[d1Order];
        cdouble *g = gdgFactors[d2Order][d3Order];

        for(int nI=0; nI<NUMISCALARS; nI++)
         if (!RetainIScalar || RetainIScalar[nI])
          Integrand[ nx*NIVD + nI*NVD + nVD ]=g[ IPairs[nI].WhichG ]*J[ IPairs[nI].WhichJ ];
      };

     if (OnInterface)
      { 
        Integrand[ ISCALAR_INDEX(nx,_EE0P_J0, 0) ] -= GEEStatic;
        Integrand[ ISCALAR_INDEX(nx,_EE2_JJ,  0) ] += GEEStatic;
        Integrand[ ISCALAR_INDEX(nx,_MM0P_J0, 0) ] -= GMMStatic;
        Integrand[ ISCALAR_INDEX(nx,_MM2_JJ,  0) ] += GMMStatic;
      };
     
     if (byqFile)
      { fprintf(byqFile,"%e ",qMag);
        fprintf(byqFile,"%e %e %e %e ",Rho[0],Rho[1],zDest,zSource);
        fprintVec(byqFile,Integrand+nx*NVD*NUMISCALARS, NVD*NUMISCALARS);
        fprintf(byqFile,"%e %e %e %e ",real(GEEStatic),imag(GEEStatic),real(GMMStatic),imag(GMMStatic));
        fprintf(byqFile,"\n");
        fflush(byqFile);
      };

   }; // for(int nx=0, ni=0; nx<XMatrix->NR; nx++)
}

/***************************************************************/
/* DerivativeDimension == 0 --> no derivatives computed        */
/* DerivativeDimension == 1 --> need Rho derivatives           */
/* DerivativeDimension == 2 --> need Rho, zD derivatives       */
/* DerivativeDimension == 3 --> need Rho, zD, zS derivatives   */
/***************************************************************/
void LayeredSubstrate::GetIScalars(cdouble Omega,
                                   HMatrix *XMatrix,
                                   cdouble *IScalars,
                                   int DerivativeDimension)
{
  UpdateCachedEpsMu(Omega);

  /*--------------------------------------------------------------*/
  /*- evaluate scalar I integrals for all evaluation points  .    */
  /*--------------------------------------------------------------*/
  qFunctionData MyData, *Data=&MyData;
  Data->XMatrix    = XMatrix;
  Data->RTwiddle   = new HMatrix(6, 4*NumInterfaces, LHM_COMPLEX);
  Data->WMatrix    = new HMatrix(4*NumInterfaces, 4*NumInterfaces, LHM_COMPLEX);
  Data->STwiddle   = new HMatrix(4*NumInterfaces, 6,               LHM_COMPLEX);
  Data->GTwiddle   = new HMatrix(6,               6,               LHM_COMPLEX);
  Data->byqFile    = WritebyqFiles ? fopen("/tmp/IScalars.log","w") : 0;

  Data->dRho       = (DerivativeDimension>=1);
  Data->dzDest     = (DerivativeDimension>=2);
  Data->dzSource   = (DerivativeDimension==3);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Data->RetainIScalar = 0;
  bool RetainIScalar[NUMISCALARS];
  char *s=getenv("SCUFF_SUBSTRATE_ISCALARS");
  if(s)
   { 
     for(int nI=0; nI<NUMISCALARS; nI++)
      RetainIScalar[nI]=false;

     char sCopy[1000];
     strncpy(sCopy, s, 1000);
     char *Tokens[NUMISCALARS];
     int NumTokens=Tokenize(sCopy, Tokens, NUMISCALARS);
     for(int nt=0; nt<NumTokens; nt++)
      { int nI;
        if (1==sscanf(Tokens[nt],"%i",&nI))
         { RetainIScalar[nI]=true;
           Log("Retaining IScalar #%i.",nI);
           printf("Retaining IScalar #%i.\n",nI);
         };
      };
     Data->RetainIScalar=RetainIScalar;
   };

  int NX   = XMatrix->NR;
  int FDim = NX*NUMISCALARS*(1<<DerivativeDimension);
  bool ThetaSymmetric=true;
  qIntegrate(Omega, qFunction_IScalars, (void *)Data, IScalars, FDim, ThetaSymmetric);
  
  if (Data->byqFile) fclose(Data->byqFile);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void IScalars2Gij(cdouble *IScalars, double rHat[2], cdouble Gij[6][6])
{
  // EE, MM parallel-parallel components
  for(int i=0; i<2; i++)
   for(int j=0; j<2; j++)
    { 
      Gij[0+i][0+j] = IScalars[_EE2_J2]*rHat[i]*rHat[j];
      Gij[3+i][3+j] = IScalars[_MM2_J2]*rHat[i]*rHat[j];
      if (i!=j) continue;
      Gij[0+i][0+j] += IScalars[_EE0P_J0] + IScalars[_EE2_JJ];
      Gij[3+i][3+j] += IScalars[_MM0P_J0] + IScalars[_MM2_JJ];
    };

  // EE, MM parallel-z components
  for(int i=0; i<2; i++)
   { Gij[0+i][0+2] = IScalars[_EE1_J1]*rHat[i];
     Gij[0+2][0+i] = IScalars[_EE1X_J1]*rHat[i];
     Gij[3+i][3+2] = IScalars[_MM1_J1]*rHat[i];
     Gij[3+2][3+i] = IScalars[_MM1X_J1]*rHat[i];
   };
 
  // EE, MM zz components
  Gij[0+2][0+2] = IScalars[_EE0Z_J0];
  Gij[3+2][3+2] = IScalars[_MM0Z_J0];

  // EM parallel-parallel components
  Gij[0+0][3+0] = +1.0*IScalars[_EM2_J2]*rHat[0]*rHat[1];
  Gij[0+0][3+1] = +1.0*(IScalars[_EM0_J0] + IScalars[_EM2_JJ]) + IScalars[_EM2_J2]*rHat[1]*rHat[1];
  Gij[0+1][3+0] = -1.0*(IScalars[_EM0_J0] + IScalars[_EM2_JJ]) - IScalars[_EM2_J2]*rHat[0]*rHat[0];
  Gij[0+1][3+1] = -1.0*IScalars[_EM2_J2]*rHat[0]*rHat[1];

  // EM parallel-z components
  Gij[0+0][3+2] = -1.0*IScalars[_EM1A_J1]*rHat[1];
  Gij[0+1][3+2] = +1.0*IScalars[_EM1A_J1]*rHat[0];
  Gij[0+2][3+0] = -1.0*IScalars[_EM1B_J1]*rHat[1];
  Gij[0+2][3+1] = +1.0*IScalars[_EM1B_J1]*rHat[0];

  // ME parallel-parallel components
  Gij[3+0][0+0] = +1.0*IScalars[_ME2_J2]*rHat[0]*rHat[1];
  Gij[3+1][0+0] = +1.0*(IScalars[_ME0_J0] + IScalars[_ME2_JJ]) + IScalars[_ME2_J2]*rHat[1]*rHat[1];
  Gij[3+0][0+1] = -1.0*(IScalars[_ME0_J0] + IScalars[_ME2_JJ]) - IScalars[_ME2_J2]*rHat[0]*rHat[0];
  Gij[3+1][0+1] = -1.0*IScalars[_ME2_J2]*rHat[0]*rHat[1];

  // ME parallel-z components
  Gij[3+2][0+0] = -1.0*IScalars[_ME1A_J1]*rHat[1];
  Gij[3+2][0+1] = +1.0*IScalars[_ME1A_J1]*rHat[0];
  Gij[3+0][0+2] = -1.0*IScalars[_ME1B_J1]*rHat[1];
  Gij[3+1][0+2] = +1.0*IScalars[_ME1B_J1]*rHat[0];

#if 0
  // ME 
  Gij[3+0][0+0] = -1.0*Gij[0+0][3+0];
  Gij[3+0][0+1] = -1.0*Gij[0+1][3+0];
  Gij[3+1][0+0] = -1.0*Gij[0+0][3+1];
  Gij[3+1][0+1] = -1.0*Gij[0+1][3+1];

  Gij[3+2][0+0] = +1.0*Gij[0+0][3+2];
  Gij[3+2][0+1] = +1.0*Gij[0+1][3+2];
  Gij[3+0][0+2] = +1.0*Gij[0+2][3+0];
  Gij[3+1][0+2] = +1.0*Gij[0+2][3+1];
#endif

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void LayeredSubstrate::GetSubstrateDGF_FastSurfaceCurrent(cdouble Omega,
                                                          HMatrix *XMatrix,
                                                          HMatrix *GMatrix)
{
  UpdateCachedEpsMu(Omega);

  int NX = XMatrix->NR;
  GetIScalars(Omega, XMatrix, GMatrix->ZM);

  for(int nx=NX-1; nx>=0; nx--)
   { 
     double Rho[2]; 
     Rho[0] = XMatrix->GetEntryD(nx,0) - XMatrix->GetEntryD(nx,3);
     Rho[1] = XMatrix->GetEntryD(nx,1) - XMatrix->GetEntryD(nx,4);
     double RhoMag = sqrt(Rho[0]*Rho[0] + Rho[1]*Rho[1]);
     double RhoHat[2]={1.0, 0.0};
     if (RhoMag>1.0e-12)
      { RhoHat[0] = Rho[0] / RhoMag;
        RhoHat[1] = Rho[1] / RhoMag;
      };

     cdouble Gij[6][6];
     cdouble *IScalars = GMatrix->ZM + nx*NUMISCALARS;
     IScalars2Gij(IScalars, RhoHat, Gij);
     for(int Mu=0; Mu<6; Mu++)
      for(int Nu=0; Nu<6; Nu++)
       GMatrix->SetEntry(6*Mu + Nu, nx, Gij[Mu][Nu]);
   };
}
