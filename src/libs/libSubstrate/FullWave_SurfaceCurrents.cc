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
 * FullWave_SurfaceCurrents.cc -- compute the Fourier-space dyadic Green's
 *                             -- function for a layered substrate
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
/* Get the substrate Green's function by evaluating the full   */
/* 2D Fourier integral with no fancy accelerations. This is    */
/* too slow for use in practical calculations but offers a     */
/* helpful sanity check for debugging, etc.                    */
/***************************************************************/
void LayeredSubstrate::GetSubstrateDGF_FullSurfaceCurrent(cdouble Omega,
                                                          HMatrix *XMatrix,
                                                          HMatrix *GMatrix)
{
#if 0
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
#endif
}

/***************************************************************/
/* Extract the 18 distinct theta-independent scalar quantities */
/* from the full 36-component Fourier-space Green's function.  */
/***************************************************************/
#define _EE0P  0
#define _EE0Z  1
#define _EE1A  2
#define _EE1B  3
#define _EE2   4
#define _EM0   5
#define _EM1A  6
#define _EM1B  7
#define _EM2   8
#define _ME0   9
#define _ME1A 10
#define _ME1B 11
#define _ME2  12
#define _MM0P 13 
#define _MM0Z 14
#define _MM1A 15
#define _MM1B 16
#define _MM2  17

#define NUMGSCALARS 18

// gTwiddleVD[0][0] = g scalar factors
// gTwiddleVD[0][1] = dg/dzDest   (if dzDest=true)
// gTwiddleVD[1][0] = dg/dzSource (if dzSource=true)
// gTwiddleVD[1][1] = d2g/dzDestdzSource (if dzDest=dzSource=true)
void GetgTwiddle_SurfaceCurrent(LayeredSubstrate *Substrate, cdouble Omega,
                 double qMag, double zDest, double zSource,
                 HMatrix *RTwiddle, HMatrix *WMatrix,
                 HMatrix *STwiddle, cdouble gTwiddleVD[2][2][NUMGSCALARS],
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

     cdouble *gScalars = gTwiddleVD[ndzDest][ndzSource];

     gScalars[_EE0P]  = GTwiddle.GetEntry(0+1,0+1);
     gScalars[_EE0Z]  = GTwiddle.GetEntry(0+2,0+2);
     gScalars[_EE1A]  = GTwiddle.GetEntry(0+0,0+2);
     gScalars[_EE2]   = GTwiddle.GetEntry(0+0,0+0) - gScalars[_EE0P];

     gScalars[_EM0]  = GTwiddle.GetEntry(0+0,3+1);
     gScalars[_EM1A] = GTwiddle.GetEntry(0+1,3+2);
     gScalars[_EM1B] = GTwiddle.GetEntry(0+2,3+1);
     gScalars[_EM2]  = -1.0*(GTwiddle.GetEntry(0+1,3+0) + gScalars[_EM0]);

     gScalars[_MM0P] = GTwiddle.GetEntry(3+1,3+1);
     gScalars[_MM0Z] = GTwiddle.GetEntry(3+2,3+2);
     gScalars[_MM1A] = GTwiddle.GetEntry(3+0,3+2);
     gScalars[_MM2]  = GTwiddle.GetEntry(3+0,3+0) - gScalars[_MM0P];

     // the following are only independent of the preceding
     // if zDest, zSource lie in different substrate layers
     // so in principle their calculation could be omitted
     // for the same-layer case
     gScalars[_EE1B] = GTwiddle.GetEntry(0+2,0+0);
     gScalars[_ME0 ] = GTwiddle.GetEntry(3+1,0+0);
     gScalars[_ME1A] = GTwiddle.GetEntry(3+2,0+1);
     gScalars[_ME1B] = GTwiddle.GetEntry(3+1,0+2);
     gScalars[_ME2 ] = -1.0*(GTwiddle.GetEntry(3+0,0+1) + gScalars[_ME0]);
     gScalars[_MM1B] = GTwiddle.GetEntry(3+2,3+0);
   };

}
