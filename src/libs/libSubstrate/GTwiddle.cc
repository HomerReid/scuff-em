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
 * GTwiddle.cc -- compute the Fourier-space dyadic Green's function
 *             -- function for a layered substrate using the method
 *             -- of induced surface currents
 *
 * homer reid  -- 3/2017-9/2017
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

/***************************************************************/
/* compute the 3x3 submatrices that enter the quadrants of the */
/* Fourier transform of the 6x6 DGF for a homogeneous medium   */
/***************************************************************/
void GetGC0Twiddle(cdouble k2, cdouble q2D[2], cdouble qz, double Sign,
                   cdouble GT[3][3], cdouble CT[3][3])
{
  cdouble qx  = q2D[0];
  cdouble qy  = q2D[1];

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
void LayeredSubstrate::GetScriptG0Twiddle(cdouble Omega, cdouble q2D[2],
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
  if (imag(qz)<0) qz*=-1.0;

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
void LayeredSubstrate::ComputeW(cdouble Omega, cdouble q2D[2], HMatrix *W)
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
/*                                                             */
/* RTwiddle: user-allocated cdouble-valued 6x4N  HMatrix       */ 
/* WMatrix : user-allocated cdouble-valued 4Nx4N HMatrix       */
/* STwiddle: user-allocated cdouble-valued 4Nx6  HMatrix       */
/* GTwiddle: user-allocated cdouble-valued 6x6  HMatrix        */
/*                                                             */
/*  where N = number of dielectric interface layers, not       */
/*  counting the ground plane (if present)                     */
/***************************************************************/
void LayeredSubstrate::GetScriptGTwiddle(cdouble Omega, cdouble q2D[2],
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
  /* loops over a and b run over the interfaces above and below the dest*/
  /* or source point.                                                   */
  /**********************************************************************/
  double TT=Secs();
  int DestLayer   = GetLayerIndex(zDest);
  int aMin        = (DestLayer==0               ? 1 : 0 );
  int aMax        = (DestLayer==NumInterfaces   ? 0 : 1 );
  cdouble ScriptG0TDest[2][6][6];
  for(int a=aMin, nlDest=DestLayer-1+a; a<=aMax; a++, nlDest++)
   GetScriptG0Twiddle(Omega, q2D, zDest, zInterface[nlDest], ScriptG0TDest[a], DestLayer, 0.0, false, dzDest, false);

  int SourceLayer  = GetLayerIndex(zSource);
  int bMin         = (SourceLayer==0             ? 1 : 0 );
  int bMax         = (SourceLayer==NumInterfaces ? 0 : 1 );
  cdouble ScriptG0TSource[2][6][6];
  for(int b=bMin, nlSource=SourceLayer-1+b; b<=bMax; b++, nlSource++)
   GetScriptG0Twiddle(Omega, q2D, zInterface[nlSource], zSource, ScriptG0TSource[b], SourceLayer, 0.0, false, false, dzSource);

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
  for(int a=aMin, nlDest=DestLayer-1+a; a<=aMax; a++, nlDest++)
   for(int i=0; i<6; i++)
    { RTwiddle->SetEntry(i, 4*nlDest+0, Sign[a]*ScriptG0TDest[a][i][0]);
      RTwiddle->SetEntry(i, 4*nlDest+1, Sign[a]*ScriptG0TDest[a][i][1]);
      RTwiddle->SetEntry(i, 4*nlDest+2, Sign[a]*ScriptG0TDest[a][i][3]);
      RTwiddle->SetEntry(i, 4*nlDest+3, Sign[a]*ScriptG0TDest[a][i][4]);
    };
  for(int b=bMin, nlSource=SourceLayer-1+b; b<=bMax; b++, nlSource++)
   for(int j=0; j<6; j++)
    { STwiddle->SetEntry(4*nlSource+0, j, -1.0*Sign[b]*ScriptG0TSource[b][0][j]);
      STwiddle->SetEntry(4*nlSource+1, j, -1.0*Sign[b]*ScriptG0TSource[b][1][j]);
      STwiddle->SetEntry(4*nlSource+2, j, -1.0*Sign[b]*ScriptG0TSource[b][3][j]);
      STwiddle->SetEntry(4*nlSource+3, j, -1.0*Sign[b]*ScriptG0TSource[b][4][j]);
    };
  WMatrix->LUSolve(STwiddle);
  RTwiddle->Multiply(STwiddle, GTwiddle);
  Times[SOLVETIME]+=Secs()-TT;
}

