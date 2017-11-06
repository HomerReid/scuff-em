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
 * GetSphericalMoments.cc  -- libscuff methods for computing the
 *                         -- spherical multipole moments induced by
 *                         -- an incident field on the scattering objects
 *
 * the algorithm used here is described in the following places:
 *
 * 1) http://homerreid.github.io/scuff-em-documentation/applications/scuff-tmatrix/scuff-tmatrix#Algorithm
 *
 * 2) http://homerreid.github.io/scuff-em-documentation/tex/scuffSpherical.pdf
 *    (see Section 10)
 *
 * homer reid              -- 3/2007  -- 11/2017
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <libhrutil.h>
#include <libSpherical.h>
#include "libscuff.h"
#include "PanelCubature.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#ifdef USE_OPENMP
 #include <omp.h>
#endif

using namespace scuff;

/***************************************************************/
/* integrand function passed to GetBFCubature2 to compute      */
/* inner products of a single RWG function with the M, N       */
/* vector spherical waves (VSWs)                               */
/***************************************************************/
typedef struct IntegrandData
 {
   double lMax;
   cdouble k;
   HMatrix *MWMatrix;
   cdouble *Workspace;
 } IntegrandData;

void VSWDotRWGIntegrand(double x[3], double b[3], double Divb,
                        void *UserData, double Weight, 
                        double *Integral)
{
  IntegrandData *Data = (IntegrandData *)UserData;
  double lMax         = Data->lMax;
  cdouble k           = Data->k;
  HMatrix *MWMatrix   = Data->MWMatrix;
  cdouble *Workspace  = Data->Workspace;

  // get spherical vector components of RWG basis function
  double r, Theta, Phi, bS[3];
  CoordinateC2S(x,&r,&Theta,&Phi);
  VectorC2S(Theta,Phi,b,bS);

  // get components of VSWs at cubature point
  // MWMatrix[i][nmw] = ith component of nmw#th VSW
  bool ConjugateRFactor=true;
  GetMWMatrix(r, Theta, Phi, k, lMax, LS_REGULAR, 
              MWMatrix, Workspace, ConjugateRFactor);
  
  cdouble *zIntegral=(cdouble *)Integral;
  for(int nmw=0; nmw<MWMatrix->NC; nmw++)
   { cdouble *W      = (cdouble *)MWMatrix->GetColumnPointer(nmw);
     zIntegral[nmw] += Weight*conj(W[0]*bS[0] + W[1]*bS[1] + W[2]*bS[2]);
   };
}

/***************************************************************/
/* get the spherical multipole moments of a single object      */
/* (surface) in a scattering geometry                          */
/***************************************************************/
HVector *GetSphericalMoments(RWGGeometry *G, int ns,
                             cdouble Omega, int lMax,
                             HVector *KN, int BFIndexOffset,
                             HVector *MomentVector)
{ 
  
  /***************************************************************/
  /* (re)allocate the MomentVector as necessary ***********************/
  /***************************************************************/
  int NumLMs = (lMax+1)*(lMax+1) - 1;
  int NumMoments = 2*NumLMs; // M-type and N-type moments for each l,m
  if ( MomentVector && (MomentVector->N != NumMoments ) )
   { Warn("wrong-size MomentVector passed to GetSphericalMoments (reallocating...)");
     delete MomentVector;
     MomentVector=0;
   };
  if ( MomentVector==0 )
   MomentVector=new HVector(NumMoments, LHM_COMPLEX);
  MomentVector->Zero();

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
#define EXTERIOR_REGION 0  // index of exterior region
  RWGSurface *S=G->Surfaces[ns];
  double Sign; 
  if(S->RegionIndices[0]==EXTERIOR_REGION)
   Sign=+1.0;
  else if (S->RegionIndices[1]==EXTERIOR_REGION)
   Sign=-1.0;
  else 
   return MomentVector; // currents on this surface do not contribute
   
  cdouble EpsRel, MuRel;
  G->RegionMPs[EXTERIOR_REGION]->GetEpsMu(Omega, &EpsRel, &MuRel);
  cdouble k2 = EpsRel*MuRel*Omega*Omega, k=sqrt(k2);

  /***************************************************************/
  /* set environment variable SCUFF_SPHERICAL_SINGLETHREADED=1   */
  /* to disable multithreading here (useful for debugging).      */
  /* set environment variable SCUFF_SPHERICAL_MOMENT_ORDER       */
  /*  = 1, 4, 7, 9, 16, 20                                       */
  /* to set order of cubature rule                               */
  /***************************************************************/
  int NT=GetNumThreads();
  char *s=getenv("SCUFF_SPHERICAL_SINGLETHREADED");
  if ( s && s[0]=='1' )
   NT=1;
  int Order=20;
  s=getenv("SCUFF_SPHERICAL_MOMENT_ORDER");
  if (s && 1==sscanf(s,"%i",&Order))
   Log("Using cubature order %i in GetSphericalMoments", Order);

  /***************************************************************/
  /* allocate per-thread storage buffers to prevent race conditions */
  /***************************************************************/
  IntegrandData *Data      = (IntegrandData *)mallocEC(NT*sizeof(Data[0]));
  HVector **PartialMoments = (HVector**)mallocEC(NT*sizeof(HVector *));
  cdouble *IBuffer         = (cdouble *)mallocEC(NT*NumMoments*sizeof(cdouble));
  Log("Computing induced spherical moments on surface %s (%i threads)",S->Label,NT);
  for(int nt=0; nt<NT; nt++)
   { Data[nt].lMax      = lMax;
     Data[nt].k         = k;
     Data[nt].MWMatrix  = new HMatrix(3, NumMoments, LHM_COMPLEX);
     Data[nt].Workspace = (cdouble *)mallocEC(7*(NumLMs+1)*sizeof(cdouble));
     PartialMoments[nt] = new HVector(NumMoments, LHM_COMPLEX);
   };

  /***************************************************************/
  /* loop to get contributions of all RWG functions              */
  /***************************************************************/
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1), num_threads(NT)
#endif
  for(int ne=0; ne<S->NumEdges; ne++)
   { 
     int nt=0;
#ifdef USE_OPENMP
     nt = omp_get_thread_num();
#endif

     // get projections of this basis function onto M, N spherical waves
     // Integral[nmw] = <W_{nmw} | b>
     IntegrandData *MyData = Data + nt;
     cdouble *Integral     = IBuffer + nt*NumMoments;
     GetBFCubature2(G, ns, ne, VSWDotRWGIntegrand, (void *)MyData,
                    2*NumMoments, Order, (double *)Integral);

     // add contributions to spherical moments
     cdouble kAlpha, nAlpha;
     G->GetKNCoefficients(KN, ns, ne, &kAlpha, &nAlpha);
     for(int nlm=0; nlm<NumLMs; nlm++)
      { 
        cdouble MdotB = Integral[2*nlm+0], NdotB = Integral[2*nlm+1];
        PartialMoments[nt]->AddEntry(2*nlm+0, -1.0*k2*(ZVAC*kAlpha*MdotB - nAlpha*NdotB));
        PartialMoments[nt]->AddEntry(2*nlm+1, -1.0*k2*(ZVAC*kAlpha*NdotB + nAlpha*MdotB));
      };
   };

  /***************************************************************/
  /* accumulate contributions of all threads                     */
  /***************************************************************/
  MomentVector->Zero();
  for(int nt=0; nt<NT; nt++)
   for(int nm=0; nm<MomentVector->N; nm++)
    MomentVector->AddEntry(nm,PartialMoments[nt]->GetEntry(nm));

  /***************************************************************/
  /* clean up temporary storage **********************************/
  /***************************************************************/
  for(int nt=0; nt<NT; nt++)
   { delete Data[nt].MWMatrix;
     free(Data[nt].Workspace);
     delete PartialMoments[nt];
   };
  free(IBuffer);
  free(PartialMoments);
  free(Data);

  return MomentVector;

}

/***************************************************************/
/* get the spherical moments for an entire geometry by summing */
/* the moments for each individual surface                     */
/***************************************************************/
HVector *GetSphericalMoments(RWGGeometry *G, cdouble Omega, int lMax,
                             HVector *KN, HVector *MomentVector)
{ 
  /***************************************************************/
  /* (re)allocate the MomentVector as necessary ***********************/
  /***************************************************************/
  int NumLMs = (lMax+1)*(lMax+1)-1;
  int NumMoments = 2*NumLMs; // M-type and N-type moments for each l,m
  if ( MomentVector && MomentVector->N != NumMoments )
   { Warn("wrong-size MomentVector passed to GetSphericalMoments (reallocating...)");
     delete MomentVector;
     MomentVector=0;
   };
  if ( MomentVector==0 )
   MomentVector=new HVector(NumMoments, LHM_COMPLEX);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  HVector *PartialMoments=new HVector(NumMoments, LHM_COMPLEX);

  MomentVector->Zero();
  for(int ns=0; ns<G->NumSurfaces; ns++)
   { RWGSurface *S = G->Surfaces[ns];
     int Offset    = G->BFIndexOffset[ns];  
     GetSphericalMoments(G, ns, Omega, lMax, KN, Offset, PartialMoments);

     for(int nm=0; nm<NumMoments; nm++)
      MomentVector->AddEntry(nm, PartialMoments->GetEntry(nm));
   };
  delete PartialMoments;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  return MomentVector;
   
}
