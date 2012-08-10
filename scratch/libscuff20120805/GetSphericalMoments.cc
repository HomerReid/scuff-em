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
 * GetSphericalMoments.cc  -- libscuff class method for computing the 
 *                         -- spherical multipole moments induced by 
 *                         -- the incident field on the scattering objects
 *
 * homer reid              -- 3/2007  -- 7/2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <libhrutil.h>
#include <libSpherical.h>
#include <libTriInt.h>

#include "libscuff.h"

#define HAVE_CONFIG_H
#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

using namespace scuff;

#define II cdouble(0,1)

/***************************************************************/
/* cubature integrand for the GetMNProjections routine         */
/***************************************************************/
typedef struct GMNPData
 {
   cdouble k; 
   int lMax;
   double *Q;
   double RWGPreFac;
   cdouble *MArray;
   cdouble *NArray;
  
 } GMNPData;

void GMNPIntegrand(double *X, void *parms, double *f)
{
  /*--------------------------------------------------------------*/
  /*- extract parameters from data structure.                     */
  /*- note: the MArray and NArray fields must point to caller-    */
  /*- allocated buffers with enough room for 3*(lMax+1)*(lMax+1)  */
  /*- cdoubles each.                                              */
  /*--------------------------------------------------------------*/
  GMNPData *Data    = (GMNPData *)parms;
  cdouble k         = Data->k;
  int lMax          = Data->lMax;
  double *Q         = Data->Q;
  double RWGPreFac  = Data->RWGPreFac;
  cdouble *MArray   = Data->MArray;
  cdouble *NArray   = Data->NArray;

  /*--------------------------------------------------------------*/
  /*- convert the eval point to spherical coordinates, then get   */
  /*- the M, N vector helmholtz solutions at this eval point.     */
  /*--------------------------------------------------------------*/
  double r, Theta, Phi;
  CoordinateC2S(X, &r, &Theta, &Phi);
  GetMNlmArray(lMax, k, r, Theta, Phi, LS_REGULAR, MArray, NArray);

  /*--------------------------------------------------------------*/
  /*- get the vector-valued RWG basis function at this eval point */
  /*--------------------------------------------------------------*/
  double fRWG[3];
  fRWG[0] = RWGPreFac*(X[0] - Q[0]); 
  fRWG[1] = RWGPreFac*(X[1] - Q[1]); 
  fRWG[2] = RWGPreFac*(X[2] - Q[2]); 

  /*--------------------------------------------------------------*/
  /*- and convert it to spherical components                      */
  /*--------------------------------------------------------------*/
  double fRWGS[3];
  VectorC2S(Theta, Phi, fRWG, fRWGS);

  /*--------------------------------------------------------------*/
  /*- loop over all elements in the M, N array to compute the     */
  /*- elements of the integrand vector                            */
  /*--------------------------------------------------------------*/
  cdouble MDotF, NDotF;
  int nf=0;
  int NumLMs=(lMax+1)*(lMax+1);
  for(int nLM=0; nLM<NumLMs; nLM++)
   { 
     MDotF = MArray[3*nLM + 0]*fRWGS[0] + MArray[3*nLM + 1]*fRWGS[1] + MArray[3*nLM + 2]*fRWGS[2];
     NDotF = NArray[3*nLM + 0]*fRWGS[0] + NArray[3*nLM + 1]*fRWGS[1] + NArray[3*nLM + 2]*fRWGS[2];

     // note: what we really want is the dot product with the complex conjugate
     // of M and N; since fRWG is real-valued we can just take the complex 
     // conjugate of the dot project with M and N. 
     f[nf++] = real( MDotF );
     f[nf++] = -imag( MDotF );
     f[nf++] = real( NDotF );
     f[nf++] = -imag( NDotF );

   };

}

/***************************************************************/
/* get the projections of a single RWG basis function onto the */
/* M and N spherical waves, up to a maximum l-value of lMax.   */
/*                                                             */
/* Workspace is a caller-allocated workspace with enough room  */
/* to store at least 10*(lMax+1)*(lMax+1) cdoubles.            */
/*                                                             */
/* The MProjection and NProjection output buffers are caller-  */
/* allocated arrays which each must have enough room to store  */
/* (lMax+1)*(lMax+1) cdoubles. On return, MProjection[Alpha]   */
/* is the projection of the given basis function onto the      */
/* M-type spherical wave with spherical wave indices l,m       */
/* such that Alpha=l^2 + l + m.                                */
/***************************************************************/
void GetMNProjections(RWGSurface *S, int ne, cdouble k, int lMax,
                      cdouble *Workspace, 
                      cdouble *MProjection, cdouble *NProjection)
{
  /* extract edge vertices */
  RWGEdge *E   = S->Edges[ne];

  double *QP   = S->Vertices + 3*(E->iQP);
  double *V1   = S->Vertices + 3*(E->iV1);
  double *V2   = S->Vertices + 3*(E->iV2);
  double *QM   = (E->iQM == -1) ? 0 : S->Vertices + 3*(E->iQM);

  double PArea = S->Panels[E->iPPanel]->Area;
  double MArea = (E->iQM == -1) ? 0.0 : S->Panels[E->iMPanel]->Area;

  /* set up data structure passed to GMNPIntegrand.           */
  GMNPData MyData, *Data=&MyData;
  Data->k=k;
  Data->lMax=lMax;

  // the Workspace field is used as follows:
  // first 3*NumLMs cdoubles: passed as the MArray scratch buffer 
  //                          required by the GMNPIntegrand routine
  //  next 3*NumLMs cdoubles: passed as the NArray scratch buffer 
  //                          required by the GMNPIntegrand routine
  //  next 2*NumLMs cdoubles: store the integrand vector resulting 
  //                          from the call to GMNPIntegrand for the 
  //                          positive triangle 
  //  next 2*NumLMs cdoubles: store the integrand vector resulting 
  //                          from the call to GMNPIntegrand for the 
  //                          negative triangle 
  // (obviously if memory were tight some of the above could be 
  //  consolidated to cut down on memory requirements)
  int NumLMs = (lMax+1)*(lMax+1);
  Data->MArray = Workspace + 0 ;
  Data->NArray = Workspace + 3*NumLMs;
  cdouble *IP  = Workspace + 6*NumLMs;
  cdouble *IM  = Workspace + 8*NumLMs;

  /* contribution of positive panel */
  Data->Q=QP;
  Data->RWGPreFac = E->Length / (2.0*PArea);
  TriIntFixed(GMNPIntegrand, 4*NumLMs, (void *)Data, QP, V1, V2, 25, (double *)IP);

  /* contribution of negative panel if present */
  if (QM)
   { Data->Q=QM;
     Data->RWGPreFac = E->Length / (2.0*MArea);
     TriIntFixed(GMNPIntegrand, 4*NumLMs, (void *)Data, V1, V2, QM, 25, (double *)IM);
   } 
  else
   {  
     memset(IM, 0, 2*NumLMs*sizeof(cdouble));
     // TODO: add line-charge contributions
   };

  for(int nLM=0; nLM<NumLMs; nLM++)
   { MProjection[nLM] = IP[2*nLM+0] - IM[2*nLM+0];
     NProjection[nLM] = IP[2*nLM+1] - IM[2*nLM+1];
   };

}

/***************************************************************/
/* thread routine for GetSphericalMoments **********************/
/***************************************************************/
typedef struct ThreadData
 { 
   int nt, NumTasks;

   RWGSurface *S;
   cdouble k;
   int lMax;

   HVector *KN;
   int BFIndexOffset;

   cdouble *Workspace;

   cdouble *PartialAVector;
   double Sign;

 } ThreadData;


void *GSM_Thread(void *data)
{ 
  ThreadData *TD          = (ThreadData *)data;
  RWGSurface *S           = TD->S;
  cdouble k               = TD->k;
  int lMax                = TD->lMax;
  HVector *KN             = TD->KN; 
  int BFIndexOffset       = TD->BFIndexOffset;
  cdouble *Workspace      = TD->Workspace;
  cdouble *PartialAVector = TD->PartialAVector;
  double Sign             = TD->Sign;
  
#ifdef USE_PTHREAD
  SetCPUAffinity(TD->nt);
#endif

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble KAlpha, NAlpha;
  int NumLMs = (lMax+1)*(lMax+1); 
  int NumMoments=2*NumLMs;
  cdouble *MProjection = Workspace + 0;
  cdouble *NProjection = Workspace + 2*NumLMs; 
  int nt=0, nbf, ne;
  memset(PartialAVector, 0, NumMoments * sizeof(cdouble));
  for(nbf=BFIndexOffset, ne=0; ne<S->NumEdges; ne++)
    { 
      nt++;
      if (nt==TD->NumTasks) nt=0;
      if (nt!=TD->nt) continue;

      LogPercent(ne, S->NumEdges);

      /***************************************************************/
      /* get the projections of the basis function onto the M and N  */
      /* spherical waves                                             */
      /***************************************************************/
      GetMNProjections(S, ne, k, lMax, Workspace, MProjection, NProjection);

      /***************************************************************/
      /* add the contributions of the electric and magnetic currents */
      /* described by this basis function to the a^{M,E} moments     */
      /***************************************************************/
      if ( S->IsPEC )
       { 
         KAlpha = Sign*KN->GetEntry( BFIndexOffset + ne );
         NAlpha = 0.0;
       }
      else 
       { KAlpha =       Sign*KN->GetEntry( BFIndexOffset + 2*ne + 0 );
         NAlpha = -Sign*ZVAC*KN->GetEntry( BFIndexOffset + 2*ne + 1 );
       };

      for(int nLM=0; nLM<NumLMs; nLM++)
       { 
         // contributions to a^M moment 
         PartialAVector[2*nLM + 0] += -k*k*( ZVAC*MProjection[nLM]*KAlpha
                                                 -NProjection[nLM]*NAlpha
                                           );

         // contributions to a^E moment
         PartialAVector[2*nLM + 1] += -k*k*( ZVAC*NProjection[nLM]*KAlpha
                                                 +MProjection[nLM]*NAlpha 
                                           );

       };
 
    };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HVector *GetSphericalMoments(RWGSurface *S, cdouble k, int lMax,
                             HVector *KN, int BFIndexOffset, 
                             HVector *AVector)
{ 
  
  /***************************************************************/
  /* (re)allocate the AVector as necessary ***********************/
  /***************************************************************/
  int NumLMs = (lMax+1)*(lMax+1);
  int NumMoments = 2*NumLMs; // a^E and a^M moments for each l,m
  if ( AVector && AVector->N != NumMoments )
   { Warn("wrong-size AVector passed to GetSphericalMoments (reallocating...)");
     AVector=0;
   };
  if ( AVector==0 )
   AVector=new HVector(NumMoments, LHM_COMPLEX);
  AVector->Zero();

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double Sign; 
  if(S->RegionIndices[0]==0)
   Sign=+1.0;
  else if (S->RegionIndices[1]==0)
   Sign=-1.0;
  else 
   return AVector; // currents on this surface do not contribute

  Log("Computing induced spherical moments on surface %s...",S->Label);

  /***************************************************************/
  /* set up thread data ******************************************/
  /***************************************************************/
  int NumThreads=GetNumThreads();
  
  int WorkspaceSize = 10*NumLMs;
  int PartialAVectorSize = 2*NumLMs;
  cdouble *WorkspaceBuffer = (cdouble *)mallocEC(NumThreads*WorkspaceSize*sizeof(cdouble));
  cdouble *PartialAVectorBuffer = (cdouble *)mallocEC(NumThreads*PartialAVectorSize*sizeof(cdouble));

  ThreadData *TDs = new ThreadData[NumThreads];
  int nt;
  for(nt=0; nt<NumThreads; nt++)
   { TDs[nt].nt             = nt;
     TDs[nt].NumTasks       = NumThreads;
     TDs[nt].S              = S;
     TDs[nt].k              = k;
     TDs[nt].lMax           = lMax;
     TDs[nt].KN             = KN;
     TDs[nt].BFIndexOffset  = BFIndexOffset;
     TDs[nt].Workspace      = WorkspaceBuffer + nt*WorkspaceSize;
     TDs[nt].PartialAVector = PartialAVectorBuffer + nt*PartialAVectorSize;
     TDs[nt].Sign           = Sign;
   };

  /***************************************************************/
  /* fire off the threads ****************************************/
  /***************************************************************/
#ifdef USE_PTHREAD
  pthread_t *Threads = new pthread_t[NumThreads];
  for(nt=0; nt<NumThreads; nt++)
   { 
     if (nt+1 == NumThreads)
      GSM_Thread((void *) &(TDs[nt]) );
     else
      pthread_create( &(Threads[nt]), 0, GSM_Thread, (void *)(&(TDs[nt])));
   };
  for(nt=0; nt<NumThreads-1; nt++)
   pthread_join(Threads[nt],0);
  delete[] Threads;
#else
#ifndef USE_OPENMP
 NumThreads=1;
#else
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(nt=0; nt<NumThreads; nt++)
   GSM_Thread((void *)&(TDs[nt]));
#endif

  /***************************************************************/
  /* accumulate contributions of all threads to the A vector     */
  /***************************************************************/
  AVector->Zero();
  cdouble *PAV;
  for(nt=0; nt<NumThreads; nt++)
   { PAV=PartialAVectorBuffer + nt*PartialAVectorSize;
     for(int nMoment=0; nMoment<NumMoments; nMoment++)
      AVector->AddEntry(nMoment, PAV[nMoment]);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  delete[] TDs;
  free(WorkspaceBuffer);
  free(PartialAVectorBuffer);

  return AVector;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HVector *GetSphericalMoments(RWGGeometry *G, cdouble k, int lMax,
                             HVector *KN, HVector *AVector)
{ 
  /***************************************************************/
  /* (re)allocate the AVector as necessary ***********************/
  /***************************************************************/
  int NumLMs = (lMax+1)*(lMax+1);
  int NumMoments = 2*NumLMs; // a^E and a^M moments for each l,m
  if ( AVector && AVector->N != NumMoments )
   { Warn("wrong-size AVector passed to GetSphericalMoments (reallocating...)");
     AVector=0;
   };
  if ( AVector==0 )
   AVector=new HVector(NumMoments, LHM_COMPLEX);
  AVector->Zero();

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  HVector *Scratch=new HVector(NumMoments, LHM_COMPLEX);

  AVector->Zero();
  for(int ns=0; ns<G->NumSurfaces; ns++)
   { 
     GetSphericalMoments(G->Surfaces[ns], k, lMax, KN, G->BFIndexOffset[ns], Scratch);

     for(int nm=0; nm<NumMoments; nm++)
      AVector->AddEntry(nm, Scratch->GetEntry(nm));
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  delete Scratch;
   
}
