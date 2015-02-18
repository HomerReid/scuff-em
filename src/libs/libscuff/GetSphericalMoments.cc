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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

using namespace scuff;

/***************************************************************/
/* get the projections of a single RWG basis function onto the */
/* M and N spherical waves, up to a maximum l-value of lMax.   */
/*                                                             */
/* MArray and NArray are caller-allocated workspace arrays with*/
/* enough room to store at least 3*(lMax+1)*(lMax+1) cdoubles  */
/* each.                                                       */
/*                                                             */
/* Workspace is a caller-allocated array with enough space to  */
/* store at least 4*(lMax+2) doubles.                          */
/*                                                             */
/* The MProjection and NProjection output buffers are caller-  */
/* allocated arrays which each must have enough room to store  */
/* (lMax+1)*(lMax+1) cdoubles. On return, MProjection[Alpha]   */
/* is the projection of the given basis function onto the      */
/* M-type spherical wave with spherical wave indices l,m       */
/* such that Alpha=l^2 + l + m.                                */
/***************************************************************/
void GetMNProjections(RWGSurface *S, int ne, cdouble k, int lMax,
                      cdouble *MArray, cdouble *NArray,
                      double *Workspace,
                      cdouble *MProjections, cdouble *NProjections)
{
  int NAlpha = (lMax+1)*(lMax+1);
  memset(MProjections, 0, NAlpha*sizeof(cdouble));
  memset(NProjections, 0, NAlpha*sizeof(cdouble));

  /***************************************************************/
  /* choose triangle cubature rule (order fixed at 20 for now)   */
  /***************************************************************/
  int NumPts;
  double *TCR = GetTCR(20, &NumPts);

  /***************************************************************/
  /* preliminary geometry setup **********************************/
  /***************************************************************/
  RWGEdge *E    = S->Edges[ne];
  double Length = E->Length;
  double *QP    = S->Vertices + 3*(E->iQP);
  double *V1    = S->Vertices + 3*(E->iV1);
  double *V2    = S->Vertices + 3*(E->iV2);
  double *QM    = (E->iQM == -1) ? 0 : S->Vertices + 3*(E->iQM);

  double AP[3], BP[3], AM[3], BM[3];
  for(int Mu=0; Mu<3; Mu++)
   { AP[Mu] = V1[Mu] - QP[Mu];
     BP[Mu] = V2[Mu] - QP[Mu];
     if (QM)
      { AM[Mu] = V1[Mu] - QM[Mu];
        BM[Mu] = V2[Mu] - QM[Mu];
      };
   };

  /***************************************************************/
  /* loop over quadrature points *********************************/
  /***************************************************************/
  for(int np=0, ncp=0; np<NumPts; np++)
   { 
     double u=TCR[ncp++];
     double v=TCR[ncp++];
     double w=TCR[ncp++];

     double XmQ[3], X[3], FS[3];
     double r, Theta, Phi;

     /*--------------------------------------------------------------*/
     /*- contribution of positive panel -----------------------------*/
     /*--------------------------------------------------------------*/
     for(int Mu=0; Mu<3; Mu++)
      { XmQ[Mu] = u*AP[Mu] + v*BP[Mu];
          X[Mu] = XmQ[Mu] + QP[Mu];
      };
     CoordinateC2S(X,&r,&Theta,&Phi);
     VectorC2S(Theta,Phi,XmQ,FS);
     GetMNlmArray(lMax, k, r, Theta, Phi, LS_REGULAR, MArray, NArray, Workspace);
     for(int Alpha=0; Alpha<NAlpha; Alpha++)
      { MProjections[Alpha] += w*(  FS[0]*MArray[3*Alpha+0] 
                                   +FS[1]*MArray[3*Alpha+1]
                                   +FS[2]*MArray[3*Alpha+2]
                                 );
        NProjections[Alpha] += w*(  FS[0]*NArray[3*Alpha+0] 
                                   +FS[1]*NArray[3*Alpha+1]
                                   +FS[2]*NArray[3*Alpha+2]
                                 );
      };

     if (QM==0) continue;

     /*--------------------------------------------------------------*/
     /*- contribution of negative panel if present ------------------*/
     /*--------------------------------------------------------------*/
     for(int Mu=0; Mu<3; Mu++)
      { XmQ[Mu] = u*AM[Mu] + v*BM[Mu];
          X[Mu] = XmQ[Mu] + QM[Mu];
      };
     CoordinateC2S(X,&r,&Theta,&Phi);
     VectorC2S(Theta,Phi,XmQ,FS);
     GetMNlmArray(lMax, k, r, Theta, Phi, LS_REGULAR, MArray, NArray, Workspace);
     for(int Alpha=0; Alpha<NAlpha; Alpha++)
      { MProjections[Alpha] -= w*(  FS[0]*MArray[3*Alpha+0] 
                                   +FS[1]*MArray[3*Alpha+1]
                                   +FS[2]*MArray[3*Alpha+2]
                                 );
        NProjections[Alpha] -= w*(  FS[0]*NArray[3*Alpha+0] 
                                   +FS[1]*NArray[3*Alpha+1]
                                   +FS[2]*NArray[3*Alpha+2]
                                 );
      };

   }; //for(int np=ncp=0; np<NumPts; np++)
  
  /***************************************************************/
  /* note: what we really want is the dot product with the       */
  /* complex conjugate of M and N; since fRWG is real-valued we  */
  /* can just take the complex conjugate of the dot project with */
  /* M and N.                                                    */
  /***************************************************************/
  for(int Alpha=0; Alpha<NAlpha; Alpha++)
   { MProjections[Alpha] = Length * conj(MProjections[Alpha]);
     NProjections[Alpha] = Length * conj(NProjections[Alpha]);
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
   cdouble *Workspace1;
   double *Workspace2;
   cdouble *PartialMomentVector;
   double Sign;

 } ThreadData;

void *GSM_Thread(void *data)
{ 
  ThreadData *TD               = (ThreadData *)data;
  RWGSurface *S                = TD->S;
  cdouble k                    = TD->k;
  int lMax                     = TD->lMax;
  HVector *KN                  = TD->KN; 
  int BFIndexOffset            = TD->BFIndexOffset;
  cdouble *Workspace1          = TD->Workspace1;
  double *Workspace2           = TD->Workspace2;
  cdouble *PartialMomentVector = TD->PartialMomentVector;
  double Sign                  = TD->Sign;
  
#ifdef USE_PTHREAD
  SetCPUAffinity(TD->nt);
#endif

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NumLMs = (lMax+1)*(lMax+1);
  cdouble *MArray       = Workspace1 + 0*NumLMs;
  cdouble *NArray       = Workspace1 + 3*NumLMs;
  cdouble *MProjections = Workspace1 + 6*NumLMs;
  cdouble *NProjections = Workspace1 + 7*NumLMs;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int nt=0;
  int NumMoments=2*NumLMs;
  memset(PartialMomentVector, 0, NumMoments * sizeof(cdouble));
  for(int ne=0; ne<S->NumEdges; ne++)
    { 
      nt++;
      if (nt==TD->NumTasks) nt=0;
      if (nt!=TD->nt) continue;

      LogPercent(ne, S->NumEdges);

      /***************************************************************/
      /* get the projections of the basis function onto the M and N  */
      /* spherical waves                                             */
      /***************************************************************/
      GetMNProjections(S, ne, k, lMax, MArray, NArray, Workspace2, 
                       MProjections, NProjections);

      /***************************************************************/
      /* add contributions of this RWG function to the M-type and    */
      /* N-type moments                                              */
      /***************************************************************/
      cdouble KAlpha, NAlpha;
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
         // contributions to M-type moment
         PartialMomentVector[2*nLM + 0] += -k*k*( ZVAC*MProjections[nLM]*KAlpha
                                                      -NProjections[nLM]*NAlpha
                                           );

         // contributions to N-type moment
         PartialMomentVector[2*nLM + 1] += -k*k*( ZVAC*NProjections[nLM]*KAlpha
                                                      +MProjections[nLM]*NAlpha
                                           );
       };
 
    };

  return 0;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HVector *GetSphericalMoments(RWGSurface *S, cdouble k, int lMax,
                             HVector *KN, int BFIndexOffset, 
                             HVector *MomentVector)
{ 
  
  /***************************************************************/
  /* (re)allocate the MomentVector as necessary ***********************/
  /***************************************************************/
  int NumLMs = (lMax+1)*(lMax+1);
  int NumMoments = 2*NumLMs; // M-type and N-type moments for each l,m
  if ( MomentVector && MomentVector->N != NumMoments )
   { Warn("wrong-size MomentVector passed to GetSphericalMoments (reallocating...)");
     MomentVector=0;
   };
  if ( MomentVector==0 )
   MomentVector=new HVector(NumMoments, LHM_COMPLEX);
  MomentVector->Zero();

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double Sign; 
  if(S->RegionIndices[0]==0)
   Sign=+1.0;
  else if (S->RegionIndices[1]==0)
   Sign=-1.0;
  else 
   return MomentVector; // currents on this surface do not contribute

  Log("Computing induced spherical moments on surface %s...",S->Label);

  /***************************************************************/
  /* allocate data structures for individual calls to ************/
  /***************************************************************/
  int NumThreads=GetNumThreads();

/***************************************************************/
/* 20140119 I have to disable multithreading due to what I     */
/* believe to be a bug in the f2c-translated version of the    */
/* AmosBessel routines                                         */
/***************************************************************/
char *s=getenv("SCUFF_SPHERICAL_SINGLETHREADED");
if ( s && s[0]=='1' )
 { Log("Using single-threading for spherical wave routines.");
   NumThreads=1;
 };
  
  int Workspace1Size = 8*NumLMs;
  int Workspace2Size = 4*(lMax+2);
  int PartialMomentVectorSize = NumMoments;
  cdouble *Workspace1Buffer = (cdouble *)mallocEC(NumThreads*Workspace1Size*sizeof(cdouble));
  double *Workspace2Buffer  = (double *)mallocEC(NumThreads*Workspace2Size*sizeof(double));
  cdouble *PartialMomentVectorBuffer = (cdouble *)mallocEC(NumThreads*PartialMomentVectorSize*sizeof(cdouble));

  ThreadData *TDs = new ThreadData[NumThreads];
  for(int nt=0; nt<NumThreads; nt++)
   { TDs[nt].nt                  = nt;
     TDs[nt].NumTasks            = NumThreads;
     TDs[nt].S                   = S;
     TDs[nt].k                   = k;
     TDs[nt].lMax                = lMax;
     TDs[nt].KN                  = KN;
     TDs[nt].BFIndexOffset       = BFIndexOffset;
     TDs[nt].Workspace1          = Workspace1Buffer + nt*Workspace1Size;
     TDs[nt].Workspace2          = Workspace2Buffer + nt*Workspace2Size;
     TDs[nt].PartialMomentVector = PartialMomentVectorBuffer + nt*PartialMomentVectorSize;
     TDs[nt].Sign                = Sign;
   }; 

  /***************************************************************/
  /* fire off the threads ****************************************/
  /***************************************************************/
   GSM_Thread((void *)&(TDs[0]));
#if 0
#ifdef USE_PTHREAD
  pthread_t *Threads = new pthread_t[NumThreads];
  for(int nt=0; nt<NumThreads; nt++)
   { 
     if (nt+1 == NumThreads)
      GSM_Thread((void *) &(TDs[nt]) );
     else
      pthread_create( &(Threads[nt]), 0, GSM_Thread, (void *)(&(TDs[nt])));
   };
  for(int nt=0; nt<NumThreads-1; nt++)
   pthread_join(Threads[nt],0);
  delete[] Threads;
#else
#ifndef USE_OPENMP
 NumThreads=1;
#else
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int nt=0; nt<NumThreads; nt++)
   GSM_Thread((void *)&(TDs[nt]));
#endif
#endif

  /***************************************************************/
  /* accumulate contributions of all threads                     */
  /***************************************************************/
  MomentVector->Zero();
  cdouble *PMV;
  for(int nt=0; nt<NumThreads; nt++)
   { PMV=PartialMomentVectorBuffer + nt*PartialMomentVectorSize;
     for(int nMoment=0; nMoment<NumMoments; nMoment++)
      MomentVector->AddEntry(nMoment, PMV[nMoment]);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  delete[] TDs;
  free(Workspace1Buffer);
  free(Workspace2Buffer);
  free(PartialMomentVectorBuffer);

  return MomentVector;

}

/***************************************************************/
/* get the spherical moments for an entire geometry by summing */
/* the moments for each individual surface                     */
/***************************************************************/
HVector *GetSphericalMoments(RWGGeometry *G, cdouble k, int lMax,
                             HVector *KN, HVector *MomentVector)
{ 
  /***************************************************************/
  /* (re)allocate the MomentVector as necessary ***********************/
  /***************************************************************/
  int NumLMs = (lMax+1)*(lMax+1);
  int NumMoments = 2*NumLMs; // M-type and N-type moments for each l,m
  if ( MomentVector && MomentVector->N != NumMoments )
   { Warn("wrong-size MomentVector passed to GetSphericalMoments (reallocating...)");
     MomentVector=0;
   };
  if ( MomentVector==0 )
   MomentVector=new HVector(NumMoments, LHM_COMPLEX);
  MomentVector->Zero();

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  HVector *Scratch=new HVector(NumMoments, LHM_COMPLEX);

  MomentVector->Zero();
  for(int ns=0; ns<G->NumSurfaces; ns++)
   { 
     GetSphericalMoments(G->Surfaces[ns], k, lMax, KN, G->BFIndexOffset[ns], Scratch);

     for(int nm=0; nm<NumMoments; nm++)
      MomentVector->AddEntry(nm, Scratch->GetEntry(nm));
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  delete Scratch;
  return MomentVector;
   
}
