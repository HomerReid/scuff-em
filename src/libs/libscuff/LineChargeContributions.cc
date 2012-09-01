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
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libhmat.h>
#include <libhrutil.h>
#include <omp.h>

#include "libscuff.h"
#include "libscuffInternals.h"

#include "cmatheval.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#ifdef USE_PTHREAD
#  include <pthread.h>
#endif

namespace scuff {

cdouble GetLineChargeInteraction(RWGSurface *S1, int ne1, RWGSurface *S2, int ne2, cdouble k);

#define II cdouble(0,1)

/***************************************************************/
/* Figure out if the RWGSurfaces with indices nsa and nsb have */
/* zero, one, or two regions in common.                        */
/***************************************************************/
int CountCommonRegions2(RWGGeometry *G, int nsa, int nsb, 
                        double Signs[2], cdouble Eps[2], cdouble Mu[2])
                                  
{
  RWGSurface  *Sa = G->Surfaces[nsa];
  RWGSurface  *Sb = G->Surfaces[nsb];

  int CommonRegions[2], NumCommonRegions=0;
  if ( Sa->RegionIndices[0] == Sb->RegionIndices[0] )
   { CommonRegions[NumCommonRegions] = Sa->RegionIndices[0];
     Signs[NumCommonRegions]=+1.0;
     NumCommonRegions++;
   }
  else if ( Sa->RegionIndices[0] == Sb->RegionIndices[1] )
   { CommonRegions[NumCommonRegions] = Sa->RegionIndices[0];
     NumCommonRegions++;
   }
  if ( Sa->RegionIndices[1] == Sb->RegionIndices[0] )
   { CommonRegions[NumCommonRegions] = Sa->RegionIndices[1];
     Signs[NumCommonRegions]=-1.0;
     NumCommonRegions++;
   }
  else if ( !Sa->IsPEC && !Sb->IsPEC && Sa->RegionIndices[1] == Sb->RegionIndices[1] )
   { CommonRegions[NumCommonRegions] = Sa->RegionIndices[1];
     Signs[NumCommonRegions]=+1.0;
     NumCommonRegions++;
   };

  if (NumCommonRegions==0)
   return 0;

  Eps[0] = G->EpsTF[ CommonRegions[0] ];
  Mu[0]  = G->MuTF[ CommonRegions[0] ];

  if (NumCommonRegions==1)
   return 1;

  Eps[1] = G->EpsTF[ CommonRegions[1] ];
  Mu[1]  = G->MuTF[ CommonRegions[1] ];

  return 2;
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct ThreadData
 { 
   RWGGeometry *G;
   cdouble Omega;
   HMatrix *M;
   int nt, NumTasks;

 } ThreadData;

/***************************************************************/
/***************************************************************/
/***************************************************************/
void *ALCCThread(void *data)
{
  /***************************************************************/
  /* local copies of fields in argument structure                */
  /***************************************************************/
  ThreadData *TD    = (ThreadData *)data;
  RWGGeometry *G    = TD->G;
  cdouble Omega     = TD->Omega;
  HMatrix *M        = TD->M;

  /***************************************************************/
  /* loop over all columns of the BEM matrix (that is, all edges */
  /* on all surfaces); for each half-RWG basis function h, add   */
  /* the contributions of the line-charge density associated with*/
  /* h to the corresponding column of the matrix.                */
  /***************************************************************/
  int ncs, nce; // num column surface, num column edge 
  int nrs, nre; // num column surface, num column edge 
  RWGSurface *CS, *RS; // column surface, row surface 
  int ColOfs, RowOfs;  // column offset, row offset 
  int NumCommonRegions;
  double Signs[2];
  cdouble EpsAB[2], MuAB[2];
  cdouble kA, PreFac1A, PreFac2A, PreFac3A;
  cdouble kB, PreFac1B, PreFac2B, PreFac3B;
  cdouble DeltaG1, DeltaG2;
  int nt=0;
  for(ncs=0; ncs<G->NumSurfaces; ncs++)
   for(CS=G->Surfaces[ncs], nce=0; nce<CS->NumEdges; nce++)
    { 
      if ( CS->Edges[nce]->iQM != 1 ) continue;
      ColOfs = G->BFIndexOffset[ncs];

      for(nrs=0; nrs<G->NumSurfaces; nrs++)
       { 
          // figure out if the two surfaces interact or not 
          NumCommonRegions=CountCommonRegions2(G,ncs, nrs, Signs, EpsAB, MuAB);
          if (NumCommonRegions==0) continue;
 
          kA=csqrt2(EpsAB[0]*MuAB[0])*Omega;
          PreFac1A = Signs[0]*II*MuAB[0]*Omega;
          PreFac2A = -Signs[0]*II*kA;
          PreFac3A = -1.0*Signs[0]*II*EpsAB[0]*Omega;

          if (NumCommonRegions==2)
           { kB=csqrt2(EpsAB[1]*MuAB[1])*Omega;
             PreFac1B = Signs[1]*II*MuAB[1]*Omega;
             PreFac2B = -Signs[1]*II*kB;
             PreFac3B = -1.0*Signs[1]*II*EpsAB[1]*Omega;
           }
          else 
           PreFac3B=0.0;

          RowOfs = G->BFIndexOffset[nrs];

          // loop over all BFs on surface #nrs to add the contributions 
          // of the given half-RWG basis function
          for(RS=G->Surfaces[nrs], nre=0; nre<RS->NumEdges; nre++)
           { 
             nt++;
             if (nt==TD->NumTasks) nt=0;
             if (nt!=TD->nt) continue;

             DeltaG1 = GetLineChargeInteraction(CS, nce, RS, nre, kA);

             if (NumCommonRegions==2)
              DeltaG2=GetLineChargeInteraction(CS, nce, RS, nre, kB);
             else
              DeltaG2=0.0;

             if ( CS->IsPEC && RS->IsPEC )
              { 
                M->AddEntry(RowOfs + nre, ColOfs + nce, PreFac1A*DeltaG1);
              }
             else if ( !(CS->IsPEC) && !(RS->IsPEC) )
              { 
                M->AddEntry(RowOfs + 2*nre,     ColOfs + 2*nce,     PreFac1A*DeltaG1 + PreFac1B*DeltaG2 );
                M->AddEntry(RowOfs + 2*nre + 1, ColOfs + 2*nce + 1, PreFac3A*DeltaG1 + PreFac3B*DeltaG2 );
              }
             else
              ErrExit("%s:%i: half-RWG functions not yet supported for PEC/non-PEC hybrids",__FILE__,__LINE__);

           };

       };

    }; // for(ncs=0; ...)

  return 0;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RWGGeometry::AddLineChargeContributionsToBEMMatrix(cdouble Omega, HMatrix *M)
{ 
  /***************************************************************/
  /* fire off threads ********************************************/
  /***************************************************************/
  int nt, NumTasks, NumThreads=GetNumThreads();

#ifdef USE_PTHREAD
  ThreadData *TDs = new ThreadData[nThread], *TD;
  pthread_t *Threads = new pthread_t[nThread];
  for(nt=0; nt<nThread; nt++)
   { 
     TD=&(TDs[nt]);
     TD->nt=nt;
     TD->NumTasks=NumThreads;
     TD->G=this;
     TD->Omega=Omega;
     TD->M=M;
     if (nt+1 == nThread)
      ALCCThread((void *)TD);
     else
      pthread_create( &(Threads[nt]), 0, ALCCThread, (void *)TD);
   }
  for(nt=0; nt<nThread-1; nt++)
   pthread_join(Threads[nt],0);
  delete[] Threads;
  delete[] TDs;

#else 
#ifndef USE_OPENMP
  NumTasks=NumThreads=1;
#else
  NumTasks=NumThreads*100;
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(nt=0; nt<NumThreads*100; nt++)
   { 
     ThreadData TD1;
     TD1.nt=nt;
     TD1.NumTasks=NumThreads*100;
     TD1.G=this;
     TD1.Omega=Omega;
     TD1.M=M;
     ALCCThread((void *)&TD1);
   };
#endif

}

/***************************************************************/
/* compute the extra contribution to the gradient of the       */
/* reduced scalar potential from the constant line charge on   */
/* the unmatched edge of a half-RWG basis function.            */
/***************************************************************/
static double x10[]=
 { 1.304673574141418e-02, 6.746831665550773e-02, 1.602952158504878e-01,
   2.833023029353764e-01, 4.255628305091844e-01, 5.744371694908156e-01,
   7.166976970646236e-01, 8.397047841495122e-01, 9.325316833444923e-01,
   9.869532642585859e-01
 };

static double w10[]=
 { 3.333567215434143e-02, 7.472567457529027e-02, 1.095431812579910e-01,
   1.346333596549959e-01, 1.477621123573765e-01, 1.477621123573765e-01,
   1.346333596549959e-01, 1.095431812579910e-01, 7.472567457529027e-02,
   3.333567215434143e-02
 };

void GetEdgeContributionToGradp(const double *X0, double *V1, double *V2,
                                cdouble K, cdouble *Gradp)
{
  int np;
  double u, w, X0mV1[3], V2mV1[3], X0mX[3], r;
  cdouble Psi;

  // we do a 10-point gauss-legendre quadrature over the edge
  int NumPts=10;
  double *QRX=x10, *QRW=w10;

  VecSub(V2, V1, V2mV1);
  VecSub(X0, V1, X0mV1);
  Gradp[0]=Gradp[1]=Gradp[2]=0.0;

  for(np=0; np<NumPts; np++)
   { 
     u=QRX[np]; w=QRW[np];

     X0mX[0] = X0mV1[0] - u*V2mV1[0];
     X0mX[1] = X0mV1[1] - u*V2mV1[1];
     X0mX[2] = X0mV1[2] - u*V2mV1[2];

     r=VecNorm(X0mX);

     Psi = (II*K - 1.0/r) * exp(II*K*r) / (4.0*M_PI*r*r);

     Gradp[0] += w*X0mX[0]*Psi;
     Gradp[1] += w*X0mX[1]*Psi;
     Gradp[2] += w*X0mX[2]*Psi;
     
   };

  cdouble PreFac = -1.0*VecNorm(V2mV1);
  Gradp[0] *= PreFac;
  Gradp[1] *= PreFac;
  Gradp[2] *= PreFac;
 
}

} // namespace scuff
