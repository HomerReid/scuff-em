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
int RWGGeometry::CountCommonRegions(int nsa, int nsb, 
                                    double Signs[2], cdouble Eps[2], cdouble Mu[2])
                                  
{
  RWGSurface  *Sa = Surfaces[nsa];
  RWGSurface  *Sb = Surfaces[nsb];

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

  Eps[0] = EpsTF[ CommonRegions[0] ];
  Mu[0]  = MuTF[ CommonRegions[0] ];

  if (NumCommonRegions==1)
   return 1;

  Eps[1] = EpsTF[ CommonRegions[1] ];
  Mu[1]  = MuTF[ CommonRegions[1] ];

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
          NumCommonRegions=G->CountCommonRegions(ncs, nrs, Signs, EpsAB, MuAB);
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

} // namespace scuff
