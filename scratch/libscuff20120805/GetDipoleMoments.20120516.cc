/*
 * GetDipoleMoments.cc  -- libscuff class methods for computing electric 
 *                      -- and magnetic dipole moments induced on 
 *                      -- scattering geometry by incident field
 *
 * homer reid    -- 3/2007  -- 11/2009
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <libhrutil.h>
#include <libTriInt.h>

#include "libscuff.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#ifdef USE_PTHREAD
#  include <pthread.h>
#endif

namespace scuff {

#define II cdouble(0,1)

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct ThreadData
 { 
   int nt, nThread;

   RWGGeometry *G;   
   double Frequency;
   int RealFreq;
   HVector *KN;
   cdouble (*PM)[6];

 } ThreadData;

static void *GetDipoleMoment_Thread(void *data)
{ 
  ThreadData *TD=(ThreadData *)data;

  RWGGeometry *G    = TD->G;
  double Frequency  = TD->Frequency;
  double RealFreq   = TD->RealFreq;
  HVector *KN       = TD->KN;

  RWGObject *O;
  int nt, no, ne, mu, Type, Offset;
  cdouble iw, KAlpha, NAlpha, PreFac;
  double QPmQM[3], V1pV2[3], QPmQMxV1pV2[3];
  double *DKN=0;
  cdouble *ZKN=0;

  if ( RealFreq )
   { ZKN=KN->ZV;
     iw=II*Frequency;
   }
  else
   { DKN=KN->DV;
     iw=-1.0*Frequency;
   };

  nt=0;
  for(no=0, O=G->Objects[0]; no<G->NumObjects; O=G->Objects[++no])
   { 
     Type=O->MP->Type;
     Offset=G->BFIndexOffset[no];

     for(ne=0; ne<O->NumEdges; ne++)
      { 
        nt++;
        if (nt==TD->nThread) nt=0;
        if (nt!=TD->nt) continue;

        /* extract weight of neth basis function from K vector */
        if ( Type==MP_PEC && RealFreq==1 )
         { KAlpha=ZKN[ Offset + ne ];
           NAlpha=0.0;
         }
        else if ( Type==MP_PEC && RealFreq==0 )
         { KAlpha=DKN[ Offset + ne ];
           NAlpha=0.0;
         }
        else if ( Type!=MP_PEC && RealFreq==1 )
         { KAlpha=ZKN[ Offset + 2*ne ];
           NAlpha=ZKN[ Offset + 2*ne + 1 ] / (-ZVAC);
         }
        else if ( Type!=MP_PEC && RealFreq==0 )
         { KAlpha=DKN[ Offset + 2*ne ];
           NAlpha=DKN[ Offset + 2*ne + 1 ] / (-ZVAC);
         };
      
        /* get dipole moments of neth RWG basis function */

        /* electric dipole moment of RWG function = (l / 3iw) * (Q^+ - Q^-) */
        /* magnetic dipole moment of RWG function = (l / 3) * (Q^+ - Q^-) x (V1+V2) */
        VecSub( O->Vertices + 3*O->Edges[ne]->iQP, O->Vertices + 3*O->Edges[ne]->iQM, QPmQM); 
        VecAdd( O->Vertices + 3*O->Edges[ne]->iV1, O->Vertices + 3*O->Edges[ne]->iV2, V1pV2); 
        VecCross( QPmQM, V1pV2, QPmQMxV1pV2 );

        PreFac= O->Edges[ne]->Length / 3.0;
        for(mu=0; mu<3; mu++)
         { TD->PM[no][mu]   += PreFac * (   KAlpha * QPmQM[mu] / iw
                                          + NAlpha * QPmQMxV1pV2[mu] );
           TD->PM[no][mu+3] += PreFac * (   KAlpha * QPmQMxV1pV2[mu]
                                          + NAlpha * QPmQM[mu] / iw );
         };

      }; // for (ne=0 ... 

    }; // for(no=0 ... 

  return 0;
 
}

/***************************************************************/
/* get electric and magnetic dipole moments of a current       */
/* distribution described by a set of RWG basis functions      */
/* populated with a vector of basis-function weights.          */
/*                                                             */
/* on entry, PM must be a two-dimensional array of cdoubles of */
/* the form PM[NO][6] where NO is the number of objects in     */
/* the geometry.                                               */
/* on return, PM[no][0..2] and PM[no][3..5] are the cartesian  */
/* components of the electric and magnetic dipole moments      */
/* induced on the #noth object.                                */
/***************************************************************/
void RWGGeometry::GetDipoleMoments(double Frequency, int RealFreq,
                                   HVector *KN, cdouble (*PM)[6], int nThread)
{ 
  int nt, no;

  if (nThread <= 0) nThread = GetNumThreads();

#ifdef USE_PTHREAD
  ThreadData *TDS = new ThreadData[nThread], *TD;
  pthread_t *Threads = new pthread_t[nThread];
#else
  ThreadData TD1;
#endif

  for(no=0; no<NumObjects; no++)
   memset(PM[no],0,6*sizeof(cdouble));

#ifdef USE_OPENMP
#pragma omp parallel for private(TD1), schedule(static,1), num_threads(nThread)
#endif
  /* fire off threads */
  for(nt=0; nt<nThread; nt++)
   { 
#ifdef USE_PTHREAD
     TD=&(TDS[nt]);
#else
     ThreadData *TD=&TD1;
#endif
     TD->nt=nt;
     TD->nThread=nThread;

     TD->G=this;
     TD->Frequency=Frequency;
     TD->RealFreq=RealFreq;
     TD->KN=KN;
     TD->PM=PM;

#ifdef USE_PTHREAD
     if (nt+1 == nThread)
       GetDipoleMoment_Thread((void *)TD);
     else
       pthread_create( &(Threads[nt]), 0, GetDipoleMoment_Thread, (void *)TD);
#else
     GetDipoleMoment_Thread((void *)TD);
#endif
   };

#ifdef USE_PTHREAD
  /* wait for threads to complete */
  for(nt=0; nt<nThread-1; nt++)
   pthread_join(Threads[nt],0);

  delete[] Threads;
  delete[] TDS;
#endif
}

} // namespace scuff
