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
 * AssembleRHSVector_Composite.cc
 *
 * homer reid  -- 10/2006 -- 6/2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>
#include <libTriInt.h>
#include <libhmat.h>

#include "libscuff.h"
#include "RWGComposite.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#ifdef USE_PTHREAD
#  include <pthread.h>
#endif

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
#define USE_OPENMP
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

namespace scuff {

/***************************************************************/
/* data structure used to pass parameters to                   */
/* InnerProductIntegrand routine.                              */
/***************************************************************/
typedef struct InnerProductIntegrandData 
 { 
   double *Q;
   double PreFac;

   IncField **PositiveIFs;
   int NPositiveIFs;
   IncField **NegativeIFs;
   int NNegativeIFs;

   int NeedHProd;

 } InnerProductIntegrandData;

/***************************************************************/
/* integrand routine passed to TriInt. this returns the dot    */
/* product of the RWG basis function at point X with the       */
/* electric field at point X, and optionally the dot product   */
/* with the magnteic field as well.                            */
/***************************************************************/
static void InnerProductIntegrand(double *X, void *opIPID, double *F)
{ 
  InnerProductIntegrandData *IPID=(InnerProductIntegrandData *)opIPID;

  /* get value of RWG basis function at X */
  double fRWG[3];
  VecSub(X,IPID->Q,fRWG);
  VecScale(fRWG,IPID->PreFac);

  /* get incident E and H fields at X */
  cdouble dEH[6], EH[6];
  memset(EH, 0, 6*sizeof(cdouble));
  int n, nif;
  for(nif=0; nif<IPID->NPositiveIFs; nif++)
   { IPID->PositiveIFs[nif]->GetFields(X,dEH);
     for(n=0; n<6; n++) 
      EH[n]+=dEH[n];
   };
  for(nif=0; nif<IPID->NNegativeIFs; nif++)
   { IPID->NegativeIFs[nif]->GetFields(X,dEH);
     for(n=0; n<6; n++) 
      EH[n]-=dEH[n];
   };
  
  /* compute dot products */
  cdouble *zF = (cdouble *)F;
  zF[0] = fRWG[0]*EH[0] + fRWG[1]*EH[1] + fRWG[2]*EH[2];
  if (IPID->NeedHProd)
   zF[1] = fRWG[0]*EH[3] + fRWG[1]*EH[4] + fRWG[2]*EH[5];
  
} 

/***************************************************************/
/* Calculate the inner product of given electric and magnteic  */
/* fields with the basis function associated with edge #ne on  */
/* the given RWGObject.                                        */
/*                                                             */
/* PositiveIFs[0..NPositiveIFs] are IncFields whose fields     */
/* contribute with a positive sign.                            */
/*                                                             */
/* NegativeIFs[0..NNegativeIFs] are IncFields whose fields     */
/* contribute with a negative sign.                            */
/*                                                             */
/* If pHProd==0, we skip the computation of the magnteic-field */
/* inner product.                                              */
/***************************************************************/
void GetInnerProducts(RWGEdge *E, RWGPanel **Panels, double *Vertices,
                      IncField **PositiveIFs, int NPositiveIFs,
                      IncField **NegativeIFs, int NNegativeIFs,
                      cdouble *pEProd, cdouble *pHProd)
{ 
  /* get edge vertices */
  //RWGEdge *E   = O->Edges[ne];
  double *QP   = Vertices + 3*(E->iQP);
  double *V1   = Vertices + 3*(E->iV1);
  double *V2   = Vertices + 3*(E->iV2);
  double PArea = Panels[E->iPPanel]->Area; 
  double MArea;

  /* set up data structure passed to InnerProductIntegrand */
  InnerProductIntegrandData MyIPID, *IPID=&MyIPID;
  IPID->PositiveIFs  = PositiveIFs;
  IPID->NPositiveIFs = NPositiveIFs;
  IPID->NegativeIFs  = NegativeIFs;
  IPID->NNegativeIFs = NNegativeIFs;
  IPID->NeedHProd    = (pHProd==NULL ? 0 : 1);
  
  int nFun = (pHProd==NULL ? 2 : 4);
  double I[4], IP[4], IM[4];

  /* integrate over positive panel */
  IPID->Q=QP;
  IPID->PreFac=E->Length / (2.0*PArea);
  TriIntFixed(InnerProductIntegrand, nFun, (void *)IPID, QP, V1, V2, 20, IP);

  /* integrate over negative panel */
  if ( E->iMPanel != -1)
   { double MArea = Panels[E->iMPanel]->Area;
     double *QM   = Vertices + 3*(E->iQM);
     IPID->Q=QM;
     IPID->PreFac=E->Length / (2.0*MArea);
     TriIntFixed(InnerProductIntegrand, nFun, (void *)IPID, V1, V2, QM, 20, IM);
   };

  /* total integral is difference between pos and neg pan integrals */
  for(int nf=0; nf<nFun; nf++)
   I[nf] = IP[nf] - IM[nf];
  
  /* return values */
  *pEProd = cdouble(I[0], I[1]);
  if (pHProd)
   *pHProd = cdouble(I[2], I[3]);

}

/***************************************************************/
/* data structure used to pass data to AssembleRHS_Thread      */
/***************************************************************/
typedef struct ThreadData
 { 
   int nt, nThread;

   RWGComposite *C;
   IncField *IF;
   int NIF;
   HVector *RHS;

 } ThreadData;

/***************************************************************/
/* AssembleRHS_Thread        ***********************************/
/***************************************************************/
static void *AssembleRHS_Thread(void *data)
{ 
  ThreadData *TD    = (ThreadData *)data;

#ifdef USE_PTHREAD
  SetCPUAffinity(TD->nt);
#endif

  /***************************************************************/
  /* extract fields from thread data structure *******************/
  /***************************************************************/
  RWGComposite *C  = TD->C;
  IncField *IFList = TD->IF;
  int NIF          = TD->NIF;
  HVector *RHS     = TD->RHS;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  IncField **PositiveIFs = new IncField *[NIF];
  IncField **NegativeIFs = new IncField *[NIF];
  int NPositiveIFs;
  int NNegativeIFs;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  PartialSurface *PS;
  RWGEdge *E;
  int nps, nte, Offset;
  cdouble EProd, HProd;
  IncField *IF;
  int nt=0;
  for(nps=0; nps<C->NumPartialSurfaces; nps++)
   { 
     PS=C->PartialSurfaces[nps];
     Offset = C->BFIndexOffset[nps];

     /*--------------------------------------------------------------*/
     /*- go through the list of IncField structures to identify the -*/
     /*- sublists of IFs that contribute with a plus sign and with  -*/
     /*- a minus sign to the RHS vector for this region.            -*/
     /*--------------------------------------------------------------*/
     if ( (C->PSSubRegions[2*nps+0]!=0) && (C->PSSubRegions[2*nps+1]!=0) ) 
      { Log("partial surface %i does not contribute to RHS; skipping",nps);
        continue;
      };
     for(NPositiveIFs=NNegativeIFs=0, IF=IFList; IF; IF=IF->Next)
      { 
        NegativeIFs[NNegativeIFs++] = IF;
      };

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     for(nte=0; nte<PS->NumTotalEdges; nte++)
      { 
        nt++;
        if (nt==TD->nThread) nt=0;
        if (nt!=TD->nt) continue;

        if ( nte<PS->NumEdges )
         E=PS->Edges[nte];
        else
         E=PS->HEdges[ nte - PS->NumEdges ];

        GetInnerProducts(E, PS->Panels, C->Vertices,
                         PositiveIFs, NPositiveIFs,
                         NegativeIFs, NNegativeIFs,
                         &EProd, &HProd );

        RHS->SetEntry(Offset + 2*nte+0, EProd / ZVAC);
        RHS->SetEntry(Offset + 2*nte+1, HProd);

      }; // for ne=...

   }; // for no=...

  delete[] PositiveIFs;
  delete[] NegativeIFs;

  return 0;
 
}

/***************************************************************/
/* Assemble the RHS vector.  ***********************************/
/***************************************************************/
HVector AssembleRHSVector_Composite(RWGComposite *C,
                                                  cdouble Omega,
                                                  IncField *IF,
                                                  HVector *RHS,
                                                  int nThread)
{ 
  RHS->Zero();
   
  if (nThread <= 0) 
   nThread = GetNumThreads();

  //int nt, NIF=UpdateIncFields(IF, Omega);
  int nt, NIF=1; // FIXME 

  ThreadData ReferenceTD;
 // ReferenceTD.G=this;
  ReferenceTD.C=C;
  ReferenceTD.IF=IF;
  ReferenceTD.NIF=NIF;
  ReferenceTD.RHS=RHS;

#ifdef USE_PTHREAD
  ThreadData *TDs = new ThreadData[nThread], *TD;
  pthread_t *Threads = new pthread_t[nThread];
  for(nt=0; nt<nThread; nt++)
   { 
     TD=&(TDs[nt]);
     memcpy(TD, &ReferenceTD, sizeof(ThreadData));
     TD->nt=nt;
     TD->nThread=nThread;

     if (nt+1 == nThread)
       AssembleRHS_Thread((void *)TD);
     else
       pthread_create( &(Threads[nt]), 0, AssembleRHS_Thread, (void *)TD);
   };
  for(nt=0; nt<nThread-1; nt++)
   pthread_join(Threads[nt],0);

  delete[] Threads;
  delete[] TDs;

#else

  int nTask;

#ifndef USE_OPENMP
  nThread=nTask=1;
#else
  nTask=nThread*100;
#pragma omp parallel for schedule(dynamic,1), num_threads(nThread)
#endif
  for(nt=0; nt<nTask; nt++)
   { 
     ThreadData TD1;
     memcpy(&TD1, &ReferenceTD, sizeof(ThreadData));
     TD1.nt=nt;
     TD1.nThread=nTask;
     AssembleRHS_Thread((void *)&TD1);
   };
#endif

  return RHS;
}

} // namespace scuff
