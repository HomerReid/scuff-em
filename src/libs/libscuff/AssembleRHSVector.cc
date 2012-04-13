/*
 * AssembleRHS.cc  -- libscuff routines for assembling RHS vectors
 *
 * homer reid      -- 10/2006 -- 11/2009
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>
#include <libTriInt.h>
#include <libhmat.h>

#include "libscuff.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#ifdef USE_PTHREAD
#  include <pthread.h>
#endif

namespace scuff {

/***************************************************************/
/* data structure used to pass parameters to                   */
/* InnerProductIntegrand routine.                              */
/***************************************************************/
typedef struct InnerProductIntegrandData 
 { 
   double *Q;
   double PreFac;
   IncField *IF;
   int PureImagFreq;
   int exterior_idx, interior_idx;
 } InnerProductIntegrandData;

/***************************************************************/
/* integrand routine passed to TriInt. this returns the dot    */
/* product of the RWG basis function at point X with the       */
/* electric field at point X, and optionally the dot product   */
/* with the magnetic field as well.                            */
/***************************************************************/
static void InnerProductIntegrand(double *X, void *opIPID, double *F)
{ 
  InnerProductIntegrandData *IPID=(InnerProductIntegrandData *)opIPID;
  double fRWG[3];
  cdouble EH[6];
  int nrv;

  /* get value of RWG basis function at X */
  VecSub(X,IPID->Q,fRWG);
  VecScale(fRWG,IPID->PreFac);

  /* get incident E and H fields at X */
  //IPID->EHFunc(X,IPID->EHFuncUD,EH, IPID->exterior_idx, IPID->interior_idx); 
  IF->GetFields(X,EH);
  
  /*--------------------------------------------------------------*/
  /*- now switch off to determine which return values are needed -*/
  /*--------------------------------------------------------------*/
  nrv=0;

  /* this return value is needed in all cases */ 
  F[nrv++]=    fRWG[0] * (real(EH[0]) )
             + fRWG[1] * (real(EH[1]) )
             + fRWG[2] * (real(EH[2]) );
  
  if ( !(IPID->PureImagFreq) ) 
   F[nrv++]=   fRWG[0] * (imag(EH[0]) )
             + fRWG[1] * (imag(EH[1]) )
             + fRWG[2] * (imag(EH[2]) );

  if (IPID->NeedHProd) 
   F[nrv++]=   fRWG[0] * (real(EH[3]) )
             + fRWG[1] * (real(EH[4]) )
             + fRWG[2] * (real(EH[5]) );

  if ( IPID->NeedHProd && !(IPID->PureImagFreq) ) 
   F[nrv++]=   fRWG[0] * (imag(EH[3]) )
             + fRWG[1] * (imag(EH[4]) )
             + fRWG[2] * (imag(EH[5]) );
  
} 

/***************************************************************/
/* Calculate the inner product of basis function #nbf          */
/* with the electric and magnetic fields described by EHFunc.  */
/*                                                             */
/* If RealFreq==1, then we need both the real and              */
/* imaginary parts of the results. Otherwise the imaginary     */
/* parts are zeroed out.                                       */
/*                                                             */
/* If pHProd==0, we don't compute the inner product with the   */
/* magnetic field.                                             */
/***************************************************************/
void RWGObject::GetInnerProducts(int nbf, IncField *IF,
                                 int PureImagFreq,
                                 cdouble *pEProd, cdouble *pHProd,
				 int exterior_index, int interior_index)
{ 
  double *QP, *V1, *V2, *QM;
  double PArea, MArea;
  int nf, nFun, nrv;
  RWGEdge *E;
  InnerProductIntegrandData MyIPID, *IPID=&MyIPID;
  double I[4], IP[4], IM[4];
  cdouble EProd, HProd;

  /* get edge vertices */
  E=Edges[nbf];
  QP=Vertices + 3*(E->iQP);
  V1=Vertices + 3*(E->iV1);
  V2=Vertices + 3*(E->iV2);
  QM=Vertices + 3*(E->iQM);
  PArea=Panels[E->iPPanel]->Area;
  MArea=Panels[E->iMPanel]->Area;

  /* set up data structure passed to InnerProductIntegrand */
  IPID->EHFunc=EHFunc;
  IPID->EHFuncUD=EHFuncUD;
  IPID->PureImagFreq=PureImagFreq;
  IPID->NeedHProd=0;
  IPID->exterior_idx = exterior_index;
  IPID->interior_idx = interior_index;
  nFun=1;
  if (pHProd)
   { IPID->NeedHProd=1;
     nFun=2;
   };
  if ( !(PureImagFreq) )
   { nFun*=2;
   };
 
  /* integrate over positive panel */
  IPID->Q=QP;
  IPID->PreFac=E->Length / (2.0*PArea);
  TriIntFixed(InnerProductIntegrand, nFun, (void *)IPID, QP, V1, V2, 20, IP);

  /* integrate over negative panel */
  IPID->Q=QM;
  IPID->PreFac=E->Length / (2.0*MArea);
  TriIntFixed(InnerProductIntegrand, nFun, (void *)IPID, V1, V2, QM, 20, IM);

  /* total integral is difference between pos and neg pan integrals */
  for(nf=0; nf<nFun; nf++)
   I[nf] = IP[nf] - IM[nf];
  
  /*--------------------------------------------------------------*/
  /*- now switch off to determine which return values are needed -*/
  /*--------------------------------------------------------------*/
  nrv=0;
  EProd=HProd=0.0;

  /* this return value is needed in all cases */ 
  real(EProd) = I[nrv++]; 

  if (!PureImagFreq)
   imag(EProd) = I[nrv++]; 
  if (pHProd!=0)
   real(HProd) = I[nrv++]; 
  if (pHProd!=0 && !PureImagFreq)
   imag(HProd) = I[nrv++]; 

  *pEProd=EProd;
  if (pHProd) *pHProd=HProd;
 
}

/***************************************************************/
/* data structure used to pass data to AssembleRHS_Thread      */
/***************************************************************/
typedef struct ThreadData
 { 
   int nt, nThread;

   RWGGeometry *G;
   IncField *IF;
   RWGObject *FieldSourceLocation;
   HVector *B;

 } ThreadData;

/***************************************************************/
/* AssembleRHS_Thread        ***********************************/
/***************************************************************/
void *AssembleRHS_Thread(void *data)
{ 
  /***************************************************************/
  /* extract fields from thread data structure *******************/
  /***************************************************************/
  ThreadData *TD=(ThreadData *)data;
  RWGGeometry *G    = TD->G;
  IncField    *IF   = TD->IF;
  RWGObject   *FieldSourceLocation = TD->FieldSourceLocation;

  /*--------------------------------------------------------------*/
  /*- EXPERIMENTAL -----------------------------------------------*/
  /*--------------------------------------------------------------*/
#if defined(_GNU_SOURCE) && defined(USE_PTHREAD)
  cpu_set_t cpuset;
  CPU_ZERO(&cpuset);
  CPU_SET(TD->nt,&cpuset);
  pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
#endif
  /*--------------------------------------------------------------*/
  /*- EXPERIMENTAL -----------------------------------------------*/
  /*--------------------------------------------------------------*/

  RWGObject *O;
  int no, ne, Offset; 
  int nt=0;
  cdouble EProd, HProd, *pHProd;
  for(no=0, O=G->Objects[0]; no<G->NumObjects; O=G->Objects[++no])
   { 
     // determine the sign of this object's contribution to the 
     // RHS vector, as follows: 
     // 1. if the exterior field sources lie in the medium exterior
     //    to this object, then Sign=-1
     // 2. if the exterior field sources lie in the medium interior
     //    to this object, then Sign=+1
     // 3. if neither of the above, then this object makes no 
     //    contribution to the RHS vector.
     if( O->ContainingObject == FieldSourceLocation ) 
      Sign = -1.0;
     else if (O == FieldSourceLocation )
      Sign = +1.0;
     else
      continue;

#if 0
     // find index of exterior object to O (-1 if none)
     int exterior_index = -1;
     if (O->ContainingObject) {
       for (int no2=0; no2 < G->NumObjects; ++no2)
	 if (O->ContainingObject == G->Objects[no2]) {
	   exterior_index = no2;
	   break;
	 }
       if (exterior_index == -1) ErrExit("invalid containing object");
     }
#endif

     // figure out whether or not we need to compute the H product 
     if ( O->MP->Type == MP_PEC )
      pHProd=0;
     else
      pHProd=&HProd;

     Offset=G->BFIndexOffset[no];
     for(ne=0; ne<O->NumEdges; ne++)
      { 
        nt++;
        if (nt==TD->nThread) nt=0;
        if (nt!=TD->nt) continue;

        O->GetInnerProducts(ne, IF, &EProd, pHProd);

        if ( O->MP->Type == MP_PEC )
         B->SetEntry(Offset + ne), Sign*EProd / ZVac;
        else 
         { B->SetEntry(Offset + 2*ne+0), Sign*EProd / ZVac;
           B->SetEntry(Offset + 2*ne+1), Sign*HProd;
         };

      }; // for ne=...

   }; // for no=...

  return 0;
 
}

/***************************************************************/
/* Assemble the RHS vector.  ***********************************/
/***************************************************************/
HVector *RWGGeometry::AssembleRHSVector(cdouble Omega, IncField *IF, 
                                        HVector *B, int nThread)
{ 
  if (B==0)
   B=AllocateRHSVector();
   
  if (nThread <= 0) 
   nThread = GetNumThreads();

#ifdef USE_PTHREAD
  ThreadData *TDS = new ThreadData[nThread], *TD;
  pthread_t *Threads = new pthread_t[nThread];
#else
  ThreadData TD1;
#endif

#ifdef USE_OPENMP
#pragma omp parallel for private(TD1), schedule(static,1), num_threads(nThread)
#endif
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
     TD->IF=IF;
     TD->FieldSourceLocation=FieldSourceLocation;
     TD->B=B;

     TD->nThread=nThread;
     
#ifdef USE_PTHREAD
     if (nt+1 == nThread)
       AssembleRHS_Thread((void *)TD);
     else
       pthread_create( &(Threads[nt]), 0, AssembleRHS_Thread, (void *)TD);
#else
     AssembleRHS_Thread((void *)TD);
#endif
   };

#ifdef USE_PTHREAD
  for(nt=0; nt<nThread-1; nt++)
   pthread_join(Threads[nt],0);

  delete[] Threads;
  delete[] TDS;
#endif
}

/***************************************************************/
/* Update the IncField ObjectIndex, Omega, Eps, and Mu values. */
/* (Only update ObjectIndex if ignoreOmega is true.)           */
/***************************************************************/

void RWGGeometry::UpdateIncFields(IncField *inc, cdouble Omega,
				  bool ignoreOmega) {
  for (IncField *i = inc; i; i = i->Next) {
    int io = 0; double X[3];
    if (!ignoreOmega) i->Omega = Omega;
    if ((i->Object && GetObjectByLabel(i->Object, &io))
     	|| (i->GetSourcePoint(X) && ((io = GetObjectIndex(X)) >= 0))
	|| ((io = i->ObjectIndex) >= 0 && i->ObjectIndex < NumObjects)) {
      if (Objects[io]->MP->Type == MP_PEC)
	i->ObjectIndex = -2; // disable sources inside PEC
      else {
	i->ObjectIndex = io;
	if (!ignoreOmega)
	  Objects[io]->MP->GetEpsMu(Omega, &i->Eps, &i->Mu);
      }
    }
    else if (!i->Object || i->ObjectIndex == -1 || io == -1) {
      i->ObjectIndex = -1;
      if (!ignoreOmega)
	ExteriorMP->GetEpsMu(Omega, &i->Eps, &i->Mu);
    }
  }
}

/***************************************************************/
/* Allocate an RHS vector of the appropriate size. *************/
/***************************************************************/
HVector *RWGGeometry::AllocateRHSVector(int PureImagFreq)
{ 
  HVector *V;

  if (PureImagFreq)
   V=new HVector(TotalBFs,LHM_REAL);
  else
   V=new HVector(TotalBFs,LHM_COMPLEX);

  return V;

} 
