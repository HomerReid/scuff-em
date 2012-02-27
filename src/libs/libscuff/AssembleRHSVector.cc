/*
 * AssembleRHS.cc  -- libscuff routines for assembling RHS vectors
 *
 * homer reid      -- 10/2006 -- 11/2009
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

#include <libhrutil.h>
#include <libTriInt.h>
#include <libhmat.h>

#include "libscuff.h"

namespace scuff {

/***************************************************************/
/* data structure used to pass parameters to                   */
/* InnerProductIntegrand routine.                              */
/***************************************************************/
typedef struct InnerProductIntegrandData 
 { 
   double *Q;
   double PreFac;
   EHFuncType EHFunc;
   void *EHFuncUD;
   int PureImagFreq;
   int NeedHProd;
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

  /* call user's incident field routine to get E and H fields at X */
  IPID->EHFunc(X,IPID->EHFuncUD,EH); 
  
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
void RWGObject::GetInnerProducts(int nbf, EHFuncType EHFunc, 
                                 void *EHFuncUD, int PureImagFreq,
                                 cdouble *pEProd, cdouble *pHProd)
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
   EHFuncType EHFunc;
   void *EHFuncUD;
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
  EHFuncType EHFunc = TD->EHFunc;
  void *EHFuncUD    = TD->EHFuncUD;

  /***************************************************************/
  /* get a pointer to the real- or complex-valued storage array  */
  /* inside the HVector                                          */
  /***************************************************************/
  double *DB;
  cdouble *ZB;
  int PureImagFreq;
  if (TD->B->RealComplex==LHM_COMPLEX)
   { ZB=TD->B->ZV;
     PureImagFreq=0;
   }
  else
   { DB=TD->B->DV; 
     PureImagFreq=1;
   };

  /*--------------------------------------------------------------*/
  /*- EXPERIMENTAL -----------------------------------------------*/
  /*--------------------------------------------------------------*/
#ifdef _GNU_SOURCE
  cpu_set_t cpuset;
  CPU_ZERO(&cpuset);
  CPU_SET(TD->nt,&cpuset);
  pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
#endif
  /*--------------------------------------------------------------*/
  /*- EXPERIMENTAL -----------------------------------------------*/
  /*--------------------------------------------------------------*/

  RWGObject *O;
  int no, ne, Type, Offset; 
  int nt=0;
  cdouble EProd, HProd, *pHProd;
  for(no=0, O=G->Objects[0]; no<G->NumObjects; O=G->Objects[++no])
   { 
     Offset=G->BFIndexOffset[no];
     Type=O->MP->Type;

     // figure out whether or not we need to compute the H product 
     if (Type==MP_PEC)
      pHProd=0;
     else
      pHProd=&HProd;

     for(ne=0; ne<O->NumEdges; ne++)
      { 
        nt++;
        if (nt==TD->nThread) nt=0;
        if (nt!=TD->nt) continue;

        O->GetInnerProducts(ne, EHFunc, EHFuncUD, PureImagFreq, &EProd, pHProd);

        /* there are four choices here based on whether we are at a general */
        /* or a pure imaginary frequency and whether the object is a        */
        /* perfect conductor (EFIE) or not (PMCHW)                          */
        if ( Type==MP_PEC && PureImagFreq==1 )
         { 
           DB[ Offset + ne ] = -1.0*(real(EProd)) / ZVAC;
         } 
        else if ( Type==MP_PEC && PureImagFreq==0 )
         { 
           ZB[ Offset + ne ] = -1.0*EProd / ZVAC;
         }
        else if ( Type!=MP_PEC && PureImagFreq==1 )
         { 
           DB[ Offset + 2*ne    ]  = -1.0*(real(EProd)) /  ZVAC;
           DB[ Offset + 2*ne + 1 ] = -1.0*(real(HProd));
         }
        else // ( Type!=MP_PEC && PureImagFreq==0 )
         { 
           ZB[ Offset + 2*ne    ]  = -1.0*EProd / ZVAC;
           ZB[ Offset + 2*ne + 1 ] = -1.0*HProd;
         };

      }; // for ne=...

   }; // for no=...

  return 0;
 
}

/***************************************************************/
/* Assemble the RHS vector.  ***********************************/
/***************************************************************/
void RWGGeometry::AssembleRHSVector(EHFuncType EHFunc, void *EHFuncUD,
                                    int nThread, HVector *B)
{ 
  int nt;

  ThreadData TDS[nThread], *TD;
  pthread_t Threads[nThread];

  if (nThread<=0)
   ErrExit("AssembleRHSVector called with nThread=%i",nThread);

  for(nt=0; nt<nThread; nt++)
   { 
     TD=&(TDS[nt]);
     TD->nt=nt;
     TD->nThread=nThread;

     TD->G=this;
     TD->EHFunc=EHFunc;
     TD->EHFuncUD=EHFuncUD;
     TD->B=B;

     TD->nThread=nThread;
     
     if (nt+1 == nThread)
       AssembleRHS_Thread((void *)TD);
     else
       pthread_create( &(Threads[nt]), 0, AssembleRHS_Thread, (void *)TD);
   };

  for(nt=0; nt<nThread-1; nt++)
   pthread_join(Threads[nt],0);

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

} // namespace scuff
