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
   int RealFreq;
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
  
  if (IPID->RealFreq) 
   F[nrv++]=   fRWG[0] * (imag(EH[0]) )
             + fRWG[1] * (imag(EH[1]) )
             + fRWG[2] * (imag(EH[2]) );

  if (IPID->NeedHProd) 
   F[nrv++]=   fRWG[0] * (real(EH[3]) )
             + fRWG[1] * (real(EH[4]) )
             + fRWG[2] * (real(EH[5]) );

  if (IPID->NeedHProd && IPID->RealFreq) 
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
                                 void *EHFuncUD, int RealFreq,
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
  IPID->RealFreq=0;
  IPID->NeedHProd=0;
  nFun=1;
  if (pHProd)
   { IPID->NeedHProd=1;
     nFun=2;
   };
  if (RealFreq)
   { IPID->RealFreq=1;
     nFun*=2;
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

  if (RealFreq)
   imag(EProd) = I[nrv++]; 
  if (pHProd!=0)
   real(HProd) = I[nrv++]; 
  if (pHProd!=0 && RealFreq)
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
   int RealFreq; /* =1 for real frequency, 0 for imag */
   HVector *B;

 } ThreadData;

/***************************************************************/
/* AssembleRHS_Thread        ***********************************/
/***************************************************************/
void *AssembleRHS_Thread(void *data)
{ 
  ThreadData *TD=(ThreadData *)data;

  RWGGeometry *G    = TD->G;
  int RealFreq      = TD->RealFreq;
  EHFuncType EHFunc = TD->EHFunc;
  void *EHFuncUD    = TD->EHFuncUD;

  int nt, Type, Offset;
  int no, ne;
  RWGObject *O;
  cdouble EProd, HProd;
  double *DB;
  cdouble *ZB;

  if (RealFreq)
   ZB=TD->B->ZV;
  else
   DB=TD->B->DV; 

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

        /* there are four choices here based on whether we are at real or */
        /* imaginary frequency and whether the object is a perfect        */
        /* conductor (EFIE) or not (PMCHW)                                */
        if ( Type==MP_PEC && RealFreq==1 )
         { 
           O->GetInnerProducts(ne, EHFunc, EHFuncUD, RealFreq, &EProd, 0);
           ZB[ Offset + ne ] = -1.0*EProd / ZVAC;
         } 
        else if ( Type==MP_PEC && RealFreq==0 )
         { 
           O->GetInnerProducts(ne, EHFunc, EHFuncUD, RealFreq, &EProd, 0);
           DB[ Offset + ne ] = -1.0*(real(EProd)) / ZVAC;
         }
        else if ( Type!=MP_PEC && RealFreq==1 )
         { 
           O->GetInnerProducts(ne , EHFunc, EHFuncUD, RealFreq, &EProd, &HProd);
           ZB[ Offset + 2*ne    ]  = -1.0*EProd / ZVAC;
           ZB[ Offset + 2*ne + 1 ] = -1.0*HProd;
         }
        else // ( Type!=MP_PEC && RealFreq==0 )
         { 
           O->GetInnerProducts(ne, EHFunc, EHFuncUD, RealFreq, &EProd, &HProd);
           DB[ Offset + 2*ne    ]  = -1.0*(real(EProd)) /  ZVAC;
           DB[ Offset + 2*ne + 1 ] = -1.0*(real(HProd));
         };

      }; // for ne=...

   }; // for no=...

  return 0;
 
}

/***************************************************************/
/* Assemble the RHS vector.  ***********************************/
/***************************************************************/
void RWGGeometry::AssembleRHSVector(EHFuncType EHFunc, void *EHFuncUD,
                                    int RealFreq, int nThread, HVector *B)
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
     TD->RealFreq=RealFreq;
     TD->B=B;

     TD->nThread=nThread;

     pthread_create( &(Threads[nt]), 0, AssembleRHS_Thread, (void *)TD);
   };

  for(nt=0; nt<nThread; nt++)
   pthread_join(Threads[nt],0);

}

/***************************************************************/
/* Allocate an RHS vector of the appropriate size. *************/
/***************************************************************/
HVector *RWGGeometry::AllocateRHSVector(int RealFreq)
{ 
  HVector *V;

  if (RealFreq)
   V=new HVector(TotalBFs,LHM_COMPLEX);
  else
   V=new HVector(TotalBFs,LHM_REAL);

  return V;

} 
