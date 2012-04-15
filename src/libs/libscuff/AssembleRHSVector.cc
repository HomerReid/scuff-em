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
   HVector *B;

 } ThreadData;

/***************************************************************/
/* AssembleRHS_Thread        ***********************************/
/***************************************************************/
void *AssembleRHS_Thread(void *data)
{ 
  ThreadData *TD    = (ThreadData *)data;

#ifdef USE_PTHREAD
  SetCPUAffinity(TD->nt);
#endif

  /***************************************************************/
  /* extract fields from thread data structure *******************/
  /***************************************************************/
  RWGGeometry *G       = TD->G;
  IncField    *IFList  = TD->IF;
  int         NIF      = TD->NIF;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  IncField *PositiveIFs = new *IncField[NIF];
  IncField *NegativeIFs = new *IncField[NIF];
  int NPositiveIFs;
  int NNegativeIFs;

  RWGObject *O;
  int no, ne, Offset, IsPEC;
  int nt=0;
  cdouble EProd, HProd;
  for(no=0, O=G->Objects[0]; no<G->NumObjects; O=G->Objects[++no])
   { 
     /*--------------------------------------------------------------*/
     /*- Go through the chain of IncField structures to identify     */
     /*- the subset of IncFields whose sources lie inside this       */
     /*- object, as well as the subset whose sources lie in the      */
     /*- region to the immediate exterior of this object. (The       */
     /*- former contribute to the RHS vector with a plus sign, while */
     /*- the latter contribute with a minus sign). If both subsets   */
     /*- are empty, this object does not contribute to the RHS.      */
     /*--------------------------------------------------------------*/
     NPositiveIFs = NNegativeIFs = 0;
     for(IF=IFList; IF!=0; IF=IF->Next)
      { 
        if (O->Index==IF->ObjectIndex)
         PositiveIFs[NPositiveIFs++] = IF;
        else if (ContainingObjectIndex==IF->ObjectIndex)
         NegativeIFs[NNegativeIFs++] = IF;
      };
     if ( (NPositiveIFs + NNegativeIFs) == 0 ) 
      continue;

     /*--------------------------------------------------------------*/
     /*- Loop over all basis functions (edges) on this object to get-*/
     /*- each BF's contribution to the RHS.                         -*/
     /*--------------------------------------------------------------*/
     Offset=G->BFIndexOffset[no];
     IsPEC=O->MP->IsPEC();
     for(ne=0; ne<O->NumEdges; ne++)
      { 
        nt++;
        if (nt==TD->nThread) nt=0;
        if (nt!=TD->nt) continue;

        O->GetInnerProducts( ne, 
                             PositiveIFs, NPositiveIFs,
                             NegativeIFs, NNegativeIFs,
                             &EProd, IsPEC ? 0 : &HProd );

        if ( IsPEC )
         { 
           RHS->SetEntry(Offset + ne), EProd / ZVac;
         }
        else 
         { RHS->SetEntry(Offset + 2*ne+0), EProd / ZVac;
           RHS->SetEntry(Offset + 2*ne+1), HProd;
         };

      }; // for ne=...

   }; // for no=...

  delete[] PositiveIFs;
  delete[] NegativeIFs;

  return 0;
 
}

/***************************************************************/
/* Assemble the RHS vector.  ***********************************/
/***************************************************************/
HVector *RWGGeometry::AssembleRHSVector(cdouble Omega, IncField *IF,
                                        HVector *B, int nThread)
{ 
  if (B==NULL)
   B=AllocateRHSVector();
   
  if (nThread <= 0) 
   nThread = GetNumThreads();

  // count the number of IncField structures in the chain
  int NIF=0;
  for (IncField *Node=IF; Node; Node=Node->Next)
   NIF++;

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
     TD->NIF=NIF;
     TD->RHS=RHS;

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

  return RHS;
}

/***************************************************************/
/* Prepare a chain of IncField structures for computations in  */
/* a given RWGGeometry at a given geometry:                    */
/*  (1) For each IncField in the chain, set the frequency to   */
/*      Omega, and set Eps and Mu to the material properties   */
/*      (at frequency Omega) of the RWGObject within which the */
/*      field sources are contained.                           */
/*  (2) Make sure the ObjectIndex field in the IncField        */
/*      structure matches the index of the object specified by */
/*      the ObjectLabel field.                                 */
/* Returns the total number of IncFields in the chain.         */
/***************************************************************/
int RWGGeometry::UpdateIncFields(IncField *IFList, cdouble Omega)
{
  cdouble ExteriorEps, ExteriorMu;
  cdouble ObjectEps, ObjectMu;
  RWGObject *O; 
  double X[3];
  int NIF;

  ExteriorMP->GetEpsMu(Omega, &ExteriorEps, &ExteriorMu);

  for (NIF=0, IncField *IF = IFList; IF; IF = IF->Next) 
   {
     NIF++;

     /*--------------------------------------------------------------*/
     /*- first get the index of the object containing the field     -*/
     /*- sources for IF                                             -*/
     /*--------------------------------------------------------------*/
     if ( IF->ObjectLabel ) 
      GetObjectByLabel(IF->ObjectLabel, &(IF->ObjectIndex) );
     else if ( IF->GetSourcePoint(X) )
      IF->ObjectIndex = GetObjectIndex(X);
     else
      IF->ObjectIndex = -1;

     /*--------------------------------------------------------------*/
     /*- now set the material properties of IF                      -*/
     /*--------------------------------------------------------------*/
     if ( IF->ObjectIndex == -1 )
      {  
        IF->SetFrequencyAndEpsMu(Omega, ExteriorEps, ExteriorMu, 0 );
      }
     else
      { 
        if ( IF->ObjectIndex<0  || IF->ObjectIndex>NumObjects )
         ErrExit("invalid object specification %i",IF->ObjectIndex);

        O=Objects[IF->ObjectIndex];
        O->MP->GetEpsMu(Omega, &ObjectEps, &ObjectMu);
      };
   };

  return NIF;

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
