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
/* with the magnetic field as well.                            */
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
/* Calculate the inner product of given electric and magnetic  */
/* fields with the basis function associated with edge #ne on  */
/* the given RWGObject.                                        */
/*                                                             */
/* PositiveIFs[0..NPositiveIFs] are IncFields whose fields     */
/* contribute with a positive sign.                            */
/*                                                             */
/* NegativeIFs[0..NNegativeIFs] are IncFields whose fields     */
/* contribute with a negative sign.                            */
/*                                                             */
/* If pHProd==0, we skip the computation of the magnetic-field */
/* inner product.                                              */
/***************************************************************/
void GetInnerProducts(RWGObject *O, int ne, 
                      IncField **PositiveIFs, int NPositiveIFs,
                      IncField **NegativeIFs, int NNegativeIFs,
                      cdouble *pEProd, cdouble *pHProd)
{ 
  /* get edge vertices */
  RWGEdge *E   = O->Edges[ne];
  double *QP   = O->Vertices + 3*(E->iQP);
  double *V1   = O->Vertices + 3*(E->iV1);
  double *V2   = O->Vertices + 3*(E->iV2);
  double *QM   = O->Vertices + 3*(E->iQM);
  double PArea = O->Panels[E->iPPanel]->Area;
  double MArea = O->Panels[E->iMPanel]->Area;

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
  IPID->Q=QM;
  IPID->PreFac=E->Length / (2.0*MArea);
  TriIntFixed(InnerProductIntegrand, nFun, (void *)IPID, V1, V2, QM, 20, IM);

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

   RWGGeometry *G;
   IncField *IF;
   int NIF;
   HVector *RHS;

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
  RWGGeometry *G   = TD->G;
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
  /* loop over all objects to get contributions to RHS vector    */
  /***************************************************************/
  RWGObject *O;
  int no, ne, Offset, IsPEC;
  int nt=0;
  int ContainingObjectIndex;
  cdouble EProd, HProd;
  IncField *IF;
  for(no=0; no<G->NumObjects; no++)
   { 
     O=G->Objects[no];
     Offset=G->BFIndexOffset[no];
     IsPEC=O->MP->IsPEC();
     if ( O->ContainingObject==NULL )
      ContainingObjectIndex=-1;
     else
      ContainingObjectIndex=O->ContainingObject->Index;

     /*--------------------------------------------------------------*/
     /*- Go through the chain of IncField structures to identify     */
     /*- the subset of IncFields whose sources lie inside this       */
     /*- object, as well as the subset whose sources lie in the      */
     /*- region to the immediate exterior of this object. (The       */
     /*- former contribute to the RHS vector with a plus sign, while */
     /*- the latter contribute with a minus sign). If both subsets   */
     /*- are empty, this object does not contribute to the RHS.      */
     /*--------------------------------------------------------------*/
     for(NPositiveIFs=NNegativeIFs=0, IF=IFList; IF; IF=IF->Next)
      { 
        if (O->Index==IF->ObjectIndex)
         PositiveIFs[NPositiveIFs++] = IF;
        else if (ContainingObjectIndex==IF->ObjectIndex)
         NegativeIFs[NNegativeIFs++] = IF;
      };
     if ( NPositiveIFs==0 && NNegativeIFs==0 )
      continue;

     /*--------------------------------------------------------------*/
     /*- Loop over all basis functions (edges) on this object to get-*/
     /*- each BF's contribution to the RHS.                         -*/
     /*--------------------------------------------------------------*/
     for(ne=0; ne<O->NumEdges; ne++)
      { 
        nt++;
        if (nt==TD->nThread) nt=0;
        if (nt!=TD->nt) continue;

        GetInnerProducts(O,ne, 
                         PositiveIFs, NPositiveIFs,
                         NegativeIFs, NNegativeIFs,
                         &EProd, IsPEC ? 0 : &HProd );

        if ( IsPEC )
         { 
           RHS->SetEntry(Offset + ne, EProd / ZVAC);
         }
        else 
         { RHS->SetEntry(Offset + 2*ne+0, EProd / ZVAC);
           RHS->SetEntry(Offset + 2*ne+1, HProd);
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
                                        HVector *RHS, int nThread)
{ 
  if (RHS==NULL)
   RHS=AllocateRHSVector();
   
  if (nThread <= 0) 
   nThread = GetNumThreads();

  int nt, NIF=UpdateIncFields(IF, Omega);

  ThreadData ReferenceTD;
  ReferenceTD.G=this;
  ReferenceTD.IF=IF;
  ReferenceTD.NIF=NIF;
  ReferenceTD.RHS=RHS;

#ifdef USE_PTHREAD
  ThreadData *TDS = new ThreadData[nThread], *TD;
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
  delete[] TDS;

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
  IncField *IF;

  ExteriorMP->GetEpsMu(Omega, &ExteriorEps, &ExteriorMu);
  
  for (NIF=0, IF=IFList; IF; NIF++, IF=IF->Next) 
   {
     /*--------------------------------------------------------------*/
     /*- first get the index of the object containing the field     -*/
     /*- sources for IF                                             -*/
     /*--------------------------------------------------------------*/
     if ( IF->GetSourcePoint(X) )
      IF->ObjectIndex = GetObjectIndex(X);
     else if ( IF->ObjectLabel )
      GetObjectByLabel(IF->ObjectLabel, &(IF->ObjectIndex) );
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
        IF->SetFrequencyAndEpsMu(Omega, ObjectEps, ObjectMu, 0 );
      };

   }; // for(NIF=0, IF=IFList ... 

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

} // namespace scuff
