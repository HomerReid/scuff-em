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
/* Prepare a chain of IncField structures for computations in  */
/* a given RWGGeometry at a given frequency and (optionally)   */
/* Bloch vector:                                               */
/*                                                             */
/*  (1) For each IncField in the chain, set the frequency to   */
/*      Omega, set Eps and Mu to the material properties       */
/*      (at frequency Omega) of the region within which the    */
/*      field sources are contained, and (for PBC geometries)  */
/*      set the Bloch vector and the lattice dimension.        */
/*                                                             */
/*  (2) Make sure the RegionIndex field in the IncField        */
/*      structure matches the index of the object specified by */
/*      the RegionLabel field (or, if the IncField implements  */
/*      the GetSourcePoint() routine, the index of the region  */
/*      containing the source point).                          */
/* Returns the total number of IncFields in the chain.         */
/***************************************************************/
int RWGGeometry::UpdateIncFields(IncField *IFList, cdouble Omega, double *kBloch)
{
  if (IFList==0) return 0;

  if ( (LBasis==0 && kBloch!=0) || (LBasis!=0 && kBloch==0) )
   ErrExit("%s:%i: internal error",__FILE__,__LINE__);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (LBasis)
   IFList->SetLattice(LBasis, true);

  /*--------------------------------------------------------------*/
  /*- make sure the cached epsilon and mu values for all regions  */
  /*- are up-to-date for the present frequency                    */
  /*--------------------------------------------------------------*/
  UpdateCachedEpsMuValues(Omega);

  /*--------------------------------------------------------------*/
  /*- run through the chain of IncField structures and set the   -*/
  /*- epsilon and mu values for each structure depending on the  -*/
  /*- region in which its sources are contained                  -*/
  /*--------------------------------------------------------------*/
  int NIF;
  IncField *IF;
  double X[3];
  for (NIF=0, IF=IFList; IF; NIF++, IF=IF->Next) 
   {
     /*--------------------------------------------------------------*/
     /*- first get the index of the object containing the field     -*/
     /*- sources for IF                                             -*/
     /*--------------------------------------------------------------*/
     if ( IF->GetSourcePoint(X) )
     {
      IF->RegionIndex = GetRegionIndex(X);
      IF->RegionLabel = strdupEC(RegionLabels[IF->RegionIndex]);
     }
     else if ( IF->RegionLabel )
      IF->RegionIndex = GetRegionByLabel(IF->RegionLabel);
     else
      IF->RegionIndex = 0; // exterior medium

     if ( IF->RegionIndex<0  || IF->RegionIndex>NumRegions )
      ErrExit("invalid region index %i",IF->RegionIndex);

     /*--------------------------------------------------------------*/
     /*- now set the material properties of IF as appropriate for   -*/
     /*- the region in question                                     -*/
     /*--------------------------------------------------------------*/
     IF->SetFrequencyAndEpsMu(Omega, EpsTF[ IF->RegionIndex ], MuTF[ IF->RegionIndex ] );
     if (kBloch)
      IF->SetkBloch(kBloch);

   }; // for(NIF=0, IF=IFList ... 

  return NIF;

}

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
      EH[n]-=dEH[n]; };
  
  /* compute dot products */
  cdouble *zF = (cdouble *)F;
  zF[0] = fRWG[0]*EH[0] + fRWG[1]*EH[1] + fRWG[2]*EH[2];
  if (IPID->NeedHProd)
   zF[1] = fRWG[0]*EH[3] + fRWG[1]*EH[4] + fRWG[2]*EH[5];
  
} 

/***************************************************************/
/* Calculate the inner product of given electric and magnetic  */
/* fields with the basis function associated with edge #ne on  */
/* the given RWGSurface.                                       */
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
void GetInnerProducts(RWGSurface *S, int ne, 
                      IncField **PositiveIFs, int NPositiveIFs,
                      IncField **NegativeIFs, int NNegativeIFs,
                      cdouble *pEProd, cdouble *pHProd)
{ 
  /* get edge vertices */
  RWGEdge *E   = S->Edges[ne];
  double *QP   = S->Vertices + 3*(E->iQP);
  double *V1   = S->Vertices + 3*(E->iV1);
  double *V2   = S->Vertices + 3*(E->iV2);
  double PArea = S->Panels[E->iPPanel]->Area;
  double *QM;
  double MArea=0.0;
  if ( E->iQM == -1 )
   QM = 0;
  else
   { QM = S->Vertices + 3*(E->iQM);
     MArea = S->Panels[E->iMPanel]->Area;
   };

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

  /* integrate over negative panel if present */
  if (QM)
   { IPID->Q=QM;
     IPID->PreFac=E->Length / (2.0*MArea);
     TriIntFixed(InnerProductIntegrand, nFun, (void *)IPID, V1, V2, QM, 20, IM);
   }
  else
   { 
     memset(IM, 0, 4*sizeof(double));
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
   int nt, NumTasks;

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
  /* loop over all surfaces to get contributions to RHS vector   */
  /***************************************************************/
  RWGSurface *S;
  int ne, Offset, IsPEC;
  int nt=0;
  cdouble EProd, HProd;
  IncField *IF;
  for(int ns=0; ns<G->NumSurfaces; ns++)
   { 
     S=G->Surfaces[ns];
     Offset=G->BFIndexOffset[ns];
     IsPEC=S->IsPEC;

     /*--------------------------------------------------------------*/
     /*- Go through the chain of IncField structures to identify     */
     /*- the subset of IncFields that contribute with a plus sign to */
     /*- the RHS vector entries for this surface, as well as those   */
     /*- that contribute with a minus sign.                          */
     /*- IncFields whose sources lie in the 'positive' region        */
     /*- associated with this surface contribute with a minus sign,  */
     /*- IncFields whose sources lie in the 'negative' region        */
     /*- associated with this surface contribute with a plus sign,   */
     /*- and all other IncFields do not contribute.                  */
     /*- The apparent sign inconsistency here arises because the     */
     /*- incident-field terms in the tangential-fields-must-be-equal */
     /*- equation are swung over to the RHS, acquiring a minus sign. */
     /*--------------------------------------------------------------*/
     for(NPositiveIFs=NNegativeIFs=0, IF=IFList; IF; IF=IF->Next)
      { 
        if (S->RegionIndices[0]==IF->RegionIndex)
         NegativeIFs[NNegativeIFs++] = IF;
        else if (S->RegionIndices[1]==IF->RegionIndex)
         PositiveIFs[NPositiveIFs++] = IF;
      };
     if ( NPositiveIFs==0 && NNegativeIFs==0 )
      continue;

     /*--------------------------------------------------------------*/
     /*- Loop over all basis functions (edges) on this object to get-*/
     /*- each BF's contribution to the RHS.                         -*/
     /*--------------------------------------------------------------*/
     for(ne=0; ne<S->NumEdges; ne++)
      { 
        nt++;
        if (nt==TD->NumTasks) nt=0;
        if (nt!=TD->nt) continue;

        GetInnerProducts(S,ne, 
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
HVector *RWGGeometry::AssembleRHSVector(cdouble Omega, double *kBloch,
                                        IncField *IF, HVector *RHS)
{ 
  if (RHS==NULL)
   RHS=AllocateRHSVector();

  RHS->Zero();
   
  int nt, NumTasks, NumThreads = GetNumThreads();
  int NIF=UpdateIncFields(IF, Omega, kBloch);

  ThreadData ReferenceTD;
  ReferenceTD.G=this;
  ReferenceTD.IF=IF;
  ReferenceTD.NIF=NIF;
  ReferenceTD.RHS=RHS;

#ifdef USE_PTHREAD
  ThreadData *TDs = new ThreadData[NumThreads], *TD;
  pthread_t *Threads = new pthread_t[NumThreads];
  for(nt=0; nt<NumThreads; nt++)
   { 
     TD=&(TDs[nt]);
     memcpy(TD, &ReferenceTD, sizeof(ThreadData));
     TD->nt=nt;
     TD->NumTasks=NumThreads;

     if (nt+1 == NumThreads)
       AssembleRHS_Thread((void *)TD);
     else
       pthread_create( &(Threads[nt]), 0, AssembleRHS_Thread, (void *)TD);
   };
  for(nt=0; nt<NumThreads-1; nt++)
   pthread_join(Threads[nt],0);

  delete[] Threads;
  delete[] TDs;

#else

#ifndef USE_OPENMP
  NumThreads=NumTasks=1;
#else
  NumTasks=NumThreads*100;
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(nt=0; nt<NumTasks; nt++)
   { 
     ThreadData TD1;
     memcpy(&TD1, &ReferenceTD, sizeof(ThreadData));
     TD1.nt=nt;
     TD1.NumTasks=NumTasks;
     AssembleRHS_Thread((void *)&TD1);
   };
#endif

  if (UseHRWGFunctions && NumMMJs>0 )
   ApplyMMJTransformation(0, RHS);

  return RHS;
}

/***************************************************************/
/* non-PBC entry point for AssembleRHSVector                   */
/***************************************************************/
HVector *RWGGeometry::AssembleRHSVector(cdouble Omega, IncField *IF, HVector *RHS)
{ return AssembleRHSVector(Omega, 0, IF, RHS); }

/***************************************************************/
/* Allocate an RHS vector of the appropriate size. *************/
/***************************************************************/
HVector *RWGGeometry::AllocateRHSVector(bool PureImagFreq)
{ 
  HVector *V;

  if (PureImagFreq && LBasis==0)
   V=new HVector(TotalBFs,LHM_REAL);
  else
   V=new HVector(TotalBFs,LHM_COMPLEX);

  return V;

} 

} // namespace scuff
