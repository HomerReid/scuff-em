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
 * EMTPFT.cc     -- libscuff class methods for computing power, force,
 *               -- and torque in classical deterministic scattering
 *               -- problems using the "energy/momentum-transfer" method
 *
 * homer reid    -- 10/2015
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <fenv.h>

#include <libhrutil.h>

#include "libscuff.h"
#include "libscuffInternals.h"
#include "PanelCubature.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#ifdef USE_OPENMP
#  include <omp.h>
#endif

#define II cdouble(0,1)

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
namespace scuff {

void CalcGC(double R[3], cdouble Omega,
            cdouble EpsR, cdouble MuR,
            cdouble GMuNu[3][3], cdouble CMuNu[3][3],
            cdouble GMuNuRho[3][3][3], cdouble CMuNuRho[3][3][3]);

void GetKNBilinears(HVector *KNVector, HMatrix *DRMatrix,
                    bool IsPECA, int KNIndexA,
                    bool IsPECB, int KNIndexB,
                    cdouble Bilinears[4]);

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct PFTIntegrandData
 {
   cdouble k;
   IncField *IF;
   double XTorque[3];

 } PFTIntegrandData;

/***************************************************************/
/***************************************************************/
/***************************************************************/
void PFTIntegrand_BFInc(double *x, PCData *PCD,
                        void *UserData, double *I)
{
  double b[3];
  b[0] = real(PCD->K[0]);
  b[1] = real(PCD->K[1]);
  b[2] = real(PCD->K[2]);

  PFTIntegrandData *PFTIData=(PFTIntegrandData *)UserData;
  IncField *IF    = PFTIData->IF;
  double *XTorque = PFTIData->XTorque;
  double XT[3];
  VecSub(x, XTorque, XT);

  // get fields and derivatives at eval point
  cdouble EH[6], *E=EH+0, *H=EH+3, dEH[3][6];
  IF->GetFields(x, EH);
  IF->GetFieldGradients(x, dEH);

  cdouble *zI = (cdouble *)I;
  memset(zI, 0, 2*NUMPFT*sizeof(cdouble));
  for(int Mu=0; Mu<3; Mu++)
   { 
     zI[0 + PFT_PABS]    += b[Mu]*EH[Mu];
     zI[0 + PFT_XFORCE]  += b[Mu]*dEH[0][Mu];
     zI[0 + PFT_YFORCE]  += b[Mu]*dEH[1][Mu];
     zI[0 + PFT_ZFORCE]  += b[Mu]*dEH[2][Mu];
     zI[0 + PFT_XTORQUE] += b[Mu]*(XT[1]*dEH[2][Mu]-XT[2]*dEH[1][Mu]);
     zI[0 + PFT_YTORQUE] += b[Mu]*(XT[2]*dEH[0][Mu]-XT[0]*dEH[2][Mu]);
     zI[0 + PFT_ZTORQUE] += b[Mu]*(XT[0]*dEH[1][Mu]-XT[1]*dEH[0][Mu]);

     int MP3=Mu+3;
     zI[NUMPFT + PFT_PABS]    += b[Mu]*EH[MP3];
     zI[NUMPFT + PFT_XFORCE]  += b[Mu]*dEH[0][MP3];
     zI[NUMPFT + PFT_YFORCE]  += b[Mu]*dEH[1][MP3];
     zI[NUMPFT + PFT_ZFORCE]  += b[Mu]*dEH[2][MP3];
     zI[NUMPFT + PFT_XTORQUE] += b[Mu]*(XT[1]*dEH[2][MP3]-XT[2]*dEH[1][MP3]);
     zI[NUMPFT + PFT_YTORQUE] += b[Mu]*(XT[2]*dEH[0][MP3]-XT[0]*dEH[2][MP3]);
     zI[NUMPFT + PFT_ZTORQUE] += b[Mu]*(XT[0]*dEH[1][MP3]-XT[1]*dEH[0][MP3]);
   };

  zI[0      + PFT_XTORQUE] += b[1]*E[2] - b[2]*E[1];
  zI[0      + PFT_YTORQUE] += b[2]*E[0] - b[0]*E[2];
  zI[0      + PFT_ZTORQUE] += b[0]*E[1] - b[1]*E[0];

  zI[NUMPFT + PFT_XTORQUE] += b[1]*H[2] - b[2]*H[1];
  zI[NUMPFT + PFT_YTORQUE] += b[2]*H[0] - b[0]*H[2];
  zI[NUMPFT + PFT_ZTORQUE] += b[0]*H[1] - b[1]*H[0];

}

/***************************************************************/
/* compute PFT integrals between an RWG basis function and an  */
/* external field                                              */
/***************************************************************/
void GetPFTIntegrals_BFInc(RWGGeometry *G, int ns, int ne,
                           IncField *IF, cdouble Omega,
                           cdouble IBFInc[2*NUMPFT])
{
  PFTIntegrandData MyPFTIData, *PFTIData=&MyPFTIData;
  PFTIData->k  = Omega;
  PFTIData->IF = IF;
  PFTIData->XTorque[0]=PFTIData->XTorque[1]=PFTIData->XTorque[2]=0.0;

  RWGSurface *S=G->Surfaces[ns];
  if (S->OTGT) S->OTGT->Apply(PFTIData->XTorque);
  if (S->GT) S->GT->Apply(PFTIData->XTorque);

  int fDim=4*NUMPFT;
  int MaxEvals=21;
  GetBFCubature(G, ns, ne,
                PFTIntegrand_BFInc, (void *)PFTIData, fDim,
                MaxEvals, 0, 0, 0, 0, (double *)IBFInc);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddIFContributionsToEMTPFT(RWGGeometry *G, HVector *KN,
                                IncField *IF, cdouble Omega,
                                HMatrix *PFTMatrix, bool Interior)
{
  if ( PFTMatrix->NR!=G->NumSurfaces || PFTMatrix->NC != NUMPFT )
   ErrExit("%s:%i: internal error", __FILE__, __LINE__);

  int NT=1;
#ifdef USE_OPENMP
  NT = GetNumThreads();
#endif
  int NS=G->NumSurfaces;
  int NQ=NUMPFT;
  double *PartialPFT=new double[NT*NS*NQ];
  memset(PartialPFT, 0, NT*NS*NQ*sizeof(double));

#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1), num_threads(NT)
#endif
  for(int neTot=0; neTot<G->TotalEdges; neTot++)
   { 
     int ns, ne, KNIndex;
     RWGSurface *S = G->ResolveEdge(neTot, &ns, &ne, &KNIndex);
  
     // TODO: this will fail if the IF chain contains multiple
     //       field sources in different regions
     if (IF->RegionIndex != S->RegionIndices[ Interior ? 1 : 0 ] )
      continue;

     cdouble KStar = conj(KN->GetEntry(KNIndex));
     cdouble NStar = 0.0;
     if (!(S->IsPEC)) NStar = -ZVAC*conj(KN->GetEntry(KNIndex+1));

     cdouble IBFInc[2*NUMPFT];
     GetPFTIntegrals_BFInc(G, ns, ne, IF, Omega, IBFInc);
     cdouble bDotE = IBFInc[0],      *bDotdE=IBFInc+PFT_XFORCE;
     cdouble bDotH = IBFInc[NUMPFT], *bDotdH=IBFInc+NUMPFT+PFT_XFORCE;

     int nt=0;
#ifdef USE_OPENMP
     nt = omp_get_thread_num();
#endif
     double *dPFT = PartialPFT + (nt*NS + ns)*NUMPFT;
     dPFT[PFT_PABS] += real(KStar*bDotE + NStar*bDotH);
     for(int Mu=0; Mu<6; Mu++)
      dPFT[PFT_XFORCE+Mu] += imag( KStar*bDotdE[Mu] + NStar*bDotdH[Mu] );
   };

   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   double PFactor = 0.5;
   double FTFactor = 0.5*TENTHIRDS/real(Omega);
   for(int nt=0; nt<NT; nt++)
    for(int ns=0; ns<NS; ns++)
     { 
       double *dPFT = PartialPFT + (nt*NS + ns)*NUMPFT;
       PFTMatrix->AddEntry(ns, PFT_PABS, PFactor*dPFT[PFT_PABS]);
       for(int Mu=0; Mu<6; Mu++)
        PFTMatrix->AddEntry(ns, PFT_XFORCE + Mu, FTFactor*dPFT[PFT_XFORCE+Mu]);
     };

  delete[] PartialPFT;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *GetEMTPFT(RWGGeometry *G, cdouble Omega, IncField *IF,
                   HVector *KNVector, HMatrix *DRMatrix,
                   HMatrix *PFTMatrix, bool Interior)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NS = G->NumSurfaces;
  if (    (PFTMatrix==0)
       || (PFTMatrix->NR != NS)
       || (PFTMatrix->NC != NUMPFT)
     )
   ErrExit("invalid PFTMatrix in GetEMTPFT");

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
      GetGCMEArgStruct MyArgs, *Args=&MyArgs;
      InitGetGCMEArgs(Args);
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int NT=1;
#ifdef USE_OPENMP 
  NT = GetNumThreads();
#endif
  int NQ     = NUMPFT;
  int NSNQ   = NS*NQ;
  int NTNSNQ = NT*NS*NQ; 
  static int DeltaPFTSize=0;
  static double *DeltaPFT=0;
  if ( DeltaPFTSize < NTNSNQ )
   { Log("(re)allocating DeltaPFT (%i,%i)",DeltaPFTSize,NTNSNQ);
     DeltaPFTSize=NTNSNQ;
     if (DeltaPFT) free(DeltaPFT);
     DeltaPFT = (double *)mallocEC(DeltaPFTSize*sizeof(double));
   };
  memset(DeltaPFT, 0, DeltaPFTSize*sizeof(double));

  /*--------------------------------------------------------------*/
  /*- multithreaded loop over all basis functions on all surfaces-*/
  /*--------------------------------------------------------------*/
  int TotalEdges  = G->TotalEdges;
#ifdef USE_OPENMP
  Log("EMT OpenMP multithreading (%i threads)",NT);
#pragma omp parallel for schedule(dynamic,1), num_threads(NT)
#endif
  for(int neaTot=0; neaTot<TotalEdges; neaTot++)
   for(int nebTot=neaTot; nebTot<TotalEdges; nebTot++)
    { 
      if (nebTot==neaTot) LogPercent(neaTot, TotalEdges, 10);

      int nsa, nea, KNIndexA;
      RWGSurface *SA = G->ResolveEdge(neaTot, &nsa, &nea, &KNIndexA);

      int nsb, neb, KNIndexB;
      RWGSurface *SB = G->ResolveEdge(nebTot, &nsb, &neb, &KNIndexB);
   
      double abSign=0.0, baSign=0.0;
      int RegionIndex;
      if (Interior)
       { RegionIndex=SA->RegionIndices[1];
         if ( SB->RegionIndices[1] == RegionIndex )
          abSign=baSign=-1.0;
         else if (SB->RegionIndices[0]==RegionIndex)
          abSign=1.0;
       }
      else
       { RegionIndex=SA->RegionIndices[0];
         if ( SB->RegionIndices[0] == RegionIndex )
          abSign=baSign=+1.0;
         else if (SB->RegionIndices[1]==RegionIndex)
          abSign=-1.0;
       };
      if ( RegionIndex==-1 || (abSign==0.0 && baSign==0.0) )
       continue;

      cdouble EpsR, MuR;
      G->RegionMPs[RegionIndex]->GetEpsMu(Omega, &EpsR, &MuR);
      cdouble k = Omega*sqrt(EpsR*MuR);

      cdouble KNB[4];
      GetKNBilinears(KNVector, DRMatrix,
                     SA->IsPEC, KNIndexA, SB->IsPEC, KNIndexB,
                     KNB);

      cdouble wu0KK  = KNB[0] * Omega * ZVAC;
      cdouble KNmNK  = KNB[1] - KNB[2];
      cdouble we0NN  = KNB[3] * Omega / ZVAC;
                     
      GetGCMEArgStruct MyArgs, *Args=&MyArgs;
      InitGetGCMEArgs(Args);
      Args->nsa = nsa;
      Args->nsb = nsb;
      Args->NumRegions = 1;
      Args->k[0] = k;
      Args->NeedGC=Args->NeedForce=Args->NeedTorque=true;
      Args->FIBBICache = (nsa==nsb) ? G->FIBBICaches[nsa] : 0;
      
      cdouble GabArray[2][NUMGCMES];
      cdouble ikCabArray[2][NUMGCMES];
      GetGCMatrixElements(G, Args, nea, neb, GabArray, ikCabArray);
      cdouble Gab     =   GabArray[0][GCME_GC];
      cdouble ikCab   = ikCabArray[0][GCME_GC];
      cdouble *dGab   =   GabArray[0] + GCME_FX;
      cdouble *dikCab = ikCabArray[0] + GCME_FX;
  
      double dPFT[NUMPFT];
      dPFT[PFT_PSCAT] = 0.5*(  real(wu0KK)*imag( MuR*Gab)
                              +real(we0NN)*imag(EpsR*Gab)
                              +imag(KNmNK)*imag(   ikCab)
                            );

      for(int Mu=0; Mu<6; Mu++)
       dPFT[PFT_XFORCE + Mu]
        = 0.5*TENTHIRDS*(  imag(wu0KK)*imag( MuR*dGab[Mu])
                          +imag(we0NN)*imag(EpsR*dGab[Mu])
                          -real(KNmNK)*imag(   dikCab[Mu])
                        ) / real(Omega);

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
dPFT[PFT_XFORCE] = 0.5*TENTHIRDS*( imag(wu0KK)*imag( MuR*dGab[2])
                                  +imag(we0NN)*imag(EpsR*dGab[2]) ) / real(Omega);
dPFT[PFT_YFORCE] = 0.5*TENTHIRDS*( -real(KNmNK)*imag(  dikCab[2]) ) / real(Omega);  

dPFT[PFT_XTORQUE] = 0.5*( real(wu0KK)*imag( MuR*Gab)
                         +real(we0NN)*imag(EpsR*Gab)
                        );

dPFT[PFT_YTORQUE] = 0.5*( imag(KNmNK)*imag(ikCab) );
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/


      int nt=0;
#ifdef USE_OPENMP
      nt=omp_get_thread_num();
#endif
      int OffsetA = nt*NSNQ + nsa*NQ;
      int OffsetB = nt*NSNQ + nsb*NQ;

       if (neaTot==nebTot)
        DeltaPFT[ OffsetA + PFT_PABS ] += abSign*dPFT[PFT_PABS];
       else if (nsa==nsb) // nebTot > neaTot but both on same surface
        { 
          DeltaPFT[ OffsetA + PFT_PABS ] += 2.0*abSign*dPFT[PFT_PABS];
          for(int Mu=0; Mu<6; Mu++)
           DeltaPFT[ OffsetA + PFT_XFORCE + Mu ] 
            += 2.0*abSign*dPFT[PFT_XFORCE + Mu];
        }
       else // nebTot > neaTot and on different objects
        { 
          DeltaPFT[ OffsetA + PFT_PABS ] += abSign*dPFT[PFT_PABS];
          DeltaPFT[ OffsetB + PFT_PABS ] += baSign*dPFT[PFT_PABS];
          for(int Mu=0; Mu<6; Mu++)
           { DeltaPFT[ OffsetA + PFT_XFORCE + Mu ] 
              += abSign*dPFT[PFT_XFORCE + Mu];
             DeltaPFT[ OffsetB + PFT_XFORCE + Mu ] 
              += baSign*dPFT[PFT_XFORCE + Mu];
           };
        };
 
    }; // end of multithreaded loop
  
  /*--------------------------------------------------------------*/
  /*- accumulate contributions of all threads                     */
  /*--------------------------------------------------------------*/
  PFTMatrix->Zero();
  for(int ns=0; ns<NS; ns++)
   for(int nq=0; nq<NQ; nq++)
    for(int nt=0; nt<NT; nt++)
     PFTMatrix->AddEntry(ns, nq, DeltaPFT[ nt*NSNQ + ns*NQ + nq ]);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (!Interior)
   for(int ns=0; ns<NS; ns++)
    PFTMatrix->SetEntry(ns, PFT_PSCAT,
                       -1.0*PFTMatrix->GetEntry(ns, PFT_PABS));

  /***************************************************************/
  /* add incident-field contributions ****************************/
  /***************************************************************/
  if (IF)
   AddIFContributionsToEMTPFT(G, KNVector, IF, Omega, PFTMatrix, Interior);
  
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int ns=0; ns<G->NumSurfaces; ns++)
   if (G->Mate[ns]==-1)
    { 
      StoreFIBBICache(G->FIBBICaches[ns], G->Surfaces[ns]->MeshFileName);
      int Hits, Misses;
      int CacheSize=GetFIBBICacheSize(G->FIBBICaches[ns],&Hits, &Misses);
      Log("EMTPFT surface %i: (%i/%i) hits/misses",ns,Hits,Misses);
    };

  return PFTMatrix;
}
  
} // namespace scuff
