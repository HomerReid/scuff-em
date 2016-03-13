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

void GetKNBilinears(HVector *KNVector, HMatrix *DRMatrix,
                    bool IsPECA, int KNIndexA,
                    bool IsPECB, int KNIndexB,
                    cdouble Bilinears[4]);

void GetReducedFields_Nearby(RWGSurface *S, const int ne,
                             const double X0[3], const cdouble k,
                             cdouble e[3], cdouble h[3],
                             cdouble de[3][3], cdouble dh[3][3],
                             bool *IncludeTerm=0);

inline cdouble VecDot(double *V1, cdouble *V2)
 { return V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2]; }

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct PFTIData
 {
   cdouble Omega, k, EpsRel, MuRel;
   IncField *IF;
   double *TorqueCenter;
   int EMTPFTIMethod;
   RWGGeometry *G;
   int nsb, neb;

 } PFTIData;

/***************************************************************/
/***************************************************************/
/***************************************************************/
void PFTIntegrand_BFBF2(double *xA, PCData *PCD,
                        void *UserData, double *I)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double bA[3], DivbA, *nHat;
  bA[0] = real(PCD->K[0]);
  bA[1] = real(PCD->K[1]);
  bA[2] = real(PCD->K[2]);
  DivbA = real(PCD->DivK);
  nHat  = PCD->nHat;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  PFTIData *Data=(PFTIData *)UserData;
  cdouble Omega   = Data->Omega;
  cdouble k       = Data->k;
  cdouble EpsRel  = Data->EpsRel;
  int EMTPFTIMethod = Data->EMTPFTIMethod;
  cdouble  MuRel  = Data->MuRel;
  RWGGeometry *G  = Data->G;
  int nsb         = Data->nsb;
  int neb         = Data->neb;
  double *XTorque = Data->TorqueCenter;

  double XT[3];
  VecSub(xA, XTorque, XT);

  double XEval[3];
  VecScaleAdd(xA, 1.0e-2, nHat, XEval);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble e[3], h[3], de[3][3], dh[3][3];
  GetReducedFields_Nearby(G->Surfaces[nsb], neb, XEval, k, e, h, de, dh);

  cdouble bxe[3], bxh[3];
  for(int Mu=0; Mu<3; Mu++)
   { int MP1=(Mu+1)%3, MP2=(Mu+2)%3;
     bxe[Mu] = bA[MP1]*e[MP2] - bA[MP2]*e[MP1];
     bxh[Mu] = bA[MP1]*h[MP2] - bA[MP2]*h[MP1];
   };

  double *QKK = I+0, *QNN=I+7, *QKNmNK=I+14;
  QKK[0]    =  imag(MuRel*VecDot(bA, e));
  QNN[0]    =  imag(EpsRel*VecDot(bA, e));
  QKNmNK[0] =  imag(VecDot(bA, h));
  if (EMTPFTIMethod==1)
   {
     for(int Mu=0; Mu<3; Mu++)
      { QKK[1+Mu]    =  imag(MuRel*VecDot(bA, de[Mu]));
        QNN[1+Mu]    =  imag(EpsRel*VecDot(bA, de[Mu]));
        QKNmNK[1+Mu] = -imag(VecDot(bA, dh[Mu]));

        QKK[4+Mu]    =  imag(MuRel*bxe[Mu]);
        QNN[4+Mu]    =  imag(EpsRel*bxe[Mu]);
        QKNmNK[4+Mu] = -imag(bxh[Mu]);
      };
   }
  else
   {
     for(int Mu=0; Mu<3; Mu++)
      { QKK[1+Mu]    =  imag(MuRel*(-DivbA*e[Mu] + bxh[Mu]));
        QNN[1+Mu]    =  imag(EpsRel*(-DivbA*e[Mu] + bxh[Mu]));
        QKNmNK[1+Mu] = -imag(DivbA*h[Mu] - k*k*bxe[Mu]);
      };

     for(int Mu=0; Mu<3; Mu++)
      { int MP1=(Mu+1)%3, MP2=(Mu+2)%3;
        QKK[4+Mu]    = XT[MP1]*QKK[1+MP2]    - XT[MP2]*QKK[1+MP1];
        QNN[4+Mu]    = XT[MP1]*QNN[1+MP2]    - XT[MP2]*QNN[1+MP1];
        QKNmNK[4+Mu] = XT[MP1]*QKNmNK[1+MP2] - XT[MP2]*QKNmNK[1+MP1];
      };

   };
}

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

  PFTIData *Data=(PFTIData *)UserData;
  IncField *IF    = Data->IF;
  double *XTorque = Data->TorqueCenter;
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
  PFTIData MyData, *Data=&MyData;
  Data->Omega        = Omega;
  Data->IF           = IF;
  Data->TorqueCenter = G->Surfaces[ns]->Origin;

  int fDim=4*NUMPFT;
  int MaxEvals=21;
  GetBFCubature(G, ns, ne, PFTIntegrand_BFInc, (void *)Data,
                fDim, MaxEvals, 0, 0, 0, 0, (double *)IBFInc);
}

/***************************************************************/
/* add incident-field contributions (i.e. extinction) to EMTPFT*/
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
/* compute PFT integrals between an RWG basis function and an  */
/* external field                                              */
/***************************************************************/
void GetPFTIntegrals_BFBF(RWGGeometry *G,
                          int nsa, int nea, int nsb, int neb,
                          cdouble Omega, cdouble k, 
                          cdouble EpsR, cdouble MuR,
                          int EMTPFTIMethod, double PFTIs[21])
{
  if (EMTPFTIMethod==SCUFF_EMTPFTI_EHDERIVATIVES)
  {
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
     PFTIs[0*7 + 0] = imag(  MuR*Omega*GabArray[0][GCME_GC] );
     PFTIs[1*7 + 0] = imag( EpsR*Omega*GabArray[0][GCME_GC] );
     PFTIs[2*7 + 0] = imag(    ikCabArray[0][GCME_GC] );
     double FTPreFac=TENTHIRDS;
     for(int Mu=0; Mu<6; Mu++)
      { PFTIs[0*7+1+Mu] = FTPreFac*imag( MuR*GabArray[0][GCME_FX + Mu]);
        PFTIs[1*7+1+Mu] = FTPreFac*imag(EpsR*GabArray[0][GCME_FX + Mu]);
        PFTIs[2*7+1+Mu] = -FTPreFac*imag(  ikCabArray[0][GCME_FX + Mu]);
      };
   }
  else 
   { 
     PFTIData MyData, *Data=&MyData;
     Data->Omega        = Omega;
     Data->k            = k;
     Data->EpsRel       = EpsR;
     Data->MuRel        = MuR;
     Data->TorqueCenter = G->Surfaces[nsa]->Origin;
     Data->EMTPFTIMethod = EMTPFTIMethod;

     int IDim     = 21;
     int MaxEvals = 78;

     Data->G   = G;
     Data->nsb = nsb;
     Data->neb = neb;
     GetBFCubature(G, nsa, nea, PFTIntegrand_BFBF2,
                       (void *)Data, IDim, MaxEvals,
                       0.0, 0.0, 0.0, 0, PFTIs);
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *GetEMTPFTMatrix(RWGGeometry *G, cdouble Omega, IncField *IF,
                         HVector *KNVector, HMatrix *DRMatrix,
                         HMatrix *PFTMatrix, bool Interior,
                         int EMTPFTIMethod)
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

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  static int init=0;
  if (init==0)
   { init=1;
     GetGCMEArgStruct MyArgs;
     InitGetGCMEArgs(&MyArgs);
   };

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
  int TotalEdges = G->TotalEdges;
  bool UseSymmetry=true;
  char *s=getenv("SCUFF_EMTPFT_SYMMETRY");
  if (s && s[0]=='0')
   { UseSymmetry=false;
     Log("Not using symmetry in EMTPFT calculation.");
   };
    
#ifdef USE_OPENMP
  Log("EMT OpenMP multithreading (%i threads)",NT);
#pragma omp parallel for schedule(dynamic,1), num_threads(NT)
#endif
  for(int neaTot=0; neaTot<TotalEdges; neaTot++)
   for(int nebTot=(UseSymmetry ? neaTot : 0); nebTot<TotalEdges; nebTot++)
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
      cdouble u0KK  = KNB[0] * ZVAC;
      cdouble KNmNK = KNB[1] - KNB[2];
      cdouble e0NN  = KNB[3] / ZVAC;
 
      double PFTIs[21];
      double *QKK=PFTIs+0, *QNN=PFTIs+7, *QKNmNK=PFTIs+14;
      GetPFTIntegrals_BFBF(G, nsa, nea, nsb, neb,
                           Omega, k, EpsR, MuR, EMTPFTIMethod, PFTIs);

      double dPFT[NUMPFT];
      dPFT[PFT_PSCAT] = 0.5*(  real(u0KK)*QKK[0]
                              +real(e0NN)*QNN[0]
                              +imag(KNmNK)*QKNmNK[0]
                            );

      for(int Mu=0; Mu<6; Mu++)
       dPFT[PFT_XFORCE + Mu]
        = -0.5*(  imag(u0KK)*QKK[1+Mu]
                 +imag(e0NN)*QNN[1+Mu]
                 +real(KNmNK)*QKNmNK[1+Mu]
               );

      int nt=0;
#ifdef USE_OPENMP
      nt=omp_get_thread_num();
#endif
      int OffsetA = nt*NSNQ + nsa*NQ;
      int OffsetB = nt*NSNQ + nsb*NQ;

      if (!UseSymmetry)
       VecPlusEquals(DeltaPFT + OffsetA, abSign, dPFT, NUMPFT);
      else
       { 
         if (neaTot==nebTot)
          DeltaPFT[ OffsetA + PFT_PSCAT ] += abSign*dPFT[PFT_PSCAT];
         else if (nsa==nsb) // nebTot > neaTot but both on same surface
          { 
            DeltaPFT[ OffsetA + PFT_PSCAT ] += 2.0*abSign*dPFT[PFT_PSCAT];
            for(int Mu=0; Mu<6; Mu++)
             DeltaPFT[ OffsetA + PFT_XFORCE + Mu ] 
              += 2.0*abSign*dPFT[PFT_XFORCE + Mu];
          }
         else // nebTot > neaTot and on different objects
          { 
            DeltaPFT[ OffsetA + PFT_PSCAT ] += abSign*dPFT[PFT_PSCAT];
            DeltaPFT[ OffsetB + PFT_PSCAT ] += baSign*dPFT[PFT_PSCAT];
            for(int Mu=0; Mu<6; Mu++)
             { DeltaPFT[ OffsetA + PFT_XFORCE + Mu ] 
                += abSign*dPFT[PFT_XFORCE + Mu];
               DeltaPFT[ OffsetB + PFT_XFORCE + Mu ] 
                += baSign*dPFT[PFT_XFORCE + Mu];
             };
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
  /* If we are using the interior expansion, then the quantity   */
  /* computed as the scattered power above is actually minus the */
  /* absorbed power.                                             */
  /* If we are using the exterior expansion, then the quantity   */
  /* computed as the scattered power wants to be subtracted from */
  /* the extinction (which will be added to the absorbed power   */
  /* when we compute the incident-field contributions below).    */
  /* So in either case we want to set the PABS slot of the matrix*/
  /* to the negative of PSCAT. In the interior-expansion case    */
  /* we want further to zero out the PSCAT slot since the        */
  /* interior expansion cannot compute the scattered power.      */
  /***************************************************************/
  for(int ns=0; ns<NS; ns++)
   { double PScat = PFTMatrix->GetEntryD(ns, PFT_PSCAT);
     PFTMatrix->SetEntry(ns, PFT_PABS, -PScat);
     if (Interior) PFTMatrix->SetEntry(ns, PFT_PSCAT, 0.0);
   };

  /***************************************************************/
  /* add incident-field contributions ****************************/
  /***************************************************************/
  if (IF)
   AddIFContributionsToEMTPFT(G, KNVector, IF, Omega, PFTMatrix, Interior);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
#if 0
  for(int ns=0; ns<G->NumSurfaces; ns++)
   if (G->Mate[ns]==-1 && G->FIBBICaches[ns])
    { 
      StoreFIBBICache(G->FIBBICaches[ns], G->Surfaces[ns]->MeshFileName);
      int Hits, Misses;
      int CacheSize=GetFIBBICacheSize(G->FIBBICaches[ns],&Hits, &Misses);
      Log("EMTPFT surface %i: (%i/%i) hits/misses",ns,Hits,Misses);
    };
#endif

  return PFTMatrix;
}
  
} // namespace scuff
