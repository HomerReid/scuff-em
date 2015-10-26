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
 * JDEPFT.cc     -- liscuff class methods for computing power, force,
 *               -- and torque in classical deterministic scattering
 *               -- problems using the "J \dot E" formalism
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

/***************************************************************/
/***************************************************************/
/***************************************************************/
RWGSurface *ResolveNE(RWGGeometry *G, int neFull,
                      int *pns, int *pne, int *pKNIndex)
{
  int ns=0, NSm1=G->NumSurfaces - 1;
  while( (ns < NSm1) && (neFull >= G->EdgeIndexOffset[ns+1]) )
   ns++;

  int ne  = neFull - G->EdgeIndexOffset[ns];
  
  int Mult    = G->Surfaces[ns]->IsPEC ? 1 : 2;
  int KNIndex = G->BFIndexOffset[ns] + Mult*ne;

  if (pns) *pns=ns;
  if (pne) *pne=ne;
  if (pKNIndex) *pKNIndex=KNIndex;
  return G->Surfaces[ns];
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetKNBilinears(HVector *KNVector, HMatrix *DMatrix,
                    bool IsPECA, int KNIndexA,
                    bool IsPECB, int KNIndexB,
                    cdouble Bilinears[4])
{
  cdouble KK, KN, NK, NN;
  if (KNVector && (IsPECA || IsPECB) )
   { 
     cdouble kAlpha =  KNVector->GetEntry(KNIndexA);
     cdouble kBeta  =  KNVector->GetEntry(KNIndexB);
     KK = conj(kAlpha) * kBeta;
     KN = NK = NN = 0.0;
   }
  else if (KNVector)
   { 
     cdouble kAlpha =       KNVector->GetEntry(KNIndexA+0);
     cdouble nAlpha = -ZVAC*KNVector->GetEntry(KNIndexA+1);
     cdouble kBeta  =       KNVector->GetEntry(KNIndexB+0);
     cdouble nBeta  = -ZVAC*KNVector->GetEntry(KNIndexB+1);

     KK = conj(kAlpha) * kBeta;
     KN = conj(kAlpha) * nBeta;
     NK = conj(nAlpha) * kBeta;
     NN = conj(nAlpha) * nBeta;
   }
  else
   {
     KK = DMatrix->GetEntry(KNIndexB+0, KNIndexA+0);
     KN = DMatrix->GetEntry(KNIndexB+1, KNIndexA+0);
     NK = DMatrix->GetEntry(KNIndexB+0, KNIndexA+1);
     NN = DMatrix->GetEntry(KNIndexB+1, KNIndexA+1);
   }

  Bilinears[0]=KK;
  Bilinears[1]=KN;
  Bilinears[2]=NK;
  Bilinears[3]=NN;
}

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
void PFTIntegrand_BFBF(double *xA, double *xB, PPCData *PPCD,
                       void *UserData, double *I)
{
  double bA[3], DivbA, bB[3], DivbB;
  DivbA = real(PPCD->DivK1);
  bA[0] = real(PPCD->K1[0]);
  bA[1] = real(PPCD->K1[1]);
  bA[2] = real(PPCD->K1[2]);
  DivbB = real(PPCD->DivK2);
  bB[0] = real(PPCD->K2[0]);
  bB[1] = real(PPCD->K2[1]);
  bB[2] = real(PPCD->K2[2]);

  PFTIntegrandData *PFTIData = (PFTIntegrandData *)UserData;
  double k                   = real(PFTIData->k);
  double *XTorque            = PFTIData->XTorque;

  double XmXT[3];
  VecSub(xA, XTorque, XmXT);

  double DotProduct    = bA[0]*bB[0] + bA[1]*bB[1] + bA[2]*bB[2];
  double ScalarProduct = DivbA*DivbB;
  double PEFIE         = DotProduct - ScalarProduct/(k*k);

  double R[3]; 
  VecSub(xA, xB, R);
  double r=VecNorm(R);
  if (r==0.0)
   { 
     memset(I, 0, 14*sizeof(double));
     I[0] = PEFIE * k/(4.0*M_PI);
     return;
   };

  cdouble G[3][3], C[3][3], dG[3][3][3], dC[3][3][3];
  cdouble EpsRel=1.0, MuRel=1.0; // FIXME for non-vacuum exterior medium
  CalcGC(R, k, EpsRel, MuRel, G, C, dG, dC);
  double ImGbB[3], ImdGbB[3][3];
  double ReCbB[3], RedCbB[3][3];
  for(int Mu=0; Mu<3; Mu++)
   { ImGbB[Mu]=ReCbB[Mu]=0.0;
     ImdGbB[0][Mu]=ImdGbB[1][Mu]=ImdGbB[2][Mu]=0.0;
     RedCbB[0][Mu]=RedCbB[1][Mu]=RedCbB[2][Mu]=0.0;

     for(int Nu=0; Nu<3; Nu++)
      { ImGbB[Mu]     += imag(G[Mu][Nu])*bB[Nu];
        ImdGbB[0][Mu] += imag(dG[Mu][Nu][0])*bB[Nu];
        ImdGbB[1][Mu] += imag(dG[Mu][Nu][1])*bB[Nu];
        ImdGbB[2][Mu] += imag(dG[Mu][Nu][2])*bB[Nu];

        ReCbB[Mu]     += real(C[Mu][Nu])*bB[Nu];
        RedCbB[0][Mu] += real(dC[Mu][Nu][0])*bB[Nu];
        RedCbB[1][Mu] += real(dC[Mu][Nu][1])*bB[Nu];
        RedCbB[2][Mu] += real(dC[Mu][Nu][2])*bB[Nu];
      };
   };

  memset(I, 0, 14*sizeof(double));
  for(int Mu=0; Mu<3; Mu++)
   { 
     I[0]      += bA[Mu]*ImGbB[Mu];

     I[1 + 0]  += bA[Mu]*ImdGbB[0][Mu];
     I[1 + 1]  += bA[Mu]*ImdGbB[1][Mu];
     I[1 + 2]  += bA[Mu]*ImdGbB[2][Mu];

     I[4 + 0]  += bA[Mu]*(XmXT[1]*ImdGbB[2][Mu]-XmXT[2]*ImdGbB[1][Mu]);
     I[4 + 1]  += bA[Mu]*(XmXT[2]*ImdGbB[0][Mu]-XmXT[0]*ImdGbB[2][Mu]);
     I[4 + 2]  += bA[Mu]*(XmXT[0]*ImdGbB[1][Mu]-XmXT[1]*ImdGbB[0][Mu]);
 
     I[7]      += bA[Mu]*ReCbB[Mu];

     I[8 + 0]  += bA[Mu]*RedCbB[0][Mu];
     I[8 + 1]  += bA[Mu]*RedCbB[1][Mu];
     I[8 + 2]  += bA[Mu]*RedCbB[2][Mu];

     I[11+ 0]  += bA[Mu]*(XmXT[1]*RedCbB[2][Mu]-XmXT[2]*RedCbB[1][Mu]);
     I[11+ 1]  += bA[Mu]*(XmXT[2]*RedCbB[0][Mu]-XmXT[0]*RedCbB[2][Mu]);
     I[11+ 2]  += bA[Mu]*(XmXT[0]*RedCbB[1][Mu]-XmXT[1]*RedCbB[0][Mu]);
 
   };

  I[4 + 0] += bA[1]*ImGbB[2] - bA[2]*ImGbB[1];
  I[4 + 1] += bA[2]*ImGbB[0] - bA[0]*ImGbB[2];
  I[4 + 2] += bA[0]*ImGbB[1] - bA[1]*ImGbB[0];

  I[11+ 0] += bA[1]*ReCbB[2] - bA[2]*ReCbB[1];
  I[11+ 1] += bA[2]*ReCbB[0] - bA[0]*ReCbB[2];
  I[11+ 2] += bA[0]*ReCbB[1] - bA[1]*ReCbB[0];

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetPFTIntegrals_BFBF(RWGGeometry *G,
                          int nsa, int nea, int nsb, int neb,
                          cdouble Omega, double IBFBF[14])
{
  PFTIntegrandData MyPFTIData, *PFTIData=&MyPFTIData;
  PFTIData->k  = Omega;
  PFTIData->XTorque[0]=PFTIData->XTorque[1]=PFTIData->XTorque[2]=0.0; 

  RWGSurface *Sa = G->Surfaces[nsa];
  RWGSurface *Sb = G->Surfaces[nsb];

  if (Sa->OTGT)
   Sa->OTGT->Apply(PFTIData->XTorque);
  if (Sa->GT)
   Sa->GT->Apply(PFTIData->XTorque);

  int fdim=14;
  int ncv = NumCommonBFVertices(Sa, nea, Sb, neb);
  int NumPts = (ncv > 0) ? 441 : 36;
  GetBFBFCubature(G, nsa, nea, nsb, neb, 0,
                  PFTIntegrand_BFBF, (void *)PFTIData, fdim,
                  NumPts, 0, 0, 0, 0, IBFBF);
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

  PFTIntegrandData *PFTIData=(PFTIntegrandData *)UserData;
  cdouble k       = PFTIData->k;
  double kAbs     = abs(k);
  IncField *IF    = PFTIData->IF;
  double *XTorque = PFTIData->XTorque;

  // get fields and derivatives at eval point by finite-differencing
  cdouble EH[6], *E=EH+0, *H=EH+3, dE[3][3], dH[3][3];
  IF->GetFields(x, EH);
  for(int Mu=0; Mu<3; Mu++)
   { 
     double xTweaked[3];
     xTweaked[0]=x[0];
     xTweaked[1]=x[1];
     xTweaked[2]=x[2];

     double Delta = (x[Mu]==0.0 ? 1.0e-4 : 1.0e-4*fabs(x[Mu]));
     if ( kAbs > 1.0 )
      Delta = fmin(Delta, 1.0e-4/kAbs);

     xTweaked[Mu] += Delta;
     cdouble EHP[6]; 
     IF->GetFields(xTweaked, EHP);
     xTweaked[Mu] -= 2.0*Delta;
     cdouble EHM[6]; 
     IF->GetFields(xTweaked, EHM);

     for(int Nu=0; Nu<3; Nu++)
      { dE[Mu][Nu] = (EHP[0+Nu]-EHM[0+Nu])/(2.0*Delta);
        dH[Mu][Nu] = (EHP[3+Nu]-EHM[3+Nu])/(2.0*Delta);
      };
   };

  double XmXT[3];
  VecSub(x, XTorque, XmXT);

  cdouble *zI = (cdouble *)I;
  memset(zI, 0, 14*sizeof(cdouble));
  for(int Mu=0; Mu<3; Mu++)
   { 
     zI[0]  += b[Mu]*E[Mu];
     zI[1]  += b[Mu]*dE[0][Mu];
     zI[2]  += b[Mu]*dE[1][Mu];
     zI[3]  += b[Mu]*dE[2][Mu];
     zI[4]  += b[Mu]*(XmXT[1]*dE[2][Mu]-XmXT[2]*dE[1][Mu]);
     zI[5]  += b[Mu]*(XmXT[2]*dE[0][Mu]-XmXT[0]*dE[2][Mu]);
     zI[6]  += b[Mu]*(XmXT[0]*dE[1][Mu]-XmXT[1]*dE[0][Mu]);

     zI[7]  += b[Mu]*H[Mu];
     zI[8]  += b[Mu]*dH[0][Mu];
     zI[9]  += b[Mu]*dH[1][Mu];
     zI[10] += b[Mu]*dH[2][Mu];
     zI[11] += b[Mu]*(XmXT[1]*dH[2][Mu]-XmXT[2]*dH[1][Mu]);
     zI[12] += b[Mu]*(XmXT[2]*dH[0][Mu]-XmXT[0]*dH[2][Mu]);
     zI[13] += b[Mu]*(XmXT[0]*dH[1][Mu]-XmXT[1]*dH[0][Mu]);
   };

  zI[4 + 0] += b[1]*E[2] - b[2]*E[1];
  zI[4 + 1] += b[2]*E[0] - b[0]*E[2];
  zI[4 + 2] += b[0]*E[1] - b[1]*E[0];

  zI[11 + 0] += b[1]*H[2] - b[2]*H[1];
  zI[11 + 1] += b[2]*H[0] - b[0]*H[2];
  zI[11 + 2] += b[0]*H[1] - b[1]*H[0];

}

/***************************************************************/
/* compute PFT integrals between an RWG basis function and an  */
/* external field                                              */
/***************************************************************/
void GetPFTIntegrals_BFInc(RWGGeometry *G, int ns, int ne,
                           IncField *IF, cdouble Omega,
                           cdouble IBFInc[14])
{
  PFTIntegrandData MyPFTIData, *PFTIData=&MyPFTIData;
  PFTIData->k  = Omega;
  PFTIData->IF = IF;
  PFTIData->XTorque[0]=PFTIData->XTorque[1]=PFTIData->XTorque[2]=0.0;

  RWGSurface *S=G->Surfaces[ns];
  if (S->OTGT) S->OTGT->Apply(PFTIData->XTorque);
  if (S->GT) S->GT->Apply(PFTIData->XTorque);

  int fDim=28;
  int MaxEvals=21;
  GetBFCubature(G, ns, ne,
                PFTIntegrand_BFInc, (void *)PFTIData, fDim,
                MaxEvals, 0, 0, 0, 0, (double *)IBFInc);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddIFContributionsToJDEPFT(RWGGeometry *G, HVector *KN,
                                IncField *IF, cdouble Omega,
                                HMatrix *PFTMatrix)
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
     RWGSurface *S = ResolveNE(G, neTot, &ns, &ne, &KNIndex);

     cdouble KStar = conj(KN->GetEntry(KNIndex)); 
     cdouble NStar = 0.0;
     if (!(S->IsPEC)) NStar = -ZVAC*conj(KN->GetEntry(KNIndex+1));

     cdouble IBFInc[14];
     GetPFTIntegrals_BFInc(G, ns, ne, IF, Omega, IBFInc);
     cdouble bDotE = IBFInc[0], *bDotdE=IBFInc+1;
     cdouble bDotH = IBFInc[7], *bDotdH=IBFInc+8;

     int nt=0;
#ifdef USE_OPENMP
     nt = omp_get_thread_num();
#endif
     double *dPFT = PartialPFT + (nt*NS + ns)*NUMPFT;
     dPFT[0] += real(KStar*bDotE + NStar*bDotH);
     for(int Mu=0; Mu<6; Mu++)
      dPFT[1+Mu] += imag( KStar*bDotdE[Mu] + NStar*bDotdH[Mu] );
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
       PFTMatrix->AddEntry(ns, PFT_PABS, PFactor*dPFT[0]);
       for(int Mu=0; Mu<6; Mu++)
        PFTMatrix->AddEntry(ns, PFT_XFORCE + Mu, FTFactor*dPFT[1+Mu]);
     };

  delete[] PartialPFT;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *GetJDEPFT(RWGGeometry *G, cdouble Omega, IncField *IF,
                   HVector *KNVector, HVector *RHSVector,
                   HMatrix *DMatrix, HMatrix *PFTMatrix)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NS = G->NumSurfaces;
  if (    (PFTMatrix==0)
       || (PFTMatrix->NR != NS)
       || (PFTMatrix->NC != NUMPFT)
     )
   ErrExit("invalid PFTMatrix in GetJDEPFT");
 
  double Z = ZVAC; // FIXME for non-vacuum exterior media

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
  /*- multithreaded loop over all basis functions in all volumes -*/
  /*--------------------------------------------------------------*/
  int TotalEdges  = G->TotalEdges;
#ifdef USE_OPENMP
  Log("JDE OpenMP multithreading (%i threads)",NT);
#pragma omp parallel for schedule(dynamic,1), num_threads(NT)
#endif
  for(int neaTot=0; neaTot<TotalEdges; neaTot++)
   for(int nebTot=neaTot; nebTot<TotalEdges; nebTot++)
    { 
      if (nebTot==neaTot) LogPercent(neaTot, TotalEdges, 10);

      int nsa, nea, KNIndexA;
      RWGSurface *SA = ResolveNE(G, neaTot, &nsa, &nea, &KNIndexA);

      int nsb, neb, KNIndexB;
      RWGSurface *SB = ResolveNE(G, nebTot, &nsb, &neb, &KNIndexB);
   
      cdouble KNB[4];
      GetKNBilinears(KNVector, DMatrix,
                     SA->IsPEC, KNIndexA, SB->IsPEC, KNIndexB,
                     KNB);
      cdouble KKpNN = Z*KNB[0] + KNB[3]/Z;
      cdouble KNmNK = KNB[1] - KNB[2];
                     
      if (KKpNN==0.0) continue;

      double IBFBF[14];
      GetPFTIntegrals_BFBF(G, nsa, nea, nsb, neb, Omega, IBFBF);
      double ImG=IBFBF[0], *ImdG=IBFBF+1;
      double ReC=IBFBF[7], *RedC=IBFBF+8;
  
      double dPFT[NUMPFT];
      dPFT[PFT_PABS] 
       = -0.5*real(Omega)*( real(KKpNN)*ImG + imag(KNmNK)*ReC );
      for(int Mu=0; Mu<6; Mu++)
       dPFT[PFT_XFORCE + Mu] 
        = -0.5*TENTHIRDS*( imag(KKpNN)*ImdG[Mu] + real(KNmNK)*RedC[Mu] );

      int nt=0;
#ifdef USE_OPENMP
      nt=omp_get_thread_num();
#endif
      int OffsetA = nt*NSNQ + nsa*NQ;
      int OffsetB = nt*NSNQ + nsb*NQ;

       if (neaTot==nebTot)
        DeltaPFT[ OffsetA + PFT_PABS ] += dPFT[PFT_PABS];
       else if (nsa==nsb) // nebTot > neaTot but both on same surface
        { 
          DeltaPFT[ OffsetA + PFT_PABS ] += 2.0*dPFT[PFT_PABS];
          for(int Mu=0; Mu<6; Mu++)
           DeltaPFT[ OffsetA + PFT_XFORCE + Mu ] 
            += 2.0*dPFT[PFT_XFORCE + Mu];
        }
       else // nebTot > neaTot and on different objects
        { 
          DeltaPFT[ OffsetA + PFT_PABS ] += dPFT[PFT_PABS];
          DeltaPFT[ OffsetB + PFT_PABS ] += dPFT[PFT_PABS];
          for(int Mu=0; Mu<6; Mu++)
           { DeltaPFT[ OffsetA + PFT_XFORCE + Mu ] 
              += dPFT[PFT_XFORCE + Mu];
             DeltaPFT[ OffsetB + PFT_XFORCE + Mu ] 
              += dPFT[PFT_XFORCE + Mu];
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
  /* add incident-field contributions ****************************/
  /***************************************************************/
  if (IF)
   AddIFContributionsToJDEPFT(G, KNVector, IF, Omega, PFTMatrix);

  /*--------------------------------------------------------------*/
  /*- if an RHS vector was specified, compute the extinction      */
  /*- (total power) and use it to compute the scattered power     */
  /*--------------------------------------------------------------*/
  if (KNVector && RHSVector)
   for (int ns=0; ns<NS; ns++)
    { int NE=G->Surfaces[ns]->NumEdges;
      bool IsPEC = G->Surfaces[ns]->IsPEC;
      int Offset = G->BFIndexOffset[ns];
      double Extinction=0.0;
      for (int ne=0, nbf=0; ne<NE; ne++)
       {
         cdouble kAlpha =   KNVector->GetEntry(Offset + nbf);
         cdouble vEAlpha = -ZVAC*RHSVector->GetEntry(Offset + nbf);
         nbf++;
         Extinction += 0.5*real( conj(kAlpha)*vEAlpha );
         if (IsPEC) continue;

         cdouble nAlpha  = -ZVAC*KNVector->GetEntry(Offset + nbf);
         cdouble vHAlpha = -1.0*RHSVector->GetEntry(Offset + nbf);
         nbf++;
         Extinction += 0.5*real( conj(nAlpha)*vHAlpha );
       };
      double PAbs = PFTMatrix->GetEntryD(ns, PFT_PABS);
      PFTMatrix->SetEntry(ns, PFT_PSCAT, Extinction-PAbs);
   };

  return PFTMatrix;

}
  
} // namespace scuff
