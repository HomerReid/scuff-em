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
void GetKNBilinears(HVector *KNVector, HMatrix *DRMatrix,
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
     KK = DRMatrix->GetEntry(KNIndexB+0, KNIndexA+0);
     KN = DRMatrix->GetEntry(KNIndexB+1, KNIndexA+0);
     NK = DRMatrix->GetEntry(KNIndexB+0, KNIndexA+1);
     NN = DRMatrix->GetEntry(KNIndexB+1, KNIndexA+1);
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
  memset(I, 0, 2*NUMPFT*sizeof(double));

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
  double r2=R[0]*R[0] + R[1]*R[1] + R[2]*R[2];
  if (r2==0.0)
   { I[0] = PEFIE * k/(4.0*M_PI);
     return;
   };

  double r   = sqrt(r2);
  double kr  = k*r, coskr=cos(kr), sinkr=sin(kr);
  double ImG = sinkr/(4.0*M_PI*r);
  double f1  = (kr*coskr - sinkr) / (4.0*M_PI*r*r*r);

  I[0 + PFT_PABS]    = PEFIE * ImG;
  I[0 + PFT_XFORCE]  = PEFIE * R[0] * f1;
  I[0 + PFT_YFORCE]  = PEFIE * R[1] * f1;
  I[0 + PFT_ZFORCE]  = PEFIE * R[2] * f1;
  I[0 + PFT_XTORQUE] = PEFIE * (XmXT[1]*R[2]-XmXT[2]*R[1]) * f1;
  I[0 + PFT_YTORQUE] = PEFIE * (XmXT[2]*R[0]-XmXT[0]*R[2]) * f1;
  I[0 + PFT_ZTORQUE] = PEFIE * (XmXT[0]*R[1]-XmXT[1]*R[0]) * f1;

  double k2  = k*k;
  double kr2 = kr*kr;
  double f2  = (kr*coskr + (kr2-1.0)*sinkr) / (4.0*M_PI*kr2*r);
  double f3  = (-3.0*kr*coskr + (3.0-kr2)*sinkr) / (4.0*M_PI*kr2*r*r*r);
  double bxb[3], bAxR[3], bAdotR=0.0, bBdotR=0.0, Rdbxb=0.0;
  for(int Mu=0; Mu<3; Mu++)
   { int MP1  =  (Mu+1)%3, MP2=(Mu+2)%3;
     bxb[Mu]  =  bA[MP1]*bB[MP2] - bA[MP2]*bB[MP1];
     bAxR[Mu] =  bA[MP1]*R[MP2]  - bA[MP2]*R[MP1];
     bAdotR   += bA[Mu]*R[Mu];
     bBdotR   += bB[Mu]*R[Mu];
     Rdbxb    += R[Mu]*bxb[Mu];
   };
 
  I[0 + PFT_XTORQUE] -=  DivbB*(bA[1]*R[2]-bA[2]*R[1])*f1/k2;
  I[0 + PFT_YTORQUE] -=  DivbB*(bA[2]*R[0]-bA[0]*R[2])*f1/k2;
  I[0 + PFT_ZTORQUE] -=  DivbB*(bA[0]*R[1]-bA[1]*R[0])*f1/k2;

  I[0 + PFT_XTORQUE] += f2*bxb[0] + f3*bBdotR*bAxR[0];
  I[0 + PFT_YTORQUE] += f2*bxb[1] + f3*bBdotR*bAxR[1];
  I[0 + PFT_ZTORQUE] += f2*bxb[2] + f3*bBdotR*bAxR[2];

  double g  = f1/k;
  double gP = f3;
  I[NUMPFT + PFT_PABS]    = g*Rdbxb;
  I[NUMPFT + PFT_XFORCE]  = g*bxb[0] + gP*Rdbxb*R[0];
  I[NUMPFT + PFT_YFORCE]  = g*bxb[1] + gP*Rdbxb*R[1];
  I[NUMPFT + PFT_ZFORCE]  = g*bxb[2] + gP*Rdbxb*R[2];
  I[NUMPFT + PFT_XTORQUE] = g*(bAdotR*bB[0] - DotProduct*R[0]);
  I[NUMPFT + PFT_YTORQUE] = g*(bAdotR*bB[1] - DotProduct*R[1]);
  I[NUMPFT + PFT_ZTORQUE] = g*(bAdotR*bB[2] - DotProduct*R[2]);

  I[NUMPFT + PFT_XTORQUE] +=         g*(XmXT[1]*bxb[2]-XmXT[2]*bxb[1])
                             +gP*Rdbxb*(XmXT[1]*  R[2]-XmXT[2]*  R[1]);
  I[NUMPFT + PFT_YTORQUE] +=         g*(XmXT[2]*bxb[0]-XmXT[0]*bxb[2])
                             +gP*Rdbxb*(XmXT[2]*  R[0]-XmXT[0]*  R[2]);
  I[NUMPFT + PFT_ZTORQUE] +=         g*(XmXT[0]*bxb[1]-XmXT[1]*bxb[0])
                             +gP*Rdbxb*(XmXT[0]*  R[1]-XmXT[1]*  R[0]);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetPFTIntegrals_BFBF(RWGGeometry *G,
                          int nsa, int nea, int nsb, int neb,
                          cdouble Omega, double IBFBF[2*NUMPFT])
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

  int fdim=2*NUMPFT;
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

  double XmXT[3];
  VecSub(x, XTorque, XmXT);

  // get fields and derivatives at eval point
  cdouble EH[6], dEH[3][6];
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
     zI[0 + PFT_XTORQUE] += b[Mu]*(XmXT[1]*dEH[2][Mu]-XmXT[2]*dEH[1][Mu]);
     zI[0 + PFT_YTORQUE] += b[Mu]*(XmXT[2]*dEH[0][Mu]-XmXT[0]*dEH[2][Mu]);
     zI[0 + PFT_ZTORQUE] += b[Mu]*(XmXT[0]*dEH[1][Mu]-XmXT[1]*dEH[0][Mu]);

     zI[NUMPFT + PFT_PABS]    += b[Mu]*EH[3+Mu];
     zI[NUMPFT + PFT_XFORCE]  += b[Mu]*dEH[0][3+Mu];
     zI[NUMPFT + PFT_YFORCE]  += b[Mu]*dEH[1][3+Mu];
     zI[NUMPFT + PFT_ZFORCE]  += b[Mu]*dEH[2][3+Mu];
     zI[NUMPFT + PFT_XTORQUE] += b[Mu]*(XmXT[1]*dEH[2][3+Mu]-XmXT[2]*dEH[1][3+Mu]);
     zI[NUMPFT + PFT_YTORQUE] += b[Mu]*(XmXT[2]*dEH[0][3+Mu]-XmXT[0]*dEH[2][3+Mu]);
     zI[NUMPFT + PFT_ZTORQUE] += b[Mu]*(XmXT[0]*dEH[1][3+Mu]-XmXT[1]*dEH[0][3+Mu]);
   };

  zI[0      + PFT_XTORQUE] += b[1]*EH[0+2] - b[2]*EH[0+1];
  zI[0      + PFT_YTORQUE] += b[2]*EH[0+0] - b[0]*EH[0+2];
  zI[0      + PFT_ZTORQUE] += b[0]*EH[0+1] - b[1]*EH[0+0];

  zI[NUMPFT + PFT_XTORQUE] += b[1]*EH[3+2] - b[2]*EH[3+1];
  zI[NUMPFT + PFT_YTORQUE] += b[2]*EH[3+0] - b[0]*EH[3+2];
  zI[NUMPFT + PFT_ZTORQUE] += b[0]*EH[3+1] - b[1]*EH[3+0];

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
                   HVector *KNVector, HVector *RHSVector,
                   HMatrix *DRMatrix, HMatrix *PFTMatrix)
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
  Log("EMT OpenMP multithreading (%i threads)",NT);
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
      GetKNBilinears(KNVector, DRMatrix,
                     SA->IsPEC, KNIndexA, SB->IsPEC, KNIndexB,
                     KNB);
      cdouble KKpNN = Z*KNB[0] + KNB[3]/Z;
      cdouble KNmNK = KNB[1] - KNB[2];
                     
      if (KKpNN==0.0) continue;

      double IBFBF[2*NUMPFT];
      GetPFTIntegrals_BFBF(G, nsa, nea, nsb, neb, Omega, IBFBF);
      double ImG=IBFBF[0],      *ImdG=IBFBF+PFT_XFORCE;
      double ReC=IBFBF[NUMPFT], *RedC=IBFBF+NUMPFT+PFT_XFORCE;
  
      double dPFT[NUMPFT];
      dPFT[PFT_PABS] 
       = -0.5*real(Omega)*( real(KKpNN)*ImG + imag(KNmNK)*ReC );
      for(int Mu=0; Mu<6; Mu++)
       dPFT[PFT_XFORCE + Mu] 
        = -0.5*TENTHIRDS*( imag(KKpNN)*ImdG[Mu] - real(KNmNK)*RedC[Mu] );

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
   AddIFContributionsToEMTPFT(G, KNVector, IF, Omega, PFTMatrix);

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
