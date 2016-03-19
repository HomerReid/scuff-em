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

#define NUMFT    6
#define NUMFTT   9
#define NUMPFTT  11

#define NUMPFTIS (3*NUMPFTT)

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

inline cdouble RCVecDot(double *V1, cdouble *V2)
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
void PFTIIntegrand_BFBF(double xA[3], double bA[3], double DivbA,
                        double xB[3], double bB[3], double DivbB,
                        void *UserData, double Weight, double *I)
{
  PFTIData *Data  = (PFTIData *)UserData;
  double Omega    = real(Data->Omega);
  double k        = real(Data->k);
  double EpsRel   = real(Data->EpsRel);
  double MuRel    = real(Data->MuRel);
  double *XTorque = Data->TorqueCenter;

  /***************************************************************/
  /* kernel factors **********************************************/
  /***************************************************************/
  double R[3];
  VecSub(xA, xB, R);
  double r2=R[0]*R[0] + R[1]*R[1] + R[2]*R[2];
  double r = sqrt(r2);
  double kr = k*r, kr2=kr*kr;
  double k2 = k*k;

  double ImPhi, ImPsi, ImZeta;
  if ( fabs(kr) < 1.0e-3 )
   { double k3=k2*k, k5=k3*k2;
     ImPhi  =  (1.0 - kr2/6.0)  * k/(4.0*M_PI);
     ImPsi  = -(1.0 - kr2/10.0) * k3/(12.0*M_PI);
     ImZeta =  (1.0 - kr2/14.0) * k5/(60.0*M_PI);
   }
  else
   { double CosKR = cos(kr), SinKR=sin(kr), r3=r*r2, r5=r3*r2;
     ImPhi  = SinKR/(4.0*M_PI*r);
     ImPsi  = (kr*CosKR - SinKR)/(4.0*M_PI*r3);
     ImZeta = (-3.0*kr*CosKR + (3.0-kr2)*SinKR)/(4.0*M_PI*r5);
   };

  /***************************************************************/
  /* polynomial factors ******************************************/
  /***************************************************************/
  double bdb  = VecDot(bA, bB);
  double DbDb = DivbA*DivbB;
  double bAdR = VecDot(bA, R);
  double bBdR = VecDot(bB, R);

  double bxb[3];    VecCross(bA, bB, bxb);
  double bAxR[3];   VecCross(bA, R, bAxR);

  double XT[3];     VecSub(xA, XTorque, XT);
  double XTxR[3];   VecCross(XT, R,   XTxR);
  double XTxbxb[3]; VecCross(XT, bxb, XTxbxb);

  double PEFIE = bdb - DbDb/(k*k);
  double PMFIE = VecDot(bxb, R);

  double *QKK    = I+0*NUMPFTT;
  double *QNN    = I+1*NUMPFTT;
  double *QKNmNK = I+2*NUMPFTT;

  QKK[PFT_PABS]    -= Weight*PEFIE*Omega*MuRel*ImPhi;
  QNN[PFT_PABS]    -= Weight*PEFIE*Omega*EpsRel*ImPhi;
  QKNmNK[PFT_PABS] -= Weight*PMFIE*ImPsi;
 
  double FTPreFac=TENTHIRDS*Weight;
  for(int Mu=0; Mu<3; Mu++)
   { 
     QKK[PFT_XFORCE + Mu]    -= FTPreFac*PEFIE*MuRel*R[Mu]*ImPsi;
     QNN[PFT_XFORCE + Mu]    -= FTPreFac*PEFIE*EpsRel*R[Mu]*ImPsi;
     QKNmNK[PFT_XFORCE + Mu] += FTPreFac*(bxb[Mu]*ImPsi + PMFIE*R[Mu]*ImZeta)/Omega;

     QKK[PFT_XTORQUE + Mu ]   -= FTPreFac*MuRel*( bxb[Mu]*(ImPhi + ImPsi/k2) + bAxR[Mu]*bBdR*ImZeta/k2);
     QNN[PFT_XTORQUE + Mu ]   -= FTPreFac*EpsRel*( bxb[Mu]*(ImPhi + ImPsi/k2) + bAxR[Mu]*bBdR*ImZeta/k2);
     QKNmNK[PFT_XTORQUE + Mu] += FTPreFac*(bAdR*bB[Mu] - R[Mu]*bdb)*ImPsi/Omega;

     QKK[PFT_XTORQUE + 3 + Mu]    -= FTPreFac*PEFIE*MuRel*XTxR[Mu]*ImPsi;
     QNN[PFT_XTORQUE + 3 + Mu]    -= FTPreFac*PEFIE*EpsRel*XTxR[Mu]*ImPsi;
     QKNmNK[PFT_XTORQUE + 3 + Mu] += FTPreFac*(XTxbxb[Mu]*ImPsi + PMFIE*XTxR[Mu]*ImZeta)/Omega;
   };

}

  
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
  cdouble Omega     = Data->Omega;
  cdouble k         = Data->k;
  cdouble EpsRel    = Data->EpsRel;
  cdouble MuRel     = Data->MuRel;
  RWGGeometry *G    = Data->G;
  int EMTPFTIMethod = Data->EMTPFTIMethod;
  int nsb           = Data->nsb;
  int neb           = Data->neb;
  double *XTorque   = Data->TorqueCenter;

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

  double *QKK    = I+0*NUMPFTT;
  double *QNN    = I+1*NUMPFTT;
  double *QKNmNK = I+2*NUMPFTT;

  QKK[0]    =  -imag(Omega*MuRel*RCVecDot(bA, e));
  QNN[0]    =  -imag(Omega*EpsRel*RCVecDot(bA, e));
  QKNmNK[0] =  +imag(RCVecDot(bA, h));
  if (EMTPFTIMethod==SCUFF_EMTPFTI_EHVALUES1)
   {
     for(int Mu=0; Mu<3; Mu++)
      { QKK[1+Mu]    =  -imag(MuRel*RCVecDot(bA, de[Mu]));
        QNN[1+Mu]    =  -imag(EpsRel*RCVecDot(bA, de[Mu]));
        QKNmNK[1+Mu] =  -imag(RCVecDot(bA, dh[Mu])/Omega);

        QKK[4+Mu]    =  -imag(MuRel*bxe[Mu]);
        QNN[4+Mu]    =  -imag(EpsRel*bxe[Mu]);
        QKNmNK[4+Mu] =  -imag(bxh[Mu]/Omega);

        int MP1 = (Mu+1)%3, MP2 = (Mu+2)%3;
        QKK[7+Mu]    =  XT[MP1]*QKK[1+MP2] - XT[MP2]*QKK[1+MP1];
        QNN[7+Mu]    =  XT[MP1]*QNN[1+MP2] - XT[MP2]*QNN[1+MP1];
        QKNmNK[7+Mu] =  XT[MP1]*QKNmNK[1+MP2] - XT[MP2]*QKNmNK[1+MP1];
      };
   }
  else
   {
     for(int Mu=0; Mu<3; Mu++)
      { QKK[1+Mu]    = -imag(MuRel*(-DivbA*e[Mu] + bxh[Mu]));
        QNN[1+Mu]    = -imag(EpsRel*(-DivbA*e[Mu] + bxh[Mu]));
        QKNmNK[1+Mu] = -imag(DivbA*h[Mu] - k*k*bxe[Mu]);
      };

     for(int Mu=0; Mu<3; Mu++)
      { 
        int MP1=(Mu+1)%3, MP2=(Mu+2)%3;

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
  memset(zI, 0, 2*NUMPFTT*sizeof(cdouble));
  for(int Mu=0; Mu<3; Mu++)
   { 
     zI[0 + PFT_PABS]    += b[Mu]*EH[Mu];
     zI[0 + PFT_XFORCE]  += b[Mu]*dEH[0][Mu];
     zI[0 + PFT_YFORCE]  += b[Mu]*dEH[1][Mu];
     zI[0 + PFT_ZFORCE]  += b[Mu]*dEH[2][Mu];
     zI[0 + PFT_XTORQUE+3] += b[Mu]*(XT[1]*dEH[2][Mu]-XT[2]*dEH[1][Mu]);
     zI[0 + PFT_YTORQUE+3] += b[Mu]*(XT[2]*dEH[0][Mu]-XT[0]*dEH[2][Mu]);
     zI[0 + PFT_ZTORQUE+3] += b[Mu]*(XT[0]*dEH[1][Mu]-XT[1]*dEH[0][Mu]);

     int MP3=Mu+3;
     zI[NUMPFTT + PFT_PABS]    += b[Mu]*EH[MP3];
     zI[NUMPFTT + PFT_XFORCE]  += b[Mu]*dEH[0][MP3];
     zI[NUMPFTT + PFT_YFORCE]  += b[Mu]*dEH[1][MP3];
     zI[NUMPFTT + PFT_ZFORCE]  += b[Mu]*dEH[2][MP3];
     zI[NUMPFTT + PFT_XTORQUE+3] += b[Mu]*(XT[1]*dEH[2][MP3]-XT[2]*dEH[1][MP3]);
     zI[NUMPFTT + PFT_YTORQUE+3] += b[Mu]*(XT[2]*dEH[0][MP3]-XT[0]*dEH[2][MP3]);
     zI[NUMPFTT + PFT_ZTORQUE+3] += b[Mu]*(XT[0]*dEH[1][MP3]-XT[1]*dEH[0][MP3]);
   };

  zI[0      + PFT_XTORQUE] = b[1]*E[2] - b[2]*E[1];
  zI[0      + PFT_YTORQUE] = b[2]*E[0] - b[0]*E[2];
  zI[0      + PFT_ZTORQUE] = b[0]*E[1] - b[1]*E[0];

  zI[NUMPFTT + PFT_XTORQUE] = b[1]*H[2] - b[2]*H[1];
  zI[NUMPFTT + PFT_YTORQUE] = b[2]*H[0] - b[0]*H[2];
  zI[NUMPFTT + PFT_ZTORQUE] = b[0]*H[1] - b[1]*H[0];

}

/***************************************************************/
/* compute PFT integrals between an RWG basis function and an  */
/* external field                                              */
/***************************************************************/
void GetPFTIntegrals_BFInc(RWGGeometry *G, int ns, int ne,
                           IncField *IF, cdouble Omega,
                           cdouble IBFInc[2*NUMPFTT])
{
  PFTIData MyData, *Data=&MyData;
  Data->Omega        = Omega;
  Data->IF           = IF;
  Data->TorqueCenter = G->Surfaces[ns]->Origin;

  int fDim=4*NUMPFTT;
  int MaxEvals=21;
  GetBFCubature(G, ns, ne, PFTIntegrand_BFInc, (void *)Data,
                fDim, MaxEvals, 0, 0, 0, 0, (double *)IBFInc);
}

/***************************************************************/
/* add incident-field contributions (i.e. extinction) to EMTPFT*/
/***************************************************************/
void AddIFContributionsToEMTPFTT(RWGGeometry *G, HVector *KN,
                                 IncField *IF, cdouble Omega,
                                 HMatrix *PFTTMatrix, bool Interior)
{
  if ( PFTTMatrix->NR!=G->NumSurfaces || PFTTMatrix->NC != NUMPFTT )
   ErrExit("%s:%i: internal error", __FILE__, __LINE__);

  int NT=1;
#ifdef USE_OPENMP
  NT = GetNumThreads();
#endif
  int NS=G->NumSurfaces;
  int NQ=NUMPFTT;
  double *PartialPFTT=new double[NT*NS*NQ];
  memset(PartialPFTT, 0, NT*NS*NQ*sizeof(double));

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

     cdouble IBFInc[2*NUMPFTT];
     GetPFTIntegrals_BFInc(G, ns, ne, IF, Omega, IBFInc);
     cdouble bDotE = IBFInc[0],       *bDotdE=IBFInc+PFT_XFORCE;
     cdouble bDotH = IBFInc[NUMPFTT], *bDotdH=IBFInc+NUMPFTT+PFT_XFORCE;

     int nt=0;
#ifdef USE_OPENMP
     nt = omp_get_thread_num();
#endif
     double *dPFTT = PartialPFTT + (nt*NS + ns)*NUMPFTT;
     dPFTT[PFT_PABS] += real(KStar*bDotE + NStar*bDotH);
     for(int Mu=0; Mu<NUMFTT; Mu++)
      dPFTT[PFT_XFORCE+Mu] += imag( KStar*bDotdE[Mu] + NStar*bDotdH[Mu] );
   };

   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   double PFactor = 0.5;
   double FTFactor = 0.5*TENTHIRDS/real(Omega);
   for(int nt=0; nt<NT; nt++)
    for(int ns=0; ns<NS; ns++)
     { 
       double *dPFTT = PartialPFTT + (nt*NS + ns)*NUMPFTT;
       PFTTMatrix->AddEntry(ns, PFT_PABS, PFactor*dPFTT[PFT_PABS]);
       for(int Mu=0; Mu<NUMFTT; Mu++)
        PFTTMatrix->AddEntry(ns, PFT_XFORCE + Mu, FTFactor*dPFTT[PFT_XFORCE+Mu]);
     };

  delete[] PartialPFTT;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddIFContributionsToEMTPFT(RWGGeometry *G, HVector *KN,
                                IncField *IF, cdouble Omega,
                                HMatrix *PFTMatrix, bool Interior)
{
  int NS = PFTMatrix->NR;

  static HMatrix *PFTTMatrix=0;
  if (PFTTMatrix==0 || PFTTMatrix->NR != NS )
   { if (PFTTMatrix) delete PFTTMatrix;
     PFTTMatrix=new HMatrix(NS, NUMPFTT);
   };

  AddIFContributionsToEMTPFTT(G, KN, IF, Omega, PFTTMatrix, Interior);
  for(int ns=0; ns<NS; ns++)
   for(int nq=0; nq<NUMPFT; nq++)
    { PFTMatrix->AddEntry(ns, nq, PFTTMatrix->GetEntry(ns, nq));
      if (PFT_XTORQUE <=nq && nq<=PFT_ZTORQUE)
       PFTMatrix->AddEntry(ns, nq, PFTTMatrix->GetEntry(ns, nq+3));
    };

}

/***************************************************************/
/* compute PFT integrals between an RWG basis function and an  */
/* external field                                              */
/***************************************************************/
void GetPFTIntegrals_BFBF(RWGGeometry *G,
                          int nsa, int nea, int nsb, int neb,
                          cdouble Omega, cdouble k, 
                          cdouble EpsR, cdouble MuR,
                          int EMTPFTIMethod, double PFTIs[NUMPFTIS])
{
  RWGSurface *SA=G->Surfaces[nsa];
  RWGSurface *SB=G->Surfaces[nsb];
  double rRel;
  int ncv=AssessBFPair(SA, nea, SB, neb, &rRel);

  if (EMTPFTIMethod==SCUFF_EMTPFTI_EHDERIVATIVES)
   {
     PFTIData MyData, *Data=&MyData;
     Data->Omega        = Omega;
     Data->k            = k;
     Data->EpsRel       = EpsR;
     Data->MuRel        = MuR;
     Data->TorqueCenter = G->Surfaces[nsa]->Origin;

     int IDim  = NUMPFTIS;
     int Order = (ncv>0) ? 9 : 4;
 
     GetBFBFCubature2(G, nsa, nea, nsb, neb,
                      PFTIIntegrand_BFBF, (void *)Data, IDim,
                      Order, PFTIs);
   }
  else if (EMTPFTIMethod==SCUFF_EMTPFTI_EHDERIVATIVES2)
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
     PFTIs[0*NUMPFTT + 0] = -imag(  MuR*Omega*GabArray[0][GCME_GC] );
     PFTIs[1*NUMPFTT + 0] = -imag( EpsR*Omega*GabArray[0][GCME_GC] );
     PFTIs[2*NUMPFTT + 0] = -imag(    ikCabArray[0][GCME_GC] );
     double FTPreFac=TENTHIRDS;
     for(int Mu=0; Mu<NUMFT; Mu++)
      { PFTIs[0*NUMPFTT+1+Mu] = -FTPreFac*imag( MuR*GabArray[0][GCME_FX + Mu]);
        PFTIs[1*NUMPFTT+1+Mu] = -FTPreFac*imag(EpsR*GabArray[0][GCME_FX + Mu]);
        PFTIs[2*NUMPFTT+1+Mu] = +FTPreFac*imag(  ikCabArray[0][GCME_FX + Mu]);
      };
     PFTIs[0*NUMPFTT+8]=PFTIs[0*NUMPFTT+9]=PFTIs[0*NUMPFTT+10]=0.0;
     PFTIs[1*NUMPFTT+8]=PFTIs[1*NUMPFTT+9]=PFTIs[1*NUMPFTT+10]=0.0;
     PFTIs[2*NUMPFTT+8]=PFTIs[2*NUMPFTT+9]=PFTIs[2*NUMPFTT+10]=0.0;
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

     int IDim     = NUMPFTIS;
     int MaxEvals = (ncv>0) ? 78 : 21;

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
  char *sss=getenv("SCUFF_EMTPFTI_METHOD");
  if (sss)
   sscanf(sss,"%i",&EMTPFTIMethod); 

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
  /* PFTTBySurface[ns] = contributions of surface #ns to PFTT    */
  /***************************************************************/
  static int NSSave=0;
  static HMatrix **PFTTBySurface=0, *IncidentPFTT=0;
  if (NSSave!=NS)
   { if (PFTTBySurface)
      { for(int ns=0; ns<NSSave; ns++)
         if (PFTTBySurface[ns]) 
          delete PFTTBySurface[ns];
        free(PFTTBySurface);
       if (IncidentPFTT)
        delete IncidentPFTT;
      };
     NSSave=NS;
     PFTTBySurface=(HMatrix **)mallocEC(NS*sizeof(HMatrix));
     for(int ns=0; ns<NS; ns++)
      PFTTBySurface[ns]=new HMatrix(NS, NUMPFTT);
     IncidentPFTT=new HMatrix(NS, NUMPFTT);
   };

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
  /*- DeltaPFT[ nt*NS*NS*NQ + nsA*NS*NQ + nsB*NQ + nq ]           */
  /*-  = contribution of surface B to PFT quantity nq on surface A*/
  /*--------------------------------------------------------------*/
  int NT=1;
#ifdef USE_OPENMP 
  NT = GetNumThreads();
#endif
  int NQ      = NUMPFTT;
  int NS2NQ   = NS*NS*NQ;
  int NTNS2NQ = NT*NS*NS*NQ; 
  static int DeltaPFTTSize=0;
  static double *DeltaPFTT=0;
  if ( DeltaPFTTSize < NTNS2NQ )
   { Log("(re)allocating DeltaPFTT (%i,%i)",DeltaPFTTSize,NTNS2NQ);
     DeltaPFTTSize=NTNS2NQ;
     if (DeltaPFTT) free(DeltaPFTT);
     DeltaPFTT = (double *)mallocEC(DeltaPFTTSize*sizeof(double));
   };
  memset(DeltaPFTT, 0, DeltaPFTTSize*sizeof(double));

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
  //TODO insert here a check for any nested objects and automatically 
  //     disable symmetry if present
    
#ifdef USE_OPENMP
  Log("EMT OpenMP multithreading (%i threads)",NT);
#pragma omp parallel for schedule(dynamic,1), num_threads(NT)
#endif
  for(int neaTot=0; neaTot<TotalEdges; neaTot++)
   for(int nebTot=(UseSymmetry ? neaTot : 0); nebTot<TotalEdges; nebTot++)
    { 
      if (nebTot==(UseSymmetry ? neaTot : 0)) 
       LogPercent(neaTot, TotalEdges, 10);

      int nsa, nea, KNIndexA;
      RWGSurface *SA = G->ResolveEdge(neaTot, &nsa, &nea, &KNIndexA);
      int RegionIndex = SA->RegionIndices[Interior ? 1 : 0];
      if (RegionIndex==-1) continue; // no interior PFT for PEC bodies

      int nsb, neb, KNIndexB;
      RWGSurface *SB = G->ResolveEdge(nebTot, &nsb, &neb, &KNIndexB);
   
      double abSign=0.0, baSign=0.0;
      if (nsa==nsb)
       abSign = Interior ? -1.0 : +1.0;
      else if (SA->RegionIndices[0] == SB->RegionIndices[0]) // A, B live in same region
       abSign = baSign = (Interior ? 0.0 : 1.0);
      else if (SA->RegionIndices[0] == SB->RegionIndices[1]) // A contained in B
       { abSign = Interior ? 0.0 : -1.0;
         baSign = Interior ? 1.0 : 0.0;
       }
      else if (SA->RegionIndices[1] == SB->RegionIndices[0]) // B contained in A
       { abSign = Interior ? 1.0 : 0.0;
         baSign = Interior ? 0.0 : -1.0;
       };

      if ( abSign==0.0 && baSign==0.0 )
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
 
      double PFTIs[NUMPFTIS];
      double *QKK    = PFTIs + 0*NUMPFTT;
      double *QNN    = PFTIs + 1*NUMPFTT;
      double *QKNmNK = PFTIs + 2*NUMPFTT;
      GetPFTIntegrals_BFBF(G, nsa, nea, nsb, neb,
                           Omega, k, EpsR, MuR, EMTPFTIMethod, PFTIs);

      double dPFTT[NUMPFTT];
      dPFTT[PFT_PABS ] = 0.0;
      dPFTT[PFT_PSCAT] = 0.5*(  real(u0KK)*QKK[PFT_PABS]
                               +real(e0NN)*QNN[PFT_PABS]
                               +imag(KNmNK)*QKNmNK[PFT_PABS]
                             );

      for(int Mu=0; Mu<NUMFTT; Mu++)
       dPFTT[PFT_XFORCE + Mu]
        = 0.5*(  imag(u0KK)*QKK[PFT_XFORCE + Mu]
                +imag(e0NN)*QNN[PFT_XFORCE + Mu]
                +real(KNmNK)*QKNmNK[PFT_XFORCE +Mu]
              );

      int nt=0;
#ifdef USE_OPENMP
      nt=omp_get_thread_num();
#endif

      int OffsetA = nt*NS2NQ + nsa*NS*NQ + nsb*NQ;
      int OffsetB = nt*NS2NQ + nsb*NS*NQ + nsa*NQ;

      if (!UseSymmetry)
       VecPlusEquals(DeltaPFTT + OffsetA, abSign, dPFTT, NUMPFTT);
      else
       { 
         if (neaTot==nebTot)
          DeltaPFTT[ OffsetA + PFT_PSCAT ] += abSign*dPFTT[PFT_PSCAT];
         else if (nsa==nsb) // nebTot > neaTot but both on same surface
          { 
            DeltaPFTT[ OffsetA + PFT_PSCAT ] += 2.0*abSign*dPFTT[PFT_PSCAT];
            VecPlusEquals(DeltaPFTT + OffsetA + PFT_XFORCE, 2.0*abSign, dPFTT+PFT_XFORCE, NUMPFTT-PFT_XFORCE);
          }
         else // nebTot > neaTot and on different objects
          { 
            DeltaPFTT[ OffsetA + PFT_PSCAT ] += abSign*dPFTT[PFT_PSCAT];
            DeltaPFTT[ OffsetB + PFT_PSCAT ] += baSign*dPFTT[PFT_PSCAT];
            VecPlusEquals(DeltaPFTT + OffsetA + PFT_XFORCE,      abSign, dPFTT+PFT_XFORCE, NUMPFTT-PFT_XFORCE);
            VecPlusEquals(DeltaPFTT + OffsetB + PFT_XFORCE, -1.0*baSign, dPFTT+PFT_XFORCE, NUMPFTT-PFT_XFORCE);
          };
       };
 
    }; // end of multithreaded loop
  
  /*--------------------------------------------------------------*/
  /*- accumulate contributions of all threads                     */
  /*--------------------------------------------------------------*/
  for(int ns=0; ns<NS; ns++)
   PFTTBySurface[ns]->Zero();
  for(int nsa=0; nsa<NS; nsa++)
   for(int nsb=0; nsb<NS; nsb++)
    for(int nq=0; nq<NQ; nq++)
     for(int nt=0; nt<NT; nt++)
      PFTTBySurface[nsb]->AddEntry(nsa, nq, DeltaPFTT[ nt*NS2NQ + nsa*NS*NQ + nsb*NQ + nq ]);

  /***************************************************************/
  /* get incident-field contributions ****************************/
  /***************************************************************/
  IncidentPFTT->Zero();
  if (IF)
   AddIFContributionsToEMTPFTT(G, KNVector, IF, Omega, IncidentPFTT, Interior);
   
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  PFTMatrix->Zero();
  for(int nsa=0; nsa<NS; nsa++)
   for(int nq=0; nq<NUMPFT; nq++)
    { 
      for(int nsb=0; nsb<NS; nsb++)
       PFTMatrix->AddEntry(nsa,nq,PFTTBySurface[nsb]->GetEntry(nsa,nq));

      PFTMatrix->AddEntry(nsa,nq,IncidentPFTT->GetEntry(nsa,nq));

      if (PFT_XTORQUE<=nq && nq<=PFT_ZTORQUE)
       { for(int nsb=0; nsb<NS; nsb++)
          PFTMatrix->AddEntry(nsa, nq, PFTTBySurface[nsb]->GetEntry(nsa,nq+3));
         PFTMatrix->AddEntry(nsa,nq,IncidentPFTT->GetEntry(nsa,nq+3));
       };
    };
  
  /***************************************************************/
  /* If we are using the interior expansion, then the quantity   */
  /* computed as the scattered power above is actually minus the */
  /* absorbed power, and we want to zero out the scattered power */
  /* since the interior expansion doesn't compute that.          */
  /***************************************************************/
  if (Interior)
   for(int ns=0; ns<NS; ns++)
    { double PScat = PFTMatrix->GetEntryD(ns, PFT_PSCAT);
      PFTMatrix->SetEntry(ns, PFT_PABS,  -PScat);
      PFTMatrix->SetEntry(ns, PFT_PSCAT, 0.0);
    };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  char *ss=getenv("SCUFF_ITEMIZE_PFT");
  if (ss && ss[0]=='1')
   { 
     static bool WrotePreamble=false;
     for(int nsa=0; nsa<NS; nsa++)
      { FILE *f=vfopen("%s.%s.%cEMTPFT%i","a",
                        GetFileBase(G->GeoFileName),
                        G->Surfaces[nsa]->Label,
                        (Interior ? 'I' : 'E'), EMTPFTIMethod);
        if (!f) continue;
        if (!WrotePreamble)
         { fprintf(f,"# %s EMTPFT contributions to surface %s (integral method %i)\n",
                     Interior ? "Interior" : "Exterior", G->Surfaces[nsa]->Label, EMTPFTIMethod);
           fprintf(f,"# columns: \n");
           fprintf(f,"# 1 frequency \n");
           fprintf(f,"# 2 destination surface label \n");
           fprintf(f,"# 03-10 PAbs, PScat, Fxyz, Txyz (total)\n");
           fprintf(f,"# 11-21 PAbs, PScat, Fxyz, T1xyz, T2xyz (extinction)\n");
           int nc=22;
           for(int nsb=0; nsb<NS; nsb++, nc+=NUMPFTT)
            fprintf(f,"# %i-%i PAbs, PScat, Fxyz, Txyz (surface %s)\n",nc,nc+NUMPFTT-1,G->Surfaces[nsb]->Label);
         };
        fprintf(f,"%e %s ",real(Omega),G->Surfaces[nsa]->Label);
        for(int nq=0; nq<NUMPFT; nq++)
         fprintf(f,"%e ",PFTMatrix->GetEntryD(nsa,nq));
        for(int nq=0; nq<NUMPFTT; nq++)
         fprintf(f,"%e ",IncidentPFTT->GetEntryD(nsa,nq));
        for(int nsb=0; nsb<NS; nsb++)
         for(int nq=0; nq<NUMPFTT; nq++)
          fprintf(f,"%e ",PFTTBySurface[nsb]->GetEntryD(nsa,nq));
        fprintf(f,"\n");
        fclose(f);
      };
     WrotePreamble=true;
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
#if 0
  if (EMTPFTIMethod==0)
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
