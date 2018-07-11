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

#define PFT_XTORQUE1 5
#define PFT_YTORQUE1 6
#define PFT_ZTORQUE1 7
#define PFT_XTORQUE2 8
#define PFT_YTORQUE2 9
#define PFT_ZTORQUE2 10
#define NUMFT    6
#define NUMFTT   9
#define NUMPFTT  11

#define NUMPFTQ 17
#define NUMPFTIS (3*NUMPFTQ)

/***************************************************************/
/***************************************************************/
/***************************************************************/
namespace scuff {

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
   double *TorqueCenterA, *TorqueCenterB;
   int EMTPFTIMethod;
   RWGGeometry *G;
   int nsb, neb;
   bool SameSurface;

 } PFTIData;

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ScatteredPFTIntegrand1(double xA[3], double bA[3], double DivbA,
                            double xB[3], double bB[3], double DivbB,
                            void *UserData, double Weight, double *I)
{
  PFTIData *Data    = (PFTIData *)UserData;
  double Omega      = real(Data->Omega);
  double k          = real(Data->k);
  double EpsRel     = real(Data->EpsRel);
  double MuRel      = real(Data->MuRel);
  double *XTorqueA  = Data->TorqueCenterA;
  double *XTorqueB  = Data->TorqueCenterB;
  bool SameSurface  = Data->SameSurface;

  /***************************************************************/
  /* kernel factors **********************************************/
  /***************************************************************/
  double R[3];
  VecSub(xA, xB, R);
  double r2=R[0]*R[0] + R[1]*R[1] + R[2]*R[2];
  double r = sqrt(r2);
  double kr = k*r, kr2=kr*kr;
  double k2 = k*k;

  cdouble Phi, Psi, Zeta;
  if (SameSurface)
   { double ImPhi, ImPsi, ImZeta;
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
     Phi=II*ImPhi;
     Psi=II*ImPsi;
     Zeta=II*ImZeta;
   }
  else
   { cdouble ikr=II*kr, ikr2=ikr*ikr;
     Phi  = exp(ikr)/(4.0*M_PI*r);
     Psi  = Phi*(ikr-1.0)/r2;
     Zeta = Phi * (ikr2 - 3.0*ikr + 3.0) / (r2*r2);
   };

  /***************************************************************/
  /* polynomial factors ******************************************/
  /***************************************************************/
  double bdb  = VecDot(bA, bB);
  double DbDb = DivbA*DivbB;
  double bAdR = VecDot(bA, R);
  double bBdR = VecDot(bB, R);

  double bxb[3];    VecCross(bA, bB, bxb);
  //double bAxR[3];   VecCross(bA, R, bAxR);

  double XTA[3];     VecSub(xA, XTorqueA, XTA);
  double XTAxR[3];   VecCross(XTA, R,   XTAxR);
  double XTAxbxb[3]; VecCross(XTA, bxb, XTAxbxb);

  double XTB[3];     VecSub(xB, XTorqueB, XTB);
  double XTBxR[3];   VecCross(XTB, R,   XTBxR);
  double XTBxbxb[3]; VecCross(XTB, bxb, XTBxbxb);

  double PEFIE = bdb - DbDb/(k*k);
  double PMFIE = VecDot(bxb, R);

  cdouble *zI     = (cdouble *)I;
  cdouble *QKK    = zI+0*NUMPFTQ;
  cdouble *QNN    = zI+1*NUMPFTQ;
  cdouble *QKNmNK = zI+2*NUMPFTQ;

  QKK[PFT_PSCAT]    += Weight*PEFIE*Omega*MuRel*Phi;
  QNN[PFT_PSCAT]    += Weight*PEFIE*Omega*EpsRel*Phi;
  QKNmNK[PFT_PSCAT] += Weight*PMFIE*Psi;
 
  double FTPreFac=TENTHIRDS*Weight;
  for(int Mu=0; Mu<3; Mu++)
   { 
     QKK[PFT_XFORCE + Mu]    += FTPreFac*PEFIE*MuRel*R[Mu]*Psi;
     QNN[PFT_XFORCE + Mu]    += FTPreFac*PEFIE*EpsRel*R[Mu]*Psi;
     QKNmNK[PFT_XFORCE + Mu] -= FTPreFac*(bxb[Mu]*Psi + PMFIE*R[Mu]*Zeta)/Omega;

  //   QKK[PFT_XTORQUE1+ Mu ]   += FTPreFac*MuRel*( bxb[Mu]*(Phi + Psi/k2)  + bAxR[Mu]*bBdR*Zeta/k2);
 //    QNN[PFT_XTORQUE1+ Mu ]   += FTPreFac*EpsRel*( bxb[Mu]*(Phi + Psi/k2) + bAxR[Mu]*bBdR*Zeta/k2);
     QKK[PFT_XTORQUE1+ Mu ]   += FTPreFac*MuRel*(bxb[Mu]*Phi);
     QNN[PFT_XTORQUE1+ Mu ]   += FTPreFac*EpsRel*(bxb[Mu]*Phi);
     QKNmNK[PFT_XTORQUE1+ Mu] -= FTPreFac*(bAdR*bB[Mu] - R[Mu]*bdb)*Psi/Omega;

     QKK[PFT_XTORQUE2 + Mu]    += FTPreFac*PEFIE*MuRel*XTAxR[Mu]*Psi;
     QNN[PFT_XTORQUE2 + Mu]    += FTPreFac*PEFIE*EpsRel*XTAxR[Mu]*Psi;
     QKNmNK[PFT_XTORQUE2 + Mu] -= FTPreFac*(XTAxbxb[Mu]*Psi + PMFIE*XTAxR[Mu]*Zeta)/Omega;

     QKK[PFT_XTORQUE1    + Mu + 6 ] -= FTPreFac*MuRel*(bxb[Mu]*Phi);
     QNN[PFT_XTORQUE1    + Mu + 6 ] -= FTPreFac*EpsRel*(bxb[Mu]*Phi);
     QKNmNK[PFT_XTORQUE1 + Mu + 6 ] += FTPreFac*(bBdR*bA[Mu] - R[Mu]*bdb)*Psi/Omega;

     QKK[PFT_XTORQUE2    + Mu + 6 ] -= FTPreFac*PEFIE*MuRel*XTBxR[Mu]*Psi;
     QNN[PFT_XTORQUE2    + Mu + 6 ] -= FTPreFac*PEFIE*EpsRel*XTBxR[Mu]*Psi;
     QKNmNK[PFT_XTORQUE2 + Mu + 6 ] += FTPreFac*(XTBxbxb[Mu]*Psi + PMFIE*XTBxR[Mu]*Zeta)/Omega;
   }
}

  
/***************************************************************/
/***************************************************************/
/***************************************************************/
#if 0
void ScatteredPFTIntegrand2(double *xA, PCData *PCD,
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
  double *XTorqueA  = Data->TorqueCenterA;
  double *XTorqueB  = Data->TorqueCenterB;

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

  cdouble *zI = (cdouble *)I;
  cdouble *QKK    = zI+0*NUMPFTT;
  cdouble *QNN    = zI+1*NUMPFTT;
  cdouble *QKNmNK = zI+2*NUMPFTT;

  QKK[PFT_PSCAT]    =  Omega*MuRel*RCVecDot(bA, e);
  QNN[PFT_PSCAT]    =  Omega*EpsRel*RCVecDot(bA, e);
  QKNmNK[PFT_PSCAT] = -RCVecDot(bA, h);
  if (EMTPFTIMethod==SCUFF_EMTPFTI_EHVALUES1)
   {
     for(int Mu=0; Mu<3; Mu++)
      { QKK[PFT_XFORCE+Mu]    =  TENTHIRDS*MuRel*RCVecDot(bA, de[Mu]);
        QNN[PFT_XFORCE+Mu]    =  TENTHIRDS*EpsRel*RCVecDot(bA, de[Mu]);
        QKNmNK[PFT_XFORCE+Mu] =  TENTHIRDS*RCVecDot(bA, dh[Mu])/Omega;

        QKK[PFT_XTORQUE+Mu]    =  TENTHIRDS*MuRel*bxe[Mu];
        QNN[PFT_XTORQUE+Mu]    =  TENTHIRDS*EpsRel*bxe[Mu];
        QKNmNK[PFT_XTORQUE+Mu] =  TENTHIRDS*bxh[Mu]/Omega;

        int MP1 = (Mu+1)%3, MP2 = (Mu+2)%3;
        QKK[PFT_XTORQUE2+Mu]    =  XT[MP1]*QKK[2+MP2] - XT[MP2]*QKK[2+MP1];
        QNN[PFT_XTORQUE2+Mu]    =  XT[MP1]*QNN[2+MP2] - XT[MP2]*QNN[2+MP1];
        QKNmNK[PFT_XTORQUE2+Mu] =  XT[MP1]*QKNmNK[2+MP2] - XT[MP2]*QKNmNK[2+MP1];
      };
   }
  else
   {
     for(int Mu=0; Mu<3; Mu++)
      { QKK[PFT_XFORCE+Mu]    = TENTHIRDS*MuRel*(-DivbA*e[Mu] + bxh[Mu]);
        QNN[PFT_XFORCE+Mu]    = TENTHIRDS*EpsRel*(-DivbA*e[Mu] + bxh[Mu]);
        QKNmNK[PFT_XFORCE+Mu] = TENTHIRDS*DivbA*h[Mu] - k*k*bxe[Mu];
      };

     for(int Mu=0; Mu<3; Mu++)
      { 
        int MP1=(Mu+1)%3, MP2=(Mu+2)%3;

        QKK[PFT_XTORQUE+Mu]    = XT[MP1]*QKK[1+MP2]    - XT[MP2]*QKK[1+MP1];
        QNN[PFT_XTORQUE+Mu]    = XT[MP1]*QNN[1+MP2]    - XT[MP2]*QNN[1+MP1];
        QKNmNK[PFT_XTORQUE+Mu] = XT[MP1]*QKNmNK[1+MP2] - XT[MP2]*QKNmNK[1+MP1];
      };

   };
}
#endif

/***************************************************************/
/* compute PFT integrals between pairs of RWG basis functions  */
/***************************************************************/
void GetScatteredPFTIntegrals(RWGGeometry *G,
                              int nsa, int nea, int nsb, int neb,
                              cdouble Omega, cdouble k, 
                              cdouble EpsR, cdouble MuR,
                              int EMTPFTIMethod, cdouble PFTIs[NUMPFTIS])
{
  RWGSurface *SA=G->Surfaces[nsa];
  RWGSurface *SB=G->Surfaces[nsb];
  double rRel;
  int ncv=AssessBFPair(SA, nea, SB, neb, &rRel);

  if (EMTPFTIMethod==SCUFF_EMTPFTI_EHDERIVATIVES1)
   {
     PFTIData MyData, *Data=&MyData;
     Data->Omega         = Omega;
     Data->k             = k;
     Data->EpsRel        = EpsR;
     Data->MuRel         = MuR;
     Data->TorqueCenterA = G->Surfaces[nsa]->Origin;
     Data->TorqueCenterB = G->Surfaces[nsb]->Origin;
     Data->SameSurface   = (nsa==nsb);

     int IDim  = 2*NUMPFTIS;
     int Order = (ncv>0) ? 9 : 4;
 
     GetBFBFCubature2(G, nsa, nea, nsb, neb,
                      ScatteredPFTIntegrand1, (void *)Data, IDim,
                      Order, (double *)PFTIs);
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

     cdouble *QKK    = PFTIs + 0*NUMPFTQ;
     cdouble *QNN    = PFTIs + 1*NUMPFTQ;
     cdouble *QKNmNK = PFTIs + 2*NUMPFTQ;
     QKK[PFT_PSCAT]    = MuR*Omega*GabArray[0][GCME_GC];
     QNN[PFT_PSCAT]    = EpsR*Omega*GabArray[0][GCME_GC];
     QKNmNK[PFT_PSCAT] = ikCabArray[0][GCME_GC];
     double FTPreFac=TENTHIRDS;
     for(int Mu=0; Mu<NUMFT; Mu++)
      { QKK[PFT_XFORCE+Mu]    = FTPreFac*MuR*GabArray[0][GCME_FX + Mu];
        QNN[PFT_XFORCE+Mu]    = FTPreFac*EpsR*GabArray[0][GCME_FX + Mu];
        QKNmNK[PFT_XFORCE+Mu] = FTPreFac*ikCabArray[0][GCME_FX + Mu]/Omega;
      }
     for(int nq=PFT_XTORQUE1+6; nq<PFT_ZTORQUE2+6; nq++)
      QKK[nq]=QNN[nq]=QKNmNK[nq]=0.0;
   }
  else 
   ErrExit("PFT integration method no longer supported");
/*
   { 
     PFTIData MyData, *Data=&MyData;
     Data->Omega        = Omega;
     Data->k            = k;
     Data->EpsRel       = EpsR;
     Data->MuRel        = MuR;
     Data->TorqueCenter = G->Surfaces[nsa]->Origin;
     Data->EMTPFTIMethod = EMTPFTIMethod;

     int IDim     = 2*NUMPFTIS;
     int MaxEvals = (ncv>0) ? 78 : 21;

     Data->G   = G;
     Data->nsb = nsb;
     Data->neb = neb;
     GetBFCubature(G, nsa, nea, ScatteredPFTIntegrand2,
                   (void *)Data, IDim, MaxEvals,
                   0.0, 0.0, 0.0, 0, (double *)PFTIs);
   };
*/

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ExtinctionPFTIntegrand(double *x, PCData *PCD,
                            void *UserData, double *I)
{
  double b[3];
  b[0] = real(PCD->K[0]);
  b[1] = real(PCD->K[1]);
  b[2] = real(PCD->K[2]);

  PFTIData *Data=(PFTIData *)UserData;
  IncField *IF    = Data->IF;
  double *XTorque = Data->TorqueCenterA;
  double XT[3];
  VecSub(x, XTorque, XT);

  // get fields and derivatives at eval point
  // dE[Mu][Nu] = \partial_\mu E_\nu 
  cdouble EH[6], dEH[3][6];
  IF->GetFields(x, EH);
  IF->GetFieldGradients(x, dEH);

  cdouble *E, *H, *dE[3], *dH[3];
  E     =  EH    + 0*3;   H     =  EH    + 1*3;
  dE[0] = dEH[0] + 0*3;   dH[0] = dEH[0] + 1*3;
  dE[1] = dEH[1] + 0*3;   dH[1] = dEH[1] + 1*3;
  dE[2] = dEH[2] + 0*3;   dH[2] = dEH[2] + 1*3;
  
  // fill in the integrand vector 
  cdouble *zI = (cdouble *)I;
  memset(zI, 0, 2*NUMPFTT*sizeof(cdouble));
  cdouble *QK = zI + 0*NUMPFTT;
  cdouble *QN = zI + 1*NUMPFTT;

  QK[PFT_PABS] = RCVecDot(b, E);
  QN[PFT_PABS] = RCVecDot(b, H);

  for(int Mu=0; Mu<3; Mu++)
   { 
     QK[PFT_XFORCE+Mu] = RCVecDot(b, dE[Mu]);
     QN[PFT_XFORCE+Mu] = RCVecDot(b, dH[Mu]);
     
     int MP1=(Mu+1)%3, MP2=(Mu+2)%3;
     QK[PFT_XTORQUE1 + Mu] = b[MP1]*E[MP2] - b[MP2]*E[MP1];
     QN[PFT_XTORQUE1 + Mu] = b[MP1]*H[MP2] - b[MP2]*H[MP1];

     for(int Nu=0; Nu<3; Nu++)
      { 
        QK[PFT_XTORQUE2 + Mu] += b[Nu]*(XT[MP1]*dE[MP2][Nu]-XT[MP2]*dE[MP1][Nu]);
        QN[PFT_XTORQUE2 + Mu] += b[Nu]*(XT[MP1]*dH[MP2][Nu]-XT[MP2]*dH[MP1][Nu]);
      };

   };

}

/***************************************************************/
/* compute PFT integrals between an RWG basis function and an  */
/* external field                                              */
/***************************************************************/
void GetExtinctionPFTIntegrals(RWGGeometry *G, int ns, int ne,
                               IncField *IF, cdouble Omega,
                               cdouble QKQN[2*NUMPFTT])
{
  PFTIData MyData, *Data = &MyData;
  Data->Omega            = Omega;
  Data->IF               = IF;
  Data->TorqueCenterA    = G->Surfaces[ns]->Origin;

  int fDim=4*NUMPFTT;
  int MaxEvals=21;
  GetBFCubature(G, ns, ne, ExtinctionPFTIntegrand, (void *)Data,
                fDim, MaxEvals, 0, 0, 0, 0, (double *)QKQN);
}

/***************************************************************/
/* add incident-field contributions (i.e. extinction) to EMTPFT*/
/***************************************************************/
void GetExtinctionPFTT(RWGGeometry *G, HVector *KN,
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
  int NTNSNQ=NT*NS*NQ;

  static int DeltaPFTTSize=0;
  static double *DeltaPFTT=0;
  if (DeltaPFTTSize < NTNSNQ)
   { DeltaPFTTSize=NTNSNQ;
     DeltaPFTT=(double *)reallocEC(DeltaPFTT, DeltaPFTTSize*sizeof(double));
   };
  memset(DeltaPFTT, 0, NTNSNQ*sizeof(double));

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

     cdouble KStar, NStar;
     G->GetKNCoefficients(KN, ns, ne, &KStar, &NStar);
     KStar=conj(KStar);
     NStar=conj(NStar);

     cdouble QKQN[2*NUMPFTT];
     GetExtinctionPFTIntegrals(G, ns, ne, IF, Omega, QKQN);
     cdouble *QK=QKQN + 0*NUMPFTT;
     cdouble *QN=QKQN + 1*NUMPFTT;

     int nt=0;
#ifdef USE_OPENMP
     nt = omp_get_thread_num();
#endif
     double *dPFTT = DeltaPFTT + (nt*NS + ns)*NUMPFTT;
     dPFTT[PFT_PABS] += real(KStar*QK[PFT_PABS] + NStar*QN[PFT_PABS]);
     for(int nq=PFT_XFORCE; nq<NUMPFTT; nq++)
      dPFTT[nq] += imag( KStar*QK[nq] + NStar*QN[nq] );
   };

   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   double PFactor = 0.5;
   double FTFactor = 0.5*TENTHIRDS/real(Omega);
   PFTTMatrix->Zero();
   for(int nt=0; nt<NT; nt++)
    for(int ns=0; ns<NS; ns++)
     { 
       double *dPFTT = DeltaPFTT + (nt*NS + ns)*NUMPFTT;
       PFTTMatrix->AddEntry(ns, PFT_PABS, PFactor*dPFTT[PFT_PABS]);
       for(int nq=PFT_XFORCE; nq<NUMPFTT; nq++)
        PFTTMatrix->AddEntry(ns, nq, FTFactor*dPFTT[nq]);
     };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *GetEMTPFTMatrix(RWGGeometry *G, cdouble Omega, IncField *IF,
                         HVector *KNVector, HMatrix *DRMatrix,
                         HMatrix *PFTMatrix, bool Interior,
                         int EMTPFTIMethod, bool Itemize)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  char *sss=getenv("SCUFF_EMTPFTI_METHOD");
  if (sss)
   { sscanf(sss,"%i",&EMTPFTIMethod);
     Log("Using EMTPFTI method %i.",EMTPFTIMethod);
   };

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
  /* ScatteredPFTT[ns] = contributions of surface #ns to         */
  /*                     scattered PFTT                          */
  /***************************************************************/
  static int NSSave=0;
  static HMatrix **ScatteredPFTT=0, *ExtinctionPFTT=0;
  if (NSSave!=NS)
   { if (ScatteredPFTT)
      { for(int ns=0; ns<NSSave; ns++)
         if (ScatteredPFTT[ns]) 
          delete ScatteredPFTT[ns];
        free(ScatteredPFTT);
       if (ExtinctionPFTT)
        delete ExtinctionPFTT;
      };
     NSSave=NS;
     ScatteredPFTT=(HMatrix **)mallocEC(NS*sizeof(HMatrix));
     for(int ns=0; ns<NS; ns++)
      ScatteredPFTT[ns]=new HMatrix(NS, NUMPFTT);
     ExtinctionPFTT=new HMatrix(NS, NUMPFTT);
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
  /*- loop over all edge pairs to get scattered PFT contributions */
  /*- of all surfaces                                             */
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
  memset(DeltaPFTT, 0, NTNS2NQ*sizeof(double));

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

  int nsaOnly=-1; CheckEnv("SCUFF_EMTPFT_NSAONLY",&nsaOnly);
  int nsbOnly=-1; CheckEnv("SCUFF_EMTPFT_NSBONLY",&nsbOnly);
  //char *EMTLog=0; CheckEnv("SCUFF_EMTPFT_LOGFILE",&EMTLog);
  //FILE *fLog = EMTLog ? fopen(EMTLog,"w") : 0;
    
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

      if ( (nsaOnly!=-1 && nsa!=nsaOnly) || (nsbOnly!=-1 && nsb!=nsbOnly) )
       continue;

      double Sign=0.0;
      if (nsa==nsb)
       Sign = Interior ? -1.0 : +1.0;
      else if (SA->RegionIndices[0] == SB->RegionIndices[0]) // A, B live in same region
       Sign = (Interior ? 0.0 : 1.0);
      else if (SA->RegionIndices[0] == SB->RegionIndices[1]) // A contained in B
       Sign = Interior ? 0.0 : -1.0;
      else if (SA->RegionIndices[1] == SB->RegionIndices[0]) // B contained in A
       Sign = Interior ? 1.0 : 0.0;

      if ( Sign==0.0 ) // B does not contribute to PFT on A
       continue;

      cdouble EpsR, MuR;
      G->RegionMPs[RegionIndex]->GetEpsMu(Omega, &EpsR, &MuR);
      cdouble k = Omega*sqrt(EpsR*MuR);

      cdouble PFTIs[NUMPFTIS];
      cdouble *QKK    = PFTIs + 0*NUMPFTQ;
      cdouble *QNN    = PFTIs + 1*NUMPFTQ;
      cdouble *QKNmNK = PFTIs + 2*NUMPFTQ;
      GetScatteredPFTIntegrals(G, nsa, nea, nsb, neb,
                               Omega, k, EpsR, MuR, EMTPFTIMethod, PFTIs);

      cdouble KNBab[4];
      GetKNBilinears(KNVector, DRMatrix,
                     SA->IsPEC, KNIndexA, SB->IsPEC, KNIndexB, KNBab);
      cdouble u0KKab  = Sign * KNBab[0] * ZVAC;
      cdouble KNmNKab = Sign * (KNBab[1] - KNBab[2]);
      cdouble e0NNab  = Sign * KNBab[3] / ZVAC;

      double dPFTT[NUMPFTT];
      dPFTT[PFT_PABS] = 0.0;

      if (nsa==nsb) // self contributions; use symmetry reduction
       { 
         dPFTT[PFT_PSCAT] =  real(u0KKab)*imag(QKK[PFT_PSCAT])
                            +real(e0NNab)*imag(QNN[PFT_PSCAT])
                            +imag(KNmNKab)*imag(QKNmNK[PFT_PSCAT]);

         for(int nq=PFT_XFORCE; nq<NUMPFTT; nq++)
          dPFTT[nq] =  imag(u0KKab)*imag(QKK[nq])
                      +imag(e0NNab)*imag(QNN[nq])
                      -real(KNmNKab)*imag(QKNmNK[nq]);
        }
      else
       { 
         dPFTT[PFT_PSCAT] = -1.0*real(  u0KKab*II*QKK[PFT_PSCAT]
                                       +e0NNab*II*QNN[PFT_PSCAT]
                                      +KNmNKab*QKNmNK[PFT_PSCAT]
                                     );

         for(int nq=PFT_XFORCE; nq<NUMPFTT; nq++)
          dPFTT[nq] = -1.0*imag(   u0KKab*II*QKK[nq]
                                  +e0NNab*II*QNN[nq]
                                 +KNmNKab*QKNmNK[nq]
                               );
       };

      int nt=0;
#ifdef USE_OPENMP
      nt=omp_get_thread_num();
#endif
      int Offset = nt*NS2NQ + nsa*NS*NQ + nsb*NQ;
      VecPlusEquals(DeltaPFTT + Offset, 0.5*Sign, dPFTT, NUMPFTT);

      if (!UseSymmetry || neaTot==nebTot)
       continue;

      for(int nq=PFT_XFORCE; nq<=PFT_ZFORCE; nq++)
       { QKK[nq]*=-1.0;
         QNN[nq]*=-1.0;
         QKNmNK[nq]*=-1.0;
       };
      for(int nq=PFT_XTORQUE1; nq<=PFT_ZTORQUE2; nq++)
       { QKK[nq]    = QKK[nq+6];
         QNN[nq]    = QNN[nq+6];
         QKNmNK[nq] = QKNmNK[nq+6];
       };

      cdouble KNBba[4];
      GetKNBilinears(KNVector, DRMatrix,
                     SB->IsPEC, KNIndexB, SA->IsPEC, KNIndexA, KNBba);
      cdouble u0KKba  = Sign * KNBba[0] * ZVAC;
      cdouble KNmNKba = Sign * (KNBba[1] - KNBba[2]);
      cdouble e0NNba  = Sign * KNBba[3] / ZVAC;

      dPFTT[PFT_PABS] = 0.0;

      if (nsa==nsb) // self contributions; use symmetry reduction
       { 
         dPFTT[PFT_PSCAT] =  real(u0KKba)*imag(QKK[PFT_PSCAT])
                            +real(e0NNba)*imag(QNN[PFT_PSCAT])
                            +imag(KNmNKba)*imag(QKNmNK[PFT_PSCAT]);

         for(int nq=PFT_XFORCE; nq<NUMPFTT; nq++)
          dPFTT[nq] =  imag(u0KKba)*imag(QKK[nq])
                      +imag(e0NNba)*imag(QNN[nq])
                      -real(KNmNKba)*imag(QKNmNK[nq]);
        }
      else
       { 
         dPFTT[PFT_PSCAT] = -1.0*real(  u0KKba*II*QKK[PFT_PSCAT]
                                       +e0NNba*II*QNN[PFT_PSCAT]
                                      +KNmNKba*QKNmNK[PFT_PSCAT]
                                     );

         for(int nq=PFT_XFORCE; nq<NUMPFTT; nq++)
          dPFTT[nq] = -1.0*imag(   u0KKba*II*QKK[nq]
                                  +e0NNba*II*QNN[nq]
                                 +KNmNKba*QKNmNK[nq]
                               );
       };

      Offset = nt*NS2NQ + nsb*NS*NQ + nsa*NQ;
      VecPlusEquals(DeltaPFTT + Offset, 0.5*Sign, dPFTT, NUMPFTT);
 
    }; // end of multithreaded loop
  
  /*--------------------------------------------------------------*/
  /*- accumulate contributions of all threads                     */
  /*--------------------------------------------------------------*/
  for(int ns=0; ns<NS; ns++)
   ScatteredPFTT[ns]->Zero();
  for(int nsa=0; nsa<NS; nsa++)
   for(int nsb=0; nsb<NS; nsb++)
    for(int nq=0; nq<NQ; nq++)
     for(int nt=0; nt<NT; nt++)
      ScatteredPFTT[nsb]->AddEntry(nsa, nq, DeltaPFTT[ nt*NS2NQ + nsa*NS*NQ + nsb*NQ + nq ]);

  /***************************************************************/
  /* get incident-field contributions ****************************/
  /***************************************************************/
  if (IF)
   GetExtinctionPFTT(G, KNVector, IF, Omega, ExtinctionPFTT, Interior);
  else
   ExtinctionPFTT->Zero();
   
  /***************************************************************/
  /* sum scattered contributions from all surfaces plus          */
  /* extinction contributions to get total PFT.                  */
  /* note that the formula for total PFT is                      */
  /* Q^{full} = Q^{extinction} - Q^{scattering}                  */
  /* so the scattering contributions enter with a minus sign     */
  /* (except for the scattered power).                           */
  /***************************************************************/
  PFTMatrix->Zero();
  for(int nsa=0; nsa<NS; nsa++)
   for(int nq=0; nq<NUMPFT; nq++)
    { 
      PFTMatrix->AddEntry(nsa,nq,ExtinctionPFTT->GetEntry(nsa,nq));

      for(int nsb=0; nsb<NS; nsb++)
       { 
         if (nq==PFT_PABS)
          PFTMatrix->AddEntry(nsa,nq,-1.0*ScatteredPFTT[nsb]->GetEntry(nsa,PFT_PSCAT));
         else if (nq==PFT_PSCAT)
          PFTMatrix->AddEntry(nsa,nq,+1.0*ScatteredPFTT[nsb]->GetEntry(nsa,PFT_PSCAT));
         else // force or torque 
          PFTMatrix->AddEntry(nsa,nq,-1.0*ScatteredPFTT[nsb]->GetEntry(nsa,nq));
       };

      if (PFT_XTORQUE<=nq && nq<=PFT_ZTORQUE)
       { PFTMatrix->AddEntry(nsa,nq,ExtinctionPFTT->GetEntry(nsa,nq+3));
         for(int nsb=0; nsb<NS; nsb++)
          PFTMatrix->AddEntry(nsa, nq, -1.0*ScatteredPFTT[nsb]->GetEntry(nsa,nq+3));
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
    PFTMatrix->SetEntry(ns, PFT_PSCAT, 0.0);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  char *ss=getenv("SCUFF_ITEMIZE_PFT");
  if (ss && ss[0]=='1')
   Itemize=true;
  if (Itemize)
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
           fprintf(f,"# 2 source (0=total, -1=extinction, +n=object #n scattered)\n");
           fprintf(f,"# 3,4      PAbs, PScat\n");
           fprintf(f,"# 5,6,7    Fx, Fy, Fz\n");
           fprintf(f,"# 8,9,10   Tx, Ty, Tz (term 1)\n");
           fprintf(f,"# 11,12,13 Tx, Ty, Tz (term 2)\n");
         };

        fprintf(f,"%e 00 ",real(Omega));
        for(int nq=0; nq<NUMPFT; nq++)
         fprintf(f,"%+e ",PFTMatrix->GetEntryD(nsa,nq));
        fprintf(f,"\n");

        if (IF)
         { fprintf(f,"%e -1 ",real(Omega));
           for(int nq=0; nq<NUMPFTT; nq++)
            fprintf(f,"%+e ",ExtinctionPFTT->GetEntryD(nsa,nq));
           fprintf(f,"\n");
         };

        for(int nsb=0; nsb<NS; nsb++)
         { fprintf(f,"%e %02i ",real(Omega),nsb+1);
            for(int nq=0; nq<NUMPFTT; nq++)
             fprintf(f,"%+e ",ScatteredPFTT[nsb]->GetEntryD(nsa,nq));
           fprintf(f,"\n");
         };
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

