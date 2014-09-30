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
 * EPPFT.cc    -- 'equivalence-principle power, force, and torque'
 *             -- calculation in scuff-EM
 *
 * homer reid  -- 9/2014
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libhmat.h>
#include <libhrutil.h>
#include <libTriInt.h>

#include "libscuff.h"
#include "libscuffInternals.h"
#include "PanelCubature.h"

#include "cmatheval.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#ifdef USE_PTHREAD
#  include <pthread.h>
#endif
#ifdef USE_OPENMP

#endif

#define II cdouble(0,1)

namespace scuff {

void CalcGC(double R[3], cdouble Omega,
            cdouble EpsR, cdouble MuR,
            cdouble GMuNu[3][3], cdouble CMuNu[3][3],
            cdouble GMuNuRho[3][3][3], cdouble CMuNuRho[3][3][3]);

void GetReducedFields_Nearby(RWGGeometry *G, int ns, int ne,
                             double X0[3],  cdouble k,
                             cdouble e[3], cdouble h[3],
                             cdouble de[3][3], cdouble dh[3][3]);

void GetReducedFields_Nearby(RWGGeometry *G, int ns, int ne,
                             double X0[3],  cdouble k,
                             cdouble e[3], cdouble h[3]);

/***************************************************************/
/* Fetch the particular matrix elements between RWG functions  */
/* that we need to compute the surface EPPFT.                  */
/*                                                             */
/* be[0]       = < b_\alpha | e_beta >                         */
/* bh[0]       = < b_\alpha | h_beta >                         */
/*                                                             */
/* divbe[0..2] = < \div b_\alpha | e_beta_i >  (i=0,1,2)       */
/* divbh[0..2] = < \div b_\alpha | h_beta_i >                  */
/*                                                             */
/* bxe[0..2]   = < b_\alpha \times e_beta >_i                  */
/* bxh[0..2]   = < b_\alpha \times h_beta >_i                  */
/*                                                             */
/***************************************************************/
void GetEPPFTMatrixElements(RWGGeometry *G,
                            int nsa, int nsb, int nea, int neb,
                            cdouble k,
                            cdouble be[1],      cdouble bh[1],
                            cdouble divbe[3],   cdouble divbh[3],
                            cdouble bxe[3],     cdouble bxh[3],
                            cdouble divbrxe[3], cdouble divbrxh[3],
                            cdouble rxbxe[3],   cdouble rxbxh[3],
                            int Order=4)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
#if 0
  int ncv=NumCommonBFVertices(G->Surfaces[nsa],nea,G->Surfaces[nsb],neb);
  if ( (ncv==0 || ForceDistant) && (!ForceNearby) )
   GetdGME_Far(G, nsa, nea, nsb, neb, k, Order, DX, GC, dG, dC);
  else
   GetdGME_Near(G, nsa, nea, nsb, neb, k, Order, DX, GC, dG, dC);
#endif

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGSurface *S=G->Surfaces[nsa];
  RWGEdge *E=S->Edges[nea];
  double *QP = S->Vertices + 3*E->iQP;
  double *V1 = S->Vertices + 3*E->iV1;
  double *V2 = S->Vertices + 3*E->iV2;
  double *QM = S->Vertices + 3*E->iQM;

  double X0[3]={0.0, 0.0, 0.0}; // torque center
  if (S->OTGT) S->OTGT->Apply(X0);
  if (S->GT) S->GT->Apply(X0);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double APlus[3], AMinus[3], B[3];
  VecSub(V1, QP, APlus);
  VecSub(V1, QM, AMinus);
  VecSub(V2, V1, B);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int NumPts;
  double *TCR=GetTCR(Order, &NumPts);
  be[0]=bh[0]=0.0;
  memset(divbe,   0, 3*sizeof(cdouble));
  memset(divbh,   0, 3*sizeof(cdouble));
  memset(bxe,     0, 3*sizeof(cdouble));
  memset(bxh,     0, 3*sizeof(cdouble));
  memset(divbrxe, 0, 3*sizeof(cdouble));
  memset(divbrxh, 0, 3*sizeof(cdouble));
  memset(rxbxe,   0, 3*sizeof(cdouble));
  memset(rxbxh,   0, 3*sizeof(cdouble));
  for(int np=0, ncp=0; np<NumPts; np++)
   { 
     /***************************************************************/
     /***************************************************************/
     /***************************************************************/
     double u=TCR[ncp++];
     double v=TCR[ncp++];
     double w=TCR[ncp++] * (E->Length);
     u+=v;

     double bPlus[3], XPlus[3], XPmX0[3];
     double bMinus[3], XMinus[3], XMmX0[3];
     double XPmX0dotb=0.0, XMmX0dotb=0.0;
     for(int Mu=0; Mu<3; Mu++)
      { 
        bPlus[Mu] = u*APlus[Mu] + v*B[Mu];
        XPlus[Mu] = bPlus[Mu] + QP[Mu];
        XPmX0[Mu] = XPlus[Mu] - X0[Mu];
        XPmX0dotb += XPmX0[Mu]*bPlus[Mu];

        bMinus[Mu] = u*AMinus[Mu] + v*B[Mu];
        XMinus[Mu] = bMinus[Mu] + QM[Mu];
        XMmX0[Mu]  = XMinus[Mu] - X0[Mu];
        XMmX0dotb += XMmX0[Mu]*bMinus[Mu];
      };

     /***************************************************************/
     /***************************************************************/
     /***************************************************************/
     cdouble ePlus[3], eMinus[3];
     cdouble hPlus[3], hMinus[3];
     GetReducedFields_Nearby(G, nsb, neb, XPlus, k, ePlus, hPlus);
     GetReducedFields_Nearby(G, nsb, neb, XMinus, k, eMinus, hMinus);

     cdouble XPmX0dote=0.0, XPmX0doth=0.0, XMmX0dote=0.0, XMmX0doth=0.0;
     for(int Mu=0; Mu<3; Mu++)
      { XPmX0dote += XPmX0[Mu]*ePlus[Mu];
        XPmX0doth += XPmX0[Mu]*hPlus[Mu];
        XMmX0dote += XMmX0[Mu]*eMinus[Mu];
        XMmX0doth += XMmX0[Mu]*hMinus[Mu];
      };
    
     for(int Mu=0; Mu<3; Mu++)
       {
         be[0] += w*(bPlus[Mu]*ePlus[Mu] - bMinus[Mu]*eMinus[Mu]); 
         bh[0] += w*(bPlus[Mu]*hPlus[Mu] - bMinus[Mu]*hMinus[Mu]);

         divbe[Mu] += 2.0*w*(ePlus[Mu] - eMinus[Mu]);
         divbh[Mu] += 2.0*w*(hPlus[Mu] - hMinus[Mu]);

         int MP1=(Mu+1)%3, MP2=(Mu+2)%3;
         bxe[Mu] += w*( (bPlus[MP1]*ePlus[MP2]   - bPlus[MP2]*ePlus[MP1])
                       -(bMinus[MP1]*eMinus[MP2] - bMinus[MP2]*eMinus[MP1])
                      );

         bxh[Mu] += w*( (bPlus[MP1]*hPlus[MP2]   - bPlus[MP2]*hPlus[MP1])
                       -(bMinus[MP1]*hMinus[MP2] - bMinus[MP2]*hMinus[MP1])
                      );

         divbrxe[Mu] += w*( (XPmX0[MP1]*ePlus[MP2]  - XPmX0[MP2]*ePlus[MP1])
                           -(XMmX0[MP1]*eMinus[MP2] - XMmX0[MP2]*eMinus[MP1])
                          );

         divbrxh[Mu] += w*( (XPmX0[MP1]*hPlus[MP2]  - XPmX0[MP2]*hPlus[MP1])
                           -(XMmX0[MP1]*hMinus[MP2] - XMmX0[MP2]*hMinus[MP1])
                          );

         // [Ax(BxC)]_mu = B_\mu (A\cdot C) - C_\mu (A\cdot B)
         rxbxe[Mu] += w*( (bPlus[Mu]*XPmX0dote  - ePlus[Mu]*XPmX0dotb)
                         -(bMinus[Mu]*XMmX0dote - eMinus[Mu]*XMmX0dotb)
                        );

         rxbxh[Mu] += w*( (bPlus[Mu]*XPmX0doth  - hPlus[Mu]*XPmX0dotb)
                         -(bMinus[Mu]*XMmX0doth - hMinus[Mu]*XMmX0dotb)
                        );
       };

   };

}

/***************************************************************/
/* Note: Either KNVector or SigmaMatrix should be non-null.    */
/***************************************************************/
void RWGGeometry::GetEPPFT(int ns, HVector *KNVector, HMatrix *SigmaMatrix,
                           cdouble Omega, double EPPFT[7])
{
  RWGSurface *S=Surfaces[ns];
  int NE = S->NumEdges;
  int Offset = BFIndexOffset[ns];
  int nrOut  = S->RegionIndices[0];
  int nrIn   = S->RegionIndices[1];
  if ( S->IsPEC || nrIn==-1 ) // EPPFT not defined for PEC bodies
   { memset(EPPFT, 0, 7*sizeof(double)); 
     return;
   };

  /*--------------------------------------------------------------*/
  /*- extract material parameters and define constant prefactors -*/
  /*--------------------------------------------------------------*/
  cdouble EpsIn, MuIn;
  RegionMPs[nrIn]->GetEpsMu(Omega, &EpsIn, &MuIn);
  cdouble k = Omega * sqrt(EpsIn*MuIn);
  cdouble ZRel = sqrt(MuIn/EpsIn);
  cdouble KZ  = k*ZVAC*ZRel;
  cdouble KOZ = k/(ZVAC*ZRel);

  cdouble PEE  = +II*KZ/2.0;
  cdouble PEM  = -1.0/2.0;
  cdouble PME  = +1.0/2.0;
  cdouble PMM  = +II*KOZ/2.0;
  cdouble FEE1 = -KZ/(2.0*Omega);
  cdouble FEE2 = +ZVAC/2.0;
  cdouble FEM1 = +1.0/(2.0*II*Omega);
  cdouble FEM2 = +II*KOZ*ZVAC/2.0;
  cdouble FME1 = -1.0/(2.0*II*Omega);
  cdouble FME2 = -II*KZ/(2.0*ZVAC);
  cdouble FMM1 = -1.0*KOZ/(2.0*Omega);
  cdouble FMM2 = +1.0/(2.0*ZVAC);

  /*--------------------------------------------------------------*/
  /*- these two lines assume the outer medium is vacuum ----------*/
  /*--------------------------------------------------------------*/
  cdouble EpsOut, MuOut;
  RegionMPs[nrOut]->GetEpsMu(Omega, &EpsOut, &MuOut);
  cdouble GammaE=(1.0/EpsIn - 1.0/EpsOut) * ZVAC;
  cdouble GammaM=(1.0/MuIn  - 1.0/MuOut) / ZVAC;
  cdouble Omega2=Omega*Omega;
  cdouble FEE3 = +0.25*GammaE/Omega2;
  cdouble FMM3 = +0.25*GammaM/Omega2;
  cdouble FEM3 = -0.25*ZVAC*GammaM/(II*Omega);
  cdouble FME3 = +0.25*GammaE/(II*Omega*ZVAC);

  int Order=4; // increase for greater accuracy in overlap integrals

  double PAbs=0.0;
  double Fx=0.0, Fy=0.0, Fz=0.0;
  double Taux=0.0, Tauy=0.0, Tauz=0.0;
#ifndef USE_OPENMP
  if (LogLevel>=SCUFF_VERBOSE2)
   Log(" no multithreading...");
#else
  int NumThreads=GetNumThreads();
  if (LogLevel>=SCUFF_VERBOSE2)
   Log(" OpenMP multithreading (%i threads)...",NumThreads);
#pragma omp parallel for schedule(dynamic,1),      \
                         num_threads(NumThreads),  \
                         reduction(+:PAbs, Fx, Fy, Fz, Taux, Tauy, Tauz)
#endif
//  for(int nea=0; nea<NE; nea++)
//   for(int neb=0; neb<NE; neb++)
   for(int neab=0; neab<NE*NE; neab++)
    { 
      int nea = neab/NE;
      int neb = neab%NE;

      /*--------------------------------------------------------------*/
      /*- Get various overlap integrals between basis function b_\alpha*/
      /*- and the fields of basis function b_\beta.                   */
      /*--------------------------------------------------------------*/
      cdouble be, bh;
      cdouble divbe[3], divbh[3], bxe[3], bxh[3];
      cdouble divbrxe[3], divbrxh[3], rxbxe[3], rxbxh[3];
      GetEPPFTMatrixElements(this, ns, ns, nea, neb, k, &be, &bh,
                             divbe, divbh, bxe, bxh,
                             divbrxe, divbrxh, rxbxe, rxbxh);

      /*--------------------------------------------------------------*/
      /*- Get various overlap integrals between basis function b_\alpha*/
      /*- and basis function b_\beta.                                 */
      /*--------------------------------------------------------------*/
      double Overlaps[20];
      S->GetOverlaps(nea, neb, Overlaps);
      double Divba_n_Divbb[3], nxba_Divbb[3]; 
      double Divba_rxn_Divbb[3], rxnxba_Divbb[3];
      Divba_n_Divbb[0]   = Overlaps[3];
      Divba_n_Divbb[1]   = Overlaps[6];
      Divba_n_Divbb[2]   = Overlaps[9];
      nxba_Divbb[0]      = Overlaps[4];
      nxba_Divbb[1]      = Overlaps[7];
      nxba_Divbb[2]      = Overlaps[10];
      Divba_rxn_Divbb[0] = Overlaps[12];
      Divba_rxn_Divbb[1] = Overlaps[15];
      Divba_rxn_Divbb[2] = Overlaps[18];
      rxnxba_Divbb[0]    = Overlaps[13];
      rxnxba_Divbb[1]    = Overlaps[16];
      rxnxba_Divbb[2]    = Overlaps[19];

      /*--------------------------------------------------------------*/
      /*- extract the surface-current coefficient either from the KN -*/
      /*- vector or the Sigma matrix.                                -*/
      /*--------------------------------------------------------------*/
      cdouble KK, KN, NK, NN;
      if (KNVector)
       { 
         cdouble kAlpha =       KNVector->GetEntry(Offset + 2*nea + 0);
         cdouble nAlpha = -ZVAC*KNVector->GetEntry(Offset + 2*nea + 1);
         cdouble kBeta  =       KNVector->GetEntry(Offset + 2*neb + 0);
         cdouble nBeta  = -ZVAC*KNVector->GetEntry(Offset + 2*neb + 1);

         KK = conj(kAlpha) * kBeta;
         KN = conj(kAlpha) * nBeta;
         NK = conj(nAlpha) * kBeta;
         NN = conj(nAlpha) * nBeta;
       }
      else
       { KK = SigmaMatrix->GetEntry(2*nea + 0, 2*neb+0);
         KN = SigmaMatrix->GetEntry(2*nea + 0, 2*neb+1);
         NK = SigmaMatrix->GetEntry(2*nea + 1, 2*neb+0);
         NN = SigmaMatrix->GetEntry(2*nea + 1, 2*neb+1);
       };

      /*--------------------------------------------------------------*/
      /*- assemble the various quantities.                            */
      /*--------------------------------------------------------------*/
      PAbs -= real( KK*PEE*be + KN*PEM*bh + NK*PME*bh + NN*PMM*be );

      Fx   -= real(   KK*(FEE1*divbe[0] + FEE2*bxh[0])
                    + KN*(FEM1*divbh[0] + FEM2*bxe[0])
                    + NK*(FME1*divbh[0] + FME2*bxe[0])
                    + NN*(FMM1*divbe[0] + FMM2*bxh[0])
                    + (FEE3*KK + FMM3*NN) * Divba_n_Divbb[0]
                    + (FEM3*KN + FME3*NK) * nxba_Divbb[0]
                  );

      Fy   -= real(   KK*(FEE1*divbe[1] + FEE2*bxh[1])
                    + KN*(FEM1*divbh[1] + FEM2*bxe[1])
                    + NK*(FME1*divbh[1] + FME2*bxe[1])
                    + NN*(FMM1*divbe[1] + FMM2*bxh[1]) 
                    + (FEE3*KK + FMM3*NN) * Divba_n_Divbb[1]
                    + (FEM3*KN + FME3*NK) * nxba_Divbb[1] 
                  );

      Fz   -= real(   KK*(FEE1*divbe[2] + FEE2*bxh[2])
                    + KN*(FEM1*divbh[2] + FEM2*bxe[2])
                    + NK*(FME1*divbh[2] + FME2*bxe[2])
                    + NN*(FMM1*divbe[2] + FMM2*bxh[2]) 
                    + (FEE3*KK + FMM3*NN) * Divba_n_Divbb[2]
                    + (FEM3*KN + FME3*NK) * nxba_Divbb[2] 
                  );

      Taux -= real(   KK*(FEE1*divbrxe[0] + FEE2*rxbxh[0])
                    + KN*(FEM1*divbrxh[0] + FEM2*rxbxe[0])
                    + NK*(FME1*divbrxh[0] + FME2*rxbxe[0])
                    + NN*(FMM1*divbrxe[0] + FMM2*rxbxh[0])
                    + (FEE3*KK + FMM3*NN) * Divba_rxn_Divbb[0]
                    + (FEM3*KN + FME3*NK) * rxnxba_Divbb[0]
                  );

      Tauy -= real(   KK*(FEE1*divbrxe[1] + FEE2*rxbxh[1])
                    + KN*(FEM1*divbrxh[1] + FEM2*rxbxe[1])
                    + NK*(FME1*divbrxh[1] + FME2*rxbxe[1])
                    + NN*(FMM1*divbrxe[1] + FMM2*rxbxh[1])
                    + (FEE3*KK + FMM3*NN) * Divba_rxn_Divbb[1]
                    + (FEM3*KN + FME3*NK) * rxnxba_Divbb[1]
                  );

      Tauz -= real(   KK*(FEE1*divbrxe[2] + FEE2*rxbxh[2])
                    + KN*(FEM1*divbrxh[2] + FEM2*rxbxe[2])
                    + NK*(FME1*divbrxh[2] + FME2*rxbxe[2])
                    + NN*(FMM1*divbrxe[2] + FMM2*rxbxh[2])
                    + (FEE3*KK + FMM3*NN) * Divba_rxn_Divbb[2]
                    + (FEM3*KN + FME3*NK) * rxnxba_Divbb[2]
                  );

#if 0
      // corrections to account for interior/exterior normal fields
      Taux -= 0.25*real( (GammaE*KK + GammaM*NN)/(Omega*Omega) )
              * Divba_nHat_Divbb[2];

      Tauy += 0.25*real( (ZVAC*GammaM*KN - GammaE*NK/ZVAC)/(II*Omega))
              * nxba_Divbb[2];

      Tauz += 0.0;
#endif

    };

  EPPFT[0] = PAbs;
  EPPFT[1] = (10.0/3.0)*Fx;
  EPPFT[2] = (10.0/3.0)*Fy;
  EPPFT[3] = (10.0/3.0)*Fz;
  EPPFT[4] = (10.0/3.0)*Taux;
  EPPFT[5] = (10.0/3.0)*Tauy;
  EPPFT[6] = (10.0/3.0)*Tauz;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
#if 0
void AssembledGMatrixBlock(RWGGeometry *G, int ns, cdouble Omega,
                           HMatrix *dG[6])
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  RWGSurface *S=G->Surfaces[ns];
  int nr = S->RegionIndices[1];
  if (nr==-1)
   { Warn("AssembledGMatrixBlock called for PEC object (ignoring)")
     return;
   };
  Log("Computing dG matrices for surface %i (%s)...",ns,S->Label);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble Eps, Mu;
  G->RegionMPs[nr]->GetEpsMu(Omega, &Eps, &Mu);
  cdouble k = sqrt(Eps*Mu)*Omega;
  cdouble Z = sqrt(Mu/Eps);
  cdouble PreFac1 = II*k*Z;
  cdouble PreFac2 = II*k;
  cdouble PreFac3 = -II*k/Z;

  /***************************************************************/
  /* loop over pairs of basis functions to assemble the matrix   */
  /***************************************************************/
  int NE = S->NumEdges;
  int nt, NumTasks, NumThreads = GetNumThreads();
#ifndef USE_OPENMP
  NumTasks=NumThreads=1;
  if (G->LogLevel>=SCUFF_VERBOSE2)
   Log(" no multithreading...");
#else
  NumTasks=NumThreads*100;
  if (NumTasks>S->NumEdges) 
  if (G->LogLevel>=SCUFF_VERBOSE2)
   Log(" OpenMP multithreading (%i threads,%i tasks)...",NumThreads,NumTasks);
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int nea=0; nea<S->NumEdges; nea++)
   for(int neb=nea; neb<S->NumEdges; neb++)
    { 
      if (nea==0 && G->LogLevel>SCUFF_VERBOSE2)
       LogPercent(nea,NE);     

      cdouble dGEntries[6], dCEntries[6];
      GetdEdgeEdgeInteractions(S, nea, neb, k, dGEntries, dCEntries);
      for(int Mu=0; Mu<6; Mu++)
       if (dG[Mu])
        { dG[Mu]->SetEntry(2*nea,   2*neb,   PreFac1*dGEntries[Mu]);
          dG[Mu]->SetEntry(2*nea,   2*neb+1, PreFac2*dCEntries[Mu]);
          dG[Mu]->SetEntry(2*nea+1, 2*neb,   PreFac2*dCEntries[Mu]);
          dG[Mu]->SetEntry(2*nea+1, 2*neb+1, PreFac3*dGEntries[Mu]);
        };
    };

}
#endif

} // namespace scuff
