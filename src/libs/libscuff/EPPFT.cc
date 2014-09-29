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
/***************************************************************/
/***************************************************************/
void GetdGME_Far(RWGGeometry *G, int nsa, int nea, int nsb, int neb,
                 cdouble k, int Order, double *DX, 
                 cdouble GC[2], cdouble dG[6], cdouble dC[6])
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGSurface *Sa=G->Surfaces[nsa];
  RWGEdge *Ea=Sa->Edges[nea];
  double *QPa = Sa->Vertices + 3*Ea->iQP;
  double *V1a = Sa->Vertices + 3*Ea->iV1;
  double *V2a = Sa->Vertices + 3*Ea->iV2;
  double *QMa = Sa->Vertices + 3*Ea->iQM;

  RWGSurface *Sb=G->Surfaces[nsb];
  RWGEdge *Eb=Sb->Edges[neb];
  double *QPb = Sa->Vertices + 3*Eb->iQP;
  double *V1b = Sa->Vertices + 3*Eb->iV1;
  double *V2b = Sa->Vertices + 3*Eb->iV2;
  double *QMb = Sa->Vertices + 3*Eb->iQM;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double AaPlus[3], AaMinus[3], Ba[3];
  VecSub(V1a, QPa, AaPlus);
  VecSub(V1a, QMa, AaMinus);
  VecSub(V2a, V1a, Ba);

  double AbPlus[3], AbMinus[3], Bb[3];
  VecSub(V1b, QPb, AbPlus);
  VecSub(V1b, QMb, AbMinus);
  VecSub(V2b, V1b, Bb);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int NumPts;
  double *TCR=GetTCR(Order, &NumPts);
  memset(GC, 0, 2*sizeof(cdouble));
  memset(dG, 0, 6*sizeof(cdouble));
  memset(dC, 0, 6*sizeof(cdouble));
  cdouble k2 = k*k;
  for(int npa=0, ncpa=0; npa<NumPts; npa++)
   { 
     /***************************************************************/
     /***************************************************************/
     /***************************************************************/
     double ua=TCR[ncpa++];
     double va=TCR[ncpa++];
     double wa=TCR[ncpa++];
     ua+=va;

     double XaPlus[3], FaPlus[3], XaMinus[3], FaMinus[3];
     for(int Mu=0; Mu<3; Mu++)
      { 
        FaPlus[Mu] = ua*AaPlus[Mu] + va*Ba[Mu];
        XaPlus[Mu] = FaPlus[Mu] + QPa[Mu];

        FaMinus[Mu] = ua*AaMinus[Mu] + va*Ba[Mu];
        XaMinus[Mu] = FaMinus[Mu] + QMa[Mu];
      };

     for(int npb=0, ncpb=0; npb<NumPts; npb++)
      { 
        /***************************************************************/
        /***************************************************************/
        /***************************************************************/
        double ub=TCR[ncpb++];
        double vb=TCR[ncpb++];
        double wb=TCR[ncpb++];
        ub+=vb;
   
        double XbPlus[3], FbPlus[3], XbMinus[3], FbMinus[3];
        for(int Mu=0; Mu<3; Mu++)
         { 
           FbPlus[Mu] = ub*AbPlus[Mu] + vb*Bb[Mu];
           XbPlus[Mu] = FbPlus[Mu] + QPb[Mu];
   
           FbMinus[Mu] = ub*AbMinus[Mu] + vb*Bb[Mu];
           XbMinus[Mu] = FbMinus[Mu] + QMb[Mu];
         };

        for(int aSign=0; aSign<=1; aSign++)
         for(int bSign=0; bSign<=1; bSign++)
          { 
            double Sign = aSign==bSign ? 1.0 : -1.0;
            double *Xa = aSign==1 ? XaPlus : XaMinus;
            double *Fa = aSign==1 ? FaPlus : FaMinus;
            double *Xb = bSign==1 ? XbPlus : XbMinus;
            double *Fb = bSign==1 ? FbPlus : FbMinus;

            double R[3];
            R[0]=Xa[0] - Xb[0] - (DX==0 ? 0.0 : DX[0]);
            R[1]=Xa[1] - Xb[1] - (DX==0 ? 0.0 : DX[1]);
            R[2]=Xa[2] - Xb[2] - (DX==0 ? 0.0 : DX[2]);
            double r2    = R[0]*R[0] + R[1]*R[1] + R[2]*R[2];
            double r     = sqrt(r2);
            cdouble ikr  = II*k*r;
            cdouble Phi  = exp(ikr) / (4.0*M_PI*r);
            cdouble Psi  = (ikr-1.0) * Phi / r2;
            cdouble Zeta = (3.0 + ikr*(-3.0 + ikr))*Phi/(r2*r2);

            cdouble dGMu[3], dGMuNu[3][3];
            for(int Mu=0; Mu<3; Mu++)
             { dGMu[Mu] = R[Mu]*Psi;
               for(int Nu=0; Nu<3; Nu++)
                dGMuNu[Mu][Nu] 
                 = (Mu==Nu ? Psi : 0.0) + R[Mu]*R[Nu]*Zeta;
             };

            for(int Nu=0; Nu<3; Nu++)
             { GC[0]
                += wa*wb*Sign*( Fa[Nu]*Fb[Nu]*Phi + 2.0*Fa[Nu]*dGMu[Nu]/k2);
               int NP1=(Nu+1)%3, NP2=(Nu+2)%3;
               GC[1]
                += wa*wb*Sign*( Fa[Nu]*(dGMu[NP1]*Fb[NP2]-dGMu[NP2]*Fb[NP1]));
             };

            for(int Mu=0; Mu<3; Mu++)
             for(int Nu=0; Nu<3; Nu++)
              { dG[Mu]
                 += wa*wb*Sign * ( Fa[Nu]*Fb[Nu]*dGMu[Mu]
                                  + 2.0*Fa[Nu]*dGMuNu[Mu][Nu]/k2
                                );

                int NP1=(Nu+1)%3, NP2=(Nu+2)%3;
                dC[Mu] 
                 += wa*wb*Sign*Fa[Nu]* (  dGMuNu[Mu][NP1]*Fb[NP2]
                                         -dGMuNu[Mu][NP2]*Fb[NP1] );
              };
#if 0
cdouble GMuNu[3][3], CMuNu[3][3], dGMuNuRho[3][3][3], dCMuNuRho[3][3][3];
CalcGC(R, k, 1.0, 1.0, GMuNu, CMuNu, dGMuNuRho, dCMuNuRho);
for(int Mu=0; Mu<3; Mu++)
 for(int Nu=0; Nu<3; Nu++)
  for(int Rho=0; Rho<3; Rho++)
   { dG[Mu]+=wa*wb*Sign*Fa[Nu]*dGMuNuRho[Nu][Rho][Mu]*Fb[Rho];
     dC[Mu]+=wa*wb*Sign*Fa[Nu]*dCMuNuRho[Nu][Rho][Mu]*Fb[Rho];
   };
#endif

          }; // for(int aSign=0; aSign<=1; aSign++) ... 

      }; // for(int npb=0, ncpb=0; npb<NumPts; npb++)

   }; // for(int npa=0, ncpa=0; npa<NumPts; npa++)

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double PreFac = Ea->Length * Eb->Length;
  cdouble MIK = -1.0*II*k;
  GC[0]*=PreFac;
  GC[1]*=PreFac / MIK;
  for(int Mu=0; Mu<6; Mu++)
   { dG[Mu] *= PreFac;
     dC[Mu] *= PreFac / MIK;
   };
   
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetdGME_Near(RWGGeometry *G, int nsa, int nea, int nsb, int neb,
                  cdouble k, int Order, double *DX,
                  cdouble GC[2], cdouble dG[6], cdouble dC[6])
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGSurface *S=G->Surfaces[nsa];
  RWGEdge *E=S->Edges[nea];
  double *QP = S->Vertices + 3*E->iQP;
  double *V1 = S->Vertices + 3*E->iV1;
  double *V2 = S->Vertices + 3*E->iV2;
  double *QM = S->Vertices + 3*E->iQM;

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
  memset(GC, 0, 2*sizeof(cdouble));
  memset(dG, 0, 6*sizeof(cdouble));
  memset(dC, 0, 6*sizeof(cdouble));
  for(int np=0, ncp=0; np<NumPts; np++)
   { 
     /***************************************************************/
     /***************************************************************/
     /***************************************************************/
     double u=TCR[ncp++];
     double v=TCR[ncp++];
     double w=TCR[ncp++];
     u+=v;

     double XPlus[3], FPlus[3], XMinus[3], FMinus[3];
     for(int Mu=0; Mu<3; Mu++)
      { 
        FPlus[Mu] = u*APlus[Mu] + v*B[Mu];
        XPlus[Mu] = FPlus[Mu] + QP[Mu] - (DX==0 ? 0.0 : DX[Mu]);

        FMinus[Mu] = u*AMinus[Mu] + v*B[Mu];
        XMinus[Mu] = FMinus[Mu] + QM[Mu] - (DX==0 ? 0.0 : DX[Mu]);
      };

     /***************************************************************/
     /***************************************************************/
     /***************************************************************/
     cdouble ePlus[3], dePlus[3][3], eMinus[3], deMinus[3][3];
     cdouble hPlus[3], dhPlus[3][3], hMinus[3], dhMinus[3][3];
     GetReducedFields_Nearby(G, nsb, neb, XPlus, k, ePlus, hPlus, dePlus, dhPlus);
     GetReducedFields_Nearby(G, nsb, neb, XMinus, k, eMinus, hMinus, deMinus, dhMinus);

     for(int Mu=0; Mu<3; Mu++)
      for(int Nu=0; Nu<3; Nu++)
       { 
         if (Mu==0)
          { GC[0] += w*(FPlus[Nu]*ePlus[Nu] - FMinus[Nu]*eMinus[Nu]);
            GC[1] += w*(FPlus[Nu]*hPlus[Nu] - FMinus[Nu]*hMinus[Nu]);
          };
         dG[Mu] += w*(FPlus[Nu]*dePlus[Mu][Nu] - FMinus[Nu]*deMinus[Mu][Nu]);
         dC[Mu] += w*(FPlus[Nu]*dhPlus[Mu][Nu] - FMinus[Nu]*dhMinus[Mu][Nu]);
       };

   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble MIK = -1.0*II*k;
  GC[0] *= E->Length;
  GC[1] *= E->Length / MIK;
  for(int Mu=0; Mu<6; Mu++)
   { dG[Mu] *= E->Length;
     dC[Mu] *= E->Length / MIK;
   };
   
}

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
                            cdouble be[1],    cdouble bh[1],
                            cdouble divbe[3], cdouble divbh[3],
                            cdouble bxe[3],   cdouble bxh[3],
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
  memset(divbe, 0, 3*sizeof(cdouble));
  memset(divbh, 0, 3*sizeof(cdouble));
  memset(bxe,   0, 3*sizeof(cdouble));
  memset(bxh,   0, 3*sizeof(cdouble));
  for(int np=0, ncp=0; np<NumPts; np++)
   { 
     /***************************************************************/
     /***************************************************************/
     /***************************************************************/
     double u=TCR[ncp++];
     double v=TCR[ncp++];
     double w=TCR[ncp++] * (E->Length);
     u+=v;

     double XPlus[3], FPlus[3], XMinus[3], FMinus[3];
     for(int Mu=0; Mu<3; Mu++)
      { 
        FPlus[Mu] = u*APlus[Mu] + v*B[Mu];
        XPlus[Mu] = FPlus[Mu] + QP[Mu];

        FMinus[Mu] = u*AMinus[Mu] + v*B[Mu];
        XMinus[Mu] = FMinus[Mu] + QM[Mu];
      };

     /***************************************************************/
     /***************************************************************/
     /***************************************************************/
     cdouble ePlus[3], eMinus[3];
     cdouble hPlus[3], hMinus[3];
     GetReducedFields_Nearby(G, nsb, neb, XPlus, k, ePlus, hPlus);
     GetReducedFields_Nearby(G, nsb, neb, XMinus, k, eMinus, hMinus);

     for(int Mu=0; Mu<3; Mu++)
       {
         be[0] += w*(FPlus[Mu]*ePlus[Mu] - FMinus[Mu]*eMinus[Mu]); 
         bh[0] += w*(FPlus[Mu]*hPlus[Mu] - FMinus[Mu]*hMinus[Mu]);

         divbe[Mu] += 2.0*w*(ePlus[Mu] - eMinus[Mu]);
         divbh[Mu] += 2.0*w*(hPlus[Mu] - hMinus[Mu]);

         int MP1=(Mu+1)%3, MP2=(Mu+2)%3;
         bxe[Mu] += w*( (FPlus[MP1]*ePlus[MP2]   - FPlus[MP2]*ePlus[MP1])
                       -(FMinus[MP1]*eMinus[MP2] - FMinus[MP2]*eMinus[MP1])
                      );

         bxh[Mu] += w*( (FPlus[MP1]*hPlus[MP2]   - FPlus[MP2]*hPlus[MP1])
                       -(FMinus[MP1]*hMinus[MP2] - FMinus[MP2]*hMinus[MP1])
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
  int nr = S->RegionIndices[1];
  if ( nr==-1 || S->IsPEC ) 
   { memset(EPPFT, 0, 7*sizeof(double));
     return;
   };

  cdouble Eps, Mu;
  RegionMPs[nr]->GetEpsMu(Omega, &Eps, &Mu);
  cdouble k = Omega * sqrt(Eps*Mu);
  cdouble ZRel = sqrt(Mu/Eps);
  cdouble KZ  = k*ZVAC*ZRel;
  cdouble KOZ = k/(ZVAC*ZRel);

  /*--------------------------------------------------------------*/
  /*- these two lines assume the outer medium is vacuum ----------*/
  /*--------------------------------------------------------------*/
  cdouble GammaE=(1.0/Eps - 1.0) * ZVAC;
  cdouble GammaM=(1.0/Mu  - 1.0) / ZVAC;

  int NE = S->NumEdges;
  int Offset = BFIndexOffset[ns];
  int Order=4;

  double PAbs=0.0;
  double Fx=0.0, Fy=0.0, Fz=0.0;
  double Taux=0.0, Tauy=0.0, Tauz=0.0;
  int NumThreads;
#ifndef USE_OPENMP
  NumThreads=1;
  if (LogLevel>=SCUFF_VERBOSE2)
   Log(" no multithreading...");
#else
  NumThreads=GetNumThreads();
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
      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      cdouble be, bh, divbe[3], divbh[3], bxe[3], bxh[3];
      GetEPPFTMatrixElements(this, ns, ns, nea, neb, k,
                             &be, &bh, divbe, divbh, bxe, bxh);

      double Overlaps[20];
      S->GetOverlaps(nea, neb, Overlaps);
      double Divba_nHat_Divbb[3], nxba_Divbb[3];
      Divba_nHat_Divbb[0] = Overlaps[3];
      Divba_nHat_Divbb[1] = Overlaps[6];
      Divba_nHat_Divbb[2] = Overlaps[9];
      nxba_Divbb[0]       = Overlaps[4];
      nxba_Divbb[1]       = Overlaps[7];
      nxba_Divbb[2]       = Overlaps[10];

      cdouble PEE, PEM, PME, PMM;
      cdouble FEE1, FEE2, FEM1, FEM2, FME1, FME2, FMM1, FMM2;
      PEE  = +II*KZ/2.0;
      PEM  = -1.0/2.0;
      PME  = +1.0/2.0;
      PMM  = +II*KOZ/2.0;
      FEE1 = -KZ/(2.0*Omega);
      FEE2 = ZVAC/2.0;
      FEM1 = 1.0/(2.0*II*Omega);
      FEM2 = II*KOZ*ZVAC/2.0;
      FME1 = -1.0/(2.0*II*Omega);
      FME2 = -II*KZ/(2.0*ZVAC);
      FMM1 = -1.0*KOZ/(2.0*Omega);
      FMM2 = 1.0/(2.0*ZVAC);

      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
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
      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      PAbs += -2.0*real( KK*PEE*be + KN*PEM*bh + NK*PME*bh + NN*PMM*be );

      Fx   += -2.0*real(  KK*(FEE1*divbe[0] + FEE2*bxh[0])
                    +KN*(FEM1*divbh[0] + FEM2*bxe[0])
                    +NK*(FME1*divbh[0] + FME2*bxe[0])
                    +NN*(FMM1*divbe[0] + FMM2*bxh[0]) );

      Fy   += -2.0*real(  KK*(FEE1*divbe[1] + FEE2*bxh[1])
                    +KN*(FEM1*divbh[1] + FEM2*bxe[1])
                    +NK*(FME1*divbh[1] + FME2*bxe[1])
                    +NN*(FMM1*divbe[1] + FMM2*bxh[1]) );

      Fz   += -2.0*real(  KK*(FEE1*divbe[2] + FEE2*bxh[2])
                    +KN*(FEM1*divbh[2] + FEM2*bxe[2])
                    +NK*(FME1*divbh[2] + FME2*bxe[2])
                    +NN*(FMM1*divbe[2] + FMM2*bxh[2]) );

      Taux += 0.5*real( (GammaE*KK + GammaM*NN)/(Omega*Omega) )
              * Divba_nHat_Divbb[2];

      Tauy += 0.5*real( (ZVAC*GammaM*KN - GammaE*NK/ZVAC)/(II*Omega))
              * nxba_Divbb[2];

      Tauz += 0.0;
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
