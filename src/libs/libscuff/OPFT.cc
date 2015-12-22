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
 * OPFT.cc     -- computation of power, force, and torque using overlap
 *             -- integrals between RWG basis functions
 *
 * homer reid  -- 5/2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <ctype.h>

#include <libhrutil.h>

#include "libscuff.h"
#include "libscuffInternals.h"

#include "cmatheval.h"

#define II cdouble(0.0,1.0)

namespace scuff {

/***************************************************************/
/* these constants identify various types of overlap           */
/* integrals; they are only used in this file.                 */
/* Note that these constants are talking about types of        */
/* overlap *integrals*, not to be confused with the various    */
/* types of overlap *matrices,* which are index by different   */
/* constants defined in libscuff.h. (The entries of the overlap*/
/* matrices are linear combinations of various types of overlap*/
/* integrals.)                                                 */
/***************************************************************/
#define OVERLAP_OVERLAP           0
#define OVERLAP_CROSS             1
#define OVERLAP_BULLET_X          2    
#define OVERLAP_NABLANABLA_X      3
#define OVERLAP_TIMESNABLA_X      4
#define OVERLAP_BULLET_Y          5  
#define OVERLAP_NABLANABLA_Y      6
#define OVERLAP_TIMESNABLA_Y      7
#define OVERLAP_BULLET_Z          8  
#define OVERLAP_NABLANABLA_Z      9
#define OVERLAP_TIMESNABLA_Z     10
#define OVERLAP_RXBULLET_X       11    
#define OVERLAP_RXNABLANABLA_X   12
#define OVERLAP_RXTIMESNABLA_X   13
#define OVERLAP_RXBULLET_Y       14  
#define OVERLAP_RXNABLANABLA_Y   15
#define OVERLAP_RXTIMESNABLA_Y   16
#define OVERLAP_RXBULLET_Z       17  
#define OVERLAP_RXNABLANABLA_Z   18
#define OVERLAP_RXTIMESNABLA_Z   19

#define NUMOVERLAPS 20

/***************************************************************/
/* this is a helper function for GetOverlaps that computes the */
/* contributions of a single panel to the overlap integrals.   */
/***************************************************************/
void AddOverlapContributions(RWGSurface *S, RWGPanel *P, int iQa, int iQb, 
                             double Sign, double LL, double *Overlaps)
{
  double *Qa   = S->Vertices + 3*P->VI[ iQa ];
  double *QaP1 = S->Vertices + 3*P->VI[ (iQa+1)%3 ];
  double *QaP2 = S->Vertices + 3*P->VI[ (iQa+2)%3 ];
  double *Qb   = S->Vertices + 3*P->VI[ iQb ];
  double *ZHat = P->ZHat;

  double L1[3], L2[3], DQ[3];
  VecSub(QaP1, Qa, L1);
  VecSub(QaP2, QaP1, L2);
  VecSub(Qa, Qb, DQ);

  double ZxL1[3], ZxL2[3], ZxDQ[3], ZxQa[3], QaxZxL1[3], QaxZxL2[3];
  VecCross(ZHat,   L1,     ZxL1);
  VecCross(ZHat,   L2,     ZxL2);
  VecCross(ZHat,   DQ,     ZxDQ);
  VecCross(ZHat,   Qa,     ZxQa);
  VecCross(Qa,   ZxL1,  QaxZxL1);
  VecCross(Qa,   ZxL2,  QaxZxL2);

  double PreFac = Sign * LL / (2.0*P->Area);

  double L1dL1 = L1[0]*L1[0] + L1[1]*L1[1] + L1[2]*L1[2];
  double L1dL2 = L1[0]*L2[0] + L1[1]*L2[1] + L1[2]*L2[2];
  double L1dDQ = L1[0]*DQ[0] + L1[1]*DQ[1] + L1[2]*DQ[2];
  double L2dL2 = L2[0]*L2[0] + L2[1]*L2[1] + L2[2]*L2[2];
  double L2dDQ = L2[0]*DQ[0] + L2[1]*DQ[1] + L2[2]*DQ[2];

  double TimesFactor = (  (2.0*L1[0]+L2[0])*ZxDQ[0]  
                        + (2.0*L1[1]+L2[1])*ZxDQ[1] 
                        + (2.0*L1[2]+L2[2])*ZxDQ[2]  ) / 6.0;

  double BulletFactor1 =  (L1dL1 + L1dL2)/4.0 + L1dDQ/3.0 + L2dL2/12.0 + L2dDQ/6.0;
  double BulletFactor2 =  (L1dL1 + L1dL2)/5.0 + L1dDQ/4.0 + L2dL2/15.0 + L2dDQ/8.0;
  double BulletFactor3 =  L1dL1/10.0 + 2.0*L1dL2/15.0 + L1dDQ/8.0 + L2dL2/20.0 + L2dDQ/12.0;
  double NablaCrossFactor =  (L1dL1 + L1dL2)/2.0 + L2dL2/6.0;

  Overlaps[0]  += PreFac * BulletFactor1;
  Overlaps[1]  += PreFac * TimesFactor;

  Overlaps[2]  += PreFac * ZHat[0] * BulletFactor1;
  Overlaps[3]  += PreFac * ZHat[0] * 2.0;
  Overlaps[4]  += PreFac * (2.0*ZxL1[0] + ZxL2[0]) / 3.0;

  Overlaps[5]  += PreFac * ZHat[1] * BulletFactor1;
  Overlaps[6]  += PreFac * ZHat[1] * 2.0;
  Overlaps[7]  += PreFac * (2.0*ZxL1[1] + ZxL2[1]) / 3.0;

  Overlaps[8]  += PreFac * ZHat[2] * BulletFactor1;
  Overlaps[9]  += PreFac * ZHat[2] * 2.0;
  Overlaps[10] += PreFac * (2.0*ZxL1[2] + ZxL2[2]) / 3.0;

  Overlaps[11] -= PreFac * (ZxQa[0]*BulletFactor1 + ZxL1[0]*BulletFactor2 + ZxL2[0]*BulletFactor3);
  Overlaps[12] -= PreFac * (2.0*ZxQa[0] + 4.0*ZxL1[0]/3.0 + 2.0*ZxL2[0]/3.0);
  Overlaps[13] += PreFac * (ZHat[0]*NablaCrossFactor + 2.0*QaxZxL1[0]/3.0 + QaxZxL2[0]/3.0);

  Overlaps[14] -= PreFac * (ZxQa[1]*BulletFactor1 + ZxL1[1]*BulletFactor2 + ZxL2[1]*BulletFactor3);
  Overlaps[15] -= PreFac * (2.0*ZxQa[1] + 4.0*ZxL1[1]/3.0 + 2.0*ZxL2[1]/3.0);
  Overlaps[16] += PreFac * (ZHat[1]*NablaCrossFactor + 2.0*QaxZxL1[1]/3.0 + QaxZxL2[1]/3.0);

  Overlaps[17] -= PreFac * (ZxQa[2]*BulletFactor1 + ZxL1[2]*BulletFactor2 + ZxL2[2]*BulletFactor3);
  Overlaps[18] -= PreFac * (2.0*ZxQa[2] + 4.0*ZxL1[2]/3.0 + 2.0*ZxL2[2]/3.0);
  Overlaps[19] += PreFac * (ZHat[2]*NablaCrossFactor + 2.0*QaxZxL1[2]/3.0 + QaxZxL2[2]/3.0);

}

/***************************************************************/
/* Get the overlap integrals between a single pair of RWG      */
/* basis functions on an RWG surface.                          */
/*                                                             */
/* entries of output array:                                    */
/*                                                             */
/*  Overlaps [0] = O^{\bullet}_{\alpha\beta}                   */
/*   = \int f_a \cdot f_b                                      */
/*                                                             */
/*  Overlaps [1] = O^{\times}_{\alpha\beta}                    */
/*   = \int f_a \cdot (nHat \times f_b)                        */
/*                                                             */
/*  Overlaps [2] = O^{x,\bullet}_{\alpha\beta}                 */
/*   = \int nHat_x f_a \cdot f_b                               */
/*                                                             */
/*  Overlaps [3] = O^{x,\nabla\nabla}_{\alpha\beta}            */
/*   = \int nHat_x (\nabla \cdot f_a) (\nabla \cdot f_b)       */
/*                                                             */
/*  Overlaps [4] = O^{x,\times\nabla}_{\alpha\beta}            */
/*   = \int (nHat \times \cdot f_a)_x (\nabla \cdot f_b)       */
/*                                                             */
/*  [5,6,7]  = like [2,3,4] but with x-->y                     */
/*  [8,9,10] = like [2,3,4] but with x-->z                     */
/*                                                             */
/*  [11-19]: like [2...10] but with an extra factor of         */
/*           (rHat x ) thrown in for torque purposes           */
/*                                                             */
/* note: for now, the origin about which torque is computed    */
/* coincides with the origin of the coordinate system in which */
/* the surface mesh was defined (i.e. the point with coords    */
/* (0,0,0) in the mesh file, as transformed by any             */
/* GTransformations that have been applied since the mesh file */
/* was read in.) if you want to compute the torque about a     */
/* different origin, a quick-and-dirty procedure is to         */
/* apply a GTransformation to the surface, compute the         */
/* overlaps, then undo the GTransformation.                    */
/***************************************************************/
void RWGSurface::GetOverlaps(int neAlpha, int neBeta, double *Overlaps)
{
  RWGEdge *EAlpha = Edges[neAlpha];
  RWGEdge *EBeta  = Edges[neBeta];

  RWGPanel *PAlphaP = Panels[EAlpha->iPPanel];
  RWGPanel *PAlphaM = (EAlpha->iMPanel == -1) ? 0 : Panels[EAlpha->iMPanel];
  int iQPAlpha = EAlpha->PIndex;
  int iQMAlpha = EAlpha->MIndex;
  int iQPBeta  = EBeta->PIndex;
  int iQMBeta  = EBeta->MIndex;

  double LL = EAlpha->Length * EBeta->Length;

  memset(Overlaps,0,NUMOVERLAPS*sizeof(double));

  if ( EAlpha->iPPanel == EBeta->iPPanel )
   AddOverlapContributions(this, PAlphaP, iQPAlpha, iQPBeta,  1.0, LL, Overlaps);
  if ( EAlpha->iPPanel == EBeta->iMPanel )
   AddOverlapContributions(this, PAlphaP, iQPAlpha, iQMBeta, -1.0, LL, Overlaps);
  if ( (EAlpha->iMPanel!=-1) && (EAlpha->iMPanel == EBeta->iPPanel ) )
   AddOverlapContributions(this, PAlphaM, iQMAlpha, iQPBeta, -1.0, LL, Overlaps);
  if ( (EAlpha->iMPanel!=-1) && (EAlpha->iMPanel == EBeta->iMPanel ) )
   AddOverlapContributions(this, PAlphaM, iQMAlpha, iQMBeta,  1.0, LL, Overlaps);
  
}

/*****************************************************************/
/* this is a simpler interface to the above routine that returns */
/* the simple overlap integral and sets *pOTimes = crossed       */
/* overlap integral if it is non-NULL                            */
/*****************************************************************/
double RWGSurface::GetOverlap(int neAlpha, int neBeta, double *pOTimes)
{
  double Overlaps[NUMOVERLAPS];
  GetOverlaps(neAlpha, neBeta, Overlaps);
  if (pOTimes) *pOTimes=Overlaps[1];
  return Overlaps[0];
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int GetOverlappingEdgeIndices(RWGSurface *S, int nea, int nebArray[5])
{
  nebArray[0] = nea;
  int Count=1;

  RWGEdge *E   = S->Edges[nea];
  RWGPanel *PP = S->Panels[ E->iPPanel ]; 
  int      iQP = E->PIndex;
  nebArray[Count] = PP->EI[ (iQP+1)%3 ];
  if (nebArray[Count] >= 0) Count++;
  nebArray[Count] = PP->EI[ (iQP+2)%3 ];
  if (nebArray[Count] >= 0) Count++;

  if ( E->iMPanel == -1 )
   return Count;

  RWGPanel *PM = S->Panels[ E->iMPanel ];
  int      iQM = E->MIndex;
  nebArray[Count] = PM->EI[ (iQM+1)%3 ];
  if (nebArray[Count] >= 0) Count++;
  nebArray[Count] = PM->EI[ (iQM+2)%3 ];
  if (nebArray[Count] >= 0) Count++;
   
  return Count;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetOPFT(RWGGeometry *G, int SurfaceIndex, cdouble Omega,
             HVector *KNVector, HVector *RHS, HMatrix *DRMatrix,
             double PFT[NUMPFT], double **ByEdge)
{
  
  if (SurfaceIndex<0 || SurfaceIndex>=G->NumSurfaces)
   { memset(PFT,0,NUMPFT*sizeof(double));
     Warn("GetOPFT called for unknown surface #i",SurfaceIndex);
     return;
   };

  RWGSurface *S=G->Surfaces[SurfaceIndex];
  int Offset = G->BFIndexOffset[SurfaceIndex];
  int NE=S->NumEdges;

  /*--------------------------------------------------------------*/
  /*- get material parameters of exterior medium -----------------*/
  /*--------------------------------------------------------------*/
  cdouble ZZ=ZVAC, k2=Omega*Omega;
  cdouble Eps, Mu;
  G->RegionMPs[S->RegionIndices[0]]->GetEpsMu(Omega, &Eps, &Mu);
  k2 *= Eps*Mu;
  ZZ *= sqrt(Mu/Eps);

  // 20151003 surface-conductivity contribution to absorbed power
  cdouble GZ = 0.0;
  if (S->SurfaceSigmaMP)
   GZ = S->SurfaceSigmaMP->GetEps(Omega);

  /*--------------------------------------------------------------*/
  /*- initialize edge-by-edge contributions to zero --------------*/
  /*--------------------------------------------------------------*/
  if (ByEdge)
   { for(int nq=0; nq<NUMPFT; nq++)
      if (ByEdge[nq])
       memset(ByEdge[nq],0,NE*sizeof(double));
   };

  /***************************************************************/
  /* loop over all interior edges #nea                           */
  /***************************************************************/
  double PAbs=0.0, Fx=0.0, Fy=0.0, Fz=0.0, Taux=0.0, Tauy=0.0, Tauz=0.0;
  for(int nea=0; nea<NE; nea++)
   { 
     /*--------------------------------------------------------------*/
     /* populate an array whose indices are the 3 or 5 edges         */
     /* that have nonzero overlaps with edge #nea, then loop over    */
     /* those edges                                                  */
     /*--------------------------------------------------------------*/
     int nebArray[5];
     int nebCount= GetOverlappingEdgeIndices(S, nea, nebArray);
     for(int nneb=0; nneb<nebCount; nneb++)
      { 
        int neb=nebArray[nneb];
        if (neb<0) continue;
        double Overlaps[20];
        S->GetOverlaps(nea, neb, Overlaps);

        /*--------------------------------------------------------------*/
        /*--------------------------------------------------------------*/
        /*--------------------------------------------------------------*/
        cdouble KK, KN, NK, NN;
        if (KNVector && S->IsPEC)
         { 
           cdouble kAlpha =       KNVector->GetEntry(Offset + nea);
           cdouble kBeta  =       KNVector->GetEntry(Offset + neb);
           KK = conj(kAlpha) * kBeta;
           KN = NK = NN = 0.0;
         }
        else if (KNVector && !(S->IsPEC) )
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
         {
           KK = DRMatrix->GetEntry(Offset+2*neb+0, Offset+2*nea+0);
           KN = DRMatrix->GetEntry(Offset+2*neb+1, Offset+2*nea+0);
           NK = DRMatrix->GetEntry(Offset+2*neb+0, Offset+2*nea+1);
           NN = DRMatrix->GetEntry(Offset+2*neb+1, Offset+2*nea+1);
         };

       /*--------------------------------------------------------------*/
       /*--------------------------------------------------------------*/
       /*--------------------------------------------------------------*/
       // power
       double dPAbs = 0.25*real( (KN-NK) * Overlaps[OVERLAP_CROSS] );

       // 20151003 surface-conductivity contribution to absorbed power
       if (GZ!=0.0)
        dPAbs += 0.5*real(KK/GZ)*Overlaps[OVERLAP_OVERLAP];

       // force, torque
       double dF[3], dTau[3];
       dF[0] = 0.25*TENTHIRDS*
               real( -(KK*ZZ + NN/ZZ)*(Overlaps[OVERLAP_BULLET_X] - Overlaps[OVERLAP_NABLANABLA_X]/k2)
                     +(NK-KN)*2.0*Overlaps[OVERLAP_TIMESNABLA_X] / (II*Omega)
                   );

       dF[1] = 0.25*TENTHIRDS*
               real( -(KK*ZZ + NN/ZZ)*(Overlaps[OVERLAP_BULLET_Y] - Overlaps[OVERLAP_NABLANABLA_Y]/k2)
                     +(NK-KN)*2.0*Overlaps[OVERLAP_TIMESNABLA_Y] / (II*Omega)
                   );

       dF[2] = 0.25*TENTHIRDS*
               real( -(KK*ZZ + NN/ZZ)*(Overlaps[OVERLAP_BULLET_Z] - Overlaps[OVERLAP_NABLANABLA_Z]/k2)
                     +(NK-KN)*2.0*Overlaps[OVERLAP_TIMESNABLA_Z] / (II*Omega)
                   );

       dTau[0] = 0.25*TENTHIRDS*
                 real( -(KK*ZZ + NN/ZZ)*(Overlaps[OVERLAP_RXBULLET_X] - Overlaps[OVERLAP_RXNABLANABLA_X]/k2)
                       +(NK-KN)*2.0*Overlaps[OVERLAP_RXTIMESNABLA_X] / (II*Omega)
                     );

       dTau[1] = 0.25*TENTHIRDS*
                 real( -(KK*ZZ + NN/ZZ)*(Overlaps[OVERLAP_RXBULLET_Y] - Overlaps[OVERLAP_RXNABLANABLA_Y]/k2)
                       +(NK-KN)*2.0*Overlaps[OVERLAP_RXTIMESNABLA_Y] / (II*Omega)
                     );

       dTau[2] = 0.25*TENTHIRDS*
                 real( -(KK*ZZ + NN/ZZ)*(Overlaps[OVERLAP_RXBULLET_Z] - Overlaps[OVERLAP_RXNABLANABLA_Z]/k2)
                       +(NK-KN)*2.0*Overlaps[OVERLAP_RXTIMESNABLA_Z] / (II*Omega)
                     );

       /*--------------------------------------------------------------*/
       /*- accumulate contributions to full sums ----------------------*/
       /*--------------------------------------------------------------*/
       PAbs += dPAbs;
       Fx   += dF[0];
       Fy   += dF[1];
       Fz   += dF[2];
       Taux += dTau[0];
       Tauy += dTau[1];
       Tauz += dTau[2];

       /*--------------------------------------------------------------*/
       /*- accumulate contributions to by-edge sums -------------------*/
       /*--------------------------------------------------------------*/
       if (ByEdge) 
        {  
          if (ByEdge[PFT_PABS])   ByEdge[PFT_PABS][nea]     += dPAbs;
          if (ByEdge[PFT_XFORCE]) ByEdge[PFT_XFORCE][nea]   += dF[0];
          if (ByEdge[PFT_YFORCE]) ByEdge[PFT_YFORCE][nea]   += dF[1];
          if (ByEdge[PFT_ZFORCE]) ByEdge[PFT_ZFORCE][nea]   += dF[2];
          if (ByEdge[PFT_XTORQUE]) ByEdge[PFT_XTORQUE][nea] += dTau[0];
          if (ByEdge[PFT_YTORQUE]) ByEdge[PFT_YTORQUE][nea] += dTau[1];
          if (ByEdge[PFT_ZTORQUE]) ByEdge[PFT_ZTORQUE][nea] += dTau[2];
        };

      } // for (int nneb=... 

   }; // for(int nea=0; nea<S->NE; nea++)

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  PFT[PFT_PABS]    = PAbs;
  PFT[PFT_PSCAT]   = 0.0;
  PFT[PFT_XFORCE]  = Fx;
  PFT[PFT_YFORCE]  = Fy;
  PFT[PFT_ZFORCE]  = Fz;
  PFT[PFT_XTORQUE] = Taux;
  PFT[PFT_YTORQUE] = Tauy;
  PFT[PFT_ZTORQUE] = Tauz;

  /*--------------------------------------------------------------*/
  /*- if an RHS vector was specified, compute the extinction      */
  /*- (total power) and use it to compute the scattered power     */
  /*--------------------------------------------------------------*/
  if (KNVector && RHS)
   { double Extinction=0.0;
     for (int ne=0, nbf=0; ne<NE; ne++)
      { 
        cdouble kAlpha =   KNVector->GetEntry(Offset + nbf);
        cdouble vEAlpha = -ZVAC*RHS->GetEntry(Offset + nbf);
        nbf++;
        Extinction += 0.5*real( conj(kAlpha)*vEAlpha );
        if (S->IsPEC) continue;

        cdouble nAlpha  = -ZVAC*KNVector->GetEntry(Offset + nbf);
        cdouble vHAlpha =       -1.0*RHS->GetEntry(Offset + nbf);
        nbf++;
        Extinction += 0.5*real( conj(nAlpha)*vHAlpha );
      };
     PFT[PFT_PSCAT] = Extinction - PFT[PFT_PABS];
   };

} // GetOPFT

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetOPFTMatrices(RWGGeometry *G,
                     int SurfaceIndex, cdouble Omega,
                     HMatrix *QPFT[NUMPFT],
                     bool NeedMatrix[NUMPFT])
{
  if (SurfaceIndex<0 || SurfaceIndex>=G->NumSurfaces)
   { Warn("GetOPFTMatrices called for unknown surface #i",SurfaceIndex);
     return;
   };

  RWGSurface *S = G->Surfaces[SurfaceIndex];
  bool IsPEC    = S->IsPEC;
  int NE        = S->NumEdges;
  int NBF       = IsPEC ? NE : 2*NE;

  /***************************************************************/
  /* for all requested quantities, double check that the user    */
  /* gave us either a matrix of the correct size or else a NULL  */
  /* pointer. In the latter case, allocate a matrix of the       */
  /* correct size.                                               */
  /***************************************************************/
  for(int nq=0; nq<NUMPFT; nq++)
   { 
     if (nq==PFT_PSCAT) // there is no overlap matrix for scattered power
      continue;
  
     if (NeedMatrix[nq])
      { 
        if ( QPFT[nq] && (     QPFT[nq]->NR!=NBF
                            || QPFT[nq]->NC!=NBF
                            || QPFT[nq]->RealComplex!=LHM_COMPLEX
                         )
           )
         { Warn("incorrect Q matrix passed to GetOPFTMatrices (reallocating...)");
           QPFT[nq]=0;
         };

        if (!QPFT[nq])
         QPFT[nq]=new HMatrix(NBF, NBF, LHM_COMPLEX);

        QPFT[nq]->Zero();
      };
   };
  HMatrix *QPAbs, *QF[3], *QT[3];
  QPAbs = (NeedMatrix[PFT_PABS] && !IsPEC) ? QPFT[PFT_PABS] : 0;
  QF[0] = NeedMatrix[PFT_XFORCE]  ? QPFT[PFT_XFORCE]  : 0;
  QF[1] = NeedMatrix[PFT_YFORCE]  ? QPFT[PFT_YFORCE]  : 0;
  QF[2] = NeedMatrix[PFT_ZFORCE]  ? QPFT[PFT_ZFORCE]  : 0;
  QT[0] = NeedMatrix[PFT_XTORQUE] ? QPFT[PFT_XTORQUE] : 0;
  QT[1] = NeedMatrix[PFT_YTORQUE] ? QPFT[PFT_YTORQUE] : 0;
  QT[2] = NeedMatrix[PFT_ZTORQUE] ? QPFT[PFT_ZTORQUE] : 0;

  /*--------------------------------------------------------------*/
  /*- get material parameters of exterior medium -----------------*/
  /*--------------------------------------------------------------*/
  cdouble ZZ=ZVAC, k2=Omega*Omega;
  cdouble EpsRel, MuRel;
  G->RegionMPs[S->RegionIndices[0]]->GetEpsMu(Omega, &EpsRel, &MuRel);
  k2 *= EpsRel*MuRel;
  ZZ *= sqrt(MuRel/EpsRel);

  // 20151003 surface-conductivity contribution to absorbed power
  cdouble GZ = 0.0;
  if (S->SurfaceSigmaMP)
   GZ = S->SurfaceSigmaMP->GetEps(Omega);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int nea=0; nea<NE; nea++)
   { 
     /*--------------------------------------------------------------*/
     /* populate an array whose indices are the 3 or 5 edges         */
     /* that have nonzero overlaps with edge #nea, then loop over    */
     /* those edges                                                  */
     /*--------------------------------------------------------------*/
     int nebArray[5];
     int nebCount= GetOverlappingEdgeIndices(S, nea, nebArray);
     for(int nneb=0; nneb<nebCount; nneb++)
      { 
        int neb=nebArray[nneb];
        double Overlaps[20];
        S->GetOverlaps(nea, neb, Overlaps);

       // absorbed power
       if (QPAbs)
        { QPAbs->SetEntry(2*nea,2*neb+1, 0.25*Overlaps[OVERLAP_CROSS]);
          QPAbs->SetEntry(2*nea+1,2*neb,-0.25*Overlaps[OVERLAP_CROSS]);
          // 20151003 surface-conductivity contribution to absorbed power
          if (GZ!=0.0)
           QPAbs->SetEntry(2*nea,2*neb,0.5*real(1.0/GZ)*Overlaps[OVERLAP_OVERLAP]);
        };

       // force, torque
       for(int Mu=0; Mu<3; Mu++)
        { 
          if (QF[Mu])
           { cdouble Term1=-0.25*TENTHIRDS*(Overlaps[OVERLAP_BULLET_X + 3*Mu] - Overlaps[OVERLAP_NABLANABLA_X + 3*Mu]/k2);
             cdouble Term2=-0.50*TENTHIRDS*Overlaps[OVERLAP_TIMESNABLA_X + 3*Mu] / (II*Omega);
             if(IsPEC)
              QF[Mu]->SetEntry(nea,neb,ZZ*Term1);
             else
              { QF[Mu]->SetEntry(2*nea+0,2*neb+0,ZZ*Term1);
                QF[Mu]->SetEntry(2*nea+0,2*neb+1,Term2);
                QF[Mu]->SetEntry(2*nea+1,2*neb+0,-Term2);
                QF[Mu]->SetEntry(2*nea+1,2*neb+1,Term1/ZZ);
              };
           };

          if (QT[Mu])
           { cdouble Term1=-0.25*TENTHIRDS*(Overlaps[OVERLAP_RXBULLET_X + 3*Mu] - Overlaps[OVERLAP_RXNABLANABLA_X + 3*Mu]/k2);
             cdouble Term2=-0.50*TENTHIRDS*Overlaps[OVERLAP_RXTIMESNABLA_X + 3*Mu] / (II*Omega);
             if(IsPEC)
              QT[Mu]->SetEntry(nea,neb,ZZ*Term1);
             else
              { QT[Mu]->SetEntry(2*nea+0,2*neb+0,ZZ*Term1);
                QT[Mu]->SetEntry(2*nea+0,2*neb+1,Term2);
                QT[Mu]->SetEntry(2*nea+1,2*neb+0,-Term2);
                QT[Mu]->SetEntry(2*nea+1,2*neb+1,Term1/ZZ);
              };
           };
        };

      } // for (int nneb=... 

   }; // for(int nea=0; nea<S->NE; nea++)

} // GetOPFTMatrices

}// namespace scuff
