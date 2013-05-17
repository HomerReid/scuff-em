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
 * Overlap.cc  -- computation of overlap integrals between RWG basis functions
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
/* contributions of a single panel to the overlap integrals    */
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

  double ZxL1[3], ZxL2[3], ZxDQ[3];
  VecCross(ZHat, L1, ZxL1);
  VecCross(ZHat, L2, ZxL2);
  VecCross(ZHat, DQ, ZxDQ);

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

  Overlaps[11] += PreFac * (-ZxL1[0]*BulletFactor2 - ZxL2[0]*BulletFactor3);
  Overlaps[12] += PreFac * (-4.0*ZxL1[0] - 2.0*ZxL2[0]) / 3.0;
  Overlaps[13] += PreFac * ZHat[0] * NablaCrossFactor;

  Overlaps[14] += PreFac * (-ZxL1[1]*BulletFactor2 - ZxL2[1]*BulletFactor3);
  Overlaps[15] += PreFac * (-4.0*ZxL1[1] - 2.0*ZxL2[1]) / 3.0;
  Overlaps[16] += PreFac * ZHat[1] * NablaCrossFactor;

  Overlaps[17] += PreFac * (-ZxL1[2]*BulletFactor2 - ZxL2[2]*BulletFactor3);
  Overlaps[18] += PreFac * (-4.0*ZxL1[2] - 2.0*ZxL2[2]) / 3.0;
  Overlaps[19] += PreFac * ZHat[2] * NablaCrossFactor;

}


/***************************************************************/
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
/*           (rHat x ) thrown in for torque purposes.          */
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

  RWGPanel *PAlphaP=Panels[EAlpha->iPPanel];
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

/*****************************************************************/
/* on entry, NeedMatrix is an array of 8 boolean flags, with     */
/* NeedMatrix[n] = 1 if the user wants overlap matrix #n.        */
/* (here 8 = SCUFF_NUM_OMATRICES).                               */
/*                                                               */
/* If NeedMatrix[n] = 1, then SArray must have at least n slots. */
/*                                                               */
/* If SArray[n] = 0 on entry, then a new SMatrix of the correct  */
/* size is allocated for that slot. If SArray[n] is nonzero but  */
/* points to an SMatrix of the incorrect size, then a new        */
/* SMatrix of the correct size is allocated, and SArray[n] is    */
/* overwritten with a pointer to this new SMatrix; the memory    */
/* allocated for the old SMatrix is not deallocated.             */
/*                                                               */
/* Omega is only referenced for the computation of force overlap */
/* matrices.                                                     */
/*                                                               */
/* ExteriorMP is only referenced for the computation of force    */
/* overlap matrices, and then only if it non-null; if ExteriorMP */
/* is null then the exterior medium is assumed to be vacuum.     */
/*****************************************************************/
void RWGSurface::GetOverlapMatrices(const bool NeedMatrix[SCUFF_NUM_OMATRICES],
                                    SMatrix *SArray[SCUFF_NUM_OMATRICES],
                                    cdouble Omega,
                                    MatProp *ExteriorMP)
{
  int NR  = NumBFs;

  /*--------------------------------------------------------------*/
  /*- the number of nonzero entries per row of the overlap matrix */
  /*- is fixed by the definition of the RWG basis function; each  */
  /*- RWG function overlaps with at most 5 RWG functions          */
  /*- (including itself), which gives 10 if we have both electric */
  /*- and magnetic currents.                                      */
  /*--------------------------------------------------------------*/
  int nnz = IsPEC ? 5 : 10;

  /*--------------------------------------------------------------*/
  /*- make sure each necessary slot of SArray points to an SMatrix*/  
  /*- of the appropriate size                                     */  
  /*--------------------------------------------------------------*/
  for(int n=0; n<SCUFF_NUM_OMATRICES; n++)
   { 
     if ( NeedMatrix[n] ) 
      { 
        if (     SArray[n] 
             && ( (SArray[n]->NR != NR) || (SArray[n]->NC != NR) )
           )
         { Warn("wrong-sized matrix passed to GetOverlapMatrices (reallocating)...");
           SArray[n]=0;
         };

        if (SArray[n]==0)
         SArray[n]=new SMatrix(NR, NR, LHM_COMPLEX);

	// TODO: avoid reallocation if shape is okay?
        SArray[n]->BeginAssembly(nnz*NR);

      };   
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble ZZ=ZVAC, K2=Omega*Omega;
  if (ExteriorMP)
   { 
     cdouble Eps, Mu;
     ExteriorMP->GetEpsMu(Omega, &Eps, &Mu);
     K2 *= Eps*Mu;
     ZZ *= sqrt(Mu/Eps);
   };
 
  /*--------------------------------------------------------------*/
  /*- sadly we use an N^2 algorithm for computing the O(N) matrix */
  /*- elements; this could be corrected at the expense of adding  */
  /*- some bookkeeping data (namely, a list within each RWGEdge   */
  /*- structure of the edges with which it overlaps), but for now */ 
  /*- the cost of this operation is negligible anyway so i don't  */
  /*- bother.                                                     */
  /*--------------------------------------------------------------*/
  int neAlpha, neBeta;
  double Overlaps[NUMOVERLAPS];
  cdouble XForce1, XForce2, YForce1, YForce2, ZForce1, ZForce2;
  cdouble XTorque1, XTorque2, YTorque1, YTorque2, ZTorque1, ZTorque2;
  for(neAlpha=0; neAlpha<NumEdges; neAlpha++)
   for(neBeta=0; neBeta<NumEdges; neBeta++)
    { 
      GetOverlaps(neAlpha, neBeta, Overlaps);
      if (Overlaps[0]==0.0) continue; 

      XForce1 = Overlaps[OVERLAP_BULLET_X] - Overlaps[OVERLAP_NABLANABLA_X]/K2;
      XForce2 = 2.0*Overlaps[OVERLAP_TIMESNABLA_X] / (II*Omega);

      YForce1 = Overlaps[OVERLAP_BULLET_Y] - Overlaps[OVERLAP_NABLANABLA_Y]/K2;
      YForce2 = 2.0*Overlaps[OVERLAP_TIMESNABLA_Y] / (II*Omega);

      ZForce1 = Overlaps[OVERLAP_BULLET_Z] - Overlaps[OVERLAP_NABLANABLA_Z]/K2;
      ZForce2 = 2.0*Overlaps[OVERLAP_TIMESNABLA_Z] / (II*Omega);

      XTorque1 = Overlaps[OVERLAP_RXBULLET_X] - Overlaps[OVERLAP_RXNABLANABLA_X]/K2;
      XTorque2 = 2.0*Overlaps[OVERLAP_RXTIMESNABLA_X] / (II*Omega);

      YTorque1 = Overlaps[OVERLAP_RXBULLET_Y] - Overlaps[OVERLAP_RXNABLANABLA_Y]/K2;
      YTorque2 = 2.0*Overlaps[OVERLAP_RXTIMESNABLA_Y] / (II*Omega);

      ZTorque1 = Overlaps[OVERLAP_RXBULLET_Z] - Overlaps[OVERLAP_RXNABLANABLA_Z]/K2;
      ZTorque2 = 2.0*Overlaps[OVERLAP_RXTIMESNABLA_Z] / (II*Omega);

      if (IsPEC)
       { 
         if ( NeedMatrix[SCUFF_OMATRIX_OVERLAP] )
          SArray[SCUFF_OMATRIX_OVERLAP]->SetEntry(neAlpha, neBeta, Overlaps[OVERLAP_OVERLAP]);

         // note in this case there is no entry in the power matrix 

         if ( NeedMatrix[SCUFF_OMATRIX_XFORCE] )
          SArray[SCUFF_OMATRIX_XFORCE]->SetEntry(neAlpha, neBeta, ZZ*XForce1);
         if ( NeedMatrix[SCUFF_OMATRIX_YFORCE] )
          SArray[SCUFF_OMATRIX_YFORCE]->SetEntry(neAlpha, neBeta, ZZ*YForce1);
         if ( NeedMatrix[SCUFF_OMATRIX_ZFORCE] )
          SArray[SCUFF_OMATRIX_ZFORCE]->SetEntry(neAlpha, neBeta, ZZ*ZForce1);
       }
      else
       {

         if ( NeedMatrix[SCUFF_OMATRIX_OVERLAP] )
          { SArray[SCUFF_OMATRIX_OVERLAP]->SetEntry(2*neAlpha+0, 2*neBeta+0, Overlaps[OVERLAP_OVERLAP]);
            SArray[SCUFF_OMATRIX_OVERLAP]->SetEntry(2*neAlpha+1, 2*neBeta+1, Overlaps[OVERLAP_OVERLAP]);
          };

         // EXPLAIN ME 
         if ( NeedMatrix[SCUFF_OMATRIX_POWER] )
          { SArray[SCUFF_OMATRIX_POWER]->SetEntry(2*neAlpha+0, 2*neBeta+1, -Overlaps[OVERLAP_CROSS]);
            SArray[SCUFF_OMATRIX_POWER]->SetEntry(2*neAlpha+1, 2*neBeta+0, Overlaps[OVERLAP_CROSS]);
          };

         if ( NeedMatrix[SCUFF_OMATRIX_XFORCE] )
          { SArray[SCUFF_OMATRIX_XFORCE]->SetEntry(2*neAlpha+0, 2*neBeta+0, ZZ*XForce1);
            SArray[SCUFF_OMATRIX_XFORCE]->SetEntry(2*neAlpha+0, 2*neBeta+1, XForce2);
            SArray[SCUFF_OMATRIX_XFORCE]->SetEntry(2*neAlpha+1, 2*neBeta+0, -XForce2);
            SArray[SCUFF_OMATRIX_XFORCE]->SetEntry(2*neAlpha+1, 2*neBeta+1, XForce1/ZZ);
          };

         if ( NeedMatrix[SCUFF_OMATRIX_YFORCE] )
          { SArray[SCUFF_OMATRIX_YFORCE]->SetEntry(2*neAlpha+0, 2*neBeta+0, ZZ*YForce1);
            SArray[SCUFF_OMATRIX_YFORCE]->SetEntry(2*neAlpha+0, 2*neBeta+1, YForce2);
            SArray[SCUFF_OMATRIX_YFORCE]->SetEntry(2*neAlpha+1, 2*neBeta+0, -YForce2);
            SArray[SCUFF_OMATRIX_YFORCE]->SetEntry(2*neAlpha+1, 2*neBeta+1, YForce1/ZZ);
          };

         if ( NeedMatrix[SCUFF_OMATRIX_ZFORCE] )
          { SArray[SCUFF_OMATRIX_ZFORCE]->SetEntry(2*neAlpha+0, 2*neBeta+0, ZZ*ZForce1);
            SArray[SCUFF_OMATRIX_ZFORCE]->SetEntry(2*neAlpha+0, 2*neBeta+1, ZForce2);
            SArray[SCUFF_OMATRIX_ZFORCE]->SetEntry(2*neAlpha+1, 2*neBeta+0, -ZForce2);
            SArray[SCUFF_OMATRIX_ZFORCE]->SetEntry(2*neAlpha+1, 2*neBeta+1, ZForce1/ZZ);
          };

         if ( NeedMatrix[SCUFF_OMATRIX_XTORQUE] )
          { SArray[SCUFF_OMATRIX_XTORQUE]->SetEntry(2*neAlpha+0, 2*neBeta+0, ZZ*XTorque1);
            SArray[SCUFF_OMATRIX_XTORQUE]->SetEntry(2*neAlpha+0, 2*neBeta+1, XTorque2);
            SArray[SCUFF_OMATRIX_XTORQUE]->SetEntry(2*neAlpha+1, 2*neBeta+0, -XTorque2);
            SArray[SCUFF_OMATRIX_XTORQUE]->SetEntry(2*neAlpha+1, 2*neBeta+1, XTorque1/ZZ);
          };

         if ( NeedMatrix[SCUFF_OMATRIX_YTORQUE] )
          { SArray[SCUFF_OMATRIX_YTORQUE]->SetEntry(2*neAlpha+0, 2*neBeta+0, ZZ*YTorque1);
            SArray[SCUFF_OMATRIX_YTORQUE]->SetEntry(2*neAlpha+0, 2*neBeta+1, YTorque2);
            SArray[SCUFF_OMATRIX_YTORQUE]->SetEntry(2*neAlpha+1, 2*neBeta+0, -YTorque2);
            SArray[SCUFF_OMATRIX_YTORQUE]->SetEntry(2*neAlpha+1, 2*neBeta+1, YTorque1/ZZ);
          };

         if ( NeedMatrix[SCUFF_OMATRIX_ZTORQUE] )
          { SArray[SCUFF_OMATRIX_ZTORQUE]->SetEntry(2*neAlpha+0, 2*neBeta+0, ZZ*ZTorque1);
            SArray[SCUFF_OMATRIX_ZTORQUE]->SetEntry(2*neAlpha+0, 2*neBeta+1, ZTorque2);
            SArray[SCUFF_OMATRIX_ZTORQUE]->SetEntry(2*neAlpha+1, 2*neBeta+0, -ZTorque2);
            SArray[SCUFF_OMATRIX_ZTORQUE]->SetEntry(2*neAlpha+1, 2*neBeta+1, ZTorque1/ZZ);
          };

       }; // if (IsPEC) ... else ... 
         
   }; // for(neAlpha...) ... for (neBeta...)

  for(int n=0; n<SCUFF_NUM_OMATRICES; n++)
   if ( NeedMatrix[n] ) 
    SArray[n]->EndAssembly();

}

/***************************************************************/
/* get power, force, and torque on a surface.                  */
/*                                                             */
/* return values:                                              */
/*  PFT[0]    = absorbed power                                 */
/*  PFT[1]    = total (absorbed + scattered) power (extinction)*/
/*  PFT[2..4] = x, y, z force                                  */
/*  PFT[5..7] = x, y, z torque                                 */
/***************************************************************/
void RWGGeometry::GetPFT2(HVector *KN, HVector *RHS, cdouble Omega,
                          int SurfaceIndex, double PFT[8])
{
  if (SurfaceIndex<0 || SurfaceIndex>=NumSurfaces)
   { memset(PFT,0,8*sizeof(double));
     Warn("GetPFT called for unknown surface #i",SurfaceIndex);
     return;  
   };
  RWGSurface *S=Surfaces[SurfaceIndex];

  /***************************************************************/
  /* compute overlap matrices.                                   */
  /* TODO: allow caller to pass caller-allocated storage for     */
  /* these matrices and for the KNTS vector below to avoid       */
  /* allocating on the fly                                       */
  /***************************************************************/
  bool NeedMatrix[SCUFF_NUM_OMATRICES];
  SMatrix *OMatrices[SCUFF_NUM_OMATRICES];
  for(int nm=0; nm<SCUFF_NUM_OMATRICES; nm++)
   { NeedMatrix[nm]=true;
     OMatrices[nm] = 0; 
   };
  S->GetOverlapMatrices(NeedMatrix, OMatrices, Omega, RegionMPs[0]);

  /***************************************************************/
  /* extract the chunk of the KN vector that is relevant for this*/
  /* surface and in the process (1) undo the SCUFF normalization */
  /* of the magnetic currents and (2) compute the total power    */
  /***************************************************************/
  int N = S->NumBFs;
  int Offset = BFIndexOffset[SurfaceIndex];
  HVector *KNTS=new HVector(N,LHM_COMPLEX); // 'KN, this surface'
  cdouble KAlpha, NAlpha, vEAlpha, vHAlpha;
  PFT[1]=0.0;
  if (S->IsPEC)
   { for(int ne=0; ne<S->NumEdges; ne++)
      { 
        KAlpha  = KN->GetEntry(Offset + ne);
        vEAlpha = RHS ? -ZVAC*RHS->GetEntry( Offset + ne ) : 0.0;

        PFT[1] += real( conj(KAlpha)*vEAlpha );

        KNTS->SetEntry(ne,  KAlpha );
      }
   }
  else //non-PEC
   { for(int ne=0; ne<S->NumEdges; ne++)
      { 
        KAlpha =       KN->GetEntry(Offset + 2*ne + 0);
        NAlpha = -ZVAC*KN->GetEntry(Offset + 2*ne + 1);

        vEAlpha = RHS ? -ZVAC*RHS->GetEntry( Offset + 2*ne + 0 ) : 0.0;
        vHAlpha = RHS ?  -1.0*RHS->GetEntry( Offset + 2*ne + 1 ) : 0.0;

        KNTS->SetEntry( 2*ne + 0,  KAlpha );
        KNTS->SetEntry( 2*ne + 1,  NAlpha );

        PFT[1] += real( conj(KAlpha)*vEAlpha + conj(NAlpha)*vHAlpha );
      };
   }
  PFT[1] *= 0.25;

  PFT[0] = 0.25*OMatrices[SCUFF_OMATRIX_POWER]->BilinearProductD(KNTS,KNTS);

  PFT[2] = 0.25*OMatrices[SCUFF_OMATRIX_XFORCE]->BilinearProductD(KNTS,KNTS);
  PFT[3] = 0.25*OMatrices[SCUFF_OMATRIX_YFORCE]->BilinearProductD(KNTS,KNTS);
  PFT[4] = 0.25*OMatrices[SCUFF_OMATRIX_ZFORCE]->BilinearProductD(KNTS,KNTS);

  PFT[5] = 0.25*OMatrices[SCUFF_OMATRIX_XTORQUE]->BilinearProductD(KNTS,KNTS);
  PFT[6] = 0.25*OMatrices[SCUFF_OMATRIX_YTORQUE]->BilinearProductD(KNTS,KNTS);
  PFT[7] = 0.25*OMatrices[SCUFF_OMATRIX_ZTORQUE]->BilinearProductD(KNTS,KNTS);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  // insert prefactors to get units right.
  // how it works: the force quantity that we just computed 
  // has units of 1 watt / c = (1 joule/s) * (1 s/nm) / 3  
  //                         = (1 nanoNewton / 3 )
  // so we multiply it by 3 to get a force in nanonewtons.
  // similarly for the torque: multiplying by 3 gives the torque
  // in nanoNewtons*microns (assuming the incident field was 
  // measured in units of volts / micron)
  for(int n=2; n<=7; n++)
   PFT[n] *= 3.0;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  delete KNTS;
  for(int nm=0; nm<SCUFF_NUM_OMATRICES; nm++)
   delete OMatrices[nm];

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RWGGeometry::GetPFT(HVector *KN, HVector *RHS, cdouble Omega,
                         int SurfaceIndex, double PFT[8])
{
  if (SurfaceIndex<0 || SurfaceIndex>=NumSurfaces)
   { Warn("invalid surface index passed to GetPFT",SurfaceIndex);
     memset(PFT, 0, 8*sizeof(double));
     return;
   };

  char *SSParmNames[4]={ const_cast<char *>("w"), const_cast<char *>("x"), 
                         const_cast<char *>("y"), const_cast<char *>("z") };
  cdouble SSParmValues[4];
  SSParmValues[0]=Omega*(MatProp::FreqUnit);

  RWGSurface *S=Surfaces[SurfaceIndex];
  int Offset=BFIndexOffset[SurfaceIndex];

  /*--------------------------------------------------------------*/
  /*- we need the material properties of the exterior medium for -*/
  /*- force and torque calculations                              -*/
  /*--------------------------------------------------------------*/
  UpdateCachedEpsMuValues(Omega);
  cdouble Eps = EpsTF[ S->RegionIndices[0] ];
  cdouble Mu  = MuTF[  S->RegionIndices[0] ];
  cdouble Z  = ZVAC*sqrt(Mu/Eps);
  cdouble K2 = Eps*Mu*Omega*Omega;

  /*--------------------------------------------------------------*/
  /*- sum contributions of all edges -----------------------------*/
  /*--------------------------------------------------------------*/
  int neAlpha, neBeta;
  cdouble KAlpha, KBeta, NAlpha, NBeta, vEAlpha, vHAlpha;
  cdouble M11, M12, M21, M22;
  double Overlaps[NUMOVERLAPS];
  memset(PFT, 0, 8*sizeof(double));
  for(neAlpha=0; neAlpha<S->NumEdges; neAlpha++)
   for(neBeta=0; neBeta<S->NumEdges; neBeta++)
    { 
      S->GetOverlaps(neAlpha, neBeta, Overlaps);
      if (Overlaps[0]==0.0) continue;

      if (S->IsPEC) 
       { KAlpha  = KN->GetEntry( Offset + neAlpha );
         KBeta   = KN->GetEntry( Offset + neBeta  );
         vEAlpha = RHS ? -ZVAC*RHS->GetEntry(Offset + neAlpha) : 0.0;
         vHAlpha = 0.0;
       }
      else
       { KAlpha  =       KN->GetEntry( Offset + 2*neAlpha + 0 );
         NAlpha  = -ZVAC*KN->GetEntry( Offset + 2*neAlpha + 1 );
         KBeta   =       KN->GetEntry( Offset + 2*neBeta  + 0 );
         NBeta   = -ZVAC*KN->GetEntry( Offset + 2*neBeta  + 1 );
         vEAlpha = RHS ? -ZVAC*RHS->GetEntry( Offset + 2*neAlpha + 0 ) : 0.0;
         vHAlpha = RHS ?  -1.0*RHS->GetEntry( Offset + 2*neAlpha + 1 ) : 0.0;
       };
  
     /*- absorbed power */
     if ( S->IsPEC && (S->SurfaceSigma || S->SurfaceSigmaMP) && Overlaps[OVERLAP_OVERLAP]!=0.0 )
      { 
        cdouble GZ;
        if (S->SurfaceSigmaMP)
         { GZ=S->SurfaceSigmaMP->GetEps(Omega);
           if (neAlpha==0 && neBeta==0) 
            Log("Running Surface-Conductivity power calculation (GZ=%s)",CD2S(GZ));
         }
        else if ( S->SurfaceSigma )
         { SSParmValues[1] = S->Edges[neAlpha]->Centroid[0];
           SSParmValues[2] = S->Edges[neAlpha]->Centroid[1];
           SSParmValues[3] = S->Edges[neAlpha]->Centroid[2];
           GZ=cevaluator_evaluate(S->SurfaceSigma, 4, SSParmNames, SSParmValues);
         };
        GZ*=ZVAC;
        PFT[0] += 0.5*real( conj(KAlpha)*KBeta * Overlaps[OVERLAP_OVERLAP] / GZ);
      }
     else
      PFT[0] += 0.5*real( conj(KAlpha)*NBeta * Overlaps[OVERLAP_CROSS] );

     /*- total power */
     if (neAlpha==neBeta)
      PFT[1] += 0.5*real( conj(KAlpha)*vEAlpha + conj(NAlpha)*vHAlpha );

     /*- X force */
     M11 = Z*(Overlaps[OVERLAP_BULLET_X] - Overlaps[OVERLAP_NABLANABLA_X]/K2); 
     M12 = +2.0*Overlaps[OVERLAP_TIMESNABLA_X]/ (II*Omega);
     M21 = -2.0*Overlaps[OVERLAP_TIMESNABLA_X]/ (II*Omega);
     M22 = (Overlaps[OVERLAP_BULLET_X] - Overlaps[OVERLAP_NABLANABLA_X]/K2)/Z;
     PFT[2] += 0.25*real(   conj(KAlpha)*M11*KBeta 
                          + conj(KAlpha)*M12*NBeta
                          + conj(NAlpha)*M21*KBeta 
                          + conj(NAlpha)*M22*NBeta );

     /*- Y force */
     M11 = Z*(Overlaps[OVERLAP_BULLET_Y] - Overlaps[OVERLAP_NABLANABLA_Y]/K2); 
     M12 = +2.0*Overlaps[OVERLAP_TIMESNABLA_Y]/ (II*Omega);
     M21 = -2.0*Overlaps[OVERLAP_TIMESNABLA_Y]/ (II*Omega);
     M22 = (Overlaps[OVERLAP_BULLET_Y] - Overlaps[OVERLAP_NABLANABLA_Y]/K2)/Z;
     PFT[3] += 0.25*real(   conj(KAlpha)*M11*KBeta 
                          + conj(KAlpha)*M12*NBeta
                          + conj(NAlpha)*M21*KBeta 
                          + conj(NAlpha)*M22*NBeta );

     /*- Z force */
     M11 = Z*(Overlaps[OVERLAP_BULLET_Z] - Overlaps[OVERLAP_NABLANABLA_Z]/K2); 
     M12 = +2.0*Overlaps[OVERLAP_TIMESNABLA_Z]/ (II*Omega);
     M21 = -2.0*Overlaps[OVERLAP_TIMESNABLA_Z]/ (II*Omega);
     M22 = (Overlaps[OVERLAP_BULLET_Z] - Overlaps[OVERLAP_NABLANABLA_Z]/K2)/Z;
     PFT[4] += 0.25*real(   conj(KAlpha)*M11*KBeta 
                          + conj(KAlpha)*M12*NBeta
                          + conj(NAlpha)*M21*KBeta 
                          + conj(NAlpha)*M22*NBeta );

    }; // for (neAlpha ... neBeta...)

  // force prefactors to get units right.
  // how it works: the force quantity that we just computed 
  // has units of 1 watt / c = (1 joule/s) * (1 s/nm) / 3  
  //                         = (1 nanoNewton / 3 )
  // so we multiply it by 3 to get a force in nanonewtons.
  PFT[2]*=3.0;
  PFT[3]*=3.0;
  PFT[4]*=3.0;

}

/***************************************************************/
/* alternative interface to GetPFT in which the caller         */
/* specifies the label of the surface instead of the index     */
/***************************************************************/
void RWGGeometry::GetPFT(HVector *KN, HVector *RHS, cdouble Omega,
                         char *SurfaceLabel, double PFT[8])
{
  /*--------------------------------------------------------------*/
  /*- find the surface in question -------------------------------*/
  /*--------------------------------------------------------------*/
  RWGSurface *S=GetSurfaceByLabel(SurfaceLabel);
  if (S)
   { 
     GetPFT(KN, RHS, Omega, S->Index, PFT);
   }
  else
   { Warn("unknown surface label %s passed to GetPFT",SurfaceLabel);
     memset(PFT, 0, 8*sizeof(double));
   };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
cdouble RWGGeometry::GetScatteredPower(HVector *KN, cdouble Omega, 
                                       int SurfaceIndex)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (SurfaceIndex<0 || SurfaceIndex>=NumSurfaces)
   { Warn("invalid surface index %i passed to GetScatteredPower",SurfaceIndex);
     return 0.0;
   };

  if (KN->RealComplex==LHM_REAL)
   return 0.0;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGSurface *S=Surfaces[SurfaceIndex];
  int NBF=S->NumBFs;
  HMatrix *M=new HMatrix(NBF, NBF, LHM_COMPLEX, LHM_SYMMETRIC);

  GetSSIArgStruct MyArgs, *Args=&MyArgs;
  InitGetSSIArgs(Args);

  Args->G=this;
  Args->Sa = Args->Sb = S;
  Args->Omega=Omega;
  Args->Symmetric=1;
  Args->B=M;

  Log("GetScatteredPower: Computing M0 matrix for surface %i (%s) ... ",SurfaceIndex,S->Label);
  
  int InteriorRegionIndex = S->RegionIndices[1];
  int SaveZeroed = RegionMPs[InteriorRegionIndex]->Zeroed;
  RegionMPs[InteriorRegionIndex]->Zeroed = 1;
  GetSurfaceSurfaceInteractions(Args);
  RegionMPs[InteriorRegionIndex]->Zeroed = SaveZeroed;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double Sign;
  cdouble PScat=0.0;
  cdouble *ZM=M->ZM;
  cdouble *ZKN=KN->ZV + BFIndexOffset[SurfaceIndex];
  // nr runs over rows, nc over columns, ne over matrix entries
  int nr, nc, ne;
  for(PScat=0.0, Sign=1.0, ne=nc=0; nc<NBF; nc++)
   for(nr=0; nr<NBF; nr++, ne++, Sign*=-1.0)
    PScat -= Sign*real( conj(ZKN[nr]) * ZM[ne] * ZKN[nc] );
 
  PScat *= 0.5*ZVAC;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  delete M;
  return PScat;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
cdouble RWGGeometry::GetScatteredPower(HVector *KN, cdouble Omega,
                                       char *SurfaceLabel)
{
  /*--------------------------------------------------------------*/
  /*- find the surface in question -------------------------------*/
  /*--------------------------------------------------------------*/
  RWGSurface *S=GetSurfaceByLabel(SurfaceLabel);
  if (S)
   { 
     return GetScatteredPower(KN, Omega, S->Index);
   }
  else
   { Warn("unknown surface label %s passed to GetScatteredPower",SurfaceLabel);
     return 0.0;
   };
}

}// namespace scuff
