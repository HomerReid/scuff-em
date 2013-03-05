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
#define OVERLAP_OVERLAP     0
#define OVERLAP_CROSS       1
#define OVERLAP_XBULLET     2  
#define OVERLAP_XNABLANABLA 3
#define OVERLAP_XTIMESNABLA 4
#define OVERLAP_YBULLET     5  
#define OVERLAP_YNABLANABLA 6
#define OVERLAP_YTIMESNABLA 7
#define OVERLAP_ZBULLET     8  
#define OVERLAP_ZNABLANABLA 9
#define OVERLAP_ZTIMESNABLA 10

#define NUMOVERLAPS 11

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

  double A[3], B[3], DQ[3];
  VecSub(QaP1, Qa, A);
  VecSub(QaP2, QaP1, B);
  VecSub(Qa, Qb, DQ);

  double ZxA[3], ZxB[3], ZxDQ[3];
  VecCross(ZHat, A, ZxA);
  VecCross(ZHat, B, ZxB);
  VecCross(ZHat, DQ, ZxDQ);

  double PreFac = Sign * LL / (2.0*P->Area);

  double Bullet=0.0, Times=0.0;
  for(int Mu=0; Mu<3; Mu++)
   { Bullet += A[Mu]*(A[Mu]/4.0 + B[Mu]/4.0 + DQ[Mu]/3.0) + B[Mu]*(B[Mu]/12.0 + DQ[Mu]/6.0);
     Times  += (A[Mu]+0.5*B[Mu])*ZxDQ[Mu]/3.0;
   };

  Overlaps[0]  += PreFac*Bullet;
  Overlaps[1]  += PreFac*Times;

  Overlaps[2]  += PreFac * ZHat[0] * Bullet;
  Overlaps[3]  += PreFac * ZHat[0] * 2.0;
  Overlaps[4]  += PreFac * (ZxA[0]/3.0 + ZxB[0]/6.0) * 2.0;

  Overlaps[5]  += PreFac * ZHat[1] * Bullet;
  Overlaps[6]  += PreFac * ZHat[1] * 2.0;
  Overlaps[7]  += PreFac * (ZxA[1]/3.0 + ZxB[1]/6.0) * 2.0;

  Overlaps[8]  += PreFac * ZHat[2] * Bullet;
  Overlaps[9]  += PreFac * ZHat[2] * 2.0;
  Overlaps[10] += PreFac * (ZxA[2]/3.0 + ZxB[2]/6.0) * 2.0;

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
/***************************************************************/
void RWGSurface::GetOverlaps(int neAlpha, int neBeta, double *Overlaps)
{
  RWGEdge *EAlpha = Edges[neAlpha];
  RWGEdge *EBeta  = Edges[neBeta];

  RWGPanel *PAlphaP=Panels[EAlpha->iPPanel];
  RWGPanel *PAlphaM = (EAlpha->iMPanel == -1) ? 0 : Panels[EAlpha->iMPanel];
  int iQPAlpha = EAlpha->PIndex;
  int iQMAlpha = EAlpha->MIndex;
  int  iQPBeta  = EBeta->PIndex;
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
/* on entry, NeedMatrix is an array of 5 boolean flags, with     */
/* NeedMatrix[n] = 1 if the user wants overlap matrix #n.        */
/* (here 5 = SCUFF_NUM_OMATRICES).                               */
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
  cdouble Z=ZVAC, K2=Omega*Omega;
  if (ExteriorMP)
   { 
     cdouble Eps, Mu;
     ExteriorMP->GetEpsMu(Omega, &Eps, &Mu);
     K2 *= Eps*Mu;
     Z *= sqrt(Mu/Eps);
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
  for(neAlpha=0; neAlpha<NumEdges; neAlpha++)
   for(neBeta=0; neBeta<NumEdges; neBeta++)
    { 
      GetOverlaps(neAlpha, neBeta, Overlaps);
      if (Overlaps[0]==0.0) continue; 

      XForce1 = Overlaps[OVERLAP_XBULLET] - Overlaps[OVERLAP_XNABLANABLA]/K2;
      XForce2 = 2.0*Overlaps[OVERLAP_XTIMESNABLA] / (II*Omega);

      YForce1 = Overlaps[OVERLAP_YBULLET] - Overlaps[OVERLAP_YNABLANABLA]/K2;
      YForce2 = 2.0*Overlaps[OVERLAP_YTIMESNABLA] / (II*Omega);

      ZForce1 = Overlaps[OVERLAP_ZBULLET] - Overlaps[OVERLAP_ZNABLANABLA]/K2;
      ZForce2 = 2.0*Overlaps[OVERLAP_ZTIMESNABLA] / (II*Omega);

      if (IsPEC)
       { 
         if ( NeedMatrix[SCUFF_OMATRIX_OVERLAP] )
          SArray[SCUFF_OMATRIX_OVERLAP]->SetEntry(neAlpha, neBeta, Overlaps[OVERLAP_OVERLAP]);

         // note in this case there is no entry in the power matrix 

         if ( NeedMatrix[SCUFF_OMATRIX_XFORCE] )
          SArray[SCUFF_OMATRIX_XFORCE]->SetEntry(neAlpha, neBeta, Z*XForce1);
         if ( NeedMatrix[SCUFF_OMATRIX_YFORCE] )
          SArray[SCUFF_OMATRIX_YFORCE]->SetEntry(neAlpha, neBeta, Z*YForce1);
         if ( NeedMatrix[SCUFF_OMATRIX_ZFORCE] )
          SArray[SCUFF_OMATRIX_ZFORCE]->SetEntry(neAlpha, neBeta, Z*ZForce1);
       }
      else
       {

double Sign=1.0;
char *str=getenv("SCUFFSIGN");
if (str && str[0]=='M')
 Sign=-1.0;
if (str)
 printf("STR was %s, SIGN=%e\n",str,Sign);
else
printf("SIGN=%e\n",Sign);

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
          { SArray[SCUFF_OMATRIX_XFORCE]->SetEntry(2*neAlpha+0, 2*neBeta+0, Z*XForce1);
            SArray[SCUFF_OMATRIX_XFORCE]->SetEntry(2*neAlpha+0, 2*neBeta+1, XForce2);
            SArray[SCUFF_OMATRIX_XFORCE]->SetEntry(2*neAlpha+1, 2*neBeta+0, Sign*XForce2);
            SArray[SCUFF_OMATRIX_XFORCE]->SetEntry(2*neAlpha+1, 2*neBeta+1, XForce1/Z);
          };

         if ( NeedMatrix[SCUFF_OMATRIX_YFORCE] )
          { SArray[SCUFF_OMATRIX_YFORCE]->SetEntry(2*neAlpha+0, 2*neBeta+0, Z*YForce1);
            SArray[SCUFF_OMATRIX_YFORCE]->SetEntry(2*neAlpha+0, 2*neBeta+1, YForce2);
            SArray[SCUFF_OMATRIX_YFORCE]->SetEntry(2*neAlpha+1, 2*neBeta+0, Sign*YForce2);
            SArray[SCUFF_OMATRIX_YFORCE]->SetEntry(2*neAlpha+1, 2*neBeta+1, YForce1/Z);
          };

         if ( NeedMatrix[SCUFF_OMATRIX_ZFORCE] )
          { SArray[SCUFF_OMATRIX_ZFORCE]->SetEntry(2*neAlpha+0, 2*neBeta+0, Z*ZForce1);
            SArray[SCUFF_OMATRIX_ZFORCE]->SetEntry(2*neAlpha+0, 2*neBeta+1, ZForce2);
            SArray[SCUFF_OMATRIX_ZFORCE]->SetEntry(2*neAlpha+1, 2*neBeta+0, Sign*ZForce2);
            SArray[SCUFF_OMATRIX_ZFORCE]->SetEntry(2*neAlpha+1, 2*neBeta+1, ZForce1/Z);
          };

       }; // if (IsPEC) ... else ... 
         
   }; // for(neAlpha...) ... for (neBeta...)

  for(int n=0; n<SCUFF_NUM_OMATRICES; n++)
   if ( NeedMatrix[n] ) 
    SArray[n]->EndAssembly();

}

/***************************************************************/
/* get power, force, and torque on a surface                   */
/***************************************************************/
void RWGGeometry::GetPFT(HVector *KN, HVector *RHS, cdouble Omega,
                         int SurfaceIndex, double PFT[8])
{
  if (SurfaceIndex<0 || SurfaceIndex>=NumSurfaces)
   { Warn("invalid surface index passed to GetPFT",SurfaceIndex);
     memset(PFT, 0, 8*sizeof(double));
     return;
   };

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
     PFT[0] += 0.5*real( conj(KAlpha)*NBeta * Overlaps[OVERLAP_CROSS] );

     /*- total power */
     if (neAlpha==neBeta)
      PFT[1] += 0.5*real( conj(KAlpha)*vEAlpha + conj(NAlpha)*vHAlpha );

     /*- X force */
     M11 = Z*(Overlaps[OVERLAP_XBULLET] - Overlaps[OVERLAP_XNABLANABLA]/K2); 
     M12 = +2.0*Overlaps[OVERLAP_XTIMESNABLA]/ (II*Omega);
     M21 = -2.0*Overlaps[OVERLAP_XTIMESNABLA]/ (II*Omega);
     M22 = (Overlaps[OVERLAP_XBULLET] - Overlaps[OVERLAP_XNABLANABLA]/K2) / Z;
     PFT[2] += 0.25*real(   conj(KAlpha)*M11*KBeta 
                          + conj(KAlpha)*M12*NBeta
                          + conj(NAlpha)*M21*KBeta 
                          + conj(NAlpha)*M22*NBeta );

     /*- Y force */
     M11 = Z*(Overlaps[OVERLAP_YBULLET] - Overlaps[OVERLAP_YNABLANABLA]/K2); 
     M12 = +2.0*Overlaps[OVERLAP_YTIMESNABLA]/ (II*Omega);
     M21 = -2.0*Overlaps[OVERLAP_YTIMESNABLA]/ (II*Omega);
     M22 = (Overlaps[OVERLAP_YBULLET] - Overlaps[OVERLAP_YNABLANABLA]/K2) / Z;
     PFT[3] += 0.25*real(   conj(KAlpha)*M11*KBeta 
                          + conj(KAlpha)*M12*NBeta
                          + conj(NAlpha)*M21*KBeta 
                          + conj(NAlpha)*M22*NBeta );

     /*- Z force */
     M11 = Z*(Overlaps[OVERLAP_ZBULLET] - Overlaps[OVERLAP_ZNABLANABLA]/K2); 
     M12 = +2.0*Overlaps[OVERLAP_ZTIMESNABLA]/ (II*Omega);
     M21 = -2.0*Overlaps[OVERLAP_ZTIMESNABLA]/ (II*Omega);
     M22 = (Overlaps[OVERLAP_ZBULLET] - Overlaps[OVERLAP_ZNABLANABLA]/K2) / Z;
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
