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
 * AssembleBEMMatrix.cc -- libscuff routines for assembling the BEM matrix.
 *                      --
 *                      -- (cf. 'libscuff Implementation and Technical
 *                      --  Details', section 8.3, 'Structure of the BEM
 *                      --  Matrix.')
 *                      --
 * homer reid           -- 10/2006 -- 10/2011
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libhmat.h>
#include <libhrutil.h>

#include "libscuff.h"
#include "libscuffInternals.h"

#include "cmatheval.h"

#ifdef USE_PTHREAD
#  include <pthread.h>
#endif
#ifdef USE_OPENMP
#  include <omp.h>
#endif

#define II cdouble(0,1)

namespace scuff {

/***************************************************************/
/* This routine adds the contents of matrix block B, which has */
/* dimensions NRxNC, to the block of M whose upper-left corner */
/* has indices (RowOffset, ColOffset).                         */
/*                                                             */
/* Each entry of B is scaled by BPF ('bloch phase factor')     */
/* before being added to the corresponding entry of M.         */
/*                                                             */
/* If SameSurface is true, then the routine additionally adds  */
/* the contents of B' times the complex conjugate of BPF to    */
/* the destination block of M.                                 */
/***************************************************************/
void StampInNeighborBlock(HMatrix *B, int NR, int NC, 
                          HMatrix *M, int RowOffset, int ColOffset,
                          cdouble BPF, bool SameSurface)
{ 
  if (SameSurface)
   { for(int nr=0; nr<NR; nr++)
      for(int nc=0; nc<NC; nc++)
       M->AddEntry(RowOffset + nr, ColOffset + nc,
                   BPF*B->GetEntry(nr,nc) + conj(BPF)*B->GetEntry(nc,nr));
   }
  else
   { for(int nr=0; nr<NR; nr++)
      for(int nc=0; nc<NC; nc++)
       M->AddEntry(RowOffset + nr, ColOffset + nc, BPF*B->GetEntry(nr,nc));
   };
}

/***************************************************************/
/* initialize GBarAB9 interpolation tables for all regions at  */
/* the present frequency and bloch wavevector                  */
/***************************************************************/
void RWGGeometry::UpdateRegionInterpolators(cdouble Omega, double *kBloch)
{
  UpdateCachedEpsMuValues(Omega);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  GBarData MyGBarData, *GBD=&MyGBarData;
  GBD->ExcludeInner9=true;
  GBD->E=-1.0;
  GBD->kBloch = kBloch;
  for(int nr=0; nr<NumRegions; nr++)
   {
     if (GBarAB9Interpolators[nr]==0) 
      continue;

     /***************************************************************/
     /* figure out whether this region has 1D or 2D periodicity     */
     /***************************************************************/  
     if ( RegionIsExtended[MAXLATTICE*nr+0] && RegionIsExtended[MAXLATTICE*nr+1] ) 
      { GBD->LDim=2;
        GBD->LBV[0]=LatticeBasisVectors[0];
        GBD->LBV[1]=LatticeBasisVectors[1];
      }
     else if ( RegionIsExtended[MAXLATTICE*nr+0] && !RegionIsExtended[MAXLATTICE*nr+1] )
      { GBD->LDim=1; 
        GBD->LBV[0]=LatticeBasisVectors[0];
      }
     else if ( !RegionIsExtended[MAXLATTICE*nr+0] && RegionIsExtended[MAXLATTICE*nr+1] )
      { GBD->LDim=1; 
        GBD->LBV[0]=LatticeBasisVectors[1];
      }
     else
      continue; // region is compact; no interpolation table needed

     /***************************************************************/
     /* check whether or not the extents of the region have changed */
     /* since we last allocated an interpolator for this region, and*/
     /* re-allocate the interpolator if so                          */
     /***************************************************************/
     double RMax[3], RMin[3], DeltaR[3];
     int NPoints[3];
     GetRegionExtents(nr, RMax, RMin, DeltaR, NPoints);
     if (    GBarAB9Interpolators[nr]==0
          || GBarAB9Interpolators[nr]->X1Min !=  -DeltaR[0]
          || GBarAB9Interpolators[nr]->X2Min !=  -DeltaR[1]
          || GBarAB9Interpolators[nr]->N1    !=  NPoints[0]
          || GBarAB9Interpolators[nr]->N2    !=  NPoints[1]
          || GBarAB9Interpolators[nr]->N3    != (1+NPoints[2]/2)
        )
      {
        Log("Region %s extents have changed (resizing interpolation table)",RegionLabels[nr]);
        if (GBarAB9Interpolators[nr])
         { delete GBarAB9Interpolators[nr];
           GBarAB9Interpolators[nr]=0;
         };
        GBarAB9Interpolators[nr]=new Interp3D( -DeltaR[0], DeltaR[0], NPoints[0],
                                               -DeltaR[1], DeltaR[1], NPoints[1],
                                                      0.0, DeltaR[2], 1 + NPoints[2]/2, 
                                                       2);
        
      };

     /***************************************************************/
     /* initialize the interpolation table                          */
     /***************************************************************/
     Log("  Initializing interpolator for region %i (%s)...",nr,RegionLabels[nr]);
     GBD->k = csqrt2(EpsTF[nr]*MuTF[nr])*Omega;
     if ( GBD->k == 0.0 ) 
      continue;
     GBarAB9Interpolators[nr]->ReInitialize(GBarVDPhi3D, (void *)GBD);

   };
}

/***************************************************************/
/* This routine computes the block of the BEM matrix that      */
/* describes the interaction between surfaces nsa and nsb.     */
/* This block is stamped into M in such a way that the upper-  */ 
/* left element of the block is at the (RowOffset, ColOffset)  */ 
/* entry of M. If GradM is non-null and GradM[Mu] is non-null  */ 
/* (Mu=0,1,2) then the X_{Mu} derivative of BEM matrix is      */ 
/* similarly stamped into GradM[Mu].                           */ 
/***************************************************************/
void RWGGeometry::AssembleBEMMatrixBlock(int nsa, int nsb,
                                         cdouble Omega, double *kBloch,
                                         HMatrix *M, HMatrix **GradM,
                                         int RowOffset, int ColOffset)
{
  bool SameSurface = (nsa==nsb);

  /***************************************************************/
  /* pre-initialize arguments for GetSurfaceSurfaceInteractions **/
  /***************************************************************/
  GetSSIArgStruct GetSSIArgs, *Args=&GetSSIArgs;
  InitGetSSIArgs(Args);

  Args->G=this;
  Args->Sa=Surfaces[nsa];
  Args->Sb=Surfaces[nsb];
  Args->Omega=Omega;

  Args->dBdTheta=0;      // FIXME angular derivatives not implemented yet
  Args->NumTorqueAxes=0;
  Args->GammaMatrix=0;

  /***************************************************************/
  /* STEP 1: compute the direct interaction of the two surfaces. */
  /* This computation is always necessary, whether or not a      */
  /* lattice is present.                                         */
  /***************************************************************/
  Args->Displacement=0;
  Args->UseAB9Kernel=false;
  Args->Symmetric = SameSurface;
  Args->Accumulate = 0;
  Args->B=M;
  Args->GradB=GradM;
  Args->RowOffset=RowOffset;
  Args->ColOffset=ColOffset;

  GetSurfaceSurfaceInteractions(Args);

  if (NumLatticeBasisVectors==0) 
   return; 

  /*********************************************************************/
  /* STEP 2: prepare to make a series of calls to                      */
  /* GetSurfaceSurfaceInteractions to get the contributions of the     */
  /* innermost 8 neighboring lattice cells                             */
  /*********************************************************************/

  // B and GradB are statically maintained matrix blocks used as 
  // temporary storage within this routine; they need to be large
  // enough to store the largest single subblock of the BEM matrix.
  // FIXME these should be fields of the RWGGeometry class, not 
  // static method variables.
  static HMatrix *B=0, *GradB[3]={0,0,0};
  if (B==0)
   { int MaxNR=Surfaces[0]->NumBFs; 
     for(int ns=1; ns<NumSurfaces; ns++)
      if (Surfaces[ns]->NumBFs > MaxNR) 
       MaxNR=Surfaces[ns]->NumBFs;
     B=new HMatrix(MaxNR, MaxNR, LHM_COMPLEX);
     if ( GradM && GradM[0] ) GradB[0]=new HMatrix(MaxNR, MaxNR, LHM_COMPLEX);
     if ( GradM && GradM[1] ) GradB[1]=new HMatrix(MaxNR, MaxNR, LHM_COMPLEX);
     if ( GradM && GradM[2] ) GradB[2]=new HMatrix(MaxNR, MaxNR, LHM_COMPLEX);
   };

  int NumCommonRegions, CommonRegionIndices[2], nr1, nr2;
  double Signs[2];
  NumCommonRegions=CountCommonRegions(Args->Sa, Args->Sb, CommonRegionIndices, Signs);
  if (NumCommonRegions==0) 
   return;
  nr1=CommonRegionIndices[0];
  nr2=NumCommonRegions==2 ? CommonRegionIndices[1] : -1;

  // L is the lattice vector through which surfaces are displaced into
  // neighboring unit cells
  double L[3];
  L[2]=0.0;

  cdouble BPF; // bloch phase factor

  int NBFA=Args->Sa->NumBFs; 
  int NBFB=Args->Sb->NumBFs;

  Args->Symmetric=0;
  Args->RowOffset=0;
  Args->ColOffset=0;
  Args->B = B;
  Args->GradB = GradB;
  Args->UseAB9Kernel = false;
  Args->Accumulate = false;
  Args->Displacement=L;

  /***************************************************************/
  /* STEP 2: compute the interaction of surface #nsa with the    */
  /* images of surface #nsb in the neighboring lattice cells.    */
  /***************************************************************/
  if (NumLatticeBasisVectors>=1)
   { 
     Log("MPZ block...");
     L[0]=LatticeBasisVectors[0][0];
     L[1]=LatticeBasisVectors[0][1];
     Args->OmitRegion1 = !RegionIsExtended[MAXLATTICE*nr1+0];
     Args->OmitRegion2 = (nr2>-1) && (!RegionIsExtended[MAXLATTICE*nr2+0]);
     GetSurfaceSurfaceInteractions(Args);
     BPF=exp( II*(kBloch[0]*L[0] + kBloch[1]*L[1]) );
     StampInNeighborBlock(B, NBFA, NBFB, M, RowOffset, ColOffset, BPF, SameSurface);
     if (GradB[0]) 
      StampInNeighborBlock(GradB[0], NBFA, NBFB, GradM[0], RowOffset, ColOffset, BPF, SameSurface);
     if (GradB[1]) 
      StampInNeighborBlock(GradB[1], NBFA, NBFB, GradM[1], RowOffset, ColOffset, BPF, SameSurface);
     if (GradB[2]) 
      StampInNeighborBlock(GradB[2], NBFA, NBFB, GradM[2], RowOffset, ColOffset, BPF, SameSurface);

     if (!SameSurface)
      { 
        Log("MMZ block...");
        L[0] = -LatticeBasisVectors[0][0];
        L[1] = -LatticeBasisVectors[0][1];
        Args->OmitRegion1 = !RegionIsExtended[MAXLATTICE*nr1+0];
        Args->OmitRegion2 = (nr2>-1) && (!RegionIsExtended[MAXLATTICE*nr2+0]);
        GetSurfaceSurfaceInteractions(Args);
        BPF=exp( II*(kBloch[0]*L[0] + kBloch[1]*L[1]) );
        StampInNeighborBlock(B, NBFA, NBFB, M, RowOffset, ColOffset, BPF, SameSurface);
        if (GradB[0]) 
         StampInNeighborBlock(GradB[0], NBFA, NBFB, GradM[0], RowOffset, ColOffset, BPF, SameSurface);
        if (GradB[1]) 
         StampInNeighborBlock(GradB[1], NBFA, NBFB, GradM[1], RowOffset, ColOffset, BPF, SameSurface);
        if (GradB[2]) 
         StampInNeighborBlock(GradB[2], NBFA, NBFB, GradM[2], RowOffset, ColOffset, BPF, SameSurface);
      };
   }

  if (NumLatticeBasisVectors==2)
   { 
     Log("MPP block...");
     L[0]=LatticeBasisVectors[0][0] + LatticeBasisVectors[1][0];
     L[1]=LatticeBasisVectors[0][1] + LatticeBasisVectors[1][1];
     Args->OmitRegion1 = !RegionIsExtended[MAXLATTICE*nr1+0] || !RegionIsExtended[MAXLATTICE*nr1+1];
     Args->OmitRegion2 = nr2>-1 && (!RegionIsExtended[MAXLATTICE*nr2+0] || !RegionIsExtended[MAXLATTICE*nr2+1]);
     GetSurfaceSurfaceInteractions(Args);
     BPF=exp( II*(kBloch[0]*L[0] + kBloch[1]*L[1]) );
     StampInNeighborBlock(B, NBFA, NBFB, M, RowOffset, ColOffset, BPF, SameSurface);
     if (GradB[0]) 
      StampInNeighborBlock(GradB[0], NBFA, NBFB, GradM[0], RowOffset, ColOffset, BPF, SameSurface);
     if (GradB[1]) 
      StampInNeighborBlock(GradB[1], NBFA, NBFB, GradM[1], RowOffset, ColOffset, BPF, SameSurface);
     if (GradB[2]) 
      StampInNeighborBlock(GradB[2], NBFA, NBFB, GradM[2], RowOffset, ColOffset, BPF, SameSurface);

     Log("MPM block...");
     L[0]=LatticeBasisVectors[0][0] - LatticeBasisVectors[1][0];
     L[1]=LatticeBasisVectors[0][1] - LatticeBasisVectors[1][1];
     Args->OmitRegion1 = !RegionIsExtended[MAXLATTICE*nr1+0] || !RegionIsExtended[MAXLATTICE*nr1+1];
     Args->OmitRegion2 = nr2>-1 && (!RegionIsExtended[MAXLATTICE*nr2+0] || !RegionIsExtended[MAXLATTICE*nr2+1]);
     GetSurfaceSurfaceInteractions(Args);
     BPF=exp( II*(kBloch[0]*L[0] + kBloch[1]*L[1]) );
     StampInNeighborBlock(B, NBFA, NBFB, M, RowOffset, ColOffset, BPF, SameSurface);
     if (GradB[0]) 
      StampInNeighborBlock(GradB[0], NBFA, NBFB, GradM[0], RowOffset, ColOffset, BPF, SameSurface);
     if (GradB[1]) 
      StampInNeighborBlock(GradB[1], NBFA, NBFB, GradM[1], RowOffset, ColOffset, BPF, SameSurface);
     if (GradB[2]) 
      StampInNeighborBlock(GradB[2], NBFA, NBFB, GradM[2], RowOffset, ColOffset, BPF, SameSurface);

     Log("MZP block...");
     L[0]=LatticeBasisVectors[1][0];
     L[1]=LatticeBasisVectors[1][1];
     Args->OmitRegion1 = !RegionIsExtended[MAXLATTICE*nr1+1];
     Args->OmitRegion2 = nr2>-1 && !RegionIsExtended[MAXLATTICE*nr2+1];
     GetSurfaceSurfaceInteractions(Args);
     BPF=exp( II*(kBloch[0]*L[0] + kBloch[1]*L[1]) );
     StampInNeighborBlock(B, NBFA, NBFB, M, RowOffset, ColOffset, BPF, SameSurface);
     if (GradB[0]) 
      StampInNeighborBlock(GradB[0], NBFA, NBFB, GradM[0], RowOffset, ColOffset, BPF, SameSurface);
     if (GradB[1]) 
      StampInNeighborBlock(GradB[1], NBFA, NBFB, GradM[1], RowOffset, ColOffset, BPF, SameSurface);
     if (GradB[2]) 
      StampInNeighborBlock(GradB[2], NBFA, NBFB, GradM[2], RowOffset, ColOffset, BPF, SameSurface);

     if (!SameSurface)
      {
        Log("MMM block...");
        L[0] = -LatticeBasisVectors[0][0] - LatticeBasisVectors[1][0];
        L[1] = -LatticeBasisVectors[0][1] - LatticeBasisVectors[1][1];
        Args->OmitRegion1 = !RegionIsExtended[MAXLATTICE*nr1+0] || !RegionIsExtended[MAXLATTICE*nr1+1];
        Args->OmitRegion2 = nr2>-1 && (!RegionIsExtended[MAXLATTICE*nr2+0] || !RegionIsExtended[MAXLATTICE*nr2+1]);
        GetSurfaceSurfaceInteractions(Args);
        BPF=exp( II*(kBloch[0]*L[0] + kBloch[1]*L[1]) );
        StampInNeighborBlock(B, NBFA, NBFB, M, RowOffset, ColOffset, BPF, SameSurface);
        if (GradB[0]) 
         StampInNeighborBlock(GradB[0], NBFA, NBFB, GradM[0], RowOffset, ColOffset, BPF, SameSurface);
        if (GradB[1]) 
         StampInNeighborBlock(GradB[1], NBFA, NBFB, GradM[1], RowOffset, ColOffset, BPF, SameSurface);
        if (GradB[2]) 
         StampInNeighborBlock(GradB[2], NBFA, NBFB, GradM[2], RowOffset, ColOffset, BPF, SameSurface);
   
        Log("MMP block...");
        L[0] = -LatticeBasisVectors[0][0] + LatticeBasisVectors[1][0];
        L[1] = -LatticeBasisVectors[0][1] + LatticeBasisVectors[1][1];
        Args->OmitRegion1 = !RegionIsExtended[MAXLATTICE*nr1+0] || !RegionIsExtended[MAXLATTICE*nr1+1];
        Args->OmitRegion2 = nr2>-1 && (!RegionIsExtended[MAXLATTICE*nr2+0] || !RegionIsExtended[MAXLATTICE*nr2+1]);
        GetSurfaceSurfaceInteractions(Args);
        BPF=exp( II*(kBloch[0]*L[0] + kBloch[1]*L[1]) );
        StampInNeighborBlock(B, NBFA, NBFB, M, RowOffset, ColOffset, BPF, SameSurface);
        if (GradB[0]) 
         StampInNeighborBlock(GradB[0], NBFA, NBFB, GradM[0], RowOffset, ColOffset, BPF, SameSurface);
        if (GradB[1]) 
         StampInNeighborBlock(GradB[1], NBFA, NBFB, GradM[1], RowOffset, ColOffset, BPF, SameSurface);
        if (GradB[2]) 
         StampInNeighborBlock(GradB[2], NBFA, NBFB, GradM[2], RowOffset, ColOffset, BPF, SameSurface);
   
        Log("MZM block...");
        L[0] = -LatticeBasisVectors[1][0];
        L[1] = -LatticeBasisVectors[1][1];
        Args->OmitRegion1 = !RegionIsExtended[MAXLATTICE*nr1+1];
        Args->OmitRegion2 = nr2>-1 && !RegionIsExtended[MAXLATTICE*nr2+1];
        GetSurfaceSurfaceInteractions(Args);
        BPF=exp( II*(kBloch[0]*L[0] + kBloch[1]*L[1]) );
        StampInNeighborBlock(B, NBFA, NBFB, M, RowOffset, ColOffset, BPF, SameSurface);
        if (GradB[0]) 
         StampInNeighborBlock(GradB[0], NBFA, NBFB, GradM[0], RowOffset, ColOffset, BPF, SameSurface);
        if (GradB[1]) 
         StampInNeighborBlock(GradB[1], NBFA, NBFB, GradM[1], RowOffset, ColOffset, BPF, SameSurface);
        if (GradB[2]) 
         StampInNeighborBlock(GradB[2], NBFA, NBFB, GradM[2], RowOffset, ColOffset, BPF, SameSurface);
   
      }; // if (!SameSurface) ... 

   };

  /***************************************************************/
  /* STEP 3: compute the interaction of surface #nsa with the    */
  /* images of surface #nsb in the outer lattice cells.          */
  /***************************************************************/
  UpdateRegionInterpolators(Omega, kBloch);

  Log("Outer cell contributions...");
  Args->Displacement = 0;
  Args->Symmetric    = SameSurface;
  Args->OmitRegion1  = false;
  Args->OmitRegion2  = false;
  Args->UseAB9Kernel = true;
  Args->Accumulate   = true;
  Args->B            = M;
  Args->GradB        = GradM;
  Args->RowOffset    = RowOffset;
  Args->ColOffset    = ColOffset;

  GetSurfaceSurfaceInteractions(Args);

}

/***************************************************************/
/* this is the actual API-exposed routine for assembling the   */
/* BEM matrix, which is pretty simple and really just calls    */
/* AssembleBEMMatrixBlock() to do all the dirty work.          */
/*                                                             */
/* If the M matrix is NULL on entry, a new HMatrix of the      */
/* appropriate size is allocated and returned. Otherwise, the  */
/* return value is M.                                          */
/***************************************************************/
HMatrix *RWGGeometry::AssembleBEMMatrix(cdouble Omega, double *kBloch, HMatrix *M)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if ( NumLatticeBasisVectors==0 && kBloch!=0 && (kBloch[0]!=0.0 || kBloch[1]!=0.0) )
   ErrExit("%s:%i: Bloch wavevector is undefined for compact geometries");
  if ( NumLatticeBasisVectors!=0 && kBloch==0 )
   ErrExit("%s:%i: Bloch wavevector must be specified for PBC geometries");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (M==NULL)
   M=AllocateBEMMatrix();
  else if ( M->NR != TotalBFs || M->NC != TotalBFs )
   { Warn("wrong-size matrix passed to AssembleBEMMatrix; reallocating...");
     M=AllocateBEMMatrix();
   };

  /***************************************************************/
  /* loop over all pairs of objects to assemble the diagonal and */
  /* above-diagonal blocks of the matrix                         */
  /***************************************************************/
  /***************************************************************/
  int nsm; // 'number of surface mate'
  for(int ns=0; ns<NumSurfaces; ns++)
   for(int nsp=ns; nsp<NumSurfaces; nsp++)
    { 
      // attempt to reuse the diagonal block of an identical previous object
      if (ns==nsp && (nsm=Mate[ns])!=-1)
       { int ThisOffset = BFIndexOffset[ns];
         int MateOffset = BFIndexOffset[nsm];
         int Dim = Surfaces[ns]->NumBFs;
         M->InsertBlock(M, ThisOffset, ThisOffset, Dim, Dim, MateOffset, MateOffset);
       }
      else
       AssembleBEMMatrixBlock(ns, nsp, Omega, kBloch, M, 0,
                              BFIndexOffset[ns], BFIndexOffset[nsp]);
    };

  /***************************************************************/
  /* if the matrix uses normal (not packed) storage, fill in its */
  /* below-diagonal blocks. note that the BEM matrix is complex  */
  /* symmetric, not hermitian, so the below-diagonals are equal  */
  /* to the above-diagonals, not to their complex conjugates.    */
  /***************************************************************/
  if (M->StorageType==LHM_NORMAL)
   { 
     for(int nr=1; nr<TotalBFs; nr++)
      for(int nc=0; nc<nr; nc++)
       M->SetEntry(nr, nc, M->GetEntry(nc, nr) );
   };

  return M;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *RWGGeometry::AssembleBEMMatrix(cdouble Omega, HMatrix *M)
{
  return AssembleBEMMatrix(Omega, 0, M); 
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *RWGGeometry::AllocateBEMMatrix(bool PureImagFreq, bool Packed)
{
  int Storage = Packed ? LHM_SYMMETRIC : LHM_NORMAL;
  if (PureImagFreq)
    return new HMatrix(TotalBFs, TotalBFs, LHM_REAL, Storage);
  else
    return new HMatrix(TotalBFs, TotalBFs, LHM_COMPLEX, Storage);
    
}

} // namespace scuff
