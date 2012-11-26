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
/* BPF = 'bloch phase factor' **********************************/
/***************************************************************/
void StampInNeighborBlock(HMatrix *M, HMatrix *B, int RowOffset, int ColOffset,
                          int NBFA, int NBFB, cdouble BPF)
{ 
  for(int nr=0; nr<NBFA; nr++)
   for(int nc=0; nc<NBFB; nc++)
    M->AddEntry(RowOffset + nr, ColOffset + nc,
                BPF*B->GetEntry(nr, nc) + conj(BPF)*B->GetEntry(nc,nr));
}

/***************************************************************/
/* initialize GBarAB9 interpolation tables for all regions at  */
/* the present frequency and bloch wavevector                  */
/***************************************************************/
void RWGGeometry::UpdateRegionInterpolators(cdouble Omega, double *kBloch)
{
  /*--------------------------------------------------------------*/
  /*- check to see if the region interpolators are already        */
  /*- initialized for the given (Omega,kBloch) value              */
  /*- note: instead of static variables within this method, we    */
  /*-       should really use class variables for these FIXME     */
  /*--------------------------------------------------------------*/
  static cdouble LastOmega;
  static double LastkBloch[2];

  if ( LastOmega==Omega && kBloch[0]==LastkBloch[0] && kBloch[1]==LastkBloch[1] )
   return;
  LastOmega=Omega;
  LastkBloch[0]=kBloch[0];
  LastkBloch[1]=kBloch[1];

  UpdateCachedEpsMuValues(Omega);

  GBarData MyGBarData, *GBD=&MyGBarData;
  GBD->kBloch = kBloch;
  GBD->ExcludeInner9=true;
  GBD->E=-1.0;
  GBD->LBV[0]=LatticeBasisVectors[0];
  GBD->LBV[1]=LatticeBasisVectors[1];
  for(int nr=0; nr<NumRegions; nr++)
   {
     if (GBarAB9Interpolators[nr]==0) 
      continue;
     Log("  Initializing interpolator for region %i (%s)...",nr,RegionLabels[nr]);
     GBD->k = csqrt2(EpsTF[nr]*MuTF[nr])*Omega;
     if ( GBD->k == 0.0 ) 
      continue;
     GBarAB9Interpolators[nr]->ReInitialize(GBarVDPhi3D, (void *)GBD);
   };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RWGGeometry::AssembleBEMMatrixBlock(int nsa, int nsb,
                                         int RowOffset, int ColOffset,
                                         cdouble Omega, double *kBloch,
                                         HMatrix *M, HMatrix **GradM)
{
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
  Args->Symmetric = (nsa==nsb) ? 1 : 0;
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

  Args->Symmetric=0;
  Args->RowOffset=Args->ColOffset=0;
  Args->B = B;
  Args->GradB = GradB;
  Args->UseAB9Kernel = false;
  Args->Accumulate = false;

  int NumCommonRegions, CommonRegionIndices[2], nr1, nr2;
  double Signs[2];
  NumCommonRegions=CountCommonRegions(Args->Sa, Args->Sb, CommonRegionIndices, Signs);
  if (NumCommonRegions==0) 
   return;
  nr1=CommonRegionIndices[0];
  nr2=NumCommonRegions==2 ? CommonRegionIndices[1] : -1;

  double Displacement[3];
  Displacement[2]=0.0;
  Args->Displacement=Displacement;

  cdouble BPF; // bloch phase factor

  int NBFA=Args->Sa->NumBFs, NBFB=Args->Sb->NumBFs;
  double *LBV[2];
  LBV[0] = LatticeBasisVectors[0];
  LBV[1] = LatticeBasisVectors[1];

  /***************************************************************/
  /* STEP 2: compute the interaction of surface #nsa with the    */
  /* images of surface #nsb in the neighboring lattice cells.    */
  /***************************************************************/
  if (NumLatticeBasisVectors>=1)
   { 
     Log("MPZ block...");
     Displacement[0]=LatticeBasisVectors[0][0];
     Displacement[1]=LatticeBasisVectors[0][1];
     Args->OmitRegion1 = !RegionIsExtended[MAXLATTICE*nr1+0];
     Args->OmitRegion2 = (nr2>-1) && (!RegionIsExtended[MAXLATTICE*nr2+0]);
     GetSurfaceSurfaceInteractions(Args);
     BPF=exp( II*(kBloch[0]*LBV[0][0] + kBloch[1]*LBV[0][1] ) );
     StampInNeighborBlock(M, B, RowOffset, ColOffset, NBFA, NBFB, BPF);
     if (GradB[0]) 
      StampInNeighborBlock(GradM[0], GradB[0], RowOffset, ColOffset, NBFA, NBFB, BPF);
     if (GradB[1]) 
      StampInNeighborBlock(GradM[1], GradB[1], RowOffset, ColOffset, NBFA, NBFB, BPF);
     if (GradB[2]) 
      StampInNeighborBlock(GradM[2], GradB[2], RowOffset, ColOffset, NBFA, NBFB, BPF);
   }

  if (NumLatticeBasisVectors==2)
   { 
     Log("MPP block...");
     Displacement[0]=LatticeBasisVectors[0][0] + LatticeBasisVectors[1][0];
     Displacement[1]=LatticeBasisVectors[0][1] + LatticeBasisVectors[1][1];
     Args->OmitRegion1 = !RegionIsExtended[MAXLATTICE*nr1+0] || !RegionIsExtended[MAXLATTICE*nr1+1];
     Args->OmitRegion2 = nr2>-1 && (!RegionIsExtended[MAXLATTICE*nr2+0] || !RegionIsExtended[MAXLATTICE*nr2+1]);
     GetSurfaceSurfaceInteractions(Args);
     BPF=exp( II*( kBloch[0]*(LBV[0][0]+LBV[1][0]) + kBloch[1]*(LBV[0][1]+LBV[1][1]) ) );
     StampInNeighborBlock(M, B, RowOffset, ColOffset, NBFA, NBFB, BPF);
     if (GradB[0]) 
      StampInNeighborBlock(GradM[0], GradB[0], RowOffset, ColOffset, NBFA, NBFB, BPF);
     if (GradB[1]) 
      StampInNeighborBlock(GradM[1], GradB[1], RowOffset, ColOffset, NBFA, NBFB, BPF);
     if (GradB[2]) 
      StampInNeighborBlock(GradM[2], GradB[2], RowOffset, ColOffset, NBFA, NBFB, BPF);

     Log("MPM block...");
     Displacement[0]=LatticeBasisVectors[0][0] - LatticeBasisVectors[1][0];
     Displacement[1]=LatticeBasisVectors[0][1] - LatticeBasisVectors[1][1];
     Args->OmitRegion1 = !RegionIsExtended[MAXLATTICE*nr1+0] || !RegionIsExtended[MAXLATTICE*nr1+1];
     Args->OmitRegion2 = nr2>-1 && (!RegionIsExtended[MAXLATTICE*nr2+0] || !RegionIsExtended[MAXLATTICE*nr2+1]);
     GetSurfaceSurfaceInteractions(Args);
     BPF=exp( II*( kBloch[0]*(LBV[0][0]-LBV[1][0]) + kBloch[1]*(LBV[0][1]-LBV[1][1]) ) );
     StampInNeighborBlock(M, B, RowOffset, ColOffset, NBFA, NBFB, BPF);
     if (GradB[0]) 
      StampInNeighborBlock(GradM[0], GradB[0], RowOffset, ColOffset, NBFA, NBFB, BPF);
     if (GradB[1]) 
      StampInNeighborBlock(GradM[1], GradB[1], RowOffset, ColOffset, NBFA, NBFB, BPF);
     if (GradB[2]) 
      StampInNeighborBlock(GradM[2], GradB[2], RowOffset, ColOffset, NBFA, NBFB, BPF);

     Log("MZP block...");
     Displacement[0]=LatticeBasisVectors[1][0];
     Displacement[1]=LatticeBasisVectors[1][1];
     Args->OmitRegion1 = !RegionIsExtended[MAXLATTICE*nr1+0];
     Args->OmitRegion2 = nr2>-1 && !RegionIsExtended[MAXLATTICE*nr2+0];
     GetSurfaceSurfaceInteractions(Args);
     BPF=exp( II*(kBloch[0]*LBV[1][0] + kBloch[1]*LBV[1][1]) ) ;
     StampInNeighborBlock(M, B, RowOffset, ColOffset, NBFA, NBFB, BPF);
     if (GradB[0]) 
      StampInNeighborBlock(GradM[0], GradB[0], RowOffset, ColOffset, NBFA, NBFB, BPF);
     if (GradB[1]) 
      StampInNeighborBlock(GradM[1], GradB[1], RowOffset, ColOffset, NBFA, NBFB, BPF);
     if (GradB[2]) 
      StampInNeighborBlock(GradM[2], GradB[2], RowOffset, ColOffset, NBFA, NBFB, BPF);
   };

  /***************************************************************/
  /* STEP 3: compute the interaction of surface #nsa with the    */
  /* images of surface #nsb in the outer lattice cells.          */
  /***************************************************************/
  UpdateRegionInterpolators(Omega, kBloch);

  Log("Outer cell contributions...");
  Args->Displacement = 0;
  Args->Symmetric    = (nsa==nsb) ? 1 : 0;
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
  if ( NumLatticeBasisVectors==0 && kBloch!=0 )
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
  for(int ns=0; ns<NumSurfaces; ns++)
   for(int nsp=ns; nsp<NumSurfaces; nsp++)
    { AssembleBEMMatrixBlock(ns, nsp, BFIndexOffset[ns], BFIndexOffset[nsp],
                             Omega, kBloch, M, 0 );
    };

  /***************************************************************/
  /* if the matrix uses normal (not packed) storage, fill in its */
  /* below-diagonal blocks. note that the BEM matrix is complex  */
  /* symmetric, not hermitian, so the below-diagonals are equal  */
  /* to the above-diagonals, not to their complex conjugates.    */
  /***************************************************************/
  if (M->StorageType==LHM_NORMAL)
   { int nr, nc;
     for(nr=1; nr<TotalBFs; nr++)
      for(nc=0; nc<nr; nc++)
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
