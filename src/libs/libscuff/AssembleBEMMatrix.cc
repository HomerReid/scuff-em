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
/* how kludgy am i? ********************************************/
/***************************************************************/
int ZeroRegionsAsNecessary(RWGGeometry *G, RWGSurface *Sa, RWGSurface *Sb, 
                           int CommonRegionIndices[2], int SaveZeroed[2])
                            
                            
{
  // see NOTE above
  int NumCommonRegions, CommonRegionIndices[2], nr1, nr2;
  bool Region1Contributes, Region2Contributes;
  int SaveZeroed1, SaveZeroed2;
  double Signs[2];
  NumCommonRegions=CountCommonRegions(Sa, Sb, CommonRegionIndices, Signs);
  if (NumCommonRegions==0)
   return true;

  nr1=CommonRegionIndices[0]; 
  Region1Contributes = (RegionIsExtended[MAXLATTICE*nr1+0] && RegionIsExtended[MAXLATTICE*nr1+1]);
  if (NumCommonRegions==2)
   { nr2=CommonRegionIndices[1];
     Region2Contributes = (RegionIsExtended[MAXLATTICE*nr2+0] && RegionIsExtended[MAXLATTICE*nr2+1]);
   }
  else
   Region2Contributes=false;

  if ( Region1Contributes==false  && Region2Contributes==false )
   return true;

  if ( Region1Contributes==false )
   { SaveZeroed1=RegionMPs[nr1]->Zeroed; 
     RegionMPs[nr1]->Zero();
   };
  if ( NumCommonRegions==2 && Region2Contributes==false )
   { SaveZeroed2=RegionMPs[nr2]->Zeroed; 
     RegionMPs[nr2]->Zero();
   };
  ZeroRegionsAsNecessary();

}

void UnZeroRegionsAsNecessary(RWGGeometry *G, int NumCommonRegions, 
                              int CommonRegionIndices[2], int SaveZeroed[2])
{ 

  if ( Region1Contributes==false )
   RegionMPs[nr1]->Zeroed=SaveZeroed1;
  if ( NumCommonRegions==2 && Region2Contributes==false )
   RegionMPs[nr2]->Zeroed=SaveZeroed2;

 // 20120924: the zeroing/unzeroing of MPs may have wreaked havoc 
 //           on the internally-stored Eps/Mu arrays, so i will 
 //           quickly restore them here. 
 //           it must be said that this zeroing/unzeroing          
 //           paradigm is fairly kludgy, and it would be nice to 
 //           redo it more elegantly at some point.
 UpdateCachedEpsMuValues(StoredOmega);
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

  Args->dBdTheta=0; // FIXME angular derivatives not implemented yet
  Args->NumTorqueAxes=0; // angular derivatives not implemented yet
  Args->GammaMatrix=0; // angular derivatives not implemented yet

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

  /***************************************************************/
  /* prepare for the next set of calls to GetSurfaceSurfaceInteractions*/
  /***************************************************************/
  Args->Symmetric 
  Args->RowOffset=Args->ColOffset=0;
  Args->B = ScratchBlock;
  Args->GradB = GradScratchBlock;
  
  //ZeroRegionsAsNecessary(G, Sa, Sb);

  /***************************************************************/
  /* STEP 2: compute the interaction of surface #nsa with the    */
  /* images of surface #nsb in the neighboring lattice cells.    */
  /***************************************************************/
  double Displacement[3];
  Displacement[3]=0.0;
  if (NumLatticeBasisVectors>=1)
   { Log("MPZ block...");
     Displacement[0]=LatticeBasisVectors[0][0];
     Displacement[1]=LatticeBasisVectors[0][1];
     Args->Symmetric=0;
   };

  //UnZeroRegionsAsNecessary(G, Sa, Sb);

  /***************************************************************/
  /* STEP 3: compute the interaction of surface #nsa with the    */
  /* images of surface #nsb in the outer lattice cells.          */
  /***************************************************************/
  
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
