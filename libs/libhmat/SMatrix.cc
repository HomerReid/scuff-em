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
 * SMatrix.cc  -- implementation of SMatrix class 
 *
 * homer reid  -- 5/2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

#include <libhrutil.h>

extern "C" {
 #include "lapack.h" 
}

#include "libhmat.h"

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/* implementation of SMatrix class methods    ------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/***************************************************************/
/* SMatrix constructor *****************************************/
/***************************************************************/
SMatrix::SMatrix(int pNR, int pNC, int pRealComplex)
{
  NR=pNR; NC=pNC;
  
  RealComplex=pRealComplex;

  nnz = nnz_alloc = 0;
  DM = 0; ZM = 0;
  ColIndices = 0;

  RowStart = (int *) mallocEC(sizeof(int) * (NR + 1));
  cur_nr = 0;
  
  ErrMsg=0;
}

/***************************************************************/
/* SMatrix constructor *****************************************/
/***************************************************************/
SMatrix::SMatrix(const char *FileName, char *MatrixName)
{ ReadFromFile(FileName, LHM_HDF5, MatrixName); }

/***************************************************************/
/* read-from-file constructor **********************************/
/***************************************************************/
void SMatrix::ReadFromFile(const char *FileName, int FileType, const char *MatrixName)
{
  NR=NC=nnz=nnz_alloc=RealComplex=cur_nr=0;
  RowStart=ColIndices=0;
  DM=0;
  ZM=0;
  ErrMsg=0;

  if (FileName==0 || MatrixName==0)
   { ErrMsg=strdup("no file/matrix name specified for matrix import");
     if (HMatrix::AbortOnIOError) 
      ErrExit(ErrMsg);
     return;
   };

  if (FileType==LHM_AUTO) 
   FileType=LHM_HDF5;

  if (FileType!=LHM_HDF5)
   { ErrMsg=strdup("only HDF5 input format supported for SMatrix");
     if (HMatrix::AbortOnIOError)
      ErrExit(ErrMsg);
     return;
   }

  ImportFromHDF5(FileName,MatrixName);

  if (ErrMsg && HMatrix::AbortOnIOError)
   ErrExit(ErrMsg);
}

/***************************************************************/
/* SMatrix destructor  *****************************************/
/***************************************************************/
SMatrix::~SMatrix()
{
  Zero();
  free(RowStart);
}


/***************************************************************/
/***************************************************************/
/***************************************************************/
void SMatrix::Zero()
{
  nnz = nnz_alloc = 0;
  free(DM); DM = 0;
  free(ZM); ZM = 0;
  memset(RowStart, 0, sizeof(int) * (NR + 1));
  free(ColIndices); ColIndices = 0;
  cur_nr = 0;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int SMatrix::GetNNZ(int nr)
{ 
  if (nr<0 || nr>=NR) return 0;
  return RowStart[nr+1] - RowStart[nr];
} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
bool SMatrix::isNZ(int nr, int nc)
{ 
  if (nr<0 || nr>=NR) return false;
  for(int n=RowStart[nr]; n<RowStart[nr+1]; n++)
   if (ColIndices[n]==nc)
    return true;
  return false;
}


/***************************************************************/
/* return pointers to an integer array and a double or cdouble */
/* array containing the column indices and the values of the   */
/* nonzero entries in row nr. the return value is the number   */
/* of nonzero entries.                                         */
/***************************************************************/
int SMatrix::GetRow(int nr, int **pColIndices, void **pNZEntries)
{ 
  if (nr<0 || nr>=NR) return 0;

  int i = RowStart[nr], nrow = RowStart[nr+1] - i;
  *pColIndices = ColIndices + i;
  if (RealComplex==LHM_REAL)
   *pNZEntries = DM + i;
  else
   *pNZEntries = ZM + i;
  return nrow;
} 

/***************************************************************/
/* return all nonzero entries in a given row. the return value */
/* is the number of nonzero entries; NZEntries must point to   */
/* an array of the correct type (real or complex) and of size  */
/* >= GetNNZ(nr).                                              */
/***************************************************************/
int SMatrix::GetRowEntries(int nr, int *pColIndices, void *NZEntries)
{ 
  if (nr<0 || nr>=NR) return 0;

  int i = RowStart[nr], nrow = RowStart[nr+1] - i;
  memcpy(pColIndices, ColIndices + i, nrow * sizeof(int));
  if (RealComplex==LHM_REAL)
    memcpy(NZEntries, DM + i, nrow * sizeof(double));
  else
    memcpy(NZEntries, ZM + i, nrow * sizeof(cdouble));
  return nrow;
} 

/***************************************************************/
/* return the given entry of the matrix                        */
/***************************************************************/
cdouble SMatrix::GetEntry(int nr, int nc)
{ 
  if (nr<0 || nr>=NR || nc<0 || nc>=NC) return 0.0;

  int iend = RowStart[nr+1];
  cdouble sum = 0.0;
  for (int i = RowStart[nr]; i < iend; ++i)
    if (ColIndices[i] == nc) {
      if (RealComplex==LHM_REAL)
	sum += DM[i];
      else
	sum += ZM[i];
    }
  return sum;
} 

/***************************************************************/
/***************************************************************/
/***************************************************************/

void SMatrix::Reallocate(int pnnz_alloc)
{
  nnz_alloc = pnnz_alloc;
  if (nnz_alloc < nnz) ErrExit("SMatrix: invalid reallocation size");
  ColIndices = (int*) reallocEC(ColIndices, sizeof(int) * nnz_alloc);
  if (RealComplex==LHM_REAL)
    DM = (double*) reallocEC(DM, sizeof(double) * nnz_alloc);
  else
    ZM = (cdouble*) reallocEC(ZM, sizeof(cdouble) * nnz_alloc);
}

// Begin matrix assembly.  As an optimization, we allow the user
// to pass an estimate of the #nonzero entries, so that we can
// pre-allocate the data structure.  If this is an underestimate
// it is okay...we grow the arrays (geometrically) as needed.
void SMatrix::BeginAssembly(int est_nnz)
{
  Zero();
  Reallocate(est_nnz < 0 ? 0 : est_nnz);
}

void SMatrix::EndAssembly()
{
  for (; cur_nr < NR-1; ++cur_nr) // insert rowstarts for any remaining rows
    RowStart[cur_nr + 2] = RowStart[cur_nr + 1];
  Reallocate(nnz); // shrink arrays to fit
}

/***************************************************************/
/***************************************************************/
/***************************************************************/

// Make an entry (if it does not exist yet) in the CSR data,
// and return the index to the data.   Can only be called during
// matrix assembly, and is fastest if done in increasing order by row.
// If force_new is true, always creates a new entry, rather
// than searching for an existing one.
int SMatrix::MakeEntry(int nr, int nc, bool force_new)
{
  if (nr<0 || nr>=NR || nc<0 || nc>=NC) ErrExit("SMatrix: invalid MakeEntry");

  // update cur_nr as needed, & create any empty rows between cur_nr and nr
  for (; cur_nr < nr; ++cur_nr)
    RowStart[cur_nr + 2] = RowStart[cur_nr + 1];

  int iend = RowStart[nr+1];

  // NOTE: searching for existing entries is slow if we have a dense row
  //       ... it makes adding a dense row cost O(NR^2)!
  //       ... hence dense rows are best added with AddEntry(..., false)
  if (!force_new)
    for (int i = RowStart[nr]; i < iend; ++i)
      if (ColIndices[i] == nc) return i; // entry already exists

  // add new entry (with zero data)
  RowStart[nr+1]++;
  if (++nnz > nnz_alloc) 
    Reallocate((nnz_alloc + 1) * 2); // allocate geometrically

  if (nr < cur_nr) { // need to insert new entry; slow if many entries to move!
    for (int i = nr+1; i <= cur_nr; ++i)
      RowStart[i+1]++;
    int nmove = nnz - (iend + 1);
    memmove(ColIndices + iend + 1, ColIndices + iend, sizeof(int) * nmove);
    if (RealComplex==LHM_REAL)
      memmove(DM + iend + 1, DM + iend, sizeof(double) * nmove);
    else
      memmove(ZM + iend + 1, ZM + iend, sizeof(cdouble) * nmove);
  }

  ColIndices[iend] = nc;
  if (RealComplex==LHM_REAL)
    DM[iend] = 0.0;
  else
    ZM[iend] = 0.0;
  return iend;
}

void SMatrix::SetEntry(int nr, int nc, cdouble Entry)
{ 
  int i = MakeEntry(nr, nc, false);
  if (RealComplex==LHM_REAL)
    DM[i] = real(Entry);
  else
    ZM[i] = Entry;
}

// Add to an existing entry.  If compress is false, we just add
// a new sparse-matrix entry without checking if there is an existing
// coefficient that we could just update ... this improves speed
// at the cost of possibly wasting space.  If we are adding a mostly
// dense row, it will be much faster (O(NR) rather than O(NR^2)) to
// call AddEntry with compress=false for each entry in the row!
void SMatrix::AddEntry(int nr, int nc, cdouble Entry, bool compress)
{
  int i = MakeEntry(nr, nc, !compress);
  if (RealComplex==LHM_REAL)
    DM[i] += real(Entry);
  else
    ZM[i] += Entry;
} 

/***************************************************************/
/* return value is <X|this|X>  *********************************/
/***************************************************************/
cdouble SMatrix::Apply(HVector *X, HVector *MX)
{ 
  if ( MX->N != NR || X->N != NC)
   ErrExit("size mismatch in SMatrix::Apply");

  cdouble XMX=0.0;
  if (RealComplex == LHM_REAL) { 
     for (int nr = 0; nr < NR; ++nr) {
       int iend = RowStart[nr+1];
       cdouble sum = 0.0;
       for (int i = RowStart[nr]; i < iend; ++i)
	 sum += DM[i] * X->GetEntry(ColIndices[i]);
       MX->SetEntry(nr, sum);
       XMX+=X->GetEntry(nr)*sum;
     }
   }
  else { 
     for (int nr = 0; nr < NR; ++nr) {
       int iend = RowStart[nr+1];
       cdouble sum = 0.0;
       for (int i = RowStart[nr]; i < iend; ++i)
	 sum += ZM[i] * X->GetEntry(ColIndices[i]);
       MX->SetEntry(nr, sum);
       XMX+=conj(X->GetEntry(nr))*sum;
     }
   }
  return XMX;
}

double SMatrix::ApplyD(HVector *X, HVector *MX)
 { return real(Apply(X,MX)); }

HVector *SMatrix::Apply(HVector *X)
{
  HVector *MX = new HVector(X->N,
			    X->RealComplex == LHM_REAL
			    && RealComplex == LHM_REAL ? LHM_REAL
			    : LHM_COMPLEX);
  Apply(X, MX);
  return MX;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SMatrix::Apply(HMatrix *X, HMatrix *MX)
{ 
  if ( (MX->NR != NR) || (X->NR != NC) || (MX->NC != X->NC) )
   ErrExit("size mismatch in SMatrix::Apply");

  if (RealComplex == LHM_REAL) { 
     for (int nc = 0; nc < X->NC; nc++) for (int nr = 0; nr < NR; ++nr) {
       int iend = RowStart[nr+1];
       cdouble sum = 0.0;
       for (int i = RowStart[nr]; i < iend; ++i)
	 sum += DM[i] * X->GetEntry(ColIndices[i], nc);
       MX->SetEntry(nr, nc, sum);
     }
   }
  else { 
     for (int nc = 0; nc < X->NC; nc++) for (int nr = 0; nr < NR; ++nr) {
       int iend = RowStart[nr+1];
       cdouble sum = 0.0;
       for (int i = RowStart[nr]; i < iend; ++i)
	 sum += ZM[i] * X->GetEntry(ColIndices[i], nc);
       MX->SetEntry(nr, nc, sum);
     }
   }
}

HMatrix *SMatrix::Apply(HMatrix *X)
{
  HMatrix *MX = new HMatrix(NR, X->NC,
			    X->RealComplex == LHM_REAL
			    && RealComplex == LHM_REAL ? LHM_REAL
			    : LHM_COMPLEX);
  Apply(X, MX);
  return MX;
}

/***************************************************************/
/* compute a vector-matrix-vector product, X' * this * Y       */
/***************************************************************/
cdouble SMatrix::BilinearProduct(HVector *X, HVector *Y)
{
  if (Y==0) Y=X;

  if ( (X->N != NR) || (Y->N != NC) )
   ErrExit("size mismatch in SMatrix::BilinearProduct: %i x (%ix%i) x %i",X->N,NR,NC,Y->N);

  cdouble Result=0.0;

  for (int nr = 0; nr < NR; ++nr) 
   { cdouble XC=conj(X->GetEntry(nr));
     int iend = RowStart[nr+1];
     for (int i = RowStart[nr]; i < iend; ++i)
      Result += XC * (RealComplex==LHM_REAL? DM[i] : ZM[i]) * Y->GetEntry(ColIndices[i]);
   };
  return Result;
}

double SMatrix::BilinearProductD(HVector *X, HVector *Y)
 { return real(BilinearProduct(X,Y)); }
