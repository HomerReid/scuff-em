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
 * HMatrix.cc  -- implementation of HMatrix class methods
 *
 * homer reid  -- 12/2009
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
/* implementation of HMatrix class methods    ------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

// initialization of static class flag
bool HMatrix::AbortOnIOError=true;

/***************************************************************/
/* HMatrix constructors that create an empty (zero) HMatrix    */
/* with known dimensions                                       */
/***************************************************************/
HMatrix::HMatrix(int NRows, int NCols, cdouble *data)
 { InitHMatrix(NRows, NCols, LHM_COMPLEX, LHM_NORMAL, (void *)data); }

HMatrix::HMatrix(int NRows, int NCols, double *data)
 { InitHMatrix(NRows, NCols, LHM_REAL, LHM_NORMAL, (void *)data); }

HMatrix::HMatrix(int NRows, int NCols, int pRealComplex, void *data)
 { InitHMatrix(NRows, NCols, pRealComplex, LHM_NORMAL, data); }

HMatrix::HMatrix(int NRows, int NCols, int pRealComplex, int pStorageType,
		 void *data)
 { InitHMatrix(NRows, NCols, pRealComplex, pStorageType, data); }

HMatrix::HMatrix(HMatrix *M, bool takedatandownership) {
  if (takedatandownership) {
    InitHMatrix(M->NR, M->NC, M->RealComplex, M->StorageType, M->RealComplex==LHM_COMPLEX ? (void*)M->ZM : (void*)M->DM);
    ownsM = true;
    M->ownsM = false;
  } else {
    InitHMatrix(M->NR, M->NC, M->RealComplex, M->StorageType);
    Copy(M);
  }
}

/***************************************************************/
/* constructor body for the above entry points *****************/
/***************************************************************/
void HMatrix::InitHMatrix(int NRows, int NCols, int pRealComplex, 
			  int pStorageType, void *data)
{ 
   RealComplex=pRealComplex;
   StorageType=pStorageType;
   ipiv=0;
   lwork=0;
   work=0;
   liwork=0;
   iwork=0;
   ErrMsg=0;

   NR=NRows;
   NC=NCols;
  
   if ( (StorageType==LHM_SYMMETRIC || StorageType==LHM_HERMITIAN) && NR!=NC )
    { StorageType=LHM_NORMAL;
      Warn("warning: packed storage not available for non-square matrices (using normal storage)");
    }
   if ( StorageType==LHM_HERMITIAN && RealComplex==LHM_REAL )
    { StorageType=LHM_SYMMETRIC;
      Warn("warning: HERMITIAN is equivalent to SYMMETRIC for real-valued matrices");
    };

   /* If the data parameter is non-NULL, it is assumed to point to a 
      user-allocated buffer of the correct size.  This is useful
      to create an HMatrix "wrapper" around an existing matrix,
      e.g. for the Python wrappers, without incurring a copy.  In this
      case, we do NOT free the data in the destructor. */

   ownsM = data == NULL;
   if (RealComplex==LHM_REAL)
    { DM=(double *)(data ? data : mallocEC(NumEntries() * sizeof(double)));
      ZM=0;
    }
   else
    { DM=0;
      ZM=(cdouble *)(data ? data : mallocEC(NumEntries() * sizeof(cdouble)));
    };

  ipiv=0; // this is only allocated when needed 
 
}

/***************************************************************/
/* HMatrix constructors that attempt to read an HMatrix from   */
/* a file (either a binary HDF5 file or an ASCII text file.    */
/*                                                             */
/* FileType may be either LHM_HDF5 or LHM_TEXT.                */
/*                                                             */
/* if FileType==LHM_HDF5, then Options should be the name of   */
/* the matrix within the HDF5 file.                            */
/*                                                             */
/* if FileType==LHM_TEXT, then Options should be a string of   */
/* zero or more of the following options:                      */
/*                                                             */
/*  --ncol xx : insist that there be a certain number of colums*/
/*  --nrow xx : insist that there be a certain number of rows  */
/*  --strict:   insist that all lines in the data file have    */
/*              the same number of entries (otherwise missing  */
/*              entries are set to zero)                       */
/*                                                             */
/* if an error occurs and HMatrix::AbortOnIOError==true, the   */
/* code will call ErrExit(). If AbortOnIOError==false, the     */
/* constructor will return an HMatrix with no data but with    */
/* a non-null value for the ErrMsg field.                      */
/***************************************************************/
HMatrix::HMatrix(const char *FileName, int FileType, const char *Options)
 { ReadFromFile(FileName, FileType, Options); }

int LHM_AUTO_FileType(const char *FileName) 
{
  int FileType = LHM_TEXT;

  if ( strstr(FileName, ".h5") == (FileName + strlen(FileName) - 3) )
   FileType = LHM_HDF5;
  else if ( strstr(FileName, ".hdf5") == (FileName + strlen(FileName) - 5) )
   FileType = LHM_HDF5;

  return FileType;
}

void HMatrix::ReadFromFile(const char *FileName, int FileType, const char *Options)
{
  DM=0;
  ZM=0;
  ipiv=0;
  lwork=0;
  work=0;
  liwork=0;
  iwork=0;
  ErrMsg=0;

  if (FileName==0)
   { ErrMsg=strdup("no filename specified for matrix import");
     if (AbortOnIOError) 
      ErrExit(ErrMsg);
     return;
   }

  if (FileType == LHM_AUTO)
    FileType = LHM_AUTO_FileType(FileName);

  switch(FileType)
   { 
     case LHM_TEXT: 
       ImportFromText(FileName,Options); 
       break;

     case LHM_HDF5: 
       ImportFromHDF5(FileName,Options);
       break;

     default:
       ErrExit("%s:%i: internal error",__FILE__,__LINE__);
       break;
   }

  if (ErrMsg && AbortOnIOError)
   ErrExit(ErrMsg);
}

/***************************************************************/
/* HMatrix constructor that copies an SMatrix.                 */
/***************************************************************/
HMatrix::HMatrix(SMatrix *S)
 { 
   NR=S->NR;
   NC=S->NC;
   
   RealComplex=S->RealComplex;
   StorageType=LHM_NORMAL;
   ipiv=0;
   lwork=0;
   work=0;
   liwork=0;
   iwork=0;
   ErrMsg=0;
   ownsM=true;

   int *RowStart = S->RowStart;
   int *ColIndices=S->ColIndices;
   double *SDM=S->DM;
   cdouble *SZM=S->ZM;

   if (RealComplex==LHM_REAL)
    { 
      ZM=0;
      DM=(double *)mallocEC(NR*NC*sizeof(double));
      // memset(DM, 0, NR*NC*sizeof(double)); done by mallocEC

      for(int nr=0; nr<NR; nr++) {
	size_t iend = RowStart[nr+1];
	for (size_t i = RowStart[nr]; i < iend; ++i)
	  SetEntry(nr, ColIndices[i], SDM[i]);
      }
    }
   else
    { 
      DM=0;
      ZM=(cdouble *)mallocEC(NR*NC*sizeof(cdouble));
      // memset(ZM, 0, NR*NC*sizeof(cdouble)); done by mallocEC

      for(int nr=0; nr<NR; nr++) {
	size_t iend = RowStart[nr+1];
	for (size_t i = RowStart[nr]; i < iend; ++i)
	  SetEntry(nr, ColIndices[i], SZM[i]);
      }
    };

 }

/***************************************************************/
/* destructor **************************************************/
/***************************************************************/
HMatrix::~HMatrix()
{
  if (ownsM) {
    if (DM) free(DM);
    if (ZM) free(ZM);
  }
  if (ipiv) free(ipiv);
  if (ErrMsg) free(ErrMsg);
  if (work) free(work);
}

/***************************************************************/
/* zero out the matrix *****************************************/
/***************************************************************/
void HMatrix::Zero()
{ 
  if (RealComplex==LHM_REAL)
   memset(DM,0,NumEntries()*sizeof(double));
  else
   memset(ZM,0,NumEntries()*sizeof(cdouble));
}

/***************************************************************/
/* zero out the NumRows x NumCols subblock of the matrix whose */
/* upper-left entry is (RowOffset, ColOffset)                  */
/***************************************************************/
void HMatrix::ZeroBlock(int RowOffset, int NumRows, int ColOffset, int NumCols)
{ 
  if (RowOffset==0 && NumRows==NR && ColOffset==0 && NumCols==NC)
   Zero();

  if ( RowOffset < 0 || (RowOffset+NumRows)>NR || ColOffset<0 || (ColOffset+NumCols)>NC )
   ErrExit("invalid call to HMatrix(%i,%i)::ZeroBlock(%i,%i,%i,%i)",NR,NC,RowOffset,NumRows,ColOffset,NumCols);

  for(int nr=RowOffset; nr<(RowOffset+NumRows); nr++)
   for(int nc=ColOffset; nc<(ColOffset+NumCols); nc++)
    SetEntry(nr,nc,0.0);
}

/***************************************************************/
/* copy data from another matrix. ******************************/
/* note that this copies the data (and makes sure the numbers  */
/* of rows and columns are the same) but does not change the   */
/* storage type of the destination matrix if that differs from */
/* the source matrix.                                          */
/***************************************************************/
void HMatrix::Copy(HMatrix *M)
{
  if ( M->NR != NR || M->NC != NC || M->RealComplex!=RealComplex )
   { Warn("%s:%i: matrix properties mismatch (copy failed)",__FILE__,__LINE__);
     return;
   }

  if ( M->StorageType==StorageType )
   { 
     if (RealComplex==LHM_REAL)
      memcpy(DM, M->DM, NumEntries()*sizeof(double));
     else 
      memcpy(ZM, M->ZM, NumEntries()*sizeof(cdouble));
   }
  else
   { for (int nr=0; nr<NR; nr++)
      for (int nc=0; nc<NC; nc++)
       SetEntry(nr,nc,M->GetEntry(nr,nc));
   }

}

/***************************************************************/
/* get the value of a matrix entry  ****************************/
/***************************************************************/
cdouble HMatrix::GetEntry(size_t nr, size_t nc)
{
  size_t Index, Flipped;

  Flipped=0;
  if (StorageType!=LHM_NORMAL && nr>nc) 
   { Index=nr; nr=nc; nc=Index; Flipped=1; }
  
  Index=(StorageType==LHM_NORMAL ? (nr + nc*NR) : (nr + nc*(nc+1)/2 ) );
  
  if ( RealComplex==LHM_REAL )
   return DM[Index];
  else if ( StorageType==LHM_HERMITIAN )
   return ( Flipped==1 ? conj(ZM[Index]) : ZM[Index] );
  else /* complex symmetric */
   return ZM[Index];

}

double HMatrix::GetEntryD(size_t nr, size_t nc)
{
  size_t Index;

  if (StorageType!=LHM_NORMAL && nr>nc) 
   { Index=nr; nr=nc; nc=Index; }
  
  Index=(StorageType==LHM_NORMAL ? (nr + nc*NR) : (nr + nc*(nc+1)/2 ) );
  
  if ( RealComplex==LHM_REAL )
   return DM[Index];
  else
   return real(ZM[Index]);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
cdouble HMatrix::GetTrace()
{
  if (NR!=NC)
   { Warn("GetTrace() called on non-square HMatrix");
     return 0.0;
   }

  cdouble Trace=0.0;
  for(int n=0; n<NR; n++)
   Trace+=GetEntry(n,n);
  return Trace;
}

double HMatrix::GetTraceD()
{ return real(GetTrace()); }

/***************************************************************/
/* set the value of a matrix entry  ****************************/
/***************************************************************/
void HMatrix::SetEntry(size_t nr, size_t nc, cdouble Entry)
{
  size_t Index, Flipped;
  
  Flipped=0;
  if (StorageType!=LHM_NORMAL && nr>nc) 
   { Index=nr; nr=nc; nc=Index; Flipped=1; }
  
  Index=(StorageType==LHM_NORMAL ? (nr + nc*NR) : (nr + nc*(nc+1)/2 ) );
  
  if (RealComplex==LHM_REAL)
   {
   DM[Index] = real(Entry);
   if (imag(Entry) != 0.0) Warn("warning: imaginary part %.2e ignored.", imag(Entry));
   }
  else if ( StorageType==LHM_HERMITIAN )
   ZM[Index] = ( Flipped==1 ? conj(Entry) : Entry );
  else /* complex symmetric */
   ZM[Index] = Entry;

}

void HMatrix::SetEntry(size_t nr, size_t nc, double Entry)
{
  size_t Index;

  if (StorageType!=LHM_NORMAL && nr>nc) 
   { Index=nr; nr=nc; nc=Index; }
  
  Index=(StorageType==LHM_NORMAL ? (nr + nc*NR) : (nr + nc*(nc+1)/2 ) );
  
  if (RealComplex==LHM_REAL)
   DM[Index] = Entry;
  else 
   ZM[Index] = Entry;

}

/***************************************************************/
/* augment a matrix entry **************************************/
/***************************************************************/
void HMatrix::AddEntry(size_t nr, size_t nc, cdouble Entry)
{
  size_t Index, Flipped;

  Flipped=0;
  if (StorageType!=LHM_NORMAL && nr>nc) 
   { Index=nr; nr=nc; nc=Index; Flipped=1; }
  
  Index=(StorageType==LHM_NORMAL ? (nr + nc*NR) : (nr + nc*(nc+1)/2 ) );
  
  if (RealComplex==LHM_REAL)
   DM[Index] += real(Entry);
  else if ( StorageType==LHM_HERMITIAN )
   ZM[Index] += ( Flipped == 1 ? conj(Entry) : Entry );
  else 
   ZM[Index] += Entry;
}

void HMatrix::AddEntry(size_t nr, size_t nc, double Entry)
{
  size_t Index;

  if (StorageType!=LHM_NORMAL && nr>nc) 
   { Index=nr; nr=nc; nc=Index; }
  
  Index=(StorageType==LHM_NORMAL ? (nr + nc*NR) : (nr + nc*(nc+1)/2 ) );
  
  if (RealComplex==LHM_REAL)
   DM[Index] += Entry;
  else 
   ZM[Index] += cdouble(Entry,0.0);
}

/***************************************************************/
/* multiply a single matrix entry by a scalar factor ***********/
/***************************************************************/
void HMatrix::ScaleEntry(size_t nr, size_t nc, cdouble ScaleFactor)
{ SetEntry(nr, nc, ScaleFactor*GetEntry(nr,nc)); }

/***************************************************************/
/* return a pointer to the head of the double-valued or        */
/* cdouble-valued array of elements for column #nc             */
/***************************************************************/
void *HMatrix::GetColumnPointer(size_t nc)
{ if (DM)
   return DM + nc*NR;
  else
   return ZM + nc*NR;
}

/***************************************************************/
/* scale the entire matrix by a scalar multiple ****************/
/***************************************************************/
void HMatrix::Scale(cdouble Alpha)
{
  size_t n, NumEntries_ = NumEntries();

  if (RealComplex==LHM_REAL)
   for(n=0; n<NumEntries_; n++)
    DM[n] *= real(Alpha);
  else
   for(n=0; n<NumEntries_; n++)
    ZM[n]*=Alpha;

}

void HMatrix::Scale(double Alpha)
{ Scale(cdouble(Alpha)); }

/***************************************************************/
/* replace the matrix with its adjoint, i.e. conjugate         */
/* transpose.                                                  */
/***************************************************************/
void HMatrix::Adjoint()
{ 
  if ( RealComplex==LHM_REAL && StorageType==LHM_NORMAL )
   { 
     // square case 
     if (NR==NC)
      { 
        for(int nr=0; nr<NR; nr++)
         for(int nc=nr+1; nc<NC; nc++)
          { double TD=DM[nr+nc*NR]; 
            DM[nr+nc*NR]=DM[nc+nr*NR]; 
            DM[nc+nr*NR]=TD; 
          }
      }
     else
      { 
        double *DM2=(double *)mallocEC(NR*NC*sizeof(double));
        for(int nr=0; nr<NR; nr++)
         for(int nc=0; nc<NC; nc++)
          DM2[nc + nr*NC] = DM[nr + nc*NR];
	memcpy(DM, DM2, NR*NC*sizeof(double));
        free(DM2);
        int temp=NR; NR=NC; NC=temp;
      };
   }
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_NORMAL )
   { 
     // handle square case in-place
     if (NR==NC)
      { 
        for(int nr=0; nr<NR; nr++)
         { 
           cdouble TZ=ZM[nr+nr*NR]; 
           ZM[nr+nr*NR]=conj(TZ);

           for(int nc=nr+1; nc<NC; nc++)
            { TZ=ZM[nr+nc*NR]; 
              ZM[nr+nc*NR]=conj(ZM[nc+nr*NR]); 
              ZM[nc+nr*NR]=conj(TZ); 
            }
         }
      }
     else
      { 
        cdouble *ZM2=(cdouble *)mallocEC(NR*NC*sizeof(cdouble));
        for(int nr=0; nr<NR; nr++)
         for(int nc=0; nc<NC; nc++)
          ZM2[nc + nr*NC] = conj(ZM[nr + nc*NR]);
	memcpy(ZM, ZM2,NR*NC*sizeof(cdouble));
        free(ZM2);
        int temp=NR; NR=NC; NC=temp;
      }
   }
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_HERMITIAN )
   { 
     /* the matrix already is its own transpose */
   } 
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_SYMMETRIC )
   { int NumEntries_ = NumEntries();
     for (int n=0; n<NumEntries_; n++)
      ZM[n] = conj(ZM[n]);
   }

}

/***************************************************************/
/* replace the matrix with its non-conjugate transpose *********/
/***************************************************************/
void HMatrix::Transpose()
{ 
  if ( RealComplex==LHM_REAL )
   Adjoint();
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_NORMAL )
   { 
     // handle square case in-place
     if (NR==NC)
      { for(int nr=0; nr<NR; nr++)
         for(int nc=nr+1; nc<NC; nc++)
          { cdouble TZ=ZM[nr+nc*NR];
            ZM[nr+nc*NR]=ZM[nc+nr*NR];
            ZM[nc+nr*NR]=TZ;
          }
      }
     else
      { cdouble *ZM2=(cdouble *)mallocEC(NR*NC*sizeof(cdouble));
        for(int nr=0; nr<NR; nr++)
         for(int nc=0; nc<NC; nc++)
          ZM2[nc + nr*NC] = ZM[nr + nc*NR];
	memcpy(ZM, ZM2, NR*NC*sizeof(cdouble));
        free(ZM2);
        int temp=NR; NR=NC; NC=temp; 
      }
   }
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_HERMITIAN )
   { 
     int NumEntries_ = NumEntries();
     for (int n=0; n<NumEntries_; n++)
      ZM[n]=conj(ZM[n]);
   } 
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_SYMMETRIC )
   { 
     /* the matrix already is its own non-conjugate transpose */
   }

}

/***************************************************************/
/* stamp the full matrix B into this matrix in such a way that */
/* the upper-left entry of B goes into slot (RowOffset,ColOffset). */
/***************************************************************/
void HMatrix::InsertBlock(HMatrix *B, int RowOffset, int ColOffset, double Scale)
{ 
  if ( ((RowOffset + B->NR) > NR) || ((ColOffset + B->NC) > NC) )
   ErrExit("InsertBlock(): block insertion exceeds matrix size");

  for (int nr=0; nr<B->NR; nr++)
   for (int nc=0; nc<B->NC; nc++)
    SetEntry(RowOffset+nr, ColOffset+nc, Scale * B->GetEntry(nr,nc) );
}

/****************************************************************/
/* like the previous routine, but only stamp in the NRBxNRC     */
/* subblock of B starting at (BRowOffset, BColOffsets)          */
/****************************************************************/
void HMatrix::InsertBlock(HMatrix *B, int RowOffset, int ColOffset,
                          int NRB, int NCB, int BRowOffset, int BColOffset, double Scale)
{
  if ( ((RowOffset + NRB) > NR) || ((ColOffset + NCB) > NC) )
   ErrExit("InsertBlock(): block insertion exceeds matrix size");
  if ( ((BRowOffset + NRB) > B->NR) || ((BColOffset + NCB) > B->NC) )
   ErrExit("InsertBlock(): block insertion exceeds block size");

  for (int nr=0; nr<NRB; nr++)
   for (int nc=0; nc<NCB; nc++)
    SetEntry(RowOffset+nr, ColOffset+nc, Scale * B->GetEntry(BRowOffset+nr,BColOffset+nc) );

}

/***************************************************************/
/* stamp the conjugate transpose of B into the matrix in       */
/* such a way that:                                            */
/* B(0,0) goes into slot (RowOffset,   ColOffset)              */
/* B(0,1) goes into slot (RowOffset+1, ColOffset)              */
/* B(1,0) goes into slot (RowOffset,   ColOffset+1)            */
/* etc.                                                        */
/***************************************************************/
void HMatrix::InsertBlockAdjoint(HMatrix *B, int RowOffset, int ColOffset, double Scale)
{ 
  if ( ((RowOffset + B->NC) > NR) || ((ColOffset + B->NR) > NC) )
   ErrExit("InsertBlockAdjoint(): block insertion exceeds matrix size");

  if (B->RealComplex==LHM_COMPLEX)
   { for (int nr=0; nr<B->NR; nr++)
      for (int nc=0; nc<B->NC; nc++)
       SetEntry(RowOffset+nc, ColOffset+nr, Scale*conj(B->GetEntry(nr,nc)) );
   }
  else
   { for (int nr=0; nr<B->NR; nr++)
      for (int nc=0; nc<B->NC; nc++)
       SetEntry(RowOffset+nc, ColOffset+nr, Scale*B->GetEntryD(nr,nc) );
   };
}

/***************************************************************/
/* stamp the non-conjugate transpose of B into the matrix in   */
/* such a way that:                                            */
/* B(0,0) goes into slot (RowOffset,   ColOffset)              */
/* B(0,1) goes into slot (RowOffset+1, ColOffset)              */
/* B(1,0) goes into slot (RowOffset,   ColOffset+1)            */
/* etc.                                                        */
/***************************************************************/
void HMatrix::InsertBlockTranspose(HMatrix *B, int RowOffset, int ColOffset, double Scale)
{ 
  if ( ((RowOffset + B->NC) > NR) || ((ColOffset + B->NR) > NC) )
   ErrExit("InsertBlockTranspose(): block insertion exceeds matrix size");

  if (B->RealComplex==LHM_COMPLEX)
   { for (int nr=0; nr<B->NR; nr++)
      for (int nc=0; nc<B->NC; nc++)
       SetEntry(RowOffset+nc, ColOffset+nr, Scale*B->GetEntry(nr,nc));
   }
  else
   { for (int nr=0; nr<B->NR; nr++)
      for (int nc=0; nc<B->NC; nc++)
       SetEntry(RowOffset+nc, ColOffset+nr, Scale*B->GetEntryD(nr,nc) );
   }
}

/***************************************************************/
/* like InsertBlock(), but entries are augmented rather than   */
/* overwritten.                                                */
/***************************************************************/
void HMatrix::AddBlock(HMatrix *B, int RowOffset, int ColOffset, cdouble ScaleFactor)
{ 
  if ( ((RowOffset + B->NR) > NR) || ((ColOffset + B->NC) > NC) )
   ErrExit("AddBlock(): block insertion exceeds matrix size");

  for (int nr=0; nr<B->NR; nr++)
   for (int nc=0; nc<B->NC; nc++)
    AddEntry(RowOffset+nr, ColOffset+nc, ScaleFactor * B->GetEntry(nr,nc) );
}

void HMatrix::AddBlockAdjoint(HMatrix *B, int RowOffset, int ColOffset, cdouble ScaleFactor)
{ 
  if ( ((RowOffset + B->NC) > NR) || ((ColOffset + B->NR) > NC) )
   ErrExit("AddBlockAdjoint(): block addition exceeds matrix size");

  if (B->RealComplex==LHM_COMPLEX)
   { for (int nr=0; nr<B->NR; nr++)
      for (int nc=0; nc<B->NC; nc++)
       AddEntry(RowOffset+nc, ColOffset+nr, ScaleFactor*conj(B->GetEntry(nr,nc)) );
   }
  else
   { for (int nr=0; nr<B->NR; nr++)
      for (int nc=0; nc<B->NC; nc++)
       AddEntry(RowOffset+nc, ColOffset+nr, ScaleFactor * B->GetEntryD(nr,nc) );
   }
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void HMatrix::AddBlock(SMatrix *B, int RowOffset, int ColOffset, cdouble ScaleFactor)
{ 
  if ( ((RowOffset + B->NR) > NR) || ((ColOffset + B->NC) > NC) )
   ErrExit("AddBlock(): block insertion exceeds matrix size");
  
  for(int nr=0; nr<B->NR; nr++)
   { int *nc;
     void *Entries;
     int NNZ = B->GetRow(nr, &nc, &Entries);
 
     if (B->RealComplex==LHM_REAL) 
      { double *dEntries=(double *)Entries, RSF=real(ScaleFactor);
        for(int nnz=0; nnz<NNZ; nnz++)
         AddEntry(RowOffset+nr, ColOffset+nc[nnz], RSF*dEntries[nnz]);
      }
     else
      { cdouble *cdEntries=(cdouble *)Entries;
        for(int nnz=0; nnz<NNZ; nnz++)
         AddEntry(RowOffset+nr, ColOffset+nc[nnz], ScaleFactor*cdEntries[nnz]);
      };
   };

}

void HMatrix::Add(HMatrix *B, cdouble ScaleFactor)
 { AddBlock(B, 0, 0, ScaleFactor); }

void HMatrix::Add(SMatrix *B, cdouble ScaleFactor)
 { AddBlock(B, 0, 0, ScaleFactor); }

/***************************************************************/
/***************************************************************/
/***************************************************************/
void HMatrix::ExtractBlock(int RowOffset, int ColOffset, HMatrix *B)
{ 
  if ( ((RowOffset + B->NR) > NR) || ((ColOffset + B->NC) > NC) )
   ErrExit("ExtractBlock(): block extraction exceeds matrix size");

  for (int nr=0; nr<B->NR; nr++)
   for (int nc=0; nc<B->NC; nc++)
    B->SetEntry(nr, nc, GetEntry(RowOffset+nr,ColOffset+nc) );
}

/***************************************************************/
/* Get the real part of the matrix  ****************************/
/***************************************************************/
HMatrix *GetReal(HMatrix *A) {
  if (A->RealComplex == LHM_REAL) return A;
  HMatrix *C;
  if (A->StorageType == LHM_SYMMETRIC || A->StorageType == LHM_HERMITIAN)
    C = new HMatrix(A->NR, A->NC, LHM_REAL, LHM_SYMMETRIC);
  else
    C = new HMatrix(A->NR, A->NC, LHM_REAL, LHM_NORMAL);
  for (size_t Index=0; Index < A->NumEntries(); Index++){
    (C->DM)[Index] = real((A->ZM)[Index]);
  }
  return C;
}

/***************************************************************/
/* Make an unpacked copy of the matrix *************************/
/***************************************************************/
HMatrix *CopyHMatrixUnpacked(HMatrix *Mpacked) {
  if (!Mpacked) return NULL;
  HMatrix *M = new HMatrix(Mpacked->NR, Mpacked->NC, Mpacked->RealComplex);
  M->Copy(Mpacked);
  return M;
}

/***************************************************************/
/* Concatenate two matrices A and B to form a new matrix C.    */
/*                                                             */
/* The 'How' argument specifies whether B is to be tacked on   */
/* to the side of A or below A, as follows:                    */
/*                                                             */
/* How=LHM_HORIZONTAL: return value = [A B]                    */
/*                                                             */
/* How=LHM_VERTICAL:   return value = [A; B]                   */
/*                                                             */
/* Note: Depending on the sizes of A and B, there may be only  */
/*       one feasible way to do the concatenation, in which    */
/*       case the How argument is silently ignored.            */
/***************************************************************/
HMatrix *Concat(HMatrix *A, HMatrix *B, int How)
{
  if ( (A->NR == B->NR) && (A->NC != B->NC) )
   How = LHM_HORIZONTAL;
  else if ( (A->NC == B->NC) && (A->NR != B->NR) )
   How = LHM_VERTICAL;
  else if ( (A->NC != B->NC) && (A->NR != B->NR) )
   ErrExit("attempt to concatenate incompatible matrices");

  // if none of the above then A->NR=B->NR and A->NC=B->NC; use caller-specified How value

  int NR = (How==LHM_HORIZONTAL) ? A->NR : (A->NR + B->NR);
  int NC = (How==LHM_VERTICAL)   ? A->NC : (A->NC + B->NC);
  int DataType = (A->RealComplex==LHM_REAL && B->RealComplex==LHM_REAL) ? LHM_REAL : LHM_COMPLEX;

  HMatrix *C = new HMatrix(NR, NC, DataType);
  C->InsertBlock(A, 0, 0);
  if (How==LHM_HORIZONTAL)
   C->InsertBlock(B, 0, A->NC);
  else
   C->InsertBlock(B, A->NR, 0);

  return C; 

}

/***************************************************************/
/* Bilinear product, X'*this*Y; no error checking or           */
/* optimization                                                */
/***************************************************************/
cdouble HMatrix::BilinearProduct(HVector *X, HVector *Y)
{
  if (Y==0) Y=X;
  cdouble Sum=0.0;
  for(int nr=0; nr<NR; nr++)
   for(int nc=0; nc<NC; nc++)
    Sum += conj(X->GetEntry(nr)) * GetEntry(nr, nc) * Y->GetEntry(nc);
  return Sum;
}

double HMatrix::BilinearProductD(HVector *X, HVector *Y)
 { return real(BilinearProduct(X,Y)); }


/***************************************************************/
/* Return a matrix of the given dimensions and data type, using*/
/* M if it is non-null and its properties match, and otherwise */
/* allocating a new matrix with the given properties (in which */
/* case M is deallocated first.)                               */
/***************************************************************/
HMatrix *CheckHMatrix(HMatrix *M, int NR, int NC, int RealComplex, const char *Context)
{
  if (M && (M->NR!=NR || M->NC!=NC || M->RealComplex!=RealComplex) )
   { Warn("incorrect matrix {%i,%i,%i}!={%i,%i,%i}%s%s (reallocating)",
           M->NR,M->NC,M->RealComplex, NR,NC,RealComplex,
           Context ? " in " : "", Context ? Context : "");
     delete M;
     M=0;
   }
  if (M==0)
   M=new HMatrix(NR, NC, RealComplex);
  return M;
}

HMatrix *CheckHMatrix(HMatrix *M, HMatrix *M0, const char *Context)
{ return CheckHMatrix(M, M0->NR, M0->NC, M0->RealComplex, Context); }
