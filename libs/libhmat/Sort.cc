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
 * Sort.cc     -- a routine for sorting an HMatrix in increasing or 
 *             -- decreasing order of one of its columns
 *
 * homer reid  -- 9/2011
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#include <algorithm>
#include <vector>

#include <libhrutil.h>
#include "libhmat.h"

typedef std::vector<int>  ivec;
typedef std::vector<char> cvec;

/***************************************************************/
/* row-sorting class passed to std::sort.                      */
/*                                                             */
/* sort type = M:    ascending by absolute magnitude           */
/*             I   : ascending by imag part                    */
/*             R, V: ascending by real part                    */
/*           = m, i, r, v: descending                          */
/***************************************************************/
struct RowSorter
{
  RowSorter(HMatrix *M, ivec SortColumns, cvec SortType);

  bool operator() (int nri, int nrj);

  HMatrix *M;
  ivec SortColumns;
  cvec SortType;
  
};

RowSorter::RowSorter(HMatrix *_M, ivec _SortColumns, cvec _SortType): M(_M), SortColumns(_SortColumns), SortType(_SortType)
{ 
  if ( SortType.size()!=SortColumns.size() )
   ErrExit("%s:%i: numbers of columns and flags disagree (%i,%i)",
            __FILE__,__LINE__,SortType.size(), SortColumns.size());

  for(unsigned nnc=0; nnc<SortColumns.size(); nnc++)
   { int nc=SortColumns[nnc];
     if ( nc<0 || nc>=M->NC )
      ErrExit("invalid column index (%i,%i) in HMatrix::Sort",nc,M->NC);
     if ( !strchr("mirvMIRV", SortType[nnc]) )
      ErrExit("unknown sort type (%c) in HMatrix::Sort",SortType[nnc]);
   }
}

bool RowSorter::operator() (int nri, int nrj)
{
  for (unsigned nnc=0, nc=SortColumns[0]; nnc<SortColumns.size(); nc=SortColumns[++nnc])
   { cdouble zvi = M->GetEntry(nri, nc), zvj = M->GetEntry(nrj, nc);
     if (zvi==zvj) continue;
     char c = tolower(SortType[nnc]);
     bool Descending = ( c == SortType[nnc] );
     bool Lesser = ( c=='m' ?  (abs(zvi)  < abs(zvj)  ) :
                     c=='i' ?  (imag(zvi) < imag(zvj) ) :
                               (real(zvi) < real(zvj) )
                   );
     return Lesser ^ Descending;
   }
  return true;
}

/***************************************************************/
/* sort type = M:    ascending by absolute magnitude           */
/*             I:    ascending by imag part                    */
/*             R, V: ascending by real part                    */
/*           = m, i, r, v: descending                          */
/***************************************************************/
HMatrix *SortHMatrix(HMatrix *M, ivec SortColumns, cvec SortType)
{
  int NR=M->NR, NC=M->NC;

  ivec RowIndices(NR);
  for(int nr=0; nr<NR; nr++) 
   RowIndices[nr]=nr;
  RowSorter MyRowSorter(M, SortColumns, SortType);
  std::sort (RowIndices.begin(), RowIndices.end(), MyRowSorter);

  HMatrix *NewM = new HMatrix(M);
  for(int nr=0; nr<NR; nr++)
   for(int nc=0; nc<NC; nc++)
    NewM->SetEntry(nr, nc, M->GetEntry(RowIndices[nr], nc));

  return NewM;
}

void HMatrix::Sort(ivec SortColumns, cvec SortType)
{ 
  if (StorageType!=LHM_NORMAL)
   ErrExit("matrix sort not yet implemented for packed matrices");
  HMatrix *NewM=SortHMatrix(this, SortColumns, SortType);
  if (RealComplex==LHM_REAL)
   memcpy(DM, NewM->DM, NR*NC*sizeof(double));
  else
   memcpy(ZM, NewM->ZM, NR*NC*sizeof(cdouble));
  delete NewM;
}

void HMatrix::Sort(ivec SortColumns)
{ Sort(SortColumns, cvec(SortColumns.size(),'V')); }


void HMatrix::Sort(int WhichColumn, const char *Options)
{
  ivec SortColumns(1, WhichColumn);
  cvec SortType(1, 'V');
  if (Options && !strcasecmp(Options,"--descending"))
   SortType[0]='v';
  Sort(SortColumns, SortType);
}

/***************************************************************/
/* sort the rows of an HMatrix by the value of one of its      */
/* columns.                                                    */
/*                                                             */
/* the matrix must have normal (non-symmetric) storage.        */
/*                                                             */
/* for complex matrices, the magnitude of the column entry     */
/* is used as the sort key.                                    */
/*                                                             */
/* the Options string may include the following:               */
/*  --ascending  (default)                                     */
/*  --descending                                               */
/***************************************************************/
#if 0
void HMatrix::Sort(int WhichColumn, const char *Options)
{
  int Ascending=1;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (StorageType!=LHM_NORMAL)
   ErrExit("matrix sort not yet implemented for packed matrices");
  if (WhichColumn<0 || WhichColumn>NC)
   ErrExit("invalid sort column (%i) passed to HMatrix::Sort",WhichColumn);
 
  /***************************************************************/
  /* process options *********************************************/
  /***************************************************************/
  if (Options && strlen(Options)>0 )
   { char *Tokens[50];
     int nt, NumTokens;
     char buffer[1000];
     strncpy(buffer,Options,1000);
     NumTokens=Tokenize(buffer,Tokens,50);
     for(nt=0; nt<NumTokens; nt++)
      { if (StrCaseCmp(Tokens[nt],"--ascending"))
         Ascending=1;
        if (StrCaseCmp(Tokens[nt],"--decending"))
         Ascending=0;
        else
         ErrExit("HMatrix::Sort: invalid option (%s)",Tokens[nt]);
      };
   };
  
  /***************************************************************/
  /* We first interchange columns so that the column on which we */
  /* are sorting is the first column, i.e. column index 0        */
  /* (this step could be avoided if there were some way of       */
  /* passing a SortColumn argument to the sorting routines above,*/
  /* but apparently the stdlibc qsort doesn't allow this?        */
  /***************************************************************/
  if (WhichColumn!=SORTCOLUMN)
   { cdouble Temp;
     for(int nr=0; nr<NR; nr++)
      { Temp=GetEntry(nr, WhichColumn);
        SetEntry(nr, WhichColumn, GetEntry(nr, SORTCOLUMN) );
        SetEntry(nr, SORTCOLUMN, Temp);
      };
   };
  
  /***************************************************************/
  /* next transpose the matrix so that the columns are stored    */
  /* contiguously in memory, not the rows                        */
  /***************************************************************/
  Transpose();

  /***************************************************************/
  /* note that 'NR' and 'NC' appear to be swapped in what        */
  /* follows; this is because we turned M into its transpose     */
  /* which interchanges these two quantities                     */
  /***************************************************************/
  if (Ascending==1 && RealComplex==LHM_REAL)
   qsort(DM, NC, NR*sizeof(double), RealAscendingSort);
  else if (Ascending==0 && RealComplex==LHM_REAL)
   qsort(DM, NC, NR*sizeof(double), RealDescendingSort);
  else if (Ascending==1 && RealComplex==LHM_COMPLEX)
   qsort(ZM, NC, NR*sizeof(cdouble), ComplexAscendingSort);
  else if (Ascending==0 && RealComplex==LHM_REAL)
   qsort(ZM, NC, NR*sizeof(cdouble), ComplexDescendingSort);
  
  /***************************************************************/
  /* now transpose back and reexchange columns                   */
  /***************************************************************/
  Transpose();
  
  if (WhichColumn!=SORTCOLUMN)
   { cdouble Temp;
     for(int nr=0; nr<NR; nr++)
      { Temp=GetEntry(nr, WhichColumn);
        SetEntry(nr, WhichColumn, GetEntry(nr, SORTCOLUMN) );
        SetEntry(nr, SORTCOLUMN, Temp);
      };
   };

} 
#endif
