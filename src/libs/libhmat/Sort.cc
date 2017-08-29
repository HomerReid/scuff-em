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

#include <libhrutil.h>
#include "libhmat.h"

#define SORTCOLUMN 0

/***************************************************************/
/***************************************************************/
/***************************************************************/
int RealAscendingSort(const void *p1, const void *p2)
{ 
  double *RP1=(double *)p1;
  double *RP2=(double *)p2;

  double Diff=RP1[SORTCOLUMN] - RP2[SORTCOLUMN];
  if (Diff>0.0) 
   return 1;
  if (Diff==0.0) 
   return 0;
  return -1;
}

int RealDescendingSort(const void *p1, const void *p2)
{ 
  double *RP1=(double *)p1;
  double *RP2=(double *)p2;

  double Diff=RP2[SORTCOLUMN] - RP1[SORTCOLUMN];
  if (Diff>0.0) 
   return 1;
  if (Diff==0.0) 
   return 0;
  return -1;
}

int ComplexAscendingSort(const void *p1, const void *p2)
{ 
  cdouble *RP1=(cdouble *)p1;
  cdouble *RP2=(cdouble *)p2;

  double Diff=abs(RP1[SORTCOLUMN]) - abs(RP2[SORTCOLUMN]);
  if (Diff>0.0) 
   return 1;
  if (Diff==0.0) 
   return 0;
  return -1;
}

int ComplexDescendingSort(const void *p1, const void *p2)
{ 
  cdouble *RP1=(cdouble *)p1;
  cdouble *RP2=(cdouble *)p2;

  double Diff=abs(RP2[SORTCOLUMN]) - abs(RP1[SORTCOLUMN]);
  if (Diff>0.0) 
   return 1;
  if (Diff==0.0) 
   return 0;
  return -1;
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

void HMatrix::Sort(int WhichColumn)
{
  Sort(WhichColumn,0);
}
