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
#include <ctype.h>

#include <libhrutil.h>
#include "libhmat.h"

#define MAXSTR 100

/***************************************************************/
/* This routine parses a string of text that selects a subset  */
/* of the indices in one dimension of an array.                */
/*                                                             */
/* Although the syntax here is reminiscent of that used in     */
/* languages like matlab or julia, an important distinction is */
/* that here indices are ZERO-BASED.                           */
/*                                                             */
/* Possible strings include                                    */
/*                                                             */
/*  :       (for all elements in the array)                    */
/*                                                             */
/*  3       (for just the 4th element)                         */
/*                                                             */
/* 0:7      (for the first 7 elements)                         */
/* 1:8      (for the 7 elements starting at the 2nd element)   */
/*                                                             */
/* 2:end    (for all but the first 2 elements)                 */
/* 0:2:end  (for every other element starting with the 1st)    */
/* 1:2:end  (for every other element starting with the 2nd)    */
/*                                                             */
/* On input, MIS is the string and End should be the number    */
/* of elements in this dimension of the array, which means that*/
/* End-1 is the number filled in for the string 'end' in MIS.  */
/*                                                             */
/* Return value:                                               */
/*  -1  --> the string was not successfully parsed.            */
/*  >=0 --> the number of elements in the subset.              */
/***************************************************************/
int ParseMatrixIndexString(const char *MIS, int End,
                           int *Start, int *Stop, int *Inc)
{
  if (!MIS)
   return false;

  /*--------------------------------------------------------------*/
  /*- Separately detect strings consisting of all white space    -*/
  /*- plus a single ':' character.                               -*/
  /*--------------------------------------------------------------*/
  int NumColons=0, NumNonWhite=0;
  for(int n=0; MIS[n]!=0; n++)
   { if (MIS[n]==':') 
      NumColons++;
     else if (!isspace(MIS[n])) 
      NumNonWhite++;
   };
  if (NumColons==1 && NumNonWhite==0)
   { *Start=0;
     *Stop=End-1;
     *Inc=1;
     return End;
   };

  /*--------------------------------------------------------------*/
  /*- more than two colons is a syntax error ---------------------*/
  /*--------------------------------------------------------------*/
  if (NumColons>2) 
   return -1;

  /*--------------------------------------------------------------*/
  /*- extract colon-separated tokens -----------------------------*/
  /*--------------------------------------------------------------*/
  char MISCopy[31];
  strncpy(MISCopy, MIS, 30);
  MISCopy[30]=0;
  char *Tokens[4];
  int NumTokens=Tokenize(MISCopy, Tokens, 4, ":");

  /*--------------------------------------------------------------*/
  /*- convert strings to integers, including handling of "end"   -*/
  /*--------------------------------------------------------------*/
  int TokenValues[3];
  for(int nt=0; nt<NumTokens; nt++)
   { if ( !strcasecmp(Tokens[nt],"end") )
      TokenValues[nt]=End-1;
     else if ( 1!=sscanf(Tokens[nt],"%i",TokenValues+nt) )
      return -1;
   };

  /*--------------------------------------------------------------*/
  /*- zero tokens is a syntax error                              -*/
  /*--------------------------------------------------------------*/
  if ( NumTokens==0 )
   return 1;

  /*--------------------------------------------------------------*/
  /*- one token is a string like "3", unless it was a string like */
  /*- "3:" which is a syntax error                                */
  /*--------------------------------------------------------------*/
  if (NumTokens==1)
   { if (NumColons!=0)
      return false;
     *Start=*Stop=TokenValues[0]; *Inc=1;
     return 1;
   };

  /*--------------------------------------------------------------*/
  /*- two tokens is something like 0:5  or 3:end  ----------------*/
  /*--------------------------------------------------------------*/
  if (NumTokens==2)
   { *Start=TokenValues[0];
     *Stop=TokenValues[1];
     *Inc=1;
     return 1 + *Stop - *Start;
   };

  /*--------------------------------------------------------------*/
  /*- three tokens is something like 1:2:7   ---------------------*/
  /*--------------------------------------------------------------*/
  if (NumTokens==3)
   { *Start=TokenValues[0];
     *Inc=TokenValues[1];
     *Stop=TokenValues[2];

     int Count=0;
     if ( *Inc > 0)
      { for (int x=*Start; x<=*Stop; x+=*Inc)
        Count++;
      }
     else if ( *Inc < 0)
      { for (int x=*Start; x>=*Stop; x+=*Inc)
        Count++;
      }
     else  // Inc==0
      return -1;

     return Count;
   };

  return false; 
  
}

/***************************************************************/
/* There are two broad classes of variants of GetEntries:      */
/*                                                             */
/*  1. variants which extract entries from an HVector, or from */
/*     a single row or column of an HMatrix, and stick them    */
/*     into C++ arrays of type double[] or cdouble[].          */
/*     (In this case, a pointer to the array is returned; if   */
/*     the user passed in a buffer of tpe ay is returned; if   */
/*                                                             */
/***************************************************************/

/***************************************************************/
/***************************************************************/
/***************************************************************/
double *HMatrix::GetEntriesD(const char *RowString, int Col, double *Entries)
{
  /***************************************************************/
  /* parse column string *****************************************/
  /***************************************************************/
  int RowStart, RowStop, RowInc;
  int Length=ParseMatrixIndexString(RowString, NR, &RowStart, &RowStop, &RowInc);
  if (Length<=0)
   ErrExit("invalid row index string %s",RowString);

  /***************************************************************/
  /* allocate output buffer if necessary *************************/
  /***************************************************************/
  if (Entries==0)
   Entries=(double *)malloc(Length*sizeof(double));

  /***************************************************************/
  /* do the extraction *******************************************/
  /***************************************************************/
  HMatrix B(1, Length, LHM_REAL, LHM_NORMAL, (void *)Entries);
  DoGetEntries(RowStart, RowStop, RowInc, Length, Col, Col, 1, 1, &B, 0, 0);

  return Entries;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
cdouble *HMatrix::GetEntries(const char *RowString, int Col, cdouble *Entries)
{
  /***************************************************************/
  /* parse column string *****************************************/
  /***************************************************************/
  int RowStart, RowStop, RowInc;
  int Length=ParseMatrixIndexString(RowString, NR, &RowStart, &RowStop, &RowInc);
  if (Length<=0)
   ErrExit("invalid row index string %s",RowString);

  /***************************************************************/
  /* allocate output buffer if necessary *************************/
  /***************************************************************/
  if (Entries==0)
   Entries=(cdouble *)malloc(Length*sizeof(cdouble));

  /***************************************************************/
  /* do the extraction *******************************************/
  /***************************************************************/
  HMatrix B(1, Length, LHM_COMPLEX, LHM_NORMAL, (void *)Entries);
  DoGetEntries(RowStart, RowStop, RowInc, Length, Col, Col, 1, 1, &B, 0, 0);

  return Entries;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
double *HMatrix::GetEntriesD(int Row, const char *ColString, double *Entries)
{
  /***************************************************************/
  /* parse column string *****************************************/
  /***************************************************************/
  int ColStart, ColStop, ColInc;
  int Length=ParseMatrixIndexString(ColString, NC, &ColStart, &ColStop, &ColInc);
  if (Length<=0)
   ErrExit("invalid column index string %s",ColString);

  /***************************************************************/
  /* allocate output buffer if necessary *************************/
  /***************************************************************/
  if (Entries==0)
   Entries=(double *)malloc(Length*sizeof(double));

  /***************************************************************/
  /* do the extraction *******************************************/
  /***************************************************************/
  HMatrix B(1, Length, LHM_REAL, LHM_NORMAL, (void *)Entries);
  DoGetEntries(Row, Row, 1, 1, ColStart, ColStop, ColInc, Length, &B, 0, 0);

  return Entries;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
cdouble *HMatrix::GetEntries(int Row, const char *ColString, cdouble *Entries)
{
  /***************************************************************/
  /* parse column string *****************************************/
  /***************************************************************/
  int ColStart, ColStop, ColInc;
  int Length=ParseMatrixIndexString(ColString, NC, &ColStart, &ColStop, &ColInc);
  if (Length<=0)
   ErrExit("invalid column index string %s",ColString);

  /***************************************************************/
  /* allocate output buffer if necessary *************************/
  /***************************************************************/
  if (Entries==0)
   Entries=(cdouble *)malloc(Length*sizeof(cdouble));

  /***************************************************************/
  /* do the extraction *******************************************/
  /***************************************************************/
  HMatrix B(1, Length, LHM_COMPLEX, LHM_NORMAL, (void *)Entries);
  DoGetEntries(Row, Row, 1, 1, ColStart, ColStop, ColInc, Length, &B, 0, 0);

  return Entries;
}

/***************************************************************/
/* RowColString is a string like "[8,2:2:4]"                   */
/***************************************************************/ 
HMatrix *HMatrix::ExtractEntries(const char *RowColString, HMatrix *B)
{ 
  /***************************************************************/
  /* extract row and column index strings from first argument ****/
  /***************************************************************/
  char RCSCopy[MAXSTR];
  strncpy(RCSCopy, RowColString, MAXSTR);

  char *Tokens[3];
  int NumTokens=Tokenize(RCSCopy,Tokens,3,"[],");
  if (NumTokens!=2)
   ErrExit("invalid array index string %s",RowColString);
 
  /***************************************************************/
  /* parse row and column index strings                          */
  /***************************************************************/
  int RowStart, RowStop, RowInc;
  int RowLen=ParseMatrixIndexString(Tokens[0], NR, &RowStart, &RowStop, &RowInc);
  if (RowLen<=0)
   ErrExit("invalid row index string %s",Tokens[0]);

  int ColStart, ColStop, ColInc;
  int ColLen=ParseMatrixIndexString(Tokens[1], NC, &ColStart, &ColStop, &ColInc);
  if (ColLen<=0)
   ErrExit("invalid column index string %s",Tokens[1]);
 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  return DoGetEntries(RowStart, RowStop, RowInc, RowLen, 
                      ColStart, ColStop, ColInc, ColLen,
                      B, 0, 0);
  
}

/***************************************************************/
/* This is the worker function that implements GetEntries().   */
/***************************************************************/
HMatrix *HMatrix::DoGetEntries(int RowStart, int RowStop, int RowInc, int RowLen,
                               int ColStart, int ColStop, int ColInc, int ColLen,
                               HMatrix *B, int RowOffset, int ColOffset)
{
  /*--------------------------------------------------------------*/
  /*- sanity checks on requested row/column swaths ---------------*/
  /*--------------------------------------------------------------*/
  if (RowStart<0 || RowStop > NR)
   ErrExit("invalid row specification (%i,%i) in GetEntries",RowStart,RowStop);
  if (ColStart<0 || ColStop > NC)
   ErrExit("invalid row specification (%i,%i) in GetEntries",ColStart,ColStop);
  if (RowLen<=0)
   ErrExit("invalid row specification (%i:%i:%i) in GetEntries",RowStart,RowInc,RowStop);
  if (ColLen<=0)
   ErrExit("invalid col specification (%i:%i:%i) in GetEntries",ColStart,ColInc,ColStop);

  /*--------------------------------------------------------------*/
  /*- (re)allocate the output matrix as necessary ----------------*/
  /*--------------------------------------------------------------*/
  if ( B && ( (B->NR < RowOffset+RowLen) || (B->NC < ColOffset+ColLen) ) )
   { Warn("wrong-size X matrix in GetEntries() (reallocating)");
     delete B;
     B=0;
   };
  if (B==0)
   B=new HMatrix(RowOffset+RowLen,ColOffset+ColLen,RealComplex);

  /*--------------------------------------------------------------*/
  /*- do the extraction ------------------------------------------*/
  /*--------------------------------------------------------------*/
  for(int nr=0; nr<RowLen; nr++)
   for(int nc=0; nc<ColLen; nc++)
    B->SetEntry(RowOffset+nr, ColOffset+nc,
                GetEntry(RowStart + nr*RowInc, ColStart + nc*ColInc)
               );

  return B;
   
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void HMatrix::SetEntriesD(const char *RowString, int Col, double *Entries, double Entry)
{
  /***************************************************************/
  /* parse column string *****************************************/
  /***************************************************************/
  int RowStart, RowStop, RowInc;
  int Length=ParseMatrixIndexString(RowString, NR, &RowStart, &RowStop, &RowInc);
  if (Length<=0)
   ErrExit("invalid row index string %s",RowString);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (Entries==0)
   { 
     DoSetEntries(RowStart, RowStop, RowInc, Length, Col, Col, 
                  1, 1, 0, Entry);
   }
  else
   { HMatrix B(1, Length, LHM_REAL, LHM_NORMAL, (void *)Entries);
     DoSetEntries(RowStart, RowStop, RowInc, Length, Col, Col, 
                  1, 1, &B);
   };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void HMatrix::SetEntries(const char *RowString, int Col, cdouble *Entries, cdouble Entry)
{
  /***************************************************************/
  /* parse column string *****************************************/
  /***************************************************************/
  int RowStart, RowStop, RowInc;
  int Length=ParseMatrixIndexString(RowString, NR, &RowStart, &RowStop, &RowInc);
  if (Length<=0)
   ErrExit("invalid row index string %s",RowString);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (Entries==0)
   { 
     DoSetEntries(RowStart, RowStop, RowInc, Length, Col, Col, 
                  1, 1, 0, Entry);
   }
  else
   { HMatrix B(1, Length, LHM_COMPLEX, LHM_NORMAL, (void *)Entries);
     DoSetEntries(RowStart, RowStop, RowInc, Length, Col, Col,
                  1, 1, &B);
   };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void HMatrix::SetEntriesD(int Row, const char *ColString, double *Entries, double Entry)
{
  /***************************************************************/
  /* parse column string *****************************************/
  /***************************************************************/
  int ColStart, ColStop, ColInc;
  int Length=ParseMatrixIndexString(ColString, NC, &ColStart, &ColStop, &ColInc);
  if (Length<=0)
   ErrExit("invalid column index string %s",ColString);

  if (Entries==0)
   { DoSetEntries(Row, Row, 1, 1, ColStart, ColStop, ColInc, 
                  Length, 0, Entry);
   }
  else
   { 
     HMatrix B(1, Length, LHM_REAL, LHM_NORMAL, (void *)Entries);
     DoSetEntries(Row, Row, 1, 1, ColStart, ColStop, ColInc, 
                  Length, &B);
   };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void HMatrix::SetEntries(int Row, const char *ColString, cdouble *Entries, cdouble Entry)
{
  /***************************************************************/
  /* parse column string *****************************************/
  /***************************************************************/
  int ColStart, ColStop, ColInc;
  int Length=ParseMatrixIndexString(ColString, NC, &ColStart, &ColStop, &ColInc);
  if (Length<=0)
   ErrExit("invalid column index string %s",ColString);

  if (Entries==0)
   { 
      DoSetEntries(Row, Row, 1, 1, ColStart, ColStop, ColInc, 
                  Length, 0, Entry);
   }
  else
   {  
     HMatrix B(1, Length, LHM_COMPLEX, LHM_NORMAL, (void *)Entries);
     DoSetEntries(Row, Row, 1, 1, ColStart, ColStop, ColInc, 
                  Length, &B);
   }
}


/***************************************************************/
/* This is the worker function that implements SetEntries().   */
/***************************************************************/
void HMatrix::DoSetEntries(int RowStart, int RowStop, int RowInc, int RowLen,
                           int ColStart, int ColStop, int ColInc, int ColLen,
                           HMatrix *B, cdouble Entry)
{
  /*--------------------------------------------------------------*/
  /*- sanity checks on requested row/column swaths ---------------*/
  /*--------------------------------------------------------------*/
  if (RowStart<0 || RowStop > NR)
   ErrExit("invalid row specification (%i,%i) in GetEntries",RowStart,RowStop);
  if (ColStart<0 || ColStop > NC)
   ErrExit("invalid row specification (%i,%i) in GetEntries",ColStart,ColStop);
  if (RowLen<=0)
   ErrExit("invalid row specification (%i:%i:%i) in GetEntries",RowStart,RowInc,RowStop);
  if (ColLen<=0)
   ErrExit("invalid col specification (%i:%i:%i) in GetEntries",ColStart,ColInc,ColStop);

  if ( !B || ( (B->NR < RowLen) || (B->NC < ColLen) ) )
   ErrExit("invalid B matrix in DoSetEntries()");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  for(int nr=0; nr<RowLen; nr++)
   for(int nc=0; nc<ColLen; nc++)
    SetEntry( RowStart + nr*RowInc, ColStart + nc*ColInc,
              B ? B->GetEntry(nr, nc) : Entry
            );

}
