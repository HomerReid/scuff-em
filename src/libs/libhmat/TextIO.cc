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
 * TextIO.cc   -- libhmat class methods for exporting/importing HMatrices 
 *             -- and HVectors to/from text files 
 *
 * homer reid  -- 12/2009
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

#include <libhrutil.h>

#include "libhmat.h"

/***************************************************************/
/***************************************************************/
/***************************************************************/
void HMatrix::ImportFromText(const char *FileName, const char *Options)
{
  FILE *f;
  char buffer[1000];

  /***************************************************************/
  /* attempt to open file ****************************************/
  /***************************************************************/
  if ( FileName==0 || !(f=fopen(FileName,"r")))
   { snprintf(buffer,1000,"could not open file %s",FileName ? FileName : "");
     ErrMsg=strdupEC(buffer);
     return;
   };

  /***************************************************************/
  /* process arguments *******************************************/
  /***************************************************************/
  char *Tokens[50];
  int nt, NumTokens, iVal;
  int nColsRequested=0, nRowsRequested=0, Strict=0;

  strncpy(buffer,Options,1000);
  NumTokens=Tokenize(buffer,Tokens,50);
  for(nt=0; nt<NumTokens; nt++)
   { 
     if ( !StrCaseCmp(Tokens[nt],"--strict") )
      Strict=1;
     else if ( !StrCaseCmp(Tokens[nt],"--ncol") )
      { if ( ++nt < NumTokens && 1==sscanf(Tokens[nt],"%i",&iVal) )
         nColsRequested=iVal;
      }
     else if ( !StrCaseCmp(Tokens[nt],"--nrow") )
      { if ( ++nt < NumTokens && 1==sscanf(Tokens[nt],"%i",&iVal) )
         nRowsRequested=iVal;
      };
   };

  /***************************************************************/
  /* first pass through file to count numbers of columns and rows*/
  /***************************************************************/
  int Line, nRows, nCols=0, nColsThisRow, AllReal;
  cdouble cdVal;

  nRows=0;
  Line=0;
  AllReal=1;
  while(fgets(buffer,1000,f))
   { 
     Line++;

     /*--------------------------------------------------------------*/
     /*- split line up into tokens ----------------------------------*/
     /*--------------------------------------------------------------*/
     NumTokens=Tokenize(buffer,Tokens,50);
     if (NumTokens==50) 
      { snprintf(buffer,1000,"%s:%i: too many columns",FileName,Line);
        ErrMsg=strdupEC(buffer);
        fclose(f);
        return;
      };

     /*--------------------------------------------------------------*/
     /*- skip blank lines and comments ------------------------------*/
     /*--------------------------------------------------------------*/
     if (NumTokens==0 || Tokens[0][0]=='#')
      continue;

     /*--------------------------------------------------------------*/
     /*- count number of tokens that can be interpreted as numbers --*/
     /*--------------------------------------------------------------*/
     nColsThisRow=0;
     for(nt=0; nt<NumTokens; nt++)
      if ( 0==S2CD(Tokens[nt],&cdVal) )
       { nColsThisRow++;
         if ( imag(cdVal)!=0.0 ) 
          AllReal=0; 
       };

     /*--------------------------------------------------------------*/
     /*- error checking:                                             */
     /*- 1. all lines (other than blank lines and comments) must     */
     /*-    have at least one column                                 */
     /*- 2. if --strict was specified, all rows must have the same   */
     /*-    number of columns                                        */
     /*- otherwise, the number of columns in the matrix will be      */
     /*- the greatest number of columns on any line.                 */
     /*--------------------------------------------------------------*/
     if (nColsThisRow==0)
      { snprintf(buffer,1000,"%s:%i: no columns found on row",FileName,Line);
        ErrMsg=strdupEC(buffer);
        fclose(f);
        return;
      };

     if (nRows==0)
      nCols=nColsThisRow;
     else if (Strict && nColsThisRow!=nCols)
      { snprintf(buffer,1000,"%s:%i: mismatch in number of data columns",FileName,Line);
        ErrMsg=strdupEC(buffer);
        fclose(f);
        return;
      }
     else
      nCols = nColsThisRow > nCols ? nColsThisRow : nCols; 

     nRows++;
   };

  /***************************************************************/
  /* if the user specified a certain number of rows and/or       */ 
  /* columns and didn't get it, that's an error                  */
  /***************************************************************/
  if (nRowsRequested && nRowsRequested!=nRows)
   { snprintf(buffer,1000,"%s: matrix has %i rows (not %i as requested)",FileName,nRows,nRowsRequested);
     ErrMsg=strdupEC(buffer);
     fclose(f);
     return;
   };
  if (nColsRequested && nColsRequested!=nCols)
   { snprintf(buffer,1000,"%s: matrix has %i columns (not %i as requested)",FileName,nCols,nColsRequested);
     ErrMsg=strdupEC(buffer);
     fclose(f);
     return;
   };

  /***************************************************************/
  /* initialize the HMatrix **************************************/
  /***************************************************************/
  InitHMatrix(nRows, nCols, AllReal ? LHM_REAL : LHM_COMPLEX);
  Zero();

  /***************************************************************/
  /* second pass to read the numbers into the matrix *************/
  /***************************************************************/
  int nr, nc;

  nr=0;
  rewind(f);
  while(fgets(buffer,1000,f))
   { 
     NumTokens=Tokenize(buffer,Tokens,50);
     if (NumTokens==0 || Tokens[0][0]=='#')
      continue;

     nc=0;
     for(nt=0; nt<NumTokens; nt++)
      if ( !S2CD(Tokens[nt],&cdVal) )
       SetEntry(nr, nc++, cdVal);
     nr++;
   };

  /***************************************************************/
  /* success *****************************************************/
  /***************************************************************/
  fclose(f);
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void HMatrix::ExportToText(const char *FileName, const char *Options)
{
  FILE *f;
  char buffer[1000];
  int nr, nc;
  cdouble Z;

  /***************************************************************/
  /* process options.                                            */
  /* for the moment the only option is "--separate," which means */
  /* that the real and imaginary parts of complex numbers are    */
  /* written as separate columns of the data file (separated by  */
  /* spaces) instead of as non-space-separated strings like 4+5I.*/
  /***************************************************************/
  char *Tokens[50];
  int nt, NumTokens;
  bool Separate=false;

  strncpy(buffer,Options,1000);
  NumTokens=Tokenize(buffer,Tokens,50);
  for(nt=0; nt<NumTokens; nt++)
   { 
     if ( !StrCaseCmp(Tokens[nt],"--separate") )
      Separate=true;
   };

  /***************************************************************/
  /* attempt to open file ****************************************/
  /***************************************************************/
  if ( !FileName || !(f=fopen(FileName,"w")))
   { snprintf(buffer,1000,"could not open file %s",FileName ? FileName : "");
     ErrMsg=strdupEC(buffer);
     return;
   };

  /***************************************************************/
  /* write the data in ascii format ******************************/
  /***************************************************************/
  if (RealComplex==LHM_REAL)
   { for(nr=0; nr<NR; nr++)
      { for(nc=0; nc<NC; nc++)
         fprintf(f,"%.15e ",GetEntryD(nr,nc));
        fprintf(f,"\n");
      };
   }
  else
   { for(nr=0; nr<NR; nr++)
      { for(nc=0; nc<NC; nc++)
         { Z=GetEntry(nr,nc);
           if (Separate)
            fprintf(f,"%.15e %.15e ",real(Z), imag(Z));
           else
            fprintf(f,"%.15e+%.15eI ",real(Z), imag(Z));
          };
        fprintf(f,"\n");
      };
   };

  fclose(f);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void HVector::ImportFromText(const char *FileName, const char *Options)
{
  FILE *f;
  char buffer[1000];

  /***************************************************************/
  /* attempt to open file ****************************************/
  /***************************************************************/
  if ( !(f=fopen(FileName,"r")))
   { snprintf(buffer,1000,"could not open file %s",FileName);
     ErrMsg=strdupEC(buffer);
     return;
   };

  /***************************************************************/
  /* process options *********************************************/
  /***************************************************************/
  char *Tokens[50];
  int nt, NumTokens, iVal;
  int nRowsRequested=0;
  int Strict=0;

  strncpy(buffer,Options,1000);
  NumTokens=Tokenize(buffer,Tokens,50);
  for(nt=0; nt<NumTokens; nt++)
   { 
     if ( !StrCaseCmp(Tokens[nt],"--strict") )
      { 
        Strict=1;
      }
     else if ( !StrCaseCmp(Tokens[nt],"--nrow") )
      { if ( ++nt < NumTokens && 1==sscanf(Tokens[nt],"%i",&iVal) )
         nRowsRequested=iVal;
      };
   };

  /***************************************************************/
  /* first pass through file to count numbers of rows            */
  /***************************************************************/
  int Line, nRows, nColsThisRow, AllReal;
  cdouble cdVal;

  nRows=0;
  Line=0;
  AllReal=1;
  while(fgets(buffer,1000,f))
   { 
     Line++;

     /*--------------------------------------------------------------*/
     /*- split line up into tokens ----------------------------------*/
     /*--------------------------------------------------------------*/
     NumTokens=Tokenize(buffer,Tokens,50);
     if (NumTokens==50) 
      { snprintf(buffer,1000,"%s:%i: too many columns",FileName,Line);
        ErrMsg=strdupEC(buffer);
        fclose(f);
        return;
      };

     /*--------------------------------------------------------------*/
     /*- skip blank lines and comments ------------------------------*/
     /*--------------------------------------------------------------*/
     if (NumTokens==0 || Tokens[0][0]=='#')
      continue;

     /*--------------------------------------------------------------*/
     /*- count number of tokens that can be interpreted as numbers --*/
     /*--------------------------------------------------------------*/
     nColsThisRow=0;
     for(nt=0; nt<NumTokens; nt++)
      if ( !S2CD(Tokens[nt],&cdVal) )
       { nColsThisRow++;
         if( imag(cdVal)!=0.0 )
          AllReal=0;
       };

     /*--------------------------------------------------------------*/
     /*- error checking:                                             */
     /*- 1. all lines (other than blank lines and comments) must     */
     /*-    have at least one column                                 */
     /*- 2. if --strict was specified, all rows must have the same   */
     /*-    number of columns, namely, 1                             */
     /*- otherwise, the first column on the line will be taken as    */
     /*- the entry for this row of the vector, and all other columns */
     /*- will be ignored                                             */
     /*--------------------------------------------------------------*/
     if (nColsThisRow==0)
      { snprintf(buffer,1000,"%s:%i: no columns found on row",FileName,Line);
        ErrMsg=strdupEC(buffer);
        fclose(f);
        return;
      }
     else if (Strict && nColsThisRow!=1)
      { snprintf(buffer,1000,"%s:%i: too many data columns",FileName,Line);
        ErrMsg=strdupEC(buffer);
        fclose(f);
        return;
      };

     nRows++;
   };

  /***************************************************************/
  /* if the user specified a certain number of rows and didn't   */ 
  /* get it, that's an error                                     */
  /***************************************************************/
  if (nRowsRequested && nRowsRequested!=nRows)
   { snprintf(buffer,1000,"%s: vector has %i rows (not %i as requested)",FileName,nRows,nRowsRequested);
     ErrMsg=strdupEC(buffer);
     fclose(f);
     return;
   };

  /***************************************************************/
  /* initialize the HVector **************************************/
  /***************************************************************/
  InitHVector(nRows, AllReal ? LHM_REAL : LHM_COMPLEX);
  Zero();

  /***************************************************************/
  /* second pass to read the numbers into the matrix *************/
  /***************************************************************/
  int n;

  n=0;
  rewind(f);
  while(fgets(buffer,1000,f))
   { 
     NumTokens=Tokenize(buffer,Tokens,50);
     if (NumTokens==0 || Tokens[0][0]=='#')
      continue;

     for(nt=0; nt<NumTokens; nt++)
      if ( !S2CD(Tokens[nt],&cdVal) )
       { SetEntry(n++,cdVal);
         break;
       };
   };
  
  /***************************************************************/
  /* success *****************************************************/
  /***************************************************************/
  fclose(f);
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void HVector::ExportToText(const char *FileName, const char *Options)
{
  FILE *f;
  char buffer[1000];
  int n;
  cdouble Z;

  (void) Options; // currently unused, silence compiler warning

  /***************************************************************/
  /* attempt to open file ****************************************/
  /***************************************************************/
  if ( !(f=fopen(FileName,"w")))
   { snprintf(buffer,1000,"could not open file %s",FileName);
     ErrMsg=strdupEC(buffer);
     return;
   };

  /***************************************************************/
  /* write the data in ascii format ******************************/
  /***************************************************************/
  if (RealComplex==LHM_REAL)
   { for(n=0; n<N; n++)
      fprintf(f,"%.15e\n ",DV[n]);
   }
  else
   { for(n=0; n<N; n++)
      fprintf(f,"%.15e+%.15eI \n",real(ZV[n]), imag(ZV[n]));
   };

  fclose(f);
}
