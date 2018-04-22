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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "libhrutil.h"
#include "libhmat.h"

/***************************************************************/
/* replacements for readline functions if necessary            */
/***************************************************************/
#ifdef HAVE_READLINE
  #include <readline/readline.h>
  #include <readline/history.h>
#else
  char *readline(const char *prompt)
   { char buffer[1000];
     printf("%s",prompt);
     (void) fgets(buffer,1000,stdin);
     return strdup(buffer);
   }
  void using_history() {}
  void read_history(const char *){}
  void add_history(const char *) {}
  void write_history(const char *) {}
#endif
/***************************************************************/
/***************************************************************/
/***************************************************************/

#define DIM 5
#define II cdouble (0.0,1.0)

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main()
{ 
  HMatrix *M=new HMatrix(DIM,DIM);

  for(int nr=0; nr<DIM; nr++)
   for(int nc=0; nc<DIM; nc++)
    { M->SetEntry(nr, nc, (double)(nr + nc*DIM) );
      printf("%+.2e %s",M->GetEntryD(nr,nc),nc==(DIM-1) ? "\n" : "");
    };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  using_history();
  read_history(0);
  for(;;)
   { 
     /*--------------------------------------------------------------*/
     /*- print prompt, fetch input string, split into tokens        -*/
     /*--------------------------------------------------------------*/
     char *p=readline("enter string: ");
     if (!p) break;
     add_history(p);
     write_history(0);
     char *Tokens[3];
     int NumTokens=Tokenize(p,Tokens,3," ");
     printf("NumTokens = %i\n",NumTokens);

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     int Row, Col;
     if (NumTokens==1)
      { 
        HMatrix *B=M->ExtractEntries(Tokens[0], 0);
        for(int nr=0; nr<B->NR; nr++)
         for(int nc=0; nc<B->NC; nc++)
          printf("%+.2e %s",B->GetEntryD(nr,nc),(nc==(B->NC-1)) ? "\n" : "");
        printf("\n\n");
      }
     else if (NumTokens==2 && 1==sscanf(Tokens[0],"%i",&Row) )
      { 
        double D[DIM];
        M->GetEntriesD(Row,Tokens[1],D);
        for(int n=0; n<DIM; n++)
         printf("%+.2e ",D[n]);
        printf("\n\n");
      }
     else if (NumTokens==2 && 1==sscanf(Tokens[1],"%i",&Col) )
      { 
        double D[DIM];
        M->GetEntriesD(Tokens[0],Col,D);
        for(int n=0; n<DIM; n++)
         printf("%+.2e ",D[n]);
        printf("\n\n");
      };

   };

}
