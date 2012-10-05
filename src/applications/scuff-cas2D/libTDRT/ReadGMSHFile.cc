/*
 * ReadGMSHFile.cc 
 *
 * homer reid 11/2008
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "libTDRT.h"
#include "libhrutil.h"

#define MAXSTR 1000

/***************************************************************/
/* Given an open file f, read lines until one is found that    */
/* contains the given search string. Returns the number of     */
/* lines read or 0 if the file ended before the string was     */
/* found. If the return value is nonzero, Line contains the    */
/* last line read from the file.                               */
/***************************************************************/
static int SkipTo(FILE *f, const char *SearchString, char *Line)
{ 
  int Lines=0;

  while( fgets(Line,MAXSTR,f) )
   { Lines++;
     if (strstr(Line,SearchString))
      return Lines;
   };
  return 0;
}

/***************************************************************/
/* Constructor helper function for reading in vertices and  ****/
/* elements for a .msh file as produced by GMSH          *******/
/***************************************************************/
void TDRTObject::ReadGMSHFile(FILE *f)
{ 
  char Line[MAXSTR];
  int nv, ns, n1, n2, nt, nTags, Type, iDum, NumElements;
  double dDum;
  int LineNum, LinesRead, strPos, strPosInc;
  char *s;
  
  LineNum=0;
 
  /***************************************************************/
  /* skip down to node definition section                        */
  /***************************************************************/
  LinesRead=SkipTo(f,"$Nodes",Line);
  if (LinesRead==0)
   ErrExit("%s: failed to find line '$Nodes'",MeshFileName);
  LineNum+=LinesRead;

  Line[0]=0; s=fgets(Line,MAXSTR,f); LineNum++;
  if (1!=sscanf(Line,"%i",&NumVertices) )
   ErrExit("%s:%i: syntax error",MeshFileName,LineNum);
  if (NumVertices<0 || NumVertices>10000000)
   ErrExit("%s:%i: too many vertices (%i)",MeshFileName,LineNum,NumVertices);

  /***************************************************************/
  /* allocate space for vertices ************************************/
  /***************************************************************/
  Vertices=(double *)malloc(2*NumVertices*sizeof(double *));

  /***************************************************************/
  /* read vertices. we ignore the Z coordinates and the node     */
  /* indices (we assume that the first node is always node 1 and */
  /* that the vertices are ordered in sequence thereafter)       */
  /***************************************************************/
  nv=0;
  for(nv=0; nv<NumVertices; nv++)
   { 
     if ( !fgets(Line,MAXSTR,f) )
      ErrExit("%s: unexpected end of file",MeshFileName);
     LineNum++;

     if ( 4!=sscanf(Line,"%i %le %le %le",
                     &iDum,Vertices + 2*nv,Vertices + 2*nv+1,&dDum) )
      ErrExit("%s:%i: syntax error",MeshFileName,LineNum);
   };

  /***************************************************************/
  /* skip down to element definition section *********************/
  /***************************************************************/
  LinesRead=SkipTo(f,"$Elements",Line);
  if (LinesRead==0)
   ErrExit("%s: failed to find line '$Elements'", MeshFileName);
  LineNum+=LinesRead;

  Line[0]=0; s=fgets(Line,MAXSTR,f); LineNum++;
  if (1!=sscanf(Line,"%i",&NumElements) )
   ErrExit("%s:%i: syntax error",MeshFileName,LineNum);

  Segments=(int *)malloc(2*NumElements*sizeof(int));

  /***************************************************************/
  /* read elements. **********************************************/
  /* note 'NumElements' is the actual number of elements (of all */
  /* types) stated to be present in the file, as distinct from   */
  /* uNumSegments' which is the number of line elements.         */
  /***************************************************************/
  ns=NumSegments=0;
  for(ns=0; ns<NumElements; ns++)
   { 
     if ( !fgets(Line,MAXSTR,f) )
      ErrExit("%s: unexpected end of file",MeshFileName);
     LineNum++;
     if ( 3!=sscanf(Line,"%i %i %i %n",&iDum,&Type,&nTags,&strPos) )
      ErrExit("%s:%i: syntax error",MeshFileName,LineNum); 

     if (Type!=1) continue; /* skip everything but lines */

     for(nt=0; nt<nTags; nt++)
      { if ( 1!=sscanf(Line+strPos,"%i %n",&iDum,&strPosInc) )
         ErrExit("%s:%i: syntax error",MeshFileName,LineNum); 
        strPos+=strPosInc;
      };
  
     if ( 2!=sscanf(Line+strPos,"%i %i",&n1,&n2))
      ErrExit("%s:%i: syntax error",MeshFileName,LineNum);

     Segments[2*NumSegments]=n1-1;
     Segments[2*NumSegments+1]=n2-1;
     NumSegments++;
   };

  /***************************************************************/
  /* ignore the rest of the file and we are done *****************/
  /***************************************************************/
  fclose(f);

}
