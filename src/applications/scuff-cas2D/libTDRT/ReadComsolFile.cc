/*
 * ReadComsolFile.cc 
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
/* Constructor helper function for reading in nodes and  *******/
/* elements for a .mphtxt file as produced by COMSOL     *******/
/***************************************************************/
void TDRTObject::ReadComsolFile(FILE *f)
{ 
  char Line[MAXSTR];
  char *s;
  int nv, ns, n1, n2, LineNum, LinesRead;
  
  LineNum=0;
 
  /***************************************************************/
  /* skip down to node definition section                        */
  /***************************************************************/
  LinesRead=SkipTo(f,"# number of mesh points",Line);
  if (LinesRead==0)
   ErrExit("%s: failed to find line '#number of mesh points'",MeshFileName);
  LineNum+=LinesRead;

  sscanf(Line,"%i ",&NumVertices);
  if (NumVertices<0 || NumVertices>10000000)
   ErrExit("%s:%i: too many vertices(%i)",MeshFileName,LineNum,NumVertices);

  LinesRead=SkipTo(f,"# Mesh point coordinates",Line);
  if (LinesRead==0)
   ErrExit("%s: failed to find line '#Mesh point coordinates'",MeshFileName);
  LineNum+=LinesRead;

  /***************************************************************/
  /* allocate space for vertices *********************************/
  /***************************************************************/
  Vertices=(double *)malloc(2*NumVertices*sizeof(double *));

  /***************************************************************/
  /* read vertices ***********************************************/
  /***************************************************************/
  for(nv=0; nv<NumVertices; nv++)
   { 
     if ( !fgets(Line,MAXSTR,f) )
      ErrExit("%s: unexpected end of file",MeshFileName);
     LineNum++;

     if ( 2!=sscanf(Line,"%le %le",Vertices + 2*nv,Vertices + 2*nv +1) )
      ErrExit("%s:%i: syntax error",MeshFileName,LineNum);
   };

  /***************************************************************/
  /* skip down to element definition section *********************/
  /***************************************************************/
  LinesRead=SkipTo(f,"2 # number of nodes per element",Line);
  if (LinesRead==0)
   ErrExit("%s: failed to find line '2 #number of nodes per element'",
                MeshFileName);
  LineNum+=LinesRead;

  s=fgets(Line,MAXSTR,f); LineNum++;
  if (1!=sscanf(Line,"%i",&NumSegments) || !strstr(Line,"# number of elements"))
   ErrExit("%s:%i: syntax error",MeshFileName,LineNum);

  s=fgets(Line,MAXSTR,f); LineNum++;
  if ( !strstr(Line,"# Elements") )
   ErrExit("%s:%i: syntax error",MeshFileName,LineNum);

  Segments=(int *)malloc(2*NumSegments*sizeof(int));

  /***************************************************************/
  /* read elements ***********************************************/
  /***************************************************************/
  ns=0;
  for(ns=0; ns<NumSegments; ns++)
   { 
     if ( !fgets(Line,MAXSTR,f) )
      ErrExit("%s: unexpected end of file",MeshFileName);
     LineNum++;
     if ( 2!=sscanf(Line,"%i %i",&n1,&n2) )
      ErrExit("%s:%i: syntax error",MeshFileName,LineNum); 
 
     Segments[2*ns]=n1;
     Segments[2*ns+1]=n2;
   };

  /***************************************************************/
  /* ignore the rest of the file and we are done *****************/
  /***************************************************************/
  fclose(f);

} 
