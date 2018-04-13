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
 * ReadComsolFile.cc 
 *
 * homer reid 6/2009
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "libscuff.h"

namespace scuff {

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
void RWGSurface::ReadComsolFile(FILE *MeshFile, char *FileName)
{ 
  char Line[MAXSTR];
  int nv, np, n1, n2, n3, LineNum, LinesRead, nConv;
  
  LineNum=0;
 
  /***************************************************************/
  /* skip down to node definition section                        */
  /***************************************************************/
  LinesRead=SkipTo(MeshFile,"# number of mesh points",Line);
  if (LinesRead==0)
   ErrExit("%s: failed to find line '#number of mesh points'",FileName);
  LineNum+=LinesRead;

  nConv=sscanf(Line,"%i ",&NumVertices);
  if (nConv!=1 || NumVertices<0 || NumVertices>10000000)
   ErrExit("%s:%i: too many vertices(%i)",FileName,LineNum,NumVertices);

  LinesRead=SkipTo(MeshFile,"# Mesh point coordinates",Line);
  if (LinesRead==0)
   ErrExit("%s: failed to find line '#Mesh point coordinates'",FileName);
  LineNum+=LinesRead;

  /***************************************************************/
  /* read vertices ***********************************************/
  /***************************************************************/
  Vertices=(double *)mallocEC(3*NumVertices*sizeof(double *));
  for(nv=0; nv<NumVertices; nv++)
   { 
     if ( !fgets(Line,MAXSTR,MeshFile) )
      ErrExit("%s: unexpected end of file",FileName);
     LineNum++;

     if ( 3!=sscanf(Line,"%le %le %le",
                    Vertices+3*nv,Vertices+3*nv+1,Vertices+3*nv+2))
      ErrExit("%s:%i: syntax error",FileName,LineNum);
   };

  /***************************************************************/
  /*- Apply geometrical transformation (if any) to all nodes.   -*/
  /***************************************************************/
  if (OTGT) OTGT->Apply(Vertices, NumVertices);

  /***************************************************************/
  /* skip down to element definition section *********************/
  /***************************************************************/
  LinesRead=SkipTo(MeshFile,"3 # number of nodes per element",Line);
  if (LinesRead==0)
   ErrExit("%s: failed to find line '3 #number of nodes per element'", FileName);
  LineNum+=LinesRead;

  if ( !fgets(Line,MAXSTR,MeshFile) )
   ErrExit("%s: unexpected end of file",FileName);
  LineNum++;
  nConv=sscanf(Line,"%i",&NumPanels);
  if (nConv!=1 || !strstr(Line,"# number of elements"))
   ErrExit("%s:%i: syntax error",FileName,LineNum);

  if ( !fgets(Line,MAXSTR,MeshFile) )
   ErrExit("%s: unexpected end of file",FileName);
  LineNum++;
  if ( !strstr(Line,"# Elements") )
   ErrExit("%s:%i: syntax error",FileName,LineNum);

  /***************************************************************/
  /* read panels    **********************************************/ 
  /***************************************************************/
  Panels=(RWGPanel **)mallocEC(NumPanels*sizeof(Panels[0]));
  for(np=0; np<NumPanels; np++)
   { 
     if ( !fgets(Line,MAXSTR,MeshFile) )
      ErrExit("%s: unexpected end of file",FileName);
     LineNum++;

     if ( 3!=sscanf(Line,"%i %i %i",&n1,&n2,&n3) )
      ErrExit("%s:%i: syntax error",FileName,LineNum); 
 
     Panels[np]=NewRWGPanel(Vertices,n1,n2,n3);
     Panels[np]->Index=np;
   };

  /***************************************************************/
  /* ignore the rest of the file and we are done *****************/
  /***************************************************************/
  fclose(MeshFile);

} 

} // namespace scuff
