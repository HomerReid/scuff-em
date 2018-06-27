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
   }
  return 0;
}

/***************************************************************/
/* Constructor helper function for reading in nodes and  *******/
/* elements for a .mphtxt file as produced by COMSOL     *******/
/***************************************************************/
char *ReadComsolFile(FILE *MeshFile, char *FileName, int MeshTag, 
                     dVec &VertexCoordinates, iVec &PanelVertexIndices)
{ 
  char Line[MAXSTR];
  int nv, np, n1, n2, n3, LineNum, LinesRead, nConv;
  LineNum=0;

  if (MeshTag!=1) 
   Warn("mesh tagging not supported for COMSOL files (ignoring mesh tag %i)",MeshTag); 
 
  /***************************************************************/
  /* skip down to node definition section                        */
  /***************************************************************/
  LinesRead=SkipTo(MeshFile,"# number of mesh points",Line);
  if (LinesRead==0)
   return vstrdup("%s: failed to find line '#number of mesh points'",FileName);
  LineNum+=LinesRead;

  int NumVertices;
  nConv=sscanf(Line,"%i ",&NumVertices);
  if (nConv!=1 || NumVertices<0 || NumVertices>10000000)
   return vstrdup("%s:%i: too many vertices(%i)",FileName,LineNum,NumVertices);

  LinesRead=SkipTo(MeshFile,"# Mesh point coordinates",Line);
  if (LinesRead==0)
   return vstrdup("%s: failed to find line '#Mesh point coordinates'",FileName);
  LineNum+=LinesRead;

  /***************************************************************/
  /* read vertices ***********************************************/
  /***************************************************************/
  VertexCoordinates.clear();
  for(nv=0; nv<NumVertices; nv++)
   { 
     if ( !fgets(Line,MAXSTR,MeshFile) )
      return vstrdup("%s: unexpected end of file",FileName);
     LineNum++;

     double V[3];
     if ( 3!=sscanf(Line,"%le %le %le",V+0, V+1, V+2) )
      return vstrdup("%s:%i: syntax error",FileName,LineNum);
     VertexCoordinates.push_back(V[0]);
     VertexCoordinates.push_back(V[1]);
     VertexCoordinates.push_back(V[2]);
   }
  if (VertexCoordinates.size() != 3*NumVertices)
   return vstrdup("%s: stated number of vertices (%i) disagrees with actual (%i)",NumVertices,VertexCoordinates.size()/3);

  /***************************************************************/
  /* skip down to element definition section *********************/
  /***************************************************************/
  LinesRead=SkipTo(MeshFile,"3 # number of nodes per element",Line);
  if (LinesRead==0)
   return vstrdup("%s: failed to find line '3 #number of nodes per element'", FileName);
  LineNum+=LinesRead;

  if ( !fgets(Line,MAXSTR,MeshFile) )
   return vstrdup("%s: unexpected end of file",FileName);
  LineNum++;
  int NumPanels;
  nConv=sscanf(Line,"%i",&NumPanels);
  if (nConv!=1 || !strstr(Line,"# number of elements"))
   return vstrdup("%s:%i: syntax error",FileName,LineNum);

  if ( !fgets(Line,MAXSTR,MeshFile) )
   return vstrdup("%s: unexpected end of file",FileName);
  LineNum++;
  if ( !strstr(Line,"# Elements") )
   return vstrdup("%s:%i: syntax error",FileName,LineNum);

  /***************************************************************/
  /* read panels    **********************************************/ 
  /***************************************************************/
  PanelVertexIndices.clear();
  for(int np=0; np<NumPanels; np++)
   { 
     if ( !fgets(Line,MAXSTR,MeshFile) )
      return vstrdup("%s: unexpected end of file",FileName);
     LineNum++;

     int VI[3];
     if ( 3!=sscanf(Line,"%i %i %i",VI+0,VI+1,VI+2) )
      return vstrdup("%s:%i: syntax error",FileName,LineNum); 
 
     PanelVertexIndices.push_back(VI[0]);
     PanelVertexIndices.push_back(VI[1]);
     PanelVertexIndices.push_back(VI[2]);
   }
  if (PanelVertexIndices.size() != 3*NumPanels)
   return vstrdup("%s: stated number of panels (%i) disagrees with actual (%i)",NumPanels,PanelVetexIndices.size()/3);

  /***************************************************************/
  /* ignore the rest of the file and we are done *****************/
  /***************************************************************/
  fclose(MeshFile);
  return 0;
} 

} // namespace scuff
