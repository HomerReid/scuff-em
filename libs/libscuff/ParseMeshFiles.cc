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
 * ReadGMSHFile.cc -- subroutine of the RWGSurface class constructor
 *
 * homer reid    -- 3/2007 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "libscuff.h"

namespace scuff {

/*************************************************************/
/* constants needed in this file only  ***********************/
/*************************************************************/
#define TYPE_TRIANGLE 2
#define TYPE_POINT    15

#define NODE_START_KEYWORD1 "$NOD"
#define NODE_START_KEYWORD2 "$Nodes"

#define NODE_END_KEYWORD1   "$ENDNOD"
#define NODE_END_KEYWORD2   "$EndNodes"

#define ELM_START_KEYWORD1  "$ELM"
#define ELM_START_KEYWORD2  "$Elements"

#define ELM_END_KEYWORD1    "$ENDELM"
#define ELM_END_KEYWORD2    "$EndElements"

#define FORMAT_LEGACY 0
#define FORMAT_NEW    1

#define MAXREFPTS 100

// vertices that are within a distance of
// PIXELSIZE of each other are considered equivalent
#define PIXELSIZE 1.0e-6

/*************************************************************/
/* Read vertices and panels from a GMSH .msh file to specify */
/* a surface.                                                */
/* If MeshTag==-1, then all panels are read.                 */
/* Otherwise, only panels on the specified physical region   */
/* are read.                                                 */
/*************************************************************/
char *ParseGMSHFile(FILE *MeshFile, char *FileName, int MeshTag,
                    dVec &VertexCoordinates, iVec &PanelVertexIndices)
{
  /*------------------------------------------------------------*/
  /* Read lines until we hit the keyword indicating the start   */
  /* of the 'nodes' section, then read the number of nodes.     */
  /*------------------------------------------------------------*/
  int LineNum=0;
  char Line[100];
  bool KeywordFound=false;
  int WhichMeshFormat;
  while(!KeywordFound)
   { if (!fgets(Line,100,MeshFile))
      return vstrdup("%s: failed to find node start keyword");
     LineNum++;
     if( !strncmp(Line,NODE_START_KEYWORD1,strlen(NODE_START_KEYWORD1)))
      { WhichMeshFormat=FORMAT_LEGACY; KeywordFound=true; }
     else if( !strncmp(Line,NODE_START_KEYWORD2,strlen(NODE_START_KEYWORD2)))
      { WhichMeshFormat=FORMAT_NEW; KeywordFound=true; }
   }
  int NumNodes;
  if ( !fgets(Line,100,MeshFile) || !(NumNodes=atoi(Line)) )
   return vstrdup("%s:%i: invalid number of nodes",FileName,LineNum);

  /*------------------------------------------------------------*/
  /*- Read in the vertices (which GMSH calls 'nodes.')          */
  /*- Note that the numbering of the vertices in GMSH does not  */
  /*- necessarily correspond to their ordering in the mesh      */
  /*- file. To remedy this situation, we construct a mapping    */
  /*- between GMSH's vertices indices and our internal vertex   */
  /*- indices, which works like this: The vertex that GMSH      */
  /*- calls 'node 3' is stored in slot GMSH2HR[3] within our    */
  /*- internal Vertices array.                                  */ 
  /*------------------------------------------------------------*/
  VertexCoordinates.clear();
  int *GMSH2HR=(int *)mallocEC(2*NumNodes*sizeof(int));
  for(int i=0; i<2*NumNodes; i++)
   GMSH2HR[i]=-1;
  for (int nn=0; nn<NumNodes; nn++)
   { if (!fgets(Line,100,MeshFile))
      return vstrdup("%s: too few nodes",FileName);
     LineNum++;
     int NodeIndex;
     double V[3];
     int nConv=sscanf(Line,"%i %le %le %le",&NodeIndex,V+0,V+1,V+2);
     if(nConv!=4)
      return vstrdup("%s:%i: invalid node specification",FileName,LineNum); 
     if (NodeIndex>2*NumNodes)
      return vstrdup("%s:%i: internal error in ReadGMSHFile",FileName,LineNum);
     VertexCoordinates.push_back(V[0]);
     VertexCoordinates.push_back(V[1]);
     VertexCoordinates.push_back(V[2]);
     GMSH2HR[NodeIndex]=nn;
   } /* for (nv=0; nv<NumVertices; nv++) */
   
  /*------------------------------------------------------------*/
  /*- 20151119 -------------------------------------------------*/
  /*------------------------------------------------------------*/
  char *s=getenv("SCUFF_PIXEL_SIZE");
  if (s)
   { double PixelSize;
     if (1!=sscanf(s,"%le",&PixelSize))
      Log("Invalid specification for SCUFF_PIXEL_SIZE (ignoring)");
     else
      Log("Rounding all vertex coordinates to be an integer multiple of %e", PixelSize);      
     for(int nn=0; nn<3*NumNodes; nn++)
      VertexCoordinates[nn] = PixelSize*round(VertexCoordinates[nn]/PixelSize);
   }
 
  /*------------------------------------------------------------*/
  /*- Eliminate any redundant vertices from the vertex list.   -*/
  /*- Silly O(N^2) algorithm used; obviously O(N) algorithm possible. */
  /*------------------------------------------------------------*/
  int NumRedundantNodes=0;
  double *Vertices = &(VertexCoordinates[0]);
  for(int i=0; i<NumNodes; i++)
   for(int j=i+1; j<NumNodes; j++)
    if ( VecDistance(Vertices+3*i, Vertices+3*j) < 1.0e-6 )
     {
       /* remap all references to my node j so that they now refer to my node i*/
       for(int jGMSH=0; jGMSH<2*NumNodes; jGMSH++)
        if (GMSH2HR[jGMSH]==j)
         GMSH2HR[jGMSH]=i;
       NumRedundantNodes++;
       //fprintf(stderr,"\n*\n* redundant nodes found!!(%i,%i)\n*\n",i,j);
     }
 
  /*------------------------------------------------------------*/
  /* Confirm that the next two lines in the file are the        */
  /* end-of-nodes-section keyword and the start-of-elements-section */
  /* keyword, then read the number of elements.                 */
  /*------------------------------------------------------------*/
  if ( !fgets(Line,100,MeshFile) )
   return vstrdup("%s: bad file format (nodes section not terminated)",FileName);
  LineNum++;

  if ( WhichMeshFormat==FORMAT_LEGACY ) 
   { if ( strncmp(Line,NODE_END_KEYWORD1,strlen(NODE_END_KEYWORD1)))
      return vstrdup("%s:%i: unexpected keyword",FileName,LineNum);
   }
  else
   { if ( strncmp(Line,NODE_END_KEYWORD2,strlen(NODE_END_KEYWORD2)))
      return vstrdup("%s:%i: unexpected keyword",FileName,LineNum);
   }

  if ( !fgets(Line,100,MeshFile) )
   return vstrdup("%s: bad file format (elements section not initiated)",FileName);
  LineNum++;

  if ( WhichMeshFormat==FORMAT_LEGACY ) 
   { if ( strncmp(Line,ELM_START_KEYWORD1,strlen(ELM_START_KEYWORD1)))
      return vstrdup("%s:%i: unexpected keyword",FileName,LineNum);
   }
  else
   { if ( strncmp(Line,ELM_START_KEYWORD2,strlen(ELM_START_KEYWORD2)))
      return vstrdup("%s:%i: unexpected keyword",FileName,LineNum);
   }

  if ( !fgets(Line,100,MeshFile) )
   return vstrdup("%s: bad file format (invalid number of elements)",FileName);
  LineNum++;
  int NumElements;
  if ( 1!=sscanf(Line,"%i",&NumElements) || NumElements<0 ) 
   return vstrdup("%s:%i: invalid number of elements",FileName,LineNum);
 
  /*------------------------------------------------------------*/
  /*- Now read each line of the elements section.               */ 
  /*------------------------------------------------------------*/
  PanelVertexIndices.clear();
  for (int ne=0; ne<NumElements; ne++)
   { 
     if (!fgets(Line,100,MeshFile))
      return vstrdup("too few elements in input file");
     LineNum++;

     int VI[3], ElType, RegPhys;
     if (WhichMeshFormat==FORMAT_LEGACY)
      { int ElNum, RegElem, NodeCnt;
        int nConv=sscanf(Line,"%i %i %i %i %i %i %i %i", &ElNum,&ElType,&RegPhys,&RegElem,&NodeCnt,VI,VI+1,VI+2);
        if (nConv<5) return vstrdup("%s:%i: invalid element specification",FileName,LineNum);
      }
     else
      { int ElNum, nTags, nRead;
        int nConv=sscanf(Line,"%i %i %i%n",&ElNum,&ElType,&nTags,&nRead);
        int bufPos=nRead;
        if (nConv<3) return vstrdup("%s:%i: invalid element specification",FileName,LineNum);

        // read the first 'tag,' which should be the physical region
        if (nTags==0)
         RegPhys=0;
        else 
         sscanf(Line+bufPos,"%i%n",&RegPhys,&nRead);
        bufPos+=nRead;

        // read any remaining tags 
        int iDummy;
        for(int nt=0; nt<nTags-1; nt++)
         { sscanf(Line+bufPos,"%i%n",&iDummy,&nRead);
           bufPos+=nRead;
         }

        // finally, read the vertex indices. 
        nConv=sscanf(Line+bufPos,"%i %i %i",VI,VI+1,VI+2);

        if (ElType==TYPE_TRIANGLE && nConv!=3)
         return vstrdup("%s:%i: invalid element specification",FileName,LineNum);
        else if (ElType==TYPE_POINT && nConv!=1)
         return vstrdup("%s:%i: invalid element specification",FileName,LineNum);
      }
   
     /* we only process elements that are triangles */
     if( ElType==TYPE_TRIANGLE && ( (MeshTag == -1) || (MeshTag==RegPhys) ) )
      { PanelVertexIndices.push_back(GMSH2HR[VI[0]]);
        PanelVertexIndices.push_back(GMSH2HR[VI[1]]);
        PanelVertexIndices.push_back(GMSH2HR[VI[2]]);
      }
   } //for (ne=0; ne<NumElements; ne++)

  free(GMSH2HR);
  fclose(MeshFile);
  return 0;
} 

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
char *ParseComsolFile(FILE *MeshFile, char *FileName, int MeshTag, 
                      dVec &VertexCoordinates, iVec &PanelVertexIndices)
{ 
  if (MeshTag!=1) 
   Warn("mesh tagging not supported for COMSOL files (ignoring mesh tag %i)",MeshTag); 
 
  /***************************************************************/
  /* skip down to node definition section                        */
  /***************************************************************/
  char Line[MAXSTR];
  int LineNum=0;
  int LinesRead=SkipTo(MeshFile,"# number of mesh points",Line);
  if (LinesRead==0)
   return vstrdup("%s: failed to find line '#number of mesh points'",FileName);
  LineNum+=LinesRead;

  int NumVertices;
  if ( 1!=sscanf(Line,"%i ",&NumVertices) || NumVertices<0) 
   return vstrdup("%s:%i: invalid number of vertices (%i)",FileName,LineNum,NumVertices);

  LinesRead=SkipTo(MeshFile,"# Mesh point coordinates",Line);
  if (LinesRead==0)
   return vstrdup("%s: failed to find line '#Mesh point coordinates'",FileName);
  LineNum+=LinesRead;

  /***************************************************************/
  /* read vertices ***********************************************/
  /***************************************************************/
  VertexCoordinates.clear();
  for(int nv=0; nv<NumVertices; nv++)
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
  if (VertexCoordinates.size() != 3*((size_t)NumVertices) )
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
  int nConv=sscanf(Line,"%i",&NumPanels);
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
  if (PanelVertexIndices.size() != 3*((size_t)NumPanels) )
   return vstrdup("%s: stated number of panels (%i) disagrees with actual (%lu)",NumPanels,PanelVertexIndices.size()/3);

  /***************************************************************/
  /* ignore the rest of the file and we are done *****************/
  /***************************************************************/
  fclose(MeshFile);
  return 0;
} 

} // namespace scuff
