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
 * RWGSurface.cc -- implementation of some methods in the RWGSurface
 *               -- class (this is the class formerly known as RWGObject)
 *
 * homer reid    -- 3/2007
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <ctype.h>

#include <libhrutil.h>

#include "libscuff.h"
#include "cmatheval.h"

namespace scuff {

#define MAXSTR 1000
#define MAXTOK 50  

/*--------------------------------------------------------------*/
/*-  RWGSurface class constructor that begins reading an open  -*/
/*-  .scuffgeo file immediately after a line like either       -*/ 
/*   OBJECT MyObjectLabel or SURFACE MySurfaceLabel.           -*/
/*-                                                            -*/
/*-  (Keyword is either "OBJECT" or "SURFACE.")                -*/
/*-                                                            -*/
/*-  On return, if an object was successfully created, its     -*/
/*-   ErrMsg field is NULL, and the file read pointer points   -*/
/*-   the line immediately following the ENDOBJECT/ENDSURFACE  -*/
/*-   line.                                                    -*/
/*-                                                            -*/
/*-  If there was an error in parsing the section, the ErrMsg  -*/
/*-   field points to an error message.                        -*/
/*--------------------------------------------------------------*/
RWGSurface::RWGSurface(FILE *f, const char *pLabel, int *LineNum, char *Keyword)
{ 
  ErrMsg=0;
  SurfaceSigma=0;
  PhysicalRegion=-1;
  MeshFileName=0;
  Label = strdup(pLabel);

  if ( !strcasecmp(Keyword, "OBJECT") )
   IsObject=1;
  else
   IsObject=0;

  /***************************************************************/
  /* read lines from the file one at a time **********************/
  /***************************************************************/
  char Line[MAXSTR], LineCopy[MAXSTR];
  char Region1[MAXSTR], Region2[MAXSTR];
  int NumTokens, TokensConsumed;
  char *Tokens[MAXTOK];
  int ReachedTheEnd=0;
  GTransformation *OTGT=0; // 'one-time geometrical transformation'
  MaterialName=0;
  RegionLabels[0]=0;
  while ( ReachedTheEnd==0 && fgets(Line, MAXSTR, f) )
   { 
     /*--------------------------------------------------------------*/
     /*- skip blank lines and comments ------------------------------*/
     /*--------------------------------------------------------------*/
     (*LineNum)++;
     strcpy(LineCopy,Line);
     NumTokens=Tokenize(LineCopy, Tokens, MAXTOK);
     if ( NumTokens==0 || Tokens[0][0]=='#' )
      continue; 

     /*--------------------------------------------------------------*/
     /*- switch off based on first token on the line ----------------*/
     /*--------------------------------------------------------------*/
     if ( !strcasecmp(Tokens[0],"MESHFILE") )
      { if (NumTokens!=2)
         { ErrMsg=strdup("MESHFILE keyword requires one argument");
           return;
         };
        MeshFileName=strdup(Tokens[1]);
      }
     else if ( !strcasecmp(Tokens[0],"PHYSICAL_REGION") )
      { if (NumTokens!=2)
         { ErrMsg=strdup("PHYSICAL_REGION keyword requires one argument");
           return;
         };
        sscanf(Tokens[1],"%i",&PhysicalRegion);
      }
     else if ( !strcasecmp(Tokens[0],"MATERIAL") )
      { 
        if (!IsObject)
         { ErrMsg=strdup("MATERIAL keyword may not be used in SURFACE...ENDSURFACE sections");
           return;
         };
        if (NumTokens!=2)
         { ErrMsg=strdup("MATERIAL keyword requires one argument");
           return;
         };
        MaterialName = strdup(Tokens[1]);
        RegionLabels[0] = strdup("EXTERIOR");
        RegionLabels[1] = strdup(Label);
        MaterialRegionsLineNum = *LineNum;
        if ( !strcasecmp(MaterialName,"PEC") )
         IsPEC=1;
      }
     else if ( !strcasecmp(Tokens[0],"REGIONS") )
      { 
        if (IsObject)
         { ErrMsg=strdup("REGIONS keyword may not be used in OBJECT...ENDOBJECT sections");
           return;
         };
        if (NumTokens==3)
         { IsPEC=0;
           RegionLabels[0] = strdup(Tokens[1]);
           RegionLabels[1] = strdup(Tokens[2]);
           MaterialRegionsLineNum = *LineNum;
         }
        else if (NumTokens==2)
         { IsPEC=1;
           RegionLabels[0] = strdup(Tokens[1]);
           RegionLabels[1] = strdup("PEC");
           MaterialRegionsLineNum = *LineNum;
         }
        else
         ErrMsg=strdup("REGIONS keyword requires one or two arguments");
      }
     else if ( !strcasecmp(Tokens[0],"DISPLACED") || !strcasecmp(Tokens[0],"ROTATED") )
      { 
        // try to parse the line as a geometrical transformation.
        // note that OTGT is used as a running GTransformation that may
        // be augmented by multiple DISPLACED ... and/or ROTATED ...
        // lines within the OBJECT...ENDOBJECT or SURFACE...ENDSURFACE section, 
        // and which is applied to the object at its birth and subsequently 
        // discarded.
        // in particular, OTGT is NOT stored as the 'GT' field inside 
        // the RWGSurface class, which is intended to be used for
        // transformations that are applied and later un-applied 
        // during the life of the object/surface.
        if (OTGT)
	 { GTransformation *OTGT2 = new GTransformation(Tokens, NumTokens, &ErrMsg, &TokensConsumed);
           if (ErrMsg)
            return;

           OTGT->Transform(OTGT2);
           delete OTGT2;
 // 20120701 why is this not the right thing to do? apparently it isn't, but i dunno why.
 //          OTGT2->Transform(OTGT);
 //          delete OTGT;
 //          OTGT=OTGT2;
	 }
        else
	 { OTGT = new GTransformation(Tokens, NumTokens, &ErrMsg, &TokensConsumed);
           if (ErrMsg)
            return;
	 };
	 
        if (TokensConsumed!=NumTokens) 
         { ErrMsg=strdup("junk at end of line");
           return;
         };
      }
     else if ( !strcasecmp(Tokens[0],"SURFACE_CONDUCTIVITY") )
      { 
        if (NumTokens<2)
         { ErrMsg=strdup("no argument specified for SURFACE_CONDUCTIVITY");
           return;
         };
;
        /* try to create a cevaluator for the user's function */
        char SigmaString[MAXSTR];
        strncpy(SigmaString,Line + strlen(Tokens[0]) + 1,MAXSTR);
        SigmaString[strlen(SigmaString)-1]=0; // remove trailing newline
        SurfaceSigma=cevaluator_create(SigmaString);
        if (SurfaceSigma==0)
         { ErrMsg=vstrdup("invalid SURFACE_CONDUCTIVITY specification");
           return;
         };
      }
     else if (   !strcasecmp(Tokens[0],"ENDOBJECT") || !strcasecmp(Tokens[0],"ENDSURFACE") )
      { 
        ReachedTheEnd=1;
      }
     else
      { ErrMsg=vstrdup("unknown keyword %s in %s section",Tokens[0],Keyword);
        return;
      };
   }; 

  if (MeshFileName==0)
   { ErrMsg=vstrdup("%s section must include a MESHFILE specification",Keyword,Tokens[0]);
     return;
   };

  if (SurfaceSigma!=0)
   Log("Surface %s has surface conductivity Sigma=%s.\n",Label,cevaluator_get_string(SurfaceSigma));

  // if we are an OBJECT and there was no MATERIAL specification,
  // or we are a SURFACE and there was no REGIONS specification, then 
  // we consider ourselves to be a PEC surface in the exterior medium.
  if ( RegionLabels[0]==0 )
   { IsPEC=1;
     RegionLabels[0] = strdup("EXTERIOR");
     RegionLabels[1] = strdup("PEC");
   };

  // ok, all checks passed, now on to the main body of the class constructor.
  InitRWGSurface(OTGT);
  
}

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- Main body of RWGSurface constructor: Create an RWGSurface   */
/*- from a mesh file describing a discretized object,           */
/*- optionally with a rotation and/or displacement applied.     */
/*-                                                             */
/*- This routine assumes that the following internal class      */
/*- fields have been initialized: MeshFileName, PhysicalRegion, */
/*- Label, SurfaceSigma, and IsPEC.                             */
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
void RWGSurface::InitRWGSurface(const GTransformation *OTGT)
{ 
  ErrMsg=0;
  kdPanels = NULL;

  /*------------------------------------------------------------*/
  /*- try to open the mesh file.                                */
  /*------------------------------------------------------------*/
  FILE *MeshFile=fopen(MeshFileName,"r");
  if (!MeshFile)
   ErrExit("could not open file %s",MeshFileName);
   
  /*------------------------------------------------------------*/
  /*- initialize simple fields ---------------------------------*/
  /*------------------------------------------------------------*/
  NumEdges=NumPanels=NumVertices=NumRefPts=0;

  /*------------------------------------------------------------*/
  /*- note: the 'OTGT' parameter to this function is distinct   */
  /*- from the 'GT' field inside the class body. the former is  */
  /*- an optional 'One-Time Geometrical Transformation' to be   */
  /*- applied to the object once at its creation. the latter    */
  /*- is designed to store a subsequent transformation that may */
  /*- be applied to the surface , and is initialized to zero.   */
  /*------------------------------------------------------------*/
  GT=0;

  /*------------------------------------------------------------*/
  /*- Switch off based on the file type to read the mesh file:  */
  /*-  1. file extension=.msh    --> ReadGMSHFile              -*/
  /*-  2. file extension=.mphtxt --> ReadComsolFile            -*/
  /*------------------------------------------------------------*/
  char *p=GetFileExtension(MeshFileName);
  if (!p)
   ErrExit("file %s: invalid extension",MeshFileName);
  else if (!strcasecmp(p,"msh"))
   ReadGMSHFile(MeshFile,MeshFileName,OTGT,PhysicalRegion);
  else if (!strcasecmp(p,"mphtxt"))
   { if ( PhysicalRegion != -1 )
      ErrExit("PHYSICAL_REGION is not yet implemented for .mphtxt files");
     ReadComsolFile(MeshFile,MeshFileName,OTGT);
   }
  else
   ErrExit("file %s: unknown extension %s",MeshFileName,p);

  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  if (NumPanels==0)
   { if ( PhysicalRegion == -1 ) 
      ErrExit("file %s: no panels found",MeshFileName);
     else
      ErrExit("file %s: no panels found on physical region %i",MeshFileName,PhysicalRegion);
   };

  /*------------------------------------------------------------*/
  /*- Now that we have put the panels in an array, go through  -*/
  /*- and fill in the Index field of each panel structure.     -*/
  /*------------------------------------------------------------*/
  int np;
  for(np=0; np<NumPanels; np++)
   Panels[np]->Index=np;
 
  /*------------------------------------------------------------*/
  /* gather necessary edge connectivity info. this is           */
  /* complicated enough to warrant its own separate routine.    */
  /*------------------------------------------------------------*/
  InitEdgeList();

  /*------------------------------------------------------------*/
  /*- By default, if we are an OBJECT we do not assign half-RWG */
  /*- basis functions to exterior edges, whereas if we are a    */
  /*- SURFACE we do.                                            */
  /*- The procedure for assigning half-RWG basis functions to   */
  /*- exterior edges is simply to append the contents of the    */
  /*- ExteriorEdges array to the end of the Edges array. The    */
  /*- content of each individual edge structure is already      */
  /*- properly initialized by the InitEdgeList() routine (with  */
  /*- the exception of the Index field.) The edges remain in    */
  /*- place within the ExteriorEdges array.                     */
  /*------------------------------------------------------------*/
  if ( !IsObject && RWGGeometry::AssignBasisFunctionsToExteriorEdges)
   { Edges = (RWGEdge **)realloc(Edges, (NumEdges + NumExteriorEdges)*sizeof(RWGEdge *));
     memcpy(Edges + NumEdges, ExteriorEdges, NumExteriorEdges * sizeof(RWGEdge *));
     for(int ne=NumEdges; ne<NumEdges + NumExteriorEdges; ne++)
      Edges[ne]->Index = ne;
     NumEdges += NumExteriorEdges;
     Log("Promoted %i exterior edges for surface %s to half-RWG basis functions.",NumExteriorEdges,Label);
   };

  /*------------------------------------------------------------*/
  /*- the number of basis functions that this object takes up   */
  /*- in the full BEM system is the number of internal edges    */
  /*- if it is a perfect electrical conductor, or 2 times the   */
  /*- number of internal edges otherwise.                       */
  /*------------------------------------------------------------*/
  NumBFs = IsPEC ? NumEdges : 2*NumEdges;
  IsClosed = (NumExteriorEdges == 0);

} 

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- Alternative RWGSurface constructor: Create an RWGSurface    */
/*- from lists of vertices and panels.                          */
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
RWGSurface::RWGSurface(double *pVertices, int pNumVertices,
                       int *PanelVertices, int pNumPanels)
{ 
  int np;

  ErrMsg=0;
  kdPanels = NULL;

  MeshFileName=strdup("ByHand.msh");
  Label=strdup("ByHand");

  NumEdges=0;
  NumVertices=pNumVertices;
  NumPanels=pNumPanels;
  RegionLabels[0]=strdup("Exterior");
  RegionLabels[1]=strdup("PEC");
  RegionIndices[0]=0;
  RegionIndices[1]=-1;
  IsPEC=1;
  IsObject=1;
  GT=0;

  Vertices=(double *)mallocEC(3*NumVertices*sizeof(double));
  memcpy(Vertices,pVertices,3*NumVertices*sizeof(double *));

  /*------------------------------------------------------------*/
  /*- add the panels -------------------------------------------*/
  /*------------------------------------------------------------*/
  Panels=(RWGPanel **)mallocEC(NumPanels*sizeof(RWGPanel *));
  for(np=0; np<NumPanels; np++)
   { Panels[np]=NewRWGPanel(Vertices,
                            PanelVertices[3*np+0],
                            PanelVertices[3*np+1],
                            PanelVertices[3*np+2]);
     Panels[np]->Index=np;
   };
 
  /*------------------------------------------------------------*/
  /* gather necessary edge connectivity info                   -*/
  /*------------------------------------------------------------*/
  InitEdgeList();

  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  if (RWGGeometry::AssignBasisFunctionsToExteriorEdges)
   { Edges = (RWGEdge **)realloc(Edges, (NumEdges + NumExteriorEdges)*sizeof(RWGEdge *));
     memcpy(Edges + NumEdges, ExteriorEdges, NumExteriorEdges * sizeof(RWGEdge *));
     for(int ne=NumEdges; ne<NumEdges + NumExteriorEdges; ne++)
      Edges[ne]->Index = ne;
     NumEdges += NumExteriorEdges;
   };
  NumBFs = IsPEC ? NumEdges : 2*NumEdges;

} 

/***************************************************************/
/* RWGSurface destructor.                                       */
/***************************************************************/
RWGSurface::~RWGSurface()
{ 
  free(Vertices);

  int np;
  for(np=0; np<NumPanels; np++)
   free(Panels[np]);
  free(Panels);

  int ne;
  for(ne=0; ne<NumEdges; ne++)
   free(Edges[ne]);
  free(Edges);

  for(ne=0; ne<NumExteriorEdges; ne++)
   free(ExteriorEdges[ne]);
  free(ExteriorEdges);

  int nbc;
  for(nbc=0; nbc<NumBCs; nbc++)
   free(BCEdges[nbc]);
  if (BCEdges) free(BCEdges);
  if (NumBCEdges) free(NumBCEdges);
  free(WhichBC);

  if (MeshFileName) free(MeshFileName);
  if (Label) free(Label);
  if (ErrMsg) free(ErrMsg);
  if (GT) delete GT;

  if (MaterialName) free(MaterialName);
  if (RegionLabels[0]) free(RegionLabels[0]);
  if (RegionLabels[1]) free(RegionLabels[1]);

  kdtri_destroy(kdPanels);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RWGSurface::Transform(const GTransformation *DeltaGT)
{ 
  /***************************************************************/
  /*- first apply the transformation to all points whose         */
  /*- coordinates we store inside the RWGSurface structure:       */
  /*- vertices, edge centroids, and panel centroids.             */
  /***************************************************************/
  /* vertices */
  DeltaGT->Apply(Vertices, NumVertices);

  /* edge centroids */
  int ne;
  for(ne=0; ne<NumEdges; ne++)
    DeltaGT->Apply(Edges[ne]->Centroid, 1);

  /***************************************************************/
  /* reinitialize geometric data on panels (which takes care of  */ 
  /* transforming the panel centroids)                           */ 
  /***************************************************************/
  int np;
  for(np=0; np<NumPanels; np++)
   InitRWGPanel(Panels[np], Vertices);

  /***************************************************************/
  /* update the internally stored GTransformation ****************/
  /***************************************************************/
  if (!GT)
    GT = new GTransformation(DeltaGT);
  else
    GT->Transform(DeltaGT);
}

void RWGSurface::Transform(char *format,...)
{
  va_list ap;
  char buffer[MAXSTR];
  va_start(ap,format);
  vsnprintf(buffer,MAXSTR,format,ap);

  GTransformation OTGT(buffer, &ErrMsg);
  if (ErrMsg)
   ErrExit(ErrMsg);
  Transform(&OTGT);
}

void RWGSurface::UnTransform()
{
  if (!GT)
   return;
 
  /***************************************************************/
  /* untransform vertices                                        */
  /***************************************************************/
  GT->UnApply(Vertices, NumVertices);

  /***************************************************************/
  /* untransform edge centroids                                  */
  /***************************************************************/
  int ne;
  for(ne=0; ne<NumEdges; ne++)
    GT->UnApply(Edges[ne]->Centroid, 1);

  /***************************************************************/
  /* reinitialize geometric data on panels (which takes care of  */ 
  /* transforming the panel centroids)                           */ 
  /***************************************************************/
  int np;
  for(np=0; np<NumPanels; np++)
   InitRWGPanel(Panels[np], Vertices);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  GT->Reset();
}

/*-----------------------------------------------------------------*/
/*- initialize geometric quantities stored within an RWGPanel      */
/*-----------------------------------------------------------------*/
void InitRWGPanel(RWGPanel *P, double *Vertices)
{
  int i;
  double *V[3]; 
  double VA[3], VB[3];

  V[0]=Vertices + 3*P->VI[0];
  V[1]=Vertices + 3*P->VI[1];
  V[2]=Vertices + 3*P->VI[2];

  /* calculate centroid */ 
  for(i=0; i<3; i++)
   P->Centroid[i]=(V[0][i] + V[1][i] + V[2][i])/3.0;

  /* calculate bounding radius */ 
  P->Radius=VecDistance(P->Centroid,V[0]);
  P->Radius=fmax(P->Radius,VecDistance(P->Centroid,V[1]));
  P->Radius=fmax(P->Radius,VecDistance(P->Centroid,V[2]));

  /* 
   * compute panel normal. notice that this computation preserves a  
   * right-hand rule: if we curl the fingers of our right hand       
   * around the triangle with the edges traversed in the order       
   * V[0]-->V[1], V[1]-->V[2], V[2]-->V[0], then our thumb points
   * in the direction of the panel normal.
   */ 
  VecSub(V[1],V[0],VA);
  VecSub(V[2],V[0],VB);
  VecCross(VA,VB,P->ZHat);

  P->Area=0.5*VecNormalize(P->ZHat);
}

/*-----------------------------------------------------------------*/
/*- Initialize and return a pointer to a new RWGPanel structure. -*/
/*-----------------------------------------------------------------*/
RWGPanel *NewRWGPanel(double *Vertices, int iV1, int iV2, int iV3)
{ 
  RWGPanel *P;

  P=(RWGPanel *)mallocEC(sizeof *P);
  P->VI[0]=iV1;
  P->VI[1]=iV2;
  P->VI[2]=iV3;

  InitRWGPanel(P, Vertices);

  return P;

} 

} // namespace scuff
