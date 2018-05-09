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
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
void RWGSurface::UpdateBoundingBox()
{ 
  RMax[0] = RMax[1] = RMax[2] = -1.0e89;
  RMin[0] = RMin[1] = RMin[2] = +1.0e89;
  double *V;
  for(int np=0; np<NumPanels; np++)
   for(int i=0; i<3; i++)
    { V = Vertices + 3*(Panels[np]->VI[i]);
      RMax[0] = fmax(RMax[0], V[0] );
      RMax[1] = fmax(RMax[1], V[1] );
      RMax[2] = fmax(RMax[2], V[2] );
      RMin[0] = fmin(RMin[0], V[0] );
      RMin[1] = fmin(RMin[1], V[1] );
      RMin[2] = fmin(RMin[2], V[2] );
    };
}

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
  SurfaceZeta=0;
  MeshTag=-1;
  MeshFileName=0;
  IsPEC=1;
  Label = strdupEC(pLabel);
  tolVecClose=0.0; // to be updated once mesh is read in
  OTGT=GT=0;

  if ( !StrCaseCmp(Keyword, "OBJECT") )
   IsObject=1;
  else
   IsObject=0;

  /***************************************************************/
  /* read lines from the file one at a time **********************/
  /***************************************************************/
  char Line[MAXSTR], LineCopy[MAXSTR];
  int NumTokens, TokensConsumed;
  char *Tokens[MAXTOK];
  int ReachedTheEnd=0;
  MaterialName=0;
  RegionLabels[0]=RegionLabels[1]=0;
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
     if ( !StrCaseCmp(Tokens[0],"MESHFILE") )
      { if (NumTokens!=2)
         { ErrMsg=strdupEC("MESHFILE keyword requires one argument");
           return;
         };
        MeshFileName=strdupEC(Tokens[1]);
      }
     else if ( !StrCaseCmp(Tokens[0],"MESHTAG") )
      { if (NumTokens!=2)
         { ErrMsg=strdupEC("MESHTAG keyword requires one argument");
           return;
         };
        sscanf(Tokens[1],"%i",&MeshTag);
      }
     else if ( !StrCaseCmp(Tokens[0],"MATERIAL") )
      { 
        if (!IsObject)
         { ErrMsg=strdupEC("MATERIAL keyword may not be used in SURFACE...ENDSURFACE sections");
           return;
         };
        if (NumTokens!=2)
         { ErrMsg=strdupEC("MATERIAL keyword requires one argument");
           return;
         };
        MaterialName = strdupEC(Tokens[1]);
        RegionLabels[0] = strdupEC("EXTERIOR");
        RegionLabels[1] = strdupEC(Label);
        MaterialRegionsLineNum = *LineNum;
        if ( !StrCaseCmp(MaterialName,"PEC") )
         IsPEC=1;
        else
         IsPEC=0;
      }
     else if ( !StrCaseCmp(Tokens[0],"REGIONS") )
      { 
        if (IsObject)
         { ErrMsg=strdupEC("REGIONS keyword may not be used in OBJECT...ENDOBJECT sections");
           return;
         };
        if (NumTokens==3)
         { IsPEC=0;
           RegionLabels[0] = strdupEC(Tokens[1]);
           RegionLabels[1] = strdupEC(Tokens[2]);
           MaterialRegionsLineNum = *LineNum;
         }
        else if (NumTokens==2)
         { IsPEC=1;
           RegionLabels[0] = strdupEC(Tokens[1]);
           RegionLabels[1] = strdupEC("PEC");
           MaterialRegionsLineNum = *LineNum;
         }
        else
         ErrMsg=strdupEC("REGIONS keyword requires one or two arguments");
      }
     else if ( !StrCaseCmp(Tokens[0],"DISPLACED") || !StrCaseCmp(Tokens[0],"ROTATED") )
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
         { ErrMsg=strdupEC("junk at end of line");
           return;
         };
      }
     else if ( !StrCaseCmp(Tokens[0],"SURFACE_IMPEDANCE") )
      { 
        if (NumTokens<2)
         { ErrMsg=strdupEC("no argument specified for SURFACE_CONDUCTIVITY");
           return;
         };

        /* try to create a cevaluator for the user's function */
        char *ZetaString=Line; 
        while (strchr(" \t\n",*ZetaString)) // skip forward to the start of the 
         ZetaString++;                        //  actual expression
        ZetaString+=strlen("SURFACE_IMPEDANCE");
        int SSLength=strlen(ZetaString);
        if (SSLength>=1 && ZetaString[SSLength-1]=='\n') 
         ZetaString[SSLength-1]=0;
        SurfaceZeta=cevaluator_create(ZetaString);
        if (SurfaceZeta==0)
         { ErrMsg=vstrdup("invalid SURFACE_IMPEDANCE %s",ZetaString);
           return;
         };

        /* these calls are required for thread safety, and they     */
        /* mandate that the w, x, y, z values be specified in that  */
        /* order in the 'values' parameter to cevaluator_evaluate() */ 
        cevaluator_set_var_index(SurfaceZeta, "w", 0);
        cevaluator_set_var_index(SurfaceZeta, "x", 1);
        cevaluator_set_var_index(SurfaceZeta, "y", 2);
        cevaluator_set_var_index(SurfaceZeta, "z", 3);

      }
     else if (   !StrCaseCmp(Tokens[0],"ENDOBJECT") || !StrCaseCmp(Tokens[0],"ENDSURFACE") )
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

  if (SurfaceZeta!=0)
   Log("Surface %s has surface impedance Zeta=%s.\n",Label,cevaluator_get_string(SurfaceZeta));

  // if we are an OBJECT and there was no MATERIAL specification,
  // or we are a SURFACE and there was no REGIONS specification, then 
  // we consider ourselves to be a PEC surface in the exterior medium.
  if ( RegionLabels[0]==0 )
   { IsPEC=1;
     RegionLabels[0] = strdupEC("EXTERIOR");
     RegionLabels[1] = strdupEC("PEC");
   };

  // ok, all checks passed, now on to the main body of the class constructor.
  ErrMsg=InitRWGSurface();
  
}

/*--------------------------------------------------------------*/
/*- This construct creates an RWGSurface that does not belong  -*/
/*- to any RWGGeometry. Useful for flux meshes, etc.           -*/
/*--------------------------------------------------------------*/
RWGSurface::RWGSurface(const char *MeshFile, int pMeshTag)
{
  MeshFileName=strdup(MeshFile);
  MeshTag=pMeshTag;
  Label=strdup(MeshFile);
  SurfaceZeta=0;
  MaterialName=0;
  RegionLabels[0]=RegionLabels[1]=0;
  IsPEC=1;
  tolVecClose=0.0; // to be updated once mesh is read in
  GT=OTGT=0;
  ErrMsg=InitRWGSurface();
}

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- Main body of RWGSurface constructor: Create an RWGSurface   */
/*- from a mesh file describing a discretized object,           */
/*- optionally with a rotation and/or displacement applied.     */
/*-                                                             */
/*- This routine assumes that the following internal class      */
/*- fields have been initialized: MeshFileName, MeshTag,        */
/*- Label, SurfaceZeta, RegionLabels, IsPEC.                    */
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
char *RWGSurface::InitRWGSurface()
{ 
  NumVertices=NumPanels=NumEdges=NumHalfRWGEdges=NumBCs=0;
  Vertices=0;
  Panels=0;
  Edges=0;
  HalfRWGEdges=0;
  BCEdges=0;
  NumBCEdges=0;
  WhichBC=0;
  ErrMsg=0;
  GT=0;
  OTGT=0;
  kdPanels = 0;

  /*------------------------------------------------------------*/
  /*- try to open the mesh file. we look in several places:     */
  /*- (a) the current working directory                         */
  /*- (b) any directories that may have been specified by       */
  /*-     MESHPATH statements in .scuffgeo files or via the     */
  /*-     SCUFF_MESH_PATH environment variable                  */
  /*------------------------------------------------------------*/
  char *WhichDir=0;
  FILE *MeshFile=fopenPath(getenv("SCUFF_MESH_PATH"),  MeshFileName, "r", &WhichDir);
#if 0
  FILE *MeshFile=fopen(MeshFileName,"r");
  if (!MeshFile)
   { for(int nmd=0; MeshFile==0 && nmd<RWGGeometry::NumMeshDirs; nmd++)
      { MeshFile=vfopen("%s/%s","r",RWGGeometry::MeshDirs[nmd],MeshFileName);
        if (MeshFile) 
         Log("Found mesh file %s/%s",RWGGeometry::MeshDirs[nmd],MeshFileName);
      }
   }
#endif
  if (!MeshFile)
   return vstrdup("could not open file %s",MeshFileName);
  if  (WhichDir)
   Log(" Opened mesh file %s/%s",WhichDir,MeshFileName);
  else
   Log(" Opened mesh file %s",MeshFileName);
   
  /*------------------------------------------------------------*/
  /*- initialize simple fields ---------------------------------*/
  /*------------------------------------------------------------*/
  NumEdges=NumPanels=NumVertices=NumRefPts=TotalStraddlers=0;
  PhasedBFCs=0;
  Origin[0]=Origin[1]=Origin[2]=0.0;
  if (OTGT) OTGT->Apply(Origin);

  /*------------------------------------------------------------*/
  /*- note: the 'OTGT' parameter to this function is distinct   */
  /*- from the 'GT' field inside the class body. the former is  */
  /*- an optional 'One-Time Geometrical Transformation' that is */
  /*- applied to the object once at its creation and is never   */
  /*- un-applied. the latter is designed to store a subsequent  */
  /*- temporary transformation that may be applied and then un- */
  /*- applied to the surface.                                   */
  /*------------------------------------------------------------*/
  GT=0;

  /*------------------------------------------------------------*/
  /*- Switch off based on the file type to read the mesh file:  */
  /*-  1. file extension=.msh    --> ReadGMSHFile              -*/
  /*-  2. file extension=.mphtxt --> ReadComsolFile            -*/
  /*------------------------------------------------------------*/
  char *p=GetFileExtension(MeshFileName);
  if (!p)
   return vstrdup("file %s: invalid extension",MeshFileName);
  else if (!StrCaseCmp(p,"msh"))
   { char *GMSHError = ReadGMSHFile(MeshFile,MeshFileName);
     if (GMSHError) return GMSHError;
   }
  else if (!StrCaseCmp(p,"mphtxt"))
   { if ( MeshTag != -1 )
      return vstrdup("MESHTAG is not yet implemented for .mphtxt files");
     ReadComsolFile(MeshFile,MeshFileName);
   }
  else
   return vstrdup("file %s: unknown extension %s",MeshFileName,p);

  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  if (NumPanels==0)
   if ( MeshTag == -1 ) 
    return vstrdup("file %s: no panels found",MeshFileName);
   else
    return vstrdup("file %s: no panels found for mesh tag %i",MeshFileName,MeshTag);

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
  ErrMsg=InitEdgeList();
  if (ErrMsg) return ErrMsg;

  /*------------------------------------------------------------*/
  /*- After InitEdgeList, the Edges array has length NumEdges   */
  /*- and contains only full RWG edges, while the HalfRWGEdges  */
  /*- array has length NumHalfRWGEdges=NumExteriorEdges and     */
  /*- contains only half-RWG edges. (If any point-like RWGPorts */
  /*- are subsequently added, the HalfRWGEdges array will grow  */
  /*- to accommodate the resulting half-RWG edges, and          */
  /*- NumHalfRWGEdges will increase, but NumExteriorEdges will  */
  /*- remain constant.)                                         */
  /*-                                                           */
  /*- By default the HalfRWGEdges are not treated as degrees of */
  /*- freedom for expanding induced surface currents. This is   */
  /*- overridden if RWGGeometry::UserHRWGFunctions==true, in    */
  /*- which case the full content of the HalfRWGEdges array is  */
  /*- on to the end of the Edges array and NumEdges (and */
  /*- NumBFs, and thus the dimension of the SIE system) grows   */
  /*- accordingly.                                              */
  /*------------------------------------------------------------*/
  if ( RWGGeometry::UseHRWGFunctions )
   { Edges = (RWGEdge **)realloc(Edges, (NumEdges + NumHalfRWGEdges)*sizeof(RWGEdge *));
     memcpy(Edges + NumEdges, HalfRWGEdges, NumHalfRWGEdges * sizeof(RWGEdge *));
     for(int ne=NumEdges; ne<NumEdges + NumHalfRWGEdges; ne++)
      Edges[ne]->Index = ne;
     NumEdges += NumHalfRWGEdges;
     Log("Promoted %i exterior edges for surface %s to half-RWG basis functions.",NumHalfRWGEdges,Label);
   }

  /*------------------------------------------------------------*/
  /*- now that we have put the edges into a list, we can go    -*/
  /*- back and fill in the EdgeIndices field in each RWGPanel. -*/
  /*- EdgeIndices[i] = index of panel edge #i within the       -*/
  /*- parent RWGSurface. (This index is negative if the edge is-*/
  /*- an exterior (half-RWG) edge and UseHRWGFunctions==false. -*/
  /*------------------------------------------------------------*/
  for(int ne=0; ne<NumEdges; ne++)
   { RWGEdge *E = Edges[ne];
     Panels[ E->iPPanel ] -> EI[ E->PIndex ] = ne;
     if ( E->iMPanel!=-1 )
      Panels[ E->iMPanel ] -> EI[ E->MIndex ] = ne;
   }
  if ( !(RWGGeometry::UseHRWGFunctions) )
   for(int ne=0; ne<NumHalfRWGEdges; ne++)
    { RWGEdge *E=HalfRWGEdges[ne];
      Panels[E->iPPanel]->EI[E->PIndex] = -(ne+1);
    }

  /*------------------------------------------------------------*/
  /*- the number of basis functions that this object takes up   */
  /*- in the full BEM system is the number of internal edges    */
  /*- if it is a perfect electrical conductor, or 2 times the   */
  /*- number of internal edges otherwise.                       */
  /*------------------------------------------------------------*/
  NumBFs = IsPEC ? NumEdges : 2*NumEdges;
  IsClosed = (NumExteriorEdges == 0);

  UpdateBoundingBox();

  /*------------------------------------------------------------*/
  /*- 20150929 initialize the kdtri by default, instead of      */
  /*-          waiting until the first call to Contains()       */
  /*-          to do so. The reason for this change is that we  */
  /*-          want to make sure we form the kdtri prior to any */
  /*-          geometric transformation that may be carried out */
  /*-          on the object; this way, the kdtri always        */
  /*-          reflects the UNTRANSFORMED object, and hence when*/
  /*-          Contains(x) is called for an object that has been*/
  /*-          transformed, we can just untransform x to make   */
  /*-          the Contain() work correctly.                    */
  /*------------------------------------------------------------*/
  InitkdPanels(false);

  return 0;
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
  ErrMsg=0;
  kdPanels = NULL;

  MeshFileName=strdupEC("ByHand.msh");
  Label=strdupEC("ByHand");

  NumEdges=0;
  NumVertices=pNumVertices;
  NumPanels=pNumPanels;
  RegionLabels[0]=strdupEC("Exterior");
  RegionLabels[1]=strdupEC("PEC");
  RegionIndices[0]=0;
  RegionIndices[1]=-1;
  IsPEC=1;
  IsObject=1;
  OTGT=GT=0;
  Origin[0]=Origin[1]=Origin[2]=0.0;

  Vertices=(double *)mallocEC(3*NumVertices*sizeof(double));
  memcpy(Vertices,pVertices,3*NumVertices*sizeof(double *));

  /*------------------------------------------------------------*/
  /*- add the panels -------------------------------------------*/
  /*------------------------------------------------------------*/
  Panels=(RWGPanel **)mallocEC(NumPanels*sizeof(RWGPanel *));
  for(int np=0; np<NumPanels; np++)
   { Panels[np]=NewRWGPanel(Vertices,
                            PanelVertices[3*np+0],
                            PanelVertices[3*np+1],
                            PanelVertices[3*np+2]);
     Panels[np]->Index=np;
   };
 
  /*------------------------------------------------------------*/
  /* gather necessary edge connectivity info                   -*/
  /*------------------------------------------------------------*/
  if ( (ErrMsg=InitEdgeList()) ) return;

  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  if (RWGGeometry::UseHRWGFunctions)
   { Edges = (RWGEdge **)realloc(Edges, (NumEdges + NumHalfRWGEdges)*sizeof(RWGEdge *));
     memcpy(Edges + NumEdges, HalfRWGEdges, NumHalfRWGEdges* sizeof(RWGEdge *));
     for(int ne=NumEdges; ne<NumEdges + NumHalfRWGEdges; ne++)
      Edges[ne]->Index = ne;
     NumEdges += NumHalfRWGEdges;
   }
  NumBFs = IsPEC ? NumEdges : 2*NumEdges;

  /*------------------------------------------------------------*/
  /*- now that we have put the edges into a list, we can go    -*/
  /*- back and fill in the EdgeIndices field in each RWGPanel. -*/
  /*- Note: EdgeIndices[i] = index of panel edge #i within the -*/
  /*- Edges[] list for the parent RWGSurface. Panel edge #i is -*/
  /*- edge opposite vertex #i.                                 -*/
  /*------------------------------------------------------------*/
  for(int ne=0; ne<NumEdges; ne++)
   { RWGEdge *E = Edges[ne];
     Panels[ E->iPPanel ] -> EI[ E->PIndex ] = ne;
     if ( E->iMPanel>=0 )
      Panels[ E->iMPanel ] -> EI[ E->MIndex ] = ne;
   }
  if ( !(RWGGeometry::UseHRWGFunctions) )
   for(int ne=0; ne<NumHalfRWGEdges; ne++)
    { RWGEdge *E=HalfRWGEdges[ne];
      Panels[E->iPPanel]->EI[E->PIndex] = -(ne+1);
    }

  /*------------------------------------------------------------*/
  /*- compute bounding box -------------------------------------*/
  /*------------------------------------------------------------*/
  RMax[0] = RMax[1] = RMax[2] = -1.0e89;
  RMin[0] = RMin[1] = RMin[2] = +1.0e89;
  double *V;
  for(int np=0; np<NumPanels; np++)
   for(int i=0; i<3; i++)
    { V = Vertices + 3*(Panels[np]->VI[i]);
      RMax[0] = fmax(RMax[0], V[0] );
      RMax[1] = fmax(RMax[1], V[1] );
      RMax[2] = fmax(RMax[2], V[2] );
      RMin[0] = fmin(RMin[0], V[0] );
      RMin[1] = fmin(RMin[1], V[1] );
      RMin[2] = fmin(RMin[2], V[2] );
    }
} 

/***************************************************************/
/* RWGSurface destructor.                                       */
/***************************************************************/
RWGSurface::~RWGSurface()
{ 
  if(Vertices) free(Vertices);

  for(int np=0; np<NumPanels; np++)
   free(Panels[np]);
  if(Panels) free(Panels);

  for(int ne=0; ne<NumEdges; ne++)
   free(Edges[ne]);
  if (Edges) free(Edges);

  for(int ne=0; ne<NumHalfRWGEdges; ne++)
   free(HalfRWGEdges[ne]);
  if (HalfRWGEdges) free(HalfRWGEdges);

  for(int nbc=0; nbc<NumBCs; nbc++)
   free(BCEdges[nbc]);
  if (BCEdges) free(BCEdges);
  if (NumBCEdges) free(NumBCEdges);
  if (WhichBC) free(WhichBC);

  if (MeshFileName) free(MeshFileName);
  if (Label) free(Label);
  if (ErrMsg) free(ErrMsg);
  if (GT) delete GT;
  if (OTGT) delete OTGT;

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
  if (DeltaGT==0)
   return;

  /***************************************************************/
  /*- first apply the transformation to all points whose         */
  /*- coordinates we store inside the RWGSurface structure:      */
  /*- vertices, edge centroids, panel centroids, and the origin  */
  /*- of coordinates.                                            */
  /***************************************************************/

  /* vertices */
  DeltaGT->Apply(Vertices, NumVertices);

  /* edge centroids */
  for(int ne=0; ne<NumEdges; ne++)
    DeltaGT->Apply(Edges[ne]->Centroid, 1);
  for(int ne=0; ne<NumHalfRWGEdges; ne++)
    DeltaGT->Apply(HalfRWGEdges[ne]->Centroid, 1);

  /* origin */
  DeltaGT->Apply(Origin);

  /***************************************************************/
  /* reinitialize geometric data on panels (which takes care of  */ 
  /* transforming the panel centroids)                           */ 
  /***************************************************************/
  int np;
  for(np=0; np<NumPanels; np++)
   InitRWGPanel(Panels[np], Vertices);

  UpdateBoundingBox();

  /***************************************************************/
  /* update the internally stored GTransformation ****************/
  /***************************************************************/
  if (!GT)
    GT = new GTransformation(DeltaGT);
  else
    GT->Transform(DeltaGT);
}

void RWGSurface::Transform(const char *format,...)
{
  va_list ap;
  char buffer[MAXSTR];
  va_start(ap,format);
  vsnprintfEC(buffer,MAXSTR,format,ap);

  GTransformation MyGT(buffer, &ErrMsg);
  if (ErrMsg)
   ErrExit(ErrMsg);
  Transform(&MyGT);
}

void RWGSurface::UnTransform()
{
  if (!GT)
   return;
 
  /***************************************************************/
  /* untransform vertices, edge centroids, origin                */
  /***************************************************************/
  GT->UnApply(Vertices, NumVertices);

  for(int ne=0; ne<NumEdges; ne++)
   GT->UnApply(Edges[ne]->Centroid, 1);
  for(int ne=0; ne<NumHalfRWGEdges; ne++)
   GT->UnApply(HalfRWGEdges[ne]->Centroid, 1);

  GT->UnApply(Origin);

  /***************************************************************/
  /* reinitialize geometric data on panels (which takes care of  */ 
  /* transforming the panel centroids)                           */ 
  /***************************************************************/
  int np;
  for(np=0; np<NumPanels; np++)
   InitRWGPanel(Panels[np], Vertices);

  UpdateBoundingBox();

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
  double *V[3]; 
  V[0]=Vertices + 3*P->VI[0];
  V[1]=Vertices + 3*P->VI[1];
  V[2]=Vertices + 3*P->VI[2];

  /* calculate centroid */ 
  for(int i=0; i<3; i++)
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
   * 
   * 20121114 update: we may subsequently need to flip the sign of 
   * the panel normal to ensure that it points into the exterior 
   * medium for the surface in question, and by that time it 
   * will be too late to reorder the panel vertices (because, 
   * among other things, RWGEdge structures that depend on them 
   * will have been defined) which means that the right-hand-rule
   * convention for the panel normal may not be preserved and must 
   * not be relied upon.
   * TODO: FIXME
   *                  
   */ 
  double VA[3], VB[3];
  VecSub(V[1],V[0],VA);
  VecSub(V[2],V[0],VB);
  VecCross(VA,VB,P->ZHat);
  if (P->ZHatFlipped)
   VecScale(P->ZHat, -1.0);

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

  P->ZHatFlipped=false;

  InitRWGPanel(P, Vertices);

  return P;

} 

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
RWGEdge *RWGSurface::GetEdgeByIndex(int ne)
{ if (ne>=0)
   { if (ne>=NumEdges) 
      ErrExit("%s:%i: internal error(%i,%i,%i)",__FILE__,__LINE__,Index,NumEdges,ne);
     return Edges[ne];
   }
  ne = -(ne+1);
  if (ne>=NumHalfRWGEdges)
   ErrExit("%s:%i: internal error(%i,%i,%i)",__FILE__,__LINE__,Index,NumHalfRWGEdges,ne);
  return HalfRWGEdges[ne];
}

} // namespace scuff
