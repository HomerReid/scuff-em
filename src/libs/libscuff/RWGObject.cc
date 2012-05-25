/*
 * RWGObject.cc -- implementation of some methods in the RWGObject
 *               -- class 
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

namespace scuff {

#define MAXSTR 1000
#define MAXTOK 50  

/*--------------------------------------------------------------*/
/*-  RWGObject class constructor that begins reading an open   -*/
/*-  file immediately after a line like OBJECT MyObjectLabel.  -*/
/*-                                                            -*/
/*-  On return, if an object was successfully created, its     -*/
/*-   ErrMsg field is NULL, and the file read pointer points   -*/
/*-   the line immediately following the ENDOBJECT line.       -*/
/*-                                                            -*/
/*-  If there was an error in parsing the OBJECT section, the  -*/
/*-   ErrMsg field points to an error message string.          -*/
/*-                                                            -*/
/*-  Note: the 'ContainingObject' field in the OBJECT          -*/
/*-        definition cannot be processed at the level of      -*/
/*-        RWGObject; instead; if that field is specified we   -*/
/*-        store its value in the 'ContainingObjetName' class  -*/
/*-        data field for subsequent processsing at the level  -*/
/*-        of RWGGeometry.                                     -*/
/*--------------------------------------------------------------*/
RWGObject::RWGObject(FILE *f, const char *pLabel, int *LineNum)
{ 
  ErrMsg=0;
  ContainingObjectLabel=0;

  /***************************************************************/
  /* read lines from the file one at a time **********************/
  /***************************************************************/
  char Line[MAXSTR];
  char MaterialName[MAXSTR];
  int NumTokens, TokensConsumed;
  char *Tokens[MAXTOK];
  int ReachedTheEnd=0;
  char *pMeshFileName=0;
  GTransformation OTGT; // 'one-time geometrical transformation'
  MaterialName[0]=0;
  while ( ReachedTheEnd==0 && fgets(Line, MAXSTR, f) )
   { 
     (*LineNum)++;
     NumTokens=Tokenize(Line, Tokens, MAXTOK);
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
        pMeshFileName=strdup(Tokens[1]);
      }
     else if ( !strcasecmp(Tokens[0],"MATERIAL") )
      { if (NumTokens!=2)
         { ErrMsg=strdup("MATERIAL keyword requires one argument");
           return;
         };
        strncpy(MaterialName, Tokens[1], MAXSTR);
      }
     else if ( !strcasecmp(Tokens[0],"INSIDE") )
      { if (NumTokens!=2)
         { ErrMsg=strdup("INSIDE keyword requires one argument");
           return;
         };
        ContainingObjectLabel=strdup(Tokens[1]);
      }
     else if ( !strcasecmp(Tokens[0],"DISPLACED") || !strcasecmp(Tokens[0],"ROTATED") )
      { 
        // try to parse the line as a geometrical transformation.
        // note that OTGT is used as a running GTransformation that may
        // be augmented by multiple DISPLACED ... and/or ROTATED ...
        // lines within the OBJECT...ENDOBJECT section, and which is 
        // applied to the object at its birth and subsequently discarded.
        // in particular, OTGT is NOT stored as the 'GT' field inside 
        // the Object class, which is intended to be used for 
        // transformations that are applied and later un-applied 
        // during the life of the object. 
	OTGT.Parse(Tokens, NumTokens, &ErrMsg, &TokensConsumed);
        if (ErrMsg)
         return;
        if (TokensConsumed!=NumTokens) 
         { ErrMsg=strdup("junk at end of line");
           return;
         };
      }
     else if ( !strcasecmp(Tokens[0],"ENDOBJECT") )
      { 
        ReachedTheEnd=1;
      }
     else
      { ErrMsg=vstrdup("unknown keyword %s in OBJECT section",Tokens[0]);
        return;
      };
   }; 

  if (pMeshFileName==0)
   ErrMsg=vstrdup("OBJECT section must include a MESHFILE specification",Tokens[0]);

  InitRWGObject(pMeshFileName, pLabel, MaterialName, &OTGT);
  
  free(pMeshFileName);

}

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- alternative constructor entry point that uses a given       */
/*- meshfile name                                               */
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
RWGObject::RWGObject(const char *pMeshFileName)
 { InitRWGObject(pMeshFileName, 0, 0, 0); }

RWGObject::RWGObject(const char *pMeshFileName, 
                     const char *pLabel,
                     const char *Material,
                     const GTransformation *OTGT)
 { InitRWGObject(pMeshFileName, pLabel, Material, OTGT); }

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- Actual body of RWGObject constructor: Create an RWGObject */
/*- from a mesh file describing a discretized object,           */
/*- optionally with a rotation and/or displacement applied.     */
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
void RWGObject::InitRWGObject(const char *pMeshFileName,
                              const char *pLabel,
                              const char *Material,
                              const GTransformation *OTGT)
{ 
  ErrMsg=0;
  kdPanels = NULL;

  /*------------------------------------------------------------*/
  /*- try to open the mesh file.                                */
  /*------------------------------------------------------------*/
  FILE *MeshFile=fopen(pMeshFileName,"r");
  if (!MeshFile)
   ErrExit("could not open file %s",pMeshFileName);
   
  /*------------------------------------------------------------*/
  /*- initialize simple fields ---------------------------------*/
  /*------------------------------------------------------------*/
  NumEdges=NumPanels=NumVertices=NumRefPts=0;
  MP=new MatProp(Material);
  ContainingObject=0;
  Label=strdup( pLabel ? pLabel : "NoLabel");
  MeshFileName=strdup(pMeshFileName);

  /*------------------------------------------------------------*/
  /*- note: the 'OTGT' parameter to this function is distinct   */
  /*- from the 'GT' field inside the class body. the former is  */
  /*- an optional 'One-Time Geometrical Transformation' to be   */
  /*- applied to the object once at its creation. the latter    */
  /*- is designed to store a subsequent transformation that may */
  /*- be applied to the object, and is initialized to zero.     */
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
   ReadGMSHFile(MeshFile,MeshFileName,OTGT);
  else if (!strcasecmp(p,"mphtxt"))
   ReadComsolFile(MeshFile,MeshFileName,OTGT);
  else
   ErrExit("file %s: unknown extension %s",MeshFileName,p);

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
  /*- the number of basis functions that this object takes up   */
  /*- in the full BEM system is the number of internal edges    */
  /*- if it is a perfect electrical conductor, or 2 times the   */
  /*- number of internal edges otherwise.                       */
  /*------------------------------------------------------------*/
  NumBFs = MP->IsPEC() ? NumEdges : 2*NumEdges;

} 

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- Alternative RWGObject constructor: Create an RWGObject    */
/*- from lists of vertices and panels.                          */
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
RWGObject::RWGObject(double *pVertices, int pNumVertices,
                     int **PanelVertexIndices, int pNumPanels)
{ 
  int np;

  ErrMsg=0;
  kdPanels = NULL;

  MeshFileName=strdup("ByHand.msh");
  Label=strdup("ByHand");

  NumEdges=0;
  NumVertices=pNumVertices;
  NumPanels=pNumPanels;
  MP=new MatProp();
  GT=0;

  Vertices=(double *)mallocEC(3*NumVertices*sizeof(double));
  memcpy(Vertices,pVertices,3*NumVertices*sizeof(double *));

  Panels=(RWGPanel **)mallocEC(NumPanels*sizeof(RWGPanel *));
  for(np=0; np<NumPanels; np++)
   { Panels[np]=NewRWGPanel(Vertices,
                             PanelVertexIndices[np][0],
                             PanelVertexIndices[np][1],
                             PanelVertexIndices[np][2]);
     Panels[np]->Index=np;
   };
 
  /*------------------------------------------------------------*/
  /* gather necessary edge connectivity info                   -*/
  /*------------------------------------------------------------*/
  InitEdgeList();

  NumBFs = MP->IsPEC() ? NumEdges : 2*NumEdges;
} 

/***************************************************************/
/* RWGObject destructor.                                       */
/***************************************************************/
RWGObject::~RWGObject()
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

  // insert code to deallocate BCEdges
  int nbc;
  for(nbc=0; nbc<NumBCs; nbc++)
   free(BCEdges[nbc]);
  if (BCEdges) free(BCEdges);
  if (NumBCEdges) free(NumBCEdges);
  free(WhichBC);

  delete MP;

  if (MeshFileName) free(MeshFileName);
  if (Label) free(Label);
  if (GT) delete GT;

  kdtri_destroy(kdPanels);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RWGObject::Transform(const GTransformation *DeltaGT)
{ 
  /***************************************************************/
  /*- first apply the transformation to all points whose         */
  /*- coordinates we store inside the RWGObject structure:       */
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

void RWGObject::Transform(char *format,...)
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

void RWGObject::UnTransform()
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
