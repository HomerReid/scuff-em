/*
 * RWGObject.cc -- implementation of some methods in the RWGObject
 *               -- class 
 *
 * homer reid    -- 3/2007 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <libhrutil.h>

#include "libscuff.h"

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
/*-  Note: even if the object is successfully read in, the     -*/
/*-        'MP' and 'ContainingObject' fields are not          -*/
/*-        initialized; instead, any user-specified values for -*/
/*-        these fields are stored in the internal class       -*/
/*-        fields 'MPName' and 'ContainingObjectName' and must -*/
/*-        be processed by the calling routine.                -*/
/*--------------------------------------------------------------*/
RWGObject::RWGObject(FILE *f, const char *pLabel, int *LineNum)
{ 
  ErrMsg=0;
  ContainingObjectLabel=0;
  MPName=0;
  MeshFileName=0;

  /***************************************************************/
  /* read lines from the file one at a time **********************/
  /***************************************************************/
  char Line[MAXSTR];
  int nt, NumTokens;
  char *p, *Tokens[MAXTOK];
  int ReachedTheEnd=0;
  GTransformation *GT=0;
  double DX[3], ZHat[3], Theta;
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
        MeshFileName=strdup(Tokens[1]);
      }
     else if ( !strcasecmp(Tokens[0],"MATERIAL") )
      { if (NumTokens!=2)
         { ErrMsg=strdup("MATERIAL keyword requires one argument");
           return;
         };
        MPName=strdup(Tokens[1]);
      }
     else if ( !strcasecmp(Tokens[0],"INSIDE") )
      { if (NumTokens!=1)
         { ErrMsg=strdup("INSIDE keyword requires one argument");
           return;
         };
        ContainingObjectLabel=strdup(Tokens[1]);
      }
     else if ( !strcasecmp(Tokens[0],"DISPLACED") )
      { 
        if (NumTokens!=4)
         { ErrMsg=strdup("DISPLACED keyword requires 3 arguments");
           return;
         };

        if (    1!=sscanf(Tokens[1],"%le",DX+0)
             || 1!=sscanf(Tokens[2],"%le",DX+1)
             || 1!=sscanf(Tokens[3],"%le",DX+2) ) 
         { ErrMsg=strdup("invalid argument to DISPLACED");
           return;
         };

        GT=CreateOrAugmentGTransformation(GT,DX);
      }
     else if ( !strcasecmp(Tokens[0],"ROTATED") )
      { 
        if (NumTokens!=5)
         { ErrMsg=strdup("ROTATED keyword requires exactly 4 arguments");
           return;
         };

        if (    1!=sscanf(Tokens[1],"%le",ZHat+0)
             || 1!=sscanf(Tokens[2],"%le",ZHat+1)
             || 1!=sscanf(Tokens[3],"%le",ZHat+2)
             || 1!=sscanf(Tokens[4],"%le",&Theta) )
         { ErrMsg=strdup("invalid argument to ROTATED");
           return;
         };

        GT=CreateOrAugmentGTransformation(GT,ZHat,Theta);

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

  /***************************************************************/
  /* note that we pass 0 for the material name here; this is     */
  /* because the material name might refer to a material that    */
  /* was defined on the fly in the .scuffgeo file, in which case */
  /* trying to assign that material would fail. instead, we      */
  /* store any user-specified material name in the MPName field  */
  /* inside the class body and leave it for whoever called this  */
  /* routine to process.                                         */
  /***************************************************************/
  InitRWGObject(MeshFileName, pLabel, 0, GT);

};

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
                     GTransformation *GT)
 { InitRWGObject(pMeshFileName, pLabel, Material, GT); }

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
                              GTransformation *OTGT)
{ 
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
  if (MeshFileName==0)
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
   RWGErrExit("file %s: invalid extension",MeshFileName);
  else if (!strcasecmp(p,"msh"))
   ReadGMSHFile(MeshFile,MeshFileName,OTGT);
  else if (!strcasecmp(p,"mphtxt"))
   ReadComsolFile(MeshFile,MeshFileName,OTGT);
  else
   RWGErrExit("file %s: unknown extension %s",MeshFileName,p);

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
  int Mu, Nu, np;

  MeshFileName=strdup("ByHand.msh");
  Label=strdup("ByHand");

  NumEdges=0;
  NumVertices=pNumVertices;
  NumPanels=pNumPanels;
  MP=new MatProp();
  GT=0;

  Vertices=(double *)RWGMalloc(3*NumVertices*sizeof(double));
  memcpy(Vertices,pVertices,3*NumVertices*sizeof(double *));

  Panels=(RWGPanel **)RWGMalloc(NumPanels*sizeof(RWGPanel *));
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
  if (GT) free(GT);
  
}

/*************************************************************/
/* Parse a "transformation" line and rotate/displace the     */
/* object accordingly.                                       */
/*                                                           */
/* The "transformation" line is a string containing one or   */
/* more of the following sections:                           */
/*                                                           */
/*  DISP x y z                    (displacement)             */
/*  ROT  Ax Ay Az Theta           (rotation)                 */
/*                                                           */
/* Rotations are through Theta degrees (NOT RADIANS) about   */
/* an axis passing through (0,0,0) and (Ax,Ay,Az).           */
/*                                                           */
/* Returns zero if the TransLine was valid, nonzero if the   */
/* transline was invalid.                                    */
/*                                                           */
/* Transformations are cumulative; two calls to Transform()  */
/* will build on each other.                                 */
/*************************************************************/
int RWGObject::Transform(const char *format, ...)
{ 
  /***************************************************************/
  /* fill out the TransLine with any printf-style arguments      */
  /***************************************************************/
  char TransLine[MAXSTR];
  va_list ap;
  va_start(ap,format);
  vsnprintf(TransLine,MAXSTR,format,ap);
  va_end(ap);

  /***************************************************************/
  /* break it up into tokens to be parsed ************************/
  /***************************************************************/
  char *Tokens[MAXTOK];
  int NumTokens=Tokenize(TransLine, Tokens, MAXTOK);

  if (NumTokens==0 || Tokens[0][0]=='#') 
   return 0; // a transformation that does nothing is considered valid

  /***************************************************************/
  /* parse the line **********************************************/
  /***************************************************************/
  int nt, nConv;
  GTransformation *DeltaGT=0; // 'new geometrical transformation'
  double DX[3], ZHat[3], Theta;
  for(nt=0; nt<NumTokens; nt++)
   { 
     /*--------------------------------------------------------------*/
     /* parse DISP element                                           */
     /*--------------------------------------------------------------*/
     if ( !strcasecmp(Tokens[nt],"DISP") )
      { 
        if ( nt+3 >= NumTokens )
         return 1;

        if (    1!=sscanf(Tokens[nt+1],"%le",DX+0)
             || 1!=sscanf(Tokens[nt+2],"%le",DX+1)
             || 1!=sscanf(Tokens[nt+3],"%le",DX+2) )
         return 1;
        
        DeltaGT=CreateOrAugmentGTransformation(DeltaGT, DX);

        nt+=3;
      }
     /*--------------------------------------------------------------*/
     /*- parse ROT element ------------------------------------------*/
     /*--------------------------------------------------------------*/
     else if ( !strcasecmp(Tokens[nt],"ROT") )
      { 
        if ( nt+4 >= NumTokens )
         return 1;

        if (    1!=sscanf(Tokens[nt+1],"%le",ZHat+0)
             || 1!=sscanf(Tokens[nt+2],"%le",ZHat+1)
             || 1!=sscanf(Tokens[nt+3],"%le",ZHat+2)
             || 1!=sscanf(Tokens[nt+4],"%le",&Theta) )
         return 1;
        
        DeltaGT=CreateOrAugmentGTransformation(DeltaGT, ZHat, Theta);

        nt+=4;
      }
     /*--------------------------------------------------------------*/
     /*- unknown keyword --------------------------------------------*/
     /*--------------------------------------------------------------*/
     else 
      return 1;

   }; // while( *p )

  /***************************************************************/
  /*- OK, now that we have constructed the transformation, we    */
  /*- need to apply it to all points whose coordinates we store  */
  /*- inside the RWGObject structure: vertices, edge centroids,  */
  /*- and panel centroids.                                       */
  /***************************************************************/
  /* vertices */
  ApplyGTransformation(DeltaGT, Vertices, NumVertices);

  /* edge centroids */
  int ne;
  for(ne=0; ne<NumEdges; ne++)
   ApplyGTransformation(DeltaGT, Edges[ne]->Centroid, 1);

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
  AugmentGTransformation(DeltaGT, GT);

  return 0;

}

void RWGObject::UnTransform()
{
  int Mu, Nu;
  int nv, ne, np;

  /***************************************************************/
  /* untransform vertices                                        */
  /***************************************************************/
  UnApplyGTransformation(GT, Vertices, NumVertices);

  /***************************************************************/
  /* untransform edge centroids                                  */
  /***************************************************************/
  for(ne=0; ne<NumEdges; ne++)
   UnApplyGTransformation(GT, Edges[ne]->Centroid, 1);

  /***************************************************************/
  /* reinitialize geometric data on panels (which takes care of  */ 
  /* transforming the panel centroids)                           */ 
  /***************************************************************/
  for(np=0; np<NumPanels; np++)
   InitRWGPanel(Panels[np], Vertices);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  ResetGTransformation(GT);

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

  P=(RWGPanel *)RWGMalloc(sizeof *P);
  P->VI[0]=iV1;
  P->VI[1]=iV2;
  P->VI[2]=iV3;

  InitRWGPanel(P, Vertices);

  return P;

} 

/***************************************************************/
/* Count how many vertices are shared in common between panels */
/* #np1 and #np2.                                              */
/* On return, Index1[0..nc] are the indices of the common      */
/* vertices in P1, and Index2[0..nc] are the indices of those  */
/* same vertices in P2, where nc is the number returned by     */
/* the routine.                                                */
/* Or, you can pass null pointers for Index1 and Index2, in    */
/* which case the routine still returns the number of common   */
/* vertices.                                                   */
/***************************************************************/
int RWGObject::CountCommonVertices(int np1, int np2, int *Index1, int *Index2)
{ 
  int i, j, nCommon;

  for(nCommon=0, i=0; i<3; i++)
   for(j=0; j<3; j++)
    if ( Panels[np1]->VI[i]==Panels[np2]->VI[j] )
     { if (Index1) Index1[nCommon]=i;
       if (Index2) Index2[nCommon]=j;
       nCommon++;
     };

  return nCommon;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
double RWGObject::GetOverlap(int neAlpha, int neBeta)
{ 
  RWGEdge *EAlpha=Edges[neAlpha], *EBeta=Edges[neBeta];
  double *V1, *V2, *QAlpha, *QBeta;
  double PreFac, Area, Term, Sum;
  int mu;

  V1=Vertices + 3*(EAlpha->iV1);
  V2=Vertices + 3*(EAlpha->iV2);
  PreFac=EAlpha->Length * EBeta->Length / (24.0);

  Sum=0.0;
  if ( EAlpha->iPPanel == EBeta->iPPanel )
   {  
      QAlpha=Vertices + 3*(EAlpha->iQP);
      QBeta =Vertices + 3*(EBeta->iQP);
      Area=Panels[EAlpha->iPPanel]->Area;

      for(Term=0.0, mu=0; mu<3; mu++)
       Term+= ( V1[mu] - QAlpha[mu] ) * ( V1[mu] + V2[mu] - 2.0*QBeta[mu] )
             +( V2[mu] - QAlpha[mu] ) * ( V2[mu] + QAlpha[mu] - 2.0*QBeta[mu] );

      Sum += PreFac * Term / Area;
   };
  if ( EAlpha->iPPanel == EBeta->iMPanel )
   {  
      QAlpha=Vertices + 3*(EAlpha->iQP);
      QBeta =Vertices + 3*(EBeta->iQM);
      Area=Panels[EAlpha->iPPanel]->Area;

      for(Term=0.0, mu=0; mu<3; mu++)
       Term+= ( V1[mu] - QAlpha[mu] ) * ( V1[mu] + V2[mu] - 2.0*QBeta[mu] )
             +( V2[mu] - QAlpha[mu] ) * ( V2[mu] + QAlpha[mu] - 2.0*QBeta[mu] );

      Sum -= PreFac * Term / Area;
   };
  if ( EAlpha->iMPanel == EBeta->iPPanel )
   {  
      QAlpha=Vertices + 3*(EAlpha->iQM);
      QBeta =Vertices + 3*(EBeta->iQP);
      Area=Panels[EAlpha->iMPanel]->Area;

      for(Term=0.0, mu=0; mu<3; mu++)
       Term+= ( V1[mu] - QAlpha[mu] ) * ( V1[mu] + V2[mu] - 2.0*QBeta[mu] )
             +( V2[mu] - QAlpha[mu] ) * ( V2[mu] + QAlpha[mu] - 2.0*QBeta[mu] );

      Sum -= PreFac * Term / Area;
   };
  if ( EAlpha->iMPanel == EBeta->iMPanel )
   {  
      QAlpha=Vertices + 3*(EAlpha->iQM);
      QBeta =Vertices + 3*(EBeta->iQM);
      Area=Panels[EAlpha->iMPanel]->Area;

      for(Term=0.0, mu=0; mu<3; mu++)
       Term+= ( V1[mu] - QAlpha[mu] ) * ( V1[mu] + V2[mu] - 2.0*QBeta[mu] )
             +( V2[mu] - QAlpha[mu] ) * ( V2[mu] + QAlpha[mu] - 2.0*QBeta[mu] );

      Sum += PreFac * Term / Area;
   };

  return Sum;

}
