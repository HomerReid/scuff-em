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
  MeshFileName=0;

  /***************************************************************/
  /* read lines from the file one at a time **********************/
  /***************************************************************/
  char Line[MAXSTR];
  char MaterialName[MAXSTR];
  int NumTokens, TokensConsumed;
  char *Tokens[MAXTOK];
  int ReachedTheEnd=0;
  GTransformation *OTGT=0; // 'one-time geometrical transformation'
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
        MeshFileName=strdup(Tokens[1]);
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
        OTGT=CreateOrAugmentGTransformation(OTGT, Tokens, NumTokens, 
                                            &ErrMsg, &TokensConsumed);
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

  InitRWGObject(MeshFileName, pLabel, MaterialName, OTGT);

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
                     GTransformation *OTGT)
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
  if (GT) free(GT);
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RWGObject::Transform(GTransformation *DeltaGT)
{ 
  /***************************************************************/
  /*- first apply the transformation to all points whose         */
  /*- coordinates we store inside the RWGObject structure:       */
  /*- vertices, edge centroids, and panel centroids.             */
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
  /* update the internally stored GTransformation ****************/
  /***************************************************************/
  GT=CreateOrAugmentGTransformation(GT, DeltaGT);

}

void RWGObject::Transform(char *format,...)
{
  va_list ap;
  char buffer[MAXSTR];
  va_start(ap,format);
  vsnprintf(buffer,MAXSTR,format,ap);

  char *ErrMsg;
  GTransformation *GT=CreateOrAugmentGTransformation(0, buffer, &ErrMsg);
  if (ErrMsg)
   ErrExit(ErrMsg);
  Transform(GT);
  free(GT);
  
}

void RWGObject::UnTransform()
{
  if (!GT)
   return;
 
  /***************************************************************/
  /* untransform vertices                                        */
  /***************************************************************/
  UnApplyGTransformation(GT, Vertices, NumVertices);

  /***************************************************************/
  /* untransform edge centroids                                  */
  /***************************************************************/
  int ne;
  for(ne=0; ne<NumEdges; ne++)
   UnApplyGTransformation(GT, Edges[ne]->Centroid, 1);

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

  P=(RWGPanel *)mallocEC(sizeof *P);
  P->VI[0]=iV1;
  P->VI[1]=iV2;
  P->VI[2]=iV3;

  InitRWGPanel(P, Vertices);

  return P;

} 

/***************************************************************/
/* Compute the overlap integral between the RWG basis functions*/
/* associated with two edges in an RWG object.                 */
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

/***************************************************************/
/* Compute all nonzero entries in the row of the overlap       */
/* matrix corresponding to basis function #neAlpha.            */
/*                                                             */
/* (The neA, neB entry in this matrix is the overlap between   */
/*  basis functions neA and neB, i.e.                          */
/*   O_{neA, neB} = \int f_{neA} \cdot f_{neB}                 */
/*  where the integral is over the intersection of the supports*/
/*  of the basis functions).                                   */
/*                                                             */
/* The return value is nnz, the number of nonzero entries in   */
/* the row. (nnz is between 1 and 5.) On return,               */
/* ColIndices[0..(nnz-1)] are the indices of the nonzero       */
/* columns, and OEntries[0..(nnz-1)] are the corresponding     */
/* entries in the overlap matrix.                              */
/*                                                             */
/* If OTimesEntries is non-NULL on entry, then on return       */
/* OTimesEntries[0..4] are the values of the corresponding     */
/* entries in the 'crossed-overlap matrix,' whose entries are  */
/* defined by                                                  */
/*                                                             */
/*   O^x_{neA, neB} = \int f_{neA} \cdot (nHat \times f_{neB}) */
/*                                                             */
/* (Note that OEntries must point to a buffer large enough to  */
/*  hold 5 doubles, as must OTimesEntries if it is non-NULL.)  */
/***************************************************************/
#if 0
int RWGObject::GetOverlaps(int neAlpha, int *ColIndices,
                           double *OEntries, double *OTimesEntries)
{ 
  /*--------------------------------------------------------------*/
  /*- precompute some geometric quantities -----------------------*/
  /*--------------------------------------------------------------*/
  RWGEdge *EAlpha=Edges[neAlpha]; 
  double *QP = Vertices + 3*(EAlpha->iQP);
  double *V1 = Vertices + 3*(EAlpha->iV1);
  double *V2 = Vertices + 3*(EAlpha->iV2);
  double *QM = Vertices + 3*(EAlpha->iQM);

  double AP[3], AM[3], B[3], APpBO2[3], AMpBO2[3];
  double Mag2APpBO2, Mag2AMpBO2, Mag2B;
  int i;

  Mag2APpBO2 = Mag2AMpBO2 = Mag2B = 0.0;
  for(i=0; i<3; i++)
   { AP[i]       = V1[i] - QP[i];
     AM[i]       = V1[i] - QM[i];
     B[i]        = V2[i] - V1[i];
     APpBO2[i]   = AP[i] + 0.5*B[i];
     APpBO2[i]   = AP[i] + 0.5*B[i];
     Mag2APpBO2 += APpBO2[i]*APpBO2[i];
     Mag2AMpBO2 += AMpBO2[i]*AMpBO2[i];
     Mag2B      += Mag2B[i]*Mag2B[i];
   };

  /*--------------------------------------------------------------*/
  /*- first nonzero element (diagonal element) -------------------*/
  /*--------------------------------------------------------------*/

  ColIndices[0]=neAlpha;

  // the diagonal element of the crossed-overlap matrix is 0
  if (OTimesEntries) 
   OTimesEntries[0]=0.0;

  // the indices of the next four nonzero columns are the edge indices
  // of the other edges in the two panels that share edge neAlpha

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
#endif

} // namespace scuff
