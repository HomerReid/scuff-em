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

#include <libTriInt.h>

#define MAXREFPTS 10

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
#if 0
RWGObject::RWGObject(FILE *f, const char *pLabel, int *LineNum)
{ 
  char Line[MAXSTR];
  int nt, nTokens;
  char *p, *Tokens[50];
    // MESHFILE 
    // MATERIAL
    // INSIDE 
    // DISPLACED
    // ROTATED
    
    MP=0;
    ContainingObject=0;
    MPName=0;
    ContainingObjectName=0;

    MPName=strdup(
    ContainingObjectName=strdup()

 };
#endif

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- various entry points to the RWGObject class constructor  --*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
RWGObject::RWGObject(const char *pMeshFileName)
 { InitRWGObject(pMeshFileName, 0, 0, 0, 0); }

RWGObject::RWGObject(const char *pMeshFileName, 
                     const char *pLabel,
                     const char *Material,
                     const char *RotFileName, 
                     double *DX)
 { InitRWGObject(pMeshFileName, pLabel, Material, RotFileName, DX); }

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
                              const char *RotFileName,
                              double *DX)
{ 
  char *p, buffer[1000];
  double RotMat[9], *pRotMat;
  FILE *MeshFile, *RotFile, *SPPIDTFile;
  int Mu, Nu, LineNum, i, np;
  char SPPIDTFileName[200];

  /*------------------------------------------------------------*/
  /*- try to open the mesh file.                                h*/
  /*- we look in a couple different places:                     */
  /*-  1. current working directory                             */
  /*-  2. $(HOME)/geomsh/                                       */
  /*------------------------------------------------------------*/
  MeshFile=fopen(pMeshFileName,"r");    /* try bare filename */
  if ( !MeshFile )
   { sprintf(buffer,"%s/geomsh/%s",getenv("HOME"),pMeshFileName);
     MeshFile=fopen(buffer,"r");
   };
  if (!MeshFile)
   RWGErrExit("could not find file %s in mesh search path",pMeshFileName);
  MeshFileName=strdup(pMeshFileName);
   
  if (pLabel)
   Label=strdup(pLabel);
  else
   Label=strdup("NoLabel");

  /*------------------------------------------------------------*/
  /*- read in the rotation file if one was specified            */
  /*------------------------------------------------------------*/
  pRotMat=0;
  if ( RotFileName && RotFileName[0]!=0 )
   { 
     /* attempt to open the file. we look in the same places as above */
     RotFile=fopen(RotFileName,"r"); 
     if ( !RotFile )  
      { sprintf(buffer,"%s/geomsh/%s",getenv("HOME"),RotFileName);
        RotFile=fopen(buffer,"r");
      };
     if (!RotFile)
      RWGErrExit("could not find file %s in mesh search path",RotFileName);

     /* read the three rows of the rotation matrix, one per line */
     LineNum=i=0;
     while( fgets(buffer,100,RotFile) )
      { 
        /* skip blank lines and comments */
        LineNum++;
        p=buffer; 
        while ( isspace(*p) )
         p++;
        if ( *p==0 || *p=='#') 
         continue;

        if(i==3) 
         RWGErrExit("too many lines in file %s",RotFileName);

        if ( 3!=sscanf(buffer,"%le %le %le",RotMat+i,RotMat+i+3,RotMat+i+6) )
         RWGErrExit("%s:%i: syntax error",RotFileName,LineNum);
        i++;

      };

      if(i!=3) 
       RWGErrExit("too few lines in file %s",RotFileName);

     fclose(RotFile);

     pRotMat=RotMat;
   };
   
  /*------------------------------------------------------------*/
  /*- initialize a new RWGObject structure ---------------------*/
  /*------------------------------------------------------------*/
  NumEdges=NumPanels=NumVertices=NumRefPts=0;
  MP=new MatProp(Material);
  SPPIDTable=0;
  ContainingObject=0;

  /* initialize the transformation to the identity transformation */
  for(Mu=0; Mu<3; Mu++)
   for(Nu=0; Nu<3; Nu++)
    MT[Mu][Nu] = (Mu==Nu) ? 1.0 : 0.0;
  memset(VT,0,3*sizeof(double));

  /*------------------------------------------------------------*/
  /*- At this point we switch off based on the file type:      -*/
  /*-  1. file extension=.msh    --> ReadGMSHFile              -*/
  /*-  2. file extension=.mphtxt --> ReadComsolFile            -*/
  /*------------------------------------------------------------*/
  p=strrchr(MeshFileName,'.');
  if (!p)
   RWGErrExit("file %s: invalid extension",MeshFileName);
  else if (!strcmp(p,".msh"))
   ReadGMSHFile(MeshFile,MeshFileName,pRotMat,DX);
  else if (!strcmp(p,".mphtxt"))
   ReadComsolFile(MeshFile,MeshFileName,pRotMat,DX);
  else
   RWGErrExit("file %s: unknown extension %s",MeshFileName,p);

  /*------------------------------------------------------------*/
  /*- Now that we have put the panels in an array, go through  -*/
  /*- and fill in the Index field of each panel structure.     -*/
  /*------------------------------------------------------------*/
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

  /* initialize the transformation to the identity transformation */
  for(Mu=0; Mu<3; Mu++)
   for(Nu=0; Nu<3; Nu++)
    MT[Mu][Nu] = (Mu==Nu) ? 1.0 : 0.0;
  memset(VT,0,3*sizeof(double));

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
/* RWGObject destructor.                                      */
/***************************************************************/
RWGObject::~RWGObject()
{ 
  int np, ne;

  for(np=0; np<NumPanels; np++)
   free(Panels[np]);
  free(Panels);

  for(ne=0; ne<NumEdges; ne++)
   free(Edges[ne]);
  free(Edges);

  free(Vertices);

  delete MP;

  if (MeshFileName)
   free(MeshFileName);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RWGObject::TransformPoint(double X[3])
{ 
  double XP[3];
  int Mu, Nu;
  
  /* XP <- MT*X */
  for(Mu=0; Mu<3; Mu++)
   for(XP[Mu]=0.0, Nu=0; Nu<3; Nu++)
    XP[Mu]+=MT[Mu][Nu]*X[Nu];

  /* X <- XP + VT */
  for(Mu=0; Mu<3; Mu++)
   X[Mu]=XP[Mu]+VT[Mu];

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RWGObject::UnTransformPoint(double X[3])
{ 
  double XP[3];
  int Mu, Nu;
  
  /* XP <- X - VT */
  for(Mu=0; Mu<3; Mu++)
   XP[Mu]=X[Mu]-VT[Mu];

  /* X <- MT^{-1}*XP */
  for(Mu=0; Mu<3; Mu++)
   for(X[Mu]=0.0, Nu=0; Nu<3; Nu++)
    X[Mu]+=MT[Nu][Mu]*XP[Nu];

}

/***************************************************************/
/* Construct the 3x3 matrix that represents a rotation of      */
/* Theta degrees (note we interpret Theta in DEGREES, NOT      */
/* RADIANS!) about the axis specified by cartesian coordinates */
/* ZHat.                                                       */
/* Algorithm:                                                  */
/*  1. Construct matrix M1 that rotates Z axis into alignment  */
/*     with ZHat.                                              */
/*  2. Construct matrix M2 that rotates through Theta about    */ 
/*     Z axis.                                                 */
/*  3. Construct matrix M=M1^{-1}*M2*M1=M1^T*M2*M1.            */  
/***************************************************************/
void ConstructRotationMatrix(double *ZHat, double Theta, double M[3][3])
{ 
  int Mu, Nu, Rho;
  double ct, st, cp, sp, CT, ST;
  double M2M1[3][3], M2[3][3], M1[3][3];

  VecNormalize(ZHat);
 
  /* construct M1 */
  ct=ZHat[2];
  st=sqrt(1.0-ct*ct);
  cp= ( st < 1.0e-8 ) ? 1.0 : ZHat[0] / st;
  sp= ( st < 1.0e-8 ) ? 0.0 : ZHat[1] / st;
  M1[0][0]=ct*cp;  M1[0][1]=ct*sp;   M1[0][2]=-st;
  M1[1][0]=-sp;    M1[1][1]=cp;      M1[1][2]=0.0;
  M1[2][0]=st*cp;  M1[2][1]=st*sp;   M1[2][2]=ct;

  /* construct M2 */
  CT=cos(Theta*M_PI/180.0);
  ST=sin(Theta*M_PI/180.0);
  M2[0][0]=CT;      M2[0][1]=-ST;     M2[0][2]=0.0;
  M2[1][0]=ST;      M2[1][1]=CT;      M2[1][2]=0.0;
  M2[2][0]=0.0;     M2[2][1]=0.0;     M2[2][2]=1.0;

  /* M2M1 <- M2*M1 */
  for(Mu=0; Mu<3; Mu++)
   for(Nu=0; Nu<3; Nu++)
    for(M2M1[Mu][Nu]=0.0, Rho=0; Rho<3; Rho++)
     M2M1[Mu][Nu] += M2[Mu][Rho] * M1[Rho][Nu];

  /* M <- M1^T * M2M1 */
  for(Mu=0; Mu<3; Mu++)
   for(Nu=0; Nu<3; Nu++)
    for(M[Mu][Nu]=0.0, Rho=0; Rho<3; Rho++)
     M[Mu][Nu] += M1[Rho][Mu]*M2M1[Rho][Nu];
  
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
/* Implementation notes:                                     */
/*                                                           */
/*  a. We keep track of transformations by maintaining a     */
/*     running matrix MT and vector VT such that XT, the     */
/*     vector of coordinates of a vertex in the object after */
/*     the transformation, is given by                       */
/*      XT = MT*X + VT                                       */
/*     where X is the vector of coordinates of the vertex in */
/*     the original geometry, i.e. as given in the .msh file.*/
/*                                                           */
/*  b. initially, we set M=unit 3x3 matrix and V=[0 0 0].    */
/*                                                           */
/*  c. as we parse the transformation line, we augment       */
/*     M and V as follows:                                   */
/*      for every DISP element: V --> V + DISP               */
/*      for every ROT  element: M --> ROT*M, V --> ROT*V     */
/*                                                           */
/*************************************************************/
int RWGObject::Transform(const char *pTransLine, ...)
{ 
  int nRead, nConv, ne, nv, np;
  char *p, TransLine[1000], Token[100], TempStr[1000];
  double D[3], Theta;
  
  int Mu, Nu, Rho;
  double MTP[3][3], MTPP[3][3];
  double VTP[3];

  va_list ap;
  va_start(ap,pTransLine);
  vsnprintf(TransLine,1000,pTransLine,ap);
  va_end(ap);

  /* skip blank lines and comments */
  p=TransLine; 
  while( isspace(*p) )
   p++;
  if ( *p==0 || *p=='#' )
   return 0;

  /*
   * parse transformation line
   */
  nRead=0;
  while( *p )
   {
     /*--------------------------------------------------------------*/
     /*- read next token off of line --------------------------------*/
     /*--------------------------------------------------------------*/
     nConv=sscanf(p,"%s%n",Token,&nRead);  
     p+=nRead;  
     if ( nConv<=0 || Token[0]=='\n' )  
      break;

     /*--------------------------------------------------------------*/
     /* parse DISP element                                           */
     /*--------------------------------------------------------------*/
     if ( !strcasecmp(Token,"DISP") )
      { 
        nConv=sscanf(p,"%le %le %le %n",D,D+1,D+2,&nRead);
        p+=nRead;

        if ( nConv!=3 )
         return 1;
        
        VecPlusEquals(VT,1.0,D);
      }
     /*--------------------------------------------------------------*/
     /*- parse ROT element ------------------------------------------*/
     /*--------------------------------------------------------------*/
     else if ( !strcasecmp(Token,"ROT") )
      { 
        nConv=sscanf(p,"%le %le %le %le %n",D,D+1,D+2,&Theta,&nRead);
        p+=nRead;

        if ( nConv!=4 )
         return 1;

        /* construct the 3x3 rotation matrix MTP for this rotation */
        ConstructRotationMatrix(D, Theta, MTP);

        /* left-multiply MT by the new matrix MTP */
        for(Mu=0; Mu<3; Mu++)
         for(Nu=0; Nu<3; Nu++)
          for(MTPP[Mu][Nu]=0.0, Rho=0; Rho<3; Rho++)
           MTPP[Mu][Nu]+=MTP[Mu][Rho]*MT[Rho][Nu];
        for(Mu=0; Mu<3; Mu++)
         for(Nu=0; Nu<3; Nu++)
          MT[Mu][Nu]=MTPP[Mu][Nu];

        /* and we also have to left-multiply VT by MTP as well */
        for(Mu=0; Mu<3; Mu++)
         for(VTP[Mu]=0.0, Nu=0; Nu<3; Nu++)
          VTP[Mu]+=MTP[Mu][Nu]*VT[Nu];

        for(Mu=0; Mu<3; Mu++)
         VT[Mu]=VTP[Mu];

      }
     /*--------------------------------------------------------------*/
     /*- unknown keyword --------------------------------------------*/
     /*--------------------------------------------------------------*/
     else 
      return 1;

   }; // while( *p )

  /*--------------------------------------------------------------*/
  /*- OK, now that we have constructed the matrix and vector that */
  /*- define the transformation, we need to apply it to all points*/
  /*- whose coordinates we store inside the RWGObject structure:  */
  /*- vertices, edge centroids, and panel centroids.              */
  /*--------------------------------------------------------------*/
  /* vertices */
  for(nv=0; nv<NumVertices; nv++)
   TransformPoint(Vertices+3*nv);

  /* edge centroids */
  for(ne=0; ne<NumEdges; ne++)
   TransformPoint(Edges[ne]->Centroid);

  /***************************************************************/
  /* reinitialize geometric data on panels (which takes care of  */ 
  /* transforming the panel centroids)                           */ 
  /***************************************************************/
  for(np=0; np<NumPanels; np++)
   InitRWGPanel(Panels[np], Vertices);

  return 0;

}

void RWGObject::UnTransform()
{
  int Mu, Nu;
  int nv, ne, np;

  /***************************************************************/
  /* untransform vertices                                        */
  /***************************************************************/
  for(nv=0; nv<NumVertices; nv++)
   UnTransformPoint(Vertices+3*nv);

  /***************************************************************/
  /* untransform edge centroids                                  */
  /***************************************************************/
  for(ne=0; ne<NumEdges; ne++)
   UnTransformPoint(Edges[ne]->Centroid);

  /***************************************************************/
  /* reinitialize geometric data on panels (which takes care of  */ 
  /* transforming the panel centroids)                           */ 
  /***************************************************************/
  for(np=0; np<NumPanels; np++)
   InitRWGPanel(Panels[np], Vertices);

  /***************************************************************/
  /* now that we are back to where we started, the transformation*/
  /* that takes us back to where we started is the identity      */
  /* transformation                                              */
  /***************************************************************/
  for(Mu=0; Mu<3; Mu++)
   for(Nu=0; Nu<3; Nu++)
    MT[Mu][Nu]=(Mu==Nu) ? 1.0 : 0.0;
  memset(VT,0,3*sizeof(double));

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
