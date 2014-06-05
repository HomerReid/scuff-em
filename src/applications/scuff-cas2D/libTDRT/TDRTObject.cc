/*
 * TDRTObject.cc -- implementation of some class methods for the TDRTObject 
 *               -- class 
 *
 * homer reid -- 11/2008
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "libhrutil.h"
#include "libMatProp.h"
#include "libTDRT.h"

#define MAXSTR 1000

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- various entry points to the TDRTObject class constructor    */
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
TDRTObject::TDRTObject(const char *pMeshFileName)
 { InitTDRTObject(pMeshFileName, 0, 0); }

TDRTObject::TDRTObject(const char *pMeshFileName, const char *pLabel, 
                       const char *pMaterial)
 { InitTDRTObject(pMeshFileName, pLabel, pMaterial); }

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- the actual body of the TDRTObject constructor: Create a    -*/
/*- TDRTObject from a mesh file describing a discretized       -*/
/*- object, with an optional label and material name.          -*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
void TDRTObject::InitTDRTObject(const char *pMeshFileName, const char *pLabel,
                                const char *Material)
{ 
  FILE *f;
  char FullPath[MAXSTR];
  char *p;
  int nv, niv, ns;
  int vi, ni;
  int *pNeighbors;

  /*------------------------------------------------------------*/
  /* try to open the mesh file.                                 */
  /*- we look in a couple different places:                     */
  /*-  1. current working directory                             */
  /*-  2. $(HOME)/geomsh/filename                               */
  /*-  3. SCUFF_MESH_PATH environment variable                  */
  /*------------------------------------------------------------*/
  MeshFileName=strdup(pMeshFileName);
  f=fopen(MeshFileName,"r");
  if ( !f )
   { sprintf(FullPath,"%s/geomsh/%s",getenv("HOME"),MeshFileName);
     f=fopen(FullPath,"r");
   };
  if ( !f )
   { char *str=getenv("SCUFF_MESH_PATH");
     if (str) 
      { sprintf(FullPath,"%s/%s",str,MeshFileName);
        f=fopen(FullPath,"r");
      };
   };
  if ( !f )
   ErrExit("could not find file %s",MeshFileName);
   
  if (pLabel)
   Label=strdup(pLabel);
  else
   Label=strdup("NoLabel");

  /*------------------------------------------------------------*/
  /*- initialize a new TDRTObject structure ----------------------*/
  /*------------------------------------------------------------*/
  NumVertices=NumSegments=NumIVs;
  Displacement[0]=Displacement[1]=0.0;
  ObjectWasRotated=0;
  Rotation=0.0;
  UnTransLine[0]=0;

  /*------------------------------------------------------------*/
  /*- initialize material properties  --------------------------*/
  /*------------------------------------------------------------*/
  MP=new MatProp(Material); // defaults to PEC if Material==0

  /*------------------------------------------------------------*/
  /* switch off based on filename extension to the appropriate  */
  /* routine for reading in vertices and segments               */
  /*------------------------------------------------------------*/
  p=strrchr(MeshFileName,'.');
  if ( p && !strcmp(p,".mphtxt") )
   ReadComsolFile(f);
  else if ( p && !strcmp(p,".msh") )
   ReadGMSHFile(f);
  else
   ErrExit("file %s: unsupported file extension",MeshFileName); 

  /*------------------------------------------------------------*/
  /* now go through the vertices and identify the 'interior     */
  /* vertices,' where an 'interior vertex' is a vertex that is  */
  /* an endpoint of exactly two line segments.                  */
  /* A vertex that is an endpoint of only one line segment is   */
  /* considered a boundary vertex and is ignored, whereas a     */
  /* vertex that is an endpoint of three or more line segments  */
  /* signifies that this topology is improper and triggers an   */
  /* abort.                                                     */
  /* After this code snippet, we have                           */
  /*  IVs[n] = index of nth interior vertex                     */
  /*  Neighbors[2*n], Neighbors[2*n+1] = indices of neighbors   */
  /*                                     of nth interior vertex */
  /*------------------------------------------------------------*/

  /* pNeighbors is a preliminary version of the Neighbors array */
  /* that we allocate here, use to construct the final version  */
  /* of the IVs and Neighbors arrays, and then deallocate.      */ 
  pNeighbors=(int *)malloc(2*NumVertices*sizeof(int));
  for(nv=0; nv<NumVertices; nv++)
   pNeighbors[2*nv]=pNeighbors[2*nv+1]=-1;

  /* pass 1: populate pNeighbors array  */
  for(ns=0; ns<NumSegments; ns++)
   { 
     vi=Segments[2*ns];                // 'vertex index' 
     ni=Segments[2*ns+1];              // 'neighbor index'
     if ( pNeighbors[2*vi] == -1 ) 
      pNeighbors[2*vi]=ni;
     else if ( pNeighbors[2*vi+1] == -1 ) 
      pNeighbors[2*vi+1]=ni;
     else
      ErrExit("%s: bad topology: vertex %i touches 3 segments (%i,%i,%i)",
               MeshFileName,vi,ni,pNeighbors[2*vi],pNeighbors[2*vi+1],ni);

     vi=Segments[2*ns+1];                // 'vertex index' 
     ni=Segments[2*ns];                  // 'neighbor index'
     if ( pNeighbors[2*vi] == -1 ) 
      pNeighbors[2*vi]=ni;
     else if ( pNeighbors[2*vi+1] == -1 ) 
      pNeighbors[2*vi+1]=ni;
     else
      ErrExit("%s: bad topology: vertex %i touches 3 segments (%i,%i,%i)",
               MeshFileName,vi,ni,pNeighbors[2*vi],pNeighbors[2*vi+1],ni);
   }; 

  /* pass 2: count vertices with exactly 2 neighbors */
  NumIVs=0;
  for(nv=0; nv<NumVertices; nv++)
   if ( pNeighbors[2*nv+1]!=-1 ) 
    NumIVs++;
  NumBFs = MP->IsPEC() ? 2*NumIVs: 4*NumIVs;

  /* pass 3: construct IVs and Neighbors tables */
  IVs=(int *)malloc(NumIVs*sizeof(int));
  Neighbors=(int *)malloc(2*NumIVs*sizeof(int));
  for(niv=nv=0; nv<NumVertices; nv++)
   { if (pNeighbors[2*nv+1]==-1) 
      continue;
     IVs[niv]=nv;
     memcpy(Neighbors+2*niv,pNeighbors+2*nv,2*sizeof(int));
     niv++;
   };
  free(pNeighbors);

  /* do a quick sanity check to look for zero-length segments */
  for(niv=0; niv<NumIVs; niv++)
   if ( VecDistance(Vertices+2*IVs[niv], Vertices+2*Neighbors[2*niv])<1.0e-12 ) 
    ErrExit("%s: bad topology: line segment between vertices %i,%i too short",
             MeshFileName, IVs[niv],Neighbors[2*niv]);

}

/***************************************************************/
/* TDRTObject class destructor  **********************************/
/***************************************************************/
TDRTObject::~TDRTObject()
{ 
  free(Vertices);
  free(Segments);
  free(IVs);
  free(Neighbors);
  free(MeshFileName);
  free(MP);
}

/***************************************************************/
/* Apply a rigid displacement to all vertices in the object ****/
/***************************************************************/
void TDRTObject::Displace(double *DX)
{ 
  int nv;

  for(nv=0; nv<NumVertices; nv++)
   VecPlusEquals(Vertices + 2*nv,1.0,DX);

  VecPlusEquals(Displacement, 1.0, DX);
}

void TDRTObject::Displace(double x, double y)
{ double DX[2];
  DX[0]=x;
  DX[1]=y;
  Displace(DX);
}

/***************************************************************/
/* Undisplace the object back to where it started out. *********/
/***************************************************************/
void TDRTObject::UnDisplace()
{ 
  Displace(-Displacement[0], -Displacement[1]);
} 

/***************************************************************/
/* Rotate the object through an angle of Theta degrees about   */
/* the Z axis.                                                 */
/***************************************************************/
void TDRTObject::Rotate(double Theta)
{ 
  int nv;
  double CT, ST;
  double *V, RV[2];

  /***************************************************************/
  /* rotate all vertices   ***************************************/
  /***************************************************************/
  CT=cos(Theta*M_PI/180.0); 
  ST=sin(Theta*M_PI/180.0);
  for(nv=0; nv<NumVertices; nv++)
   { 
     V=Vertices + 2*nv;
     RV[0]=CT*V[0] - ST*V[1];
     RV[1]=ST*V[0] + CT*V[1];
     memcpy(V,RV,2*sizeof(double));
   };

  Rotation+=Theta;
  ObjectWasRotated=1;
  
}

/***************************************************************/
/* Unrotate the object back to where it started out.           */
/***************************************************************/
void TDRTObject::UnRotate()
{ Rotate(-Rotation); }

/*************************************************************/
/* Parse a "transformation" line and rotate/displace the     */
/* object accordingly.                                       */
/*                                                           */
/* The "transformation" line is a string containing one or   */
/* more of the following sections:                           */
/*                                                           */
/*  DISP x y                 (displacement)                  */
/*  ROT  Theta               (rotation)                      */
/*                                                           */
/* Rotations are through Theta degrees (NOT RADIANS) about   */
/* the Z axis.                                              */
/*                                                           */
/* Returns zero if the TransLine was valid, nonzero if the   */
/* transline was invalid.                                    */
/*                                                           */
/* Note: As we proceed through this routine, we construct an */
/* 'UnTransLine' that will bring about the inverse of the    */
/* given transformation, i.e. that will restore the object   */
/* to its position prior to application of this transform.   */ 
/*************************************************************/
int TDRTObject::Transform(char *TransLine)
{ 
  int nRead, nConv;
  char *p, Token[100], TempStr[1000];
  double D[2], Theta;

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
  ObjectWasRotated=0;
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
     /* parse DISP element */
     /*--------------------------------------------------------------*/
     if ( !StrCaseCmp(Token,"DISP") )
      { 
        nConv=sscanf(p,"%le %le %n",D,D+1,&nRead);
        p+=nRead;

        if ( nConv!=2 )
         return 1;
        
        Displace(D);
        sprintf(TempStr,"DISP %.15e %.15e %s",-D[0],-D[1],UnTransLine);
        strcpy(UnTransLine,TempStr);
      }
     /*--------------------------------------------------------------*/
     /*- parse ROT element ------------------------------------------*/
     /*--------------------------------------------------------------*/
     else if ( !StrCaseCmp(Token,"ROT") )
      { 
        nConv=sscanf(p,"%le %n",&Theta,&nRead);
        p+=nRead;

        if ( nConv!=1 )
         return 1;
        
        Rotate(Theta);
        sprintf(TempStr,"ROT %.15e %s",-Theta,UnTransLine);
        strcpy(UnTransLine,TempStr);
      }
     else
       return 1;

   }; // while( *p )

  return 0;

}

void TDRTObject::UnTransform()
{
  char UTLCopy[1000];
  
  if (UnTransLine[0]==0)
   return;

  strcpy(UTLCopy,UnTransLine);
  Transform(UTLCopy);
  UnTransLine[0]=0;

}
