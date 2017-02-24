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
 * RWGPorts.cc -- libRWG module for handling of 'ports' in 
 *             -- RF/microwave structures
 *                                         
 * homer reid  -- 9/2011                                        
 *                                         
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <libhrutil.h>
#include <libhmat.h>
#include <libscuff.h>
#include <libscuffInternals.h>
#include <libSGJC.h>

#include "RWGPorts.h"

#define ABSTOL 1.0e-8
#define RELTOL 1.0e-4
#define FREQ2OMEGA (2.0*M_PI/300.0)
#define II cdouble(0.0,1.0)

#define MAXEDGES 100           // max # of edges in (each half of) a port
#define MAXPOLYGONVERTICES 30

#ifdef USE_OPENMP
#  include <omp.h>
#endif

using namespace scuff;

/***************************************************************/
/* slightly tricky: the panels associated with what we call the*/
/* 'positive edge' of an RWGPort are actually the *negative*   */
/* of the two panels that define an RWG basis function         */
/* straddling the port (in the sense that on these panels the  */
/* current flows *toward* the Q vertex and is sunk there, as   */
/* opposed to being sourced from there). this confounds the    */
/* usual sense of positive and negative weights for the        */
/* contributions of the two panels that make up the RWG basis  */
/* function.                                                   */  
/*                                                             */  
/* to be precise: the panels indexed by the 'PPanelIndices'    */  
/* array in the Port structure are the panels on the positive  */  
/* side of the port, which makes them the negative of the two  */  
/* panels in the RWG basis function to which they belong.      */  
/***************************************************************/
RWGPort *CreatePort(RWGSurface *PSurface, int NumPEdges, int *PEdgeIndices, 
                    RWGSurface *MSurface, int NumMEdges, int *MEdgeIndices)
{
  RWGPort *P=(RWGPort *)malloc(sizeof *P);
  RWGEdge *E;

  P->PSurface       = PSurface;
  P->NumPEdges     = NumPEdges;
  P->PPanelIndices = (int *)malloc(NumPEdges*sizeof(int));
  P->PPaneliQs     = (int *)malloc(NumPEdges*sizeof(int));
  P->PLengths      = (double *)malloc(NumPEdges*sizeof(double));

  P->MSurface       = MSurface;
  P->NumMEdges     = NumMEdges;
  P->MPanelIndices = (int *)malloc(NumMEdges*sizeof(int));
  P->MPaneliQs     = (int *)malloc(NumMEdges*sizeof(int));
  P->MLengths      = (double *)malloc(NumMEdges*sizeof(double));
  
  int ne, EI;

  /***************************************************************/
  /* gather info on the RWG edges that define the 'positive' side*/
  /* of the port                                                 */
  /***************************************************************/
  P->PPerimeter=0.0;
  for(ne=0; ne<NumPEdges; ne++)
   { 
     EI = PEdgeIndices[ne];
     if ( EI<0 || EI>PSurface->NumExteriorEdges )
      ErrExit("object %s(%s): invalid exterior edge index %i",
               PSurface->Label,PSurface->MeshFileName,EI);
     E=PSurface->ExteriorEdges[EI];
     P->PPanelIndices[ne]=E->iPPanel;
     P->PPaneliQs[ne]=E->PIndex;
     P->PLengths[ne]=E->Length;
     P->PPerimeter+=E->Length;

     if (ne==0) 
      memcpy(P->PRefPoint, E->Centroid, 3*sizeof(double));
   };

  /***************************************************************/
  /* gather info on the RWG edges that define the 'negative' side*/
  /* of the port                                                 */
  /***************************************************************/
  P->MPerimeter=0.0;
  for(ne=0; ne<NumMEdges; ne++)
   {    
     EI=MEdgeIndices[ne]; // 'minus edge index'
     if ( EI<0 || EI>MSurface->NumExteriorEdges )
      ErrExit("object %s(%s): invalid exterior edge index %i",
               MSurface->Label,MSurface->MeshFileName,EI);
     E=MSurface->ExteriorEdges[EI];
     P->MPanelIndices[ne]=E->iPPanel;
     P->MPaneliQs[ne]=E->PIndex;
     P->MLengths[ne]=E->Length;
     P->MPerimeter+=E->Length;

     if (ne==0) 
      memcpy(P->MRefPoint, E->Centroid, 3*sizeof(double));
   };

  return P;

} 

/***************************************************************/
/* compute the normal to the triangle defined by three points  */
/***************************************************************/
void GetNormal(double *V1, double *V2, double *V3, double *Z)
{ 
  double A[3], B[3];
  VecSub(V2, V1, A);
  VecSub(V3, V2, B);
  VecCross(A, B, Z);
} 

/***************************************************************/
/* return the cosine of the angle between two 3-vectors ********/
/***************************************************************/
double CosAngle(double *L1, double *L2)
{
  double L1L1 = L1[0]*L1[0] + L1[1]*L1[1] + L1[2]*L1[2];
  double L2L2 = L2[0]*L2[0] + L2[1]*L2[1] + L2[2]*L2[2];
  double L1L2 = L1[0]*L2[0] + L1[1]*L2[1] + L1[2]*L2[2];
  return L1L2 / sqrt(L1L1*L2L2);
}

/***************************************************************/
/* return 1 if X lies on the line segment connecting V1 and V2 */
/*  (more specifically: if the shortest distance from X to the */
/*   line segment is <1e-6 * the length of the line segment)   */
/***************************************************************/
bool PointOnLineSegment(double *X, double *V1, double *V2)
{
  double A[3], B[3];

  VecSub( X, V1, A);
  VecSub(V2, V1, B);
  double A2  =   A[0]*A[0] + A[1]*A[1] + A[2]*A[2];
  double B2  =   B[0]*B[0] + B[1]*B[1] + B[2]*B[2];

  if ( A2<1.0e-12*B2 || fabs(A2-B2)<1.0e-12*B2 )
   return true;
  if ( A2 > B2 )
   return false;

  // the cosine of the angle between A and B is AdB/sqrt(A2*B2)
  // define the point to lie on the segment if the shortest
  // distance from the point to the segment is < 1e-6*length of segment
  // sqrt(A2) * sin(angle) < 1e-6*sqrt(B2)
  // A2 * (1-AdB^2/A2B2) < 1e-12*B2
  // (A2*B2 - AdB^2)  < 1e-12*B2*B2

  double AdB =   A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
  if (AdB<0.0) 
   return false;
  if ( (A2*B2-AdB*AdB) > 1.0e-12*B2*B2 ) 
   return false;

  return true;
 
}

/***************************************************************/
/* return true if the point X lies within the polygon defined  */
/* by vertices V. (If V has only two vertices, i.e. it is a    */
/* line segment, this amounts to testing whether X lies on     */
/* that line segment.)                                         */
/*                                                             */
/* X[0..2] = X,Y,Z coordinates of evaluation point             */
/*                                                             */
/* V[0..2] = X,Y,Z coordinates of polygon vertex 1             */
/* V[3..5] = X,Y,Z coordinates of polygon vertex 2             */
/* etc.                                                        */
/*                                                             */
/* V must have (at least) 3*NumVertices entries.               */
/*                                                             */
/* All vertices should be coplanar; this is not checked.       */
/*                                                             */
/* If X coincides with one of the vertices, then the return    */
/* value is true.                                              */
/***************************************************************/
bool PointInPolygon(double *X, double *V, int NumVertices)
{
  if (NumVertices==0) 
   return false; 
  else if (NumVertices==1)
   return VecEqualFloat(X,V);
  else if (NumVertices==2)
   return PointOnLineSegment(X, V+0, V+3);

  /*--------------------------------------------------------------*/
  /* 3 or more vertices; check that the normals to the triangles  */
  /* formed by X and all sequential pairs of vertices are parallel*/
  /*--------------------------------------------------------------*/
  double FirstZ[3], ThisZ[3], AD;
  GetNormal(V+0, V+3, X, FirstZ);
  for(int nv=1; nv<NumVertices; nv++)
   { GetNormal(V + 3*nv, V+3*((nv+1)%NumVertices), X, ThisZ);
     AD=fabs(1.0-CosAngle(ThisZ,FirstZ)); // 'angle deviation'
     if (AD>1.0e-12)
      return false;
   };
  return true; 

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int FindEdgesInPolygon(RWGSurface *S, 
                       double *PolygonVertices, int NumVertices,
                       int *NewEdgeIndices)
{
  RWGEdge *E;
  double *V1, *V2;
  int Count=0;
  for(int nei=0; nei<S->NumExteriorEdges; nei++)
   { 
     E=S->ExteriorEdges[nei];
     V1 = S->Vertices + 3*E->iV1;
     V2 = S->Vertices + 3*E->iV2;

     if (    PointInPolygon(V1, PolygonVertices, NumVertices)
          && PointInPolygon(V2, PolygonVertices, NumVertices)  
        )
      { if (Count==MAXEDGES) return -1;
        NewEdgeIndices[Count++]=nei;
      };
   };

  return Count;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AutoSelectRefPoints(RWGPort *Port)
{ 
  RWGSurface *PS = Port->PSurface, *MS = Port->MSurface;
  RWGPanel *PPanel, *MPanel;
  int iPPanel, iMPanel;
  int PiQ, MiQ;
  double *V1P, *V2P, *V1M, *V2M, VMidP[3], VMidM[3];
  double ThisDistance, MinDistance=0.0;

  for(int nppe=0; nppe<Port->NumPEdges; nppe++)
   { 
     iPPanel = Port->PPanelIndices[nppe];
     PPanel  = PS->Panels[iPPanel];
     PiQ     = Port->PPaneliQs[nppe]; 
     V1P     = PS->Vertices + 3*PPanel->VI[ (PiQ+1)%3 ];
     V2P     = PS->Vertices + 3*PPanel->VI[ (PiQ+2)%3 ];
     VMidP[0]  = 0.5*(V1P[0] + V2P[0]);
     VMidP[1]  = 0.5*(V1P[1] + V2P[1]);
     VMidP[2]  = 0.5*(V1P[2] + V2P[2]);

     for(int nmpe=0; nmpe<Port->NumMEdges; nmpe++)
      { 
        iMPanel = Port->MPanelIndices[nmpe];
        MPanel  = MS->Panels[iMPanel];
        MiQ     = Port->MPaneliQs[nmpe]; 
        V1M     = MS->Vertices + 3*MPanel->VI[ (MiQ+1)%3 ];
        V2M     = MS->Vertices + 3*MPanel->VI[ (MiQ+2)%3 ];
        VMidM[0]  = 0.5*(V1M[0] + V2M[0]);
        VMidM[1]  = 0.5*(V1M[1] + V2M[1]);
        VMidM[2]  = 0.5*(V1M[2] + V2M[2]);
    
        if (nppe==0 && nmpe==0)
         { MinDistance = VecDistance(VMidP, VMidM);
           memcpy(Port->PRefPoint, VMidP, 3*sizeof(double));
           memcpy(Port->MRefPoint, VMidM, 3*sizeof(double));
         }
        else
         { ThisDistance = VecDistance(VMidP, VMidM);
           if (ThisDistance < MinDistance) 
            { MinDistance = ThisDistance;
              memcpy(Port->PRefPoint, VMidP, 3*sizeof(double));
              memcpy(Port->MRefPoint, VMidM, 3*sizeof(double));
            };
         };
      };

    };

} 

/***************************************************************/
/* port file syntax example                                    */
/*  PORT                                                       */
/*   POBJECT   FirstObjectLabel                                */
/*   MOBJECT   SecondObjectLabel                               */
/*   PSURFACE  FirstSurfaceLabel                               */
/*   MSURFACE  SecondSurfaceLabel                              */
/*   PEDGES pe1 pe2 ... peN                                    */
/*   MEDGES me1 me2 ... meN                                    */
/*   PPOLYGON A1 A2 A3 B1 B2 B3 ... Z1 Z2 Z3                   */
/*   MPOLYGON A1 A2 A3 B1 B2 B3 ... Z1 Z2 Z3                   */
/*   PREFPOINT x1 x2 x3 <optional>                             */
/*   MREFPOINT x1 x2 x3 <optional>                             */
/*  ENDPORT                                                    */
/***************************************************************/
RWGPort **ParsePortFile(RWGGeometry *G, 
                        const char *PortFileName, 
                        int *pNumPorts)
{
  /***************************************************************/
  /* try to open the file ****************************************/
  /***************************************************************/
  FILE *f=fopen(PortFileName,"r");
  if (f==0)
   ErrExit("could not open file %s",PortFileName);
  
  /***************************************************************/
  /* read through the file and parse lines one-at-a-time *********/
  /***************************************************************/
  RWGPort **PortArray=0;
  int NumPorts=0;
  int NumPEdges=0, NumMEdges=0;
  int PEIndices[MAXEDGES], MEIndices[MAXEDGES];
  double PRefPoint[3], MRefPoint[3];
  RWGSurface *PSurface=0, *MSurface=0;

  bool IsPositive;

  int NumNewEdges;
  int NewEdgeIndices[MAXEDGES];

  int NumPolygonVertices;
  double PolygonVertices[3*MAXPOLYGONVERTICES];

  RWGSurface *WhichSurface;

  char buffer[1000];
  char *Tokens[MAXEDGES+2];
  int nt, NumTokens;
  int LineNum=0;
  int InPortSection=0;
  int PRefPointSpecified=0;
  int MRefPointSpecified=0;
  while( fgets(buffer, 1000, f) )
   { 
     /*--------------------------------------------------------------*/
     /*- skip blank lines and tokens --------------------------------*/
     /*--------------------------------------------------------------*/
     LineNum++;
     NumTokens=Tokenize(buffer, Tokens, MAXEDGES+2); 
     if ( NumTokens==0 || Tokens[0][0]=='#')
      continue;

     /*--------------------------------------------------------------*/
     /*- parse lines separately depending on whether or not we are  -*/
     /*- in a 'PORT...ENDPORT' section                              -*/
     /*--------------------------------------------------------------*/
     if ( !InPortSection ) 
      { if ( !StrCaseCmp(Tokens[0],"PORT") )
         { InPortSection=1;
           NumPEdges=NumMEdges=0;
           PRefPointSpecified=MRefPointSpecified=0;
           PSurface=MSurface=G->Surfaces[0];
         }
        else
         ErrExit("%s:%i: syntax error",PortFileName,LineNum);
      }
     else // !InPortSection
      {
        if ( !StrCaseCmp(Tokens[0],"ENDPORT") )
         { 
           if (NumPEdges==0)
            ErrExit("%s:%i: no edges specified or detected for positive port",PortFileName,LineNum);

           if (NumMEdges==0)
            ErrExit("%s:%i: no edges specified or detected for negative port",PortFileName,LineNum);

           
           PortArray = (RWGPort **)realloc( PortArray, (NumPorts+1)*sizeof(PortArray[0]) );
           PortArray[NumPorts] = CreatePort(PSurface, NumPEdges, PEIndices,
                                            MSurface, NumMEdges, MEIndices);

           AutoSelectRefPoints(PortArray[NumPorts]);


           if (PRefPointSpecified) 
            memcpy(PortArray[NumPorts]->PRefPoint, PRefPoint, 3*sizeof(double));
           if (MRefPointSpecified)
            memcpy(PortArray[NumPorts]->MRefPoint, MRefPoint, 3*sizeof(double));

           NumPorts++;
           InPortSection=0;
         }
        else if (    !StrCaseCmp(Tokens[0],"POBJECT") || !StrCaseCmp(Tokens[0],"PSURFACE") 
                  || !StrCaseCmp(Tokens[0],"MOBJECT") || !StrCaseCmp(Tokens[0],"MSURFACE") 
                )
         {
           /*--------------------------------------------------------------*/
           /*- POBJECT/PSURFACE or MOBJECT/MSURFACE specifies the label of */
           /*- the RWGSurface on which the port edges lie                  */
           /*--------------------------------------------------------------*/
           if( NumTokens!=2 )
            ErrExit("%s:%i: syntax error",PortFileName,LineNum);

           WhichSurface=G->GetSurfaceByLabel(Tokens[1]);
           if( !WhichSurface )
            ErrExit("%s:%i: could not find object/surface %s in geometry %s",PortFileName,LineNum,Tokens[1],G->GeoFileName);

           IsPositive = (Tokens[0][0]=='P') || (Tokens[0][0]=='p');
           if (IsPositive)
            PSurface=WhichSurface;
           else
            MSurface=WhichSurface;
         }
        else if ( !StrCaseCmp(Tokens[0],"PEDGES") || !StrCaseCmp(Tokens[0],"MEDGES") )
         { 
           /*--------------------------------------------------------------*/
           /*- PEDGES / MEDGES specifies a list of exterior edge indices   */
           /*- for the positive / negative port                            */
           /*--------------------------------------------------------------*/
           NumNewEdges = NumTokens-1;
           if ( NumNewEdges > MAXEDGES )
            ErrExit("%s:%i: too many edges",PortFileName,LineNum);

           for(nt=1; nt<NumTokens; nt++)
            if (1!=sscanf(Tokens[nt],"%i",NewEdgeIndices+(nt-1)) )
             ErrExit("%s:%i: syntax error %s",PortFileName,LineNum,Tokens[nt]);

           IsPositive = (Tokens[0][0]=='P') || (Tokens[0][0]=='p');
           if (IsPositive)
            { if ( (NumPEdges+NumNewEdges) > MAXEDGES )
               ErrExit("%s:%i: too many edges",PortFileName,LineNum);
              memcpy( &(PEIndices[NumPEdges]), NewEdgeIndices, NumNewEdges*sizeof(int));
              NumPEdges+=NumNewEdges;
            }
           else
            { if ( (NumMEdges+NumNewEdges) > MAXEDGES )
               ErrExit("%s:%i: too many edges",PortFileName,LineNum);
              memcpy( &(MEIndices[NumMEdges]), NewEdgeIndices, NumNewEdges*sizeof(int));
              NumMEdges+=NumNewEdges;
            };
         }
        else if ( !StrCaseCmp(Tokens[0],"PPOLYGON") || !StrCaseCmp(Tokens[0],"MPOLYGON") )
         { 
           /*--------------------------------------------------------------*/
           /*- PPOLYGON/MPOLYGON specifies a list of polygon vertices;     */
           /*- exterior RWG edges that lie within the specified polygon    */
           /*- are added to the positive/negative port                     */
           /*--------------------------------------------------------------*/
           if ( ((NumTokens-1) % 3) != 0 )
            ErrExit("%s:%i: number of arguments to %s must be a multiple of 3",PortFileName,LineNum,Tokens[0]);
           NumPolygonVertices = (NumTokens-1)/3;
           if ( NumPolygonVertices < 2 )
            ErrExit("%s:%i: %s requires at least 2 vertices",PortFileName,LineNum,Tokens[0]);
           if ( NumPolygonVertices > MAXPOLYGONVERTICES )
            ErrExit("%s:%i: too many vertices for %s",PortFileName,LineNum,Tokens[0]);

           for(nt=1; nt<NumTokens; nt++)
            sscanf(Tokens[nt],"%le",PolygonVertices+(nt-1));

           IsPositive = (Tokens[0][0]=='P') || (Tokens[0][0]=='p');
           NumNewEdges=FindEdgesInPolygon(IsPositive ? PSurface : MSurface,
                                          PolygonVertices, NumPolygonVertices, 
                                          NewEdgeIndices);
           if (NumNewEdges==-1)
            ErrExit("%s:%i: too many edges",PortFileName,LineNum);

           if (IsPositive)
            { if ( (NumPEdges + NumNewEdges) > MAXEDGES )
               ErrExit("%s:%i: too many edges",PortFileName,LineNum);
              memcpy( &(PEIndices[NumPEdges]), NewEdgeIndices, NumNewEdges*sizeof(int));
              NumPEdges += NumNewEdges;
            }
           else 
            { if ( (NumMEdges + NumNewEdges) > MAXEDGES )
               ErrExit("%s:%i: too many edges",PortFileName,LineNum);
              memcpy( &(MEIndices[NumMEdges]), NewEdgeIndices, NumNewEdges*sizeof(int));
              NumMEdges += NumNewEdges;
            };

           Log("Port %i, %s edge: found %i edges in polygon: ", 
                NumPorts+1, IsPositive ? "positive" : "negative", NumNewEdges);
           for(int nne=0; nne<NumNewEdges; nne++)
            LogC("%i ",NewEdgeIndices[nne]);
         }
        else if ( !StrCaseCmp(Tokens[0],"PREFPOINT") || !StrCaseCmp(Tokens[0],"MREFPOINT") )
         {
           /*--------------------------------------------------------------*/
           /*- PREFPOINT/MREFPOINT specifies the cartesian coordinates of  */
           /*- the positive/negative port reference point                  */
           /*--------------------------------------------------------------*/
           IsPositive = (Tokens[0][0]=='P') || (Tokens[0][0]=='p');
           double *Target;
           if (IsPositive)
            { PRefPointSpecified=1;
              Target = PRefPoint;
            }
           else
            { MRefPointSpecified=1;
              Target = MRefPoint;
            };
           if(    NumTokens!=4
               || 1!=sscanf(Tokens[1],"%le",Target+0)
               || 1!=sscanf(Tokens[2],"%le",Target+1)
               || 1!=sscanf(Tokens[3],"%le",Target+2)
             )
            ErrExit("%s:%i: syntax error",PortFileName,LineNum);
         }
       else
        ErrExit("%s:%i: unknown keyword %s",PortFileName,LineNum,Tokens[0]);

      }; // if (InPortSection ...

   };
  fclose(f);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  *pNumPorts=NumPorts;
  return PortArray;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddPortContributionsToRHS(RWGGeometry *G, 
                               RWGPort **Ports, int NumPorts, cdouble *PortCurrents, 
                               cdouble Omega, HVector *KN)
{
  int ns, ne, npe, BFIndex, nPort;
  RWGPort *Port;
  RWGEdge *E;
  cdouble IK=II*Omega;
 
  RWGSurface *SourceSurface;
  int SourcePanelIndex;
  int SourcePaneliQ;
  double SourceLength;
  cdouble PortCurrent, Weight;

  RWGSurface *DestSurface;
  int DestPPanelIndex, DestMPanelIndex;
  int DestPPaneliQ, DestMPaneliQ;
  double DestLength;

  cdouble EFieldIntegral;

  GetPPIArgStruct GPPIArgsBuffer, *GPPIArgs=&GPPIArgsBuffer;
  InitGetPPIArgs(GPPIArgs);
  GPPIArgs->k=Omega; // this assumes the exterior medium is vacuum

  /***************************************************************/
  /* insert a check to make sure all objects are PEC             */
  /***************************************************************/

  /***************************************************************/
  /* fill in (actually augment) entries of the RHS vector one-by-*/
  /* one                                                         */
  /***************************************************************/
  int NumThreads=GetNumThreads();
//#ifdef USE_OPENMP
//#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
//#endif
  for(BFIndex=0, ns=0; ns<G->NumSurfaces; ns++)
   { 
     DestSurface=G->Surfaces[ns];
     for(ne=0; ne<DestSurface->NumEdges; ne++, BFIndex++)
      { 
        E=DestSurface->Edges[ne];

        DestPPanelIndex  = E->iPPanel;
        DestPPaneliQ     = E->PIndex;

        DestMPanelIndex  = E->iMPanel;
        DestMPaneliQ     = E->MIndex;

        DestLength       = E->Length;

        /***************************************************************/
        /* get contribution of all ports to this element of the vector.*/
        /* note: at the conclusion of this loop, we have that          */
        /* EFieldIntegral is the inner product of basis function #ne   */
        /* on object O with the total E-field produced by all driven   */
        /* ports. the corresponding contribution to the RHS vector is  */
        /* then simply -EFieldIntegral.                                */
        /***************************************************************/
        for(EFieldIntegral=0.0, nPort=0; nPort<NumPorts; nPort++)
         { 
           PortCurrent=PortCurrents[nPort];
           if (PortCurrent==0.0) continue;

           Port=Ports[nPort];

           /***************************************************************/
           /* get contributions of panels on the positive side of the port*/
           /***************************************************************/
           SourceSurface=Port->PSurface;
           Weight=PortCurrent/(Port->PPerimeter);
           for(npe=0; npe<Port->NumPEdges; npe++) // npe = 'num port edge'
            { 
              SourcePanelIndex  = Port->PPanelIndices[npe];
              SourcePaneliQ     = Port->PPaneliQs[npe];
              SourceLength      = Port->PLengths[npe];

              GPPIArgs->Sa  = SourceSurface;
              GPPIArgs->npa = SourcePanelIndex;
              GPPIArgs->iQa = SourcePaneliQ;

              GPPIArgs->Sb  = DestSurface;
              GPPIArgs->npb = DestPPanelIndex;
              GPPIArgs->iQb = DestPPaneliQ;
              GetPanelPanelInteractions(GPPIArgs);
              EFieldIntegral -= Weight*SourceLength*DestLength*IK*GPPIArgs->H[0];

              GPPIArgs->Sb  = DestSurface;
              GPPIArgs->npb = DestMPanelIndex;
              GPPIArgs->iQb = DestMPaneliQ;
              GetPanelPanelInteractions(GPPIArgs);
              EFieldIntegral += Weight*SourceLength*DestLength*IK*GPPIArgs->H[0];

            }; // for ne=0; ne<Port->NumPEdges; ne++)

           /***************************************************************/
           /* get contributions of panels on the negative side of the port*/
           /***************************************************************/
           SourceSurface=Port->MSurface;
           Weight=PortCurrent/(Port->MPerimeter);
           for(npe=0; npe<Port->NumMEdges; npe++)
            { 
              SourcePanelIndex  = Port->MPanelIndices[npe];
              SourcePaneliQ     = Port->MPaneliQs[npe];
              SourceLength      = Port->MLengths[npe];

              GPPIArgs->Sa  = SourceSurface;
              GPPIArgs->npa = SourcePanelIndex;
              GPPIArgs->iQa = SourcePaneliQ;

              GPPIArgs->Sb  = DestSurface;
              GPPIArgs->npb = DestPPanelIndex;
              GPPIArgs->iQb = DestPPaneliQ;
              GetPanelPanelInteractions(GPPIArgs);
              EFieldIntegral += Weight*SourceLength*DestLength*IK*GPPIArgs->H[0];

              GPPIArgs->Sb  = DestSurface;
              GPPIArgs->npb = DestMPanelIndex;
              GPPIArgs->iQb = DestMPaneliQ;
              GetPanelPanelInteractions(GPPIArgs);
              EFieldIntegral -= Weight*SourceLength*DestLength*IK*GPPIArgs->H[0];

            }; // for ne=0; ne<Port->NumMEdges; ne++)
         }; // for(nPort=0; nPort<NumPorts; nPort++)

        KN->AddEntry(BFIndex, -1.0*EFieldIntegral );

      }; // for(ne=0; ne<S->NumEdges; ne++)
   }; // for(ns=0=; ns<G->NumSurfaces; ns++)

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddPortContributionsToPSD(RWGGeometry *G, 
                               RWGPort **Ports, int NumPorts, cdouble *PortCurrents,
                               cdouble Omega, HMatrix *PSD)
{ 
  cdouble IW = II*Omega;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (PSD==0 || PSD->NR!=G->TotalPanels || PSD->NC!=12 || PSD->RealComplex!=LHM_COMPLEX)
   { if (PSD) 
      Warn("invalid PSD matrix passed to AddPortContributionsToPSD (skipping)");
     return;
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  RWGPort *Port;
  RWGSurface *S;
  RWGPanel *Panel;
  int PanelIndex, PaneliQ, PIOffset;
  double Length;
  double *Q, XmQ[3];
  cdouble Weight, PreFac;
  for(int nPort=0; nPort<NumPorts; nPort++)
   { 
     if (PortCurrents[nPort]==0.0) continue;
     Port=Ports[nPort];

     /*--------------------------------------------------------------*/
     /*- contributions of positive port edges -----------------------*/
     /*--------------------------------------------------------------*/
     S=Port->PSurface;
     PIOffset=G->PanelIndexOffset[S->Index]; 
     Weight=PortCurrents[nPort]/(Port->PPerimeter);
     for(int npe=0; npe<Port->NumPEdges; npe++) // npe = 'num port edge'
      { 
        PanelIndex = Port->PPanelIndices[npe];
        PaneliQ    = Port->PPaneliQs[npe];
        Length     = Port->PLengths[npe];

        Panel      = S->Panels[PanelIndex];
        Q          = S->Vertices + 3*PaneliQ;

        // prefactor is *negative* for positive port edges 
        PreFac     = -1.0*Length * Weight / (2.0*Panel->Area);
        VecSub(Panel->Centroid, Q, XmQ);

        PSD->AddEntry( PIOffset + PanelIndex,  4, 2.0*PreFac/ IW);  // rho 
        PSD->AddEntry( PIOffset + PanelIndex,  5, PreFac*XmQ[0] );  // K_x
        PSD->AddEntry( PIOffset + PanelIndex,  6, PreFac*XmQ[1] );  // K_y
        PSD->AddEntry( PIOffset + PanelIndex,  7, PreFac*XmQ[2] );  // K_z
      }; // for(npe=0; npe<Port->NumPEdges; npe++) // npe = 'num port edge'

     /*--------------------------------------------------------------*/
     /*- contributions of negative port edges -----------------------*/
     /*--------------------------------------------------------------*/
     S=Port->MSurface;
     PIOffset=G->PanelIndexOffset[S->Index]; 
     Weight=PortCurrents[nPort]/(Port->MPerimeter);
     for(int npe=0; npe<Port->NumPEdges; npe++) // npe = 'num port edge'
      { 
        PanelIndex = Port->MPanelIndices[npe];
        PaneliQ    = Port->MPaneliQs[npe];
        Length     = Port->MLengths[npe];

        Panel      = S->Panels[PanelIndex];
        Q          = S->Vertices + 3*PaneliQ;

        // prefactor is *positive* for negative port edges 
        PreFac     = +1.0*Length * Weight / (2.0*Panel->Area);
        VecSub(Panel->Centroid, Q, XmQ);

        PSD->AddEntry( PIOffset + PanelIndex,  4, 2.0*PreFac/ IW);  // rho 
        PSD->AddEntry( PIOffset + PanelIndex,  5, PreFac*XmQ[0] );  // K_x
        PSD->AddEntry( PIOffset + PanelIndex,  6, PreFac*XmQ[1] );  // K_y
        PSD->AddEntry( PIOffset + PanelIndex,  7, PreFac*XmQ[2] );  // K_z
      }; // for(npe=0; npe<Port->NumMEdges; npe++) // npe = 'num port edge'

   }; // for(int nPort=0; nPort<NumPorts; nPort++)
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void PlotPortsInGNUPLOT(const char *GPFileName, RWGPort **Ports, int NumPorts)
{
  FILE *f=fopen(GPFileName,"w");

  RWGSurface *S;
  RWGPort *Port;
  int nPort, nPanel, PanelIndex, PaneliQ;
  double *V1, *V2;
  for(nPort=0; nPort<NumPorts; nPort++)
   { 
     Port=Ports[nPort];

     fprintf(f,"%e %e %e \n", Port->PRefPoint[0], Port->PRefPoint[1], Port->PRefPoint[2]);
     fprintf(f,"\n\n");
     fprintf(f,"%e %e %e \n\n\n", Port->MRefPoint[0], Port->MRefPoint[1], Port->MRefPoint[2]);
     fprintf(f,"\n\n");

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     S=Port->PSurface;
     for(nPanel=0; nPanel<Port->NumPEdges; nPanel++)
      { 
        PanelIndex = Port->PPanelIndices[nPanel]; 
        PaneliQ    = Port->PPaneliQs[nPanel]; 
        V1 = S->Vertices + 3*S->Panels[PanelIndex]->VI[(PaneliQ+1)%3];
        V2 = S->Vertices + 3*S->Panels[PanelIndex]->VI[(PaneliQ+2)%3];
        fprintf(f,"%e %e %e \n",V1[0],V1[1],V1[2]);
        fprintf(f,"%e %e %e \n",V2[0],V2[1],V2[2]);
        fprintf(f,"\n\n");
      };

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     S=Port->MSurface;
     for(nPanel=0; nPanel<Port->NumMEdges; nPanel++)
      { 
        PanelIndex = Port->MPanelIndices[nPanel]; 
        PaneliQ    = Port->MPaneliQs[nPanel]; 
        V1 = S->Vertices + 3*S->Panels[PanelIndex]->VI[(PaneliQ+1)%3];
        V2 = S->Vertices + 3*S->Panels[PanelIndex]->VI[(PaneliQ+2)%3];
        fprintf(f,"%e %e %e \n",V1[0],V1[1],V1[2]);
        fprintf(f,"%e %e %e \n",V2[0],V2[1],V2[2]);
        fprintf(f,"\n\n");
      };

   };

  fclose(f);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void PlotPortsInGMSH(RWGPort **Ports, int NumPorts, const char *format, ...)
{
  /***************************************************************/
  /* open the file ***********************************************/
  /***************************************************************/
  va_list ap;
  char FileName[1000];
  va_start(ap,format);
  vsnprintfEC(FileName,997,format,ap);
  va_end(ap);

  FILE *f=fopen(FileName,"w");
  if (!f)
   { Warn("could not open file %s (skipping port plot)",FileName);
     return;
   };
  fprintf(f,"View.LineWidth = 5;\n");
  fprintf(f,"View.LineType  = 1;\n");
  fprintf(f,"View.CustomMax = %i;\n",+(NumPorts+1));
  fprintf(f,"View.CustomMin = %i;\n",-(NumPorts+1));
  fprintf(f,"View.RangeType = 2;\n");
  fprintf(f,"View.ShowScale = 0;\n");

  /***************************************************************/
  /* loop over all ports on all surfaces *************************/
  /***************************************************************/
  for(int nPort=0; nPort<NumPorts; nPort++)
   for(int Polarity=0; Polarity<2; Polarity++)
    { 
      fprintf(f,"View \"Port %i %s terminal\" {\n",nPort+1, Polarity ? "positive" : "negative");

      RWGPort *Port     = Ports[nPort];
      double *RefPoint  = Polarity ? Port->PRefPoint     : Port->MRefPoint;
      RWGSurface *S     = Polarity ? Port->PSurface      : Port->MSurface;
      int NumEdges      = Polarity ? Port->NumPEdges     : Port->NumMEdges;
      int *PanelIndices = Polarity ? Port->PPanelIndices : Port->MPanelIndices;
      int *PaneliQs     = Polarity ? Port->PPaneliQs     : Port->MPaneliQs;
      int Value         = (Polarity ? 1 : -1 ) * (nPort+1);

      /*--------------------------------------------------------------*/
      /*- scalar points for ref points                               -*/
      /*--------------------------------------------------------------*/
      fprintf(f,"SP(%e,%e,%e) {%i};\n", RefPoint[0],RefPoint[1],RefPoint[2],Value);

      /*--------------------------------------------------------------*/ 
      /*- scalar lines for port edges                                 */
      /*--------------------------------------------------------------*/
      for(int nPanel=0; nPanel<NumEdges; nPanel++)
       { 
         int PanelIndex = PanelIndices[nPanel];
         int PaneliQ    = PaneliQs[nPanel]; 
         double *V1 = S->Vertices + 3*S->Panels[PanelIndex]->VI[(PaneliQ+1)%3];
         double *V2 = S->Vertices + 3*S->Panels[PanelIndex]->VI[(PaneliQ+2)%3];
         fprintf(f,"SL(%e,%e,%e,%e,%e,%e) {%i,%i};\n",
                    V1[0],V1[1],V1[2],V2[0],V2[1],V2[2],Value,Value);
         fprintf(f,"\n\n");
       };

      fprintf(f,"};\n");

   };

  fclose(f);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void DrawGMSHCircle(const char *PPFile, const char *Name, 
                    double *X0, double Theta, double Phi, double Radius)
{ 
  FILE *f=fopen(PPFile,"a");
  fprintf(f,"View \"%s\" {\n",Name);

  double CT, ST, CP, SP, M[3][3];
  CT=cos(Theta);
  ST=sin(Theta);
  CP=cos(Phi);
  SP=sin(Phi);
  M[0][0]=CT*CP;     M[0][1]=CT*SP;    M[0][2]=-ST;
  M[1][0]=-SP;       M[1][1]=CP;       M[1][2]=0.0;
  M[2][0]=ST*CP;     M[2][1]=ST*SP;    M[2][2]=CT;

  double Psi, X1P[3], X1[3], X2P[3], X2[3];
  int Mu, Nu;
  for(Psi=0; Psi<2.0*M_PI; Psi+=M_PI/50.0)
   { 
     X1P[0]=Radius*cos(Psi);
     X1P[1]=Radius*sin(Psi);
     X1P[2]=0.0;

     X2P[0]=Radius*cos(Psi+M_PI/50.0);
     X2P[1]=Radius*sin(Psi+M_PI/50.0);
     X2P[2]=0.0;

     memcpy(X1,X0,3*sizeof(double));
     memcpy(X2,X0,3*sizeof(double));
     for(Mu=0; Mu<3; Mu++)
      for(Nu=0; Nu<3; Nu++)
       { X1[Mu]+=M[Nu][Mu]*X1P[Nu];
         X2[Mu]+=M[Nu][Mu]*X2P[Nu];
       };
     fprintf(f,"SL(%e,%e,%e,%e,%e,%e) {%e,%e};\n",
                X1[0],X1[1],X1[2],X2[0],X2[1],X2[2], 0.0, 0.0);
    };

  fprintf(f,"};\n");
  fclose(f);
 
}
