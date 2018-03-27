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

#include <vector>
using namespace std;

#include "RWGPorts.h"

#ifdef USE_OPENMP
#  include <omp.h>
#endif

using namespace scuff;

#define II cdouble(0.0,1.0)
static const double Signs[2]={1.0,-1.0};
static const char *PolStr[2]={"POSITIVE", "NEGATIVE"};

void GetRzMinMax(RWGGeometry *G, HMatrix *XMatrix, double RhoMinMax[2], double zMinMax[2]);

/***************************************************************/
/***************************************************************/
/***************************************************************/
int GetMOIMatrixElement(RWGGeometry *G, LayeredSubstrate *Substrate,
                        int nsa, int nea, int nsb, int neb,
                        cdouble Omega, cdouble *ME,
                        int Order, bool Subtract=true, cdouble *Terms=0);

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
bool PointInPolygon(double *X, dVec VVector)
{ 
  int NumVertices = VVector.size() / 3;
  double *V = &(VVector[0]);
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
iVec FindEdgesInPolygon(RWGSurface *S, dVec PolygonVertices)
{
  iVec NewEdges;
  double *Vertices = S->Vertices;
  for(int ne=0; ne<S->NumExteriorEdges; ne++)
   { 
     RWGEdge *E = S->ExteriorEdges[ne];
     double *V1 = S->Vertices + 3*E->iV1;
     double *V2 = S->Vertices + 3*E->iV2;

bool V1In = PointInPolygon(V1, PolygonVertices);
bool V2In = PointInPolygon(V2, PolygonVertices);

     if ( PointInPolygon(V1, PolygonVertices) && PointInPolygon(V2, PolygonVertices) )
      NewEdges.push_back(ne);
   }
  return NewEdges;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AutoSelectRefPoints(RWGGeometry *G, RWGPort *Port)
{ 
  double MinDistance=1.0e100;
  for(int ieP=0; ieP<Port->PortEdges[_PLUS].size(); ieP++)
   for(int ieM=0; ieM<Port->PortEdges[_MINUS].size(); ieM++)
    { RWGPortEdge *PPE = Port->PortEdges[_PLUS][ieP];
      RWGPortEdge *MPE = Port->PortEdges[_MINUS][ieM];
      RWGEdge *EP      = G->Surfaces[PPE->ns]->GetEdgeByIndex(PPE->ne);
      RWGEdge *EM      = G->Surfaces[MPE->ns]->GetEdgeByIndex(MPE->ne);
      double Distance  = VecDistance(EP->Centroid, EM->Centroid);
      if (Distance<MinDistance)
       { MinDistance=Distance;
         memcpy(Port->RefPoint[_PLUS],EP->Centroid,3*sizeof(double));
         memcpy(Port->RefPoint[_MINUS],EM->Centroid,3*sizeof(double));
       }
    }
} 

/***************************************************************/
/* port file syntax example                                    */
/*  PORT                                                       */
/*   POBJECT   FirstObjectLabel                                */
/*   MOBJECT   SecondObjectLabel                               */
/*   PSURFACE  FirstSurfaceLabel                               */
/*   MSURFACE  SecondSurfaceLabel                              */
/*   PEDGES    pe1 pe2 ... peN                                 */
/*   MEDGES    me1 me2 ... meN                                 */
/*   PPOLYGON  A1 A2 A3 B1 B2 B3 ... Z1 Z2 Z3                  */
/*   MPOLYGON  A1 A2 A3 B1 B2 B3 ... Z1 Z2 Z3                  */
/*   PREFPOINT x1 x2 x3 <optional>                             */
/*   MREFPOINT x1 x2 x3 <optional>                             */
/*  ENDPORT                                                    */
/***************************************************************/
RWGPortList *ParsePortFile(RWGGeometry *G, const char *PortFileName)
{
  /***************************************************************/
  /* try to open the file ****************************************/
  /***************************************************************/
  FILE *f=fopen(PortFileName,"r");
  if (f==0)
   ErrExit("could not open file %s",PortFileName);
  
  RWGPortList *PortList = new RWGPortList;

  /***************************************************************/
  /* read through the file and parse lines one-at-a-time *********/
  /***************************************************************/
  RWGPort *CurrentPort=0;
  RWGSurface *CurrentSurface[2]={0,0};
  char Line[1000];
  int LineNum=0;
  while( fgets(Line, 1000, f) )
   { 
     /*--------------------------------------------------------------*/
     /*- skip blank lines and tokens --------------------------------*/
     /*--------------------------------------------------------------*/
     LineNum++;
     sVec Tokens=Tokenize(Line);
     int NumTokens=Tokens.size();
     if ( NumTokens==0 || Tokens[0][0]=='#')
      continue;

     /*--------------------------------------------------------------*/
     /*- parse lines                                                -*/
     /*--------------------------------------------------------------*/
     if ( !StrCaseCmp(Tokens[0],"PORT") )
      { if (CurrentPort!=0) ErrExit("%s:%i: syntax error (missing ENDPORT?)",PortFileName,LineNum);
        CurrentPort = new RWGPort;
        CurrentPort->RefPoint[0][0]=CurrentPort->RefPoint[1][0]=HUGE_VAL; // to indicate unspecified
      }
     else if ( !StrCaseCmp(Tokens[0],"ENDPORT") )
      { PortList->Ports.push_back(CurrentPort);
        CurrentPort=0;
        CurrentSurface[0]=CurrentSurface[1]=0;
      }
     else if ( !StrCaseCmp(Tokens[0]+1,"OBJECT") || !StrCaseCmp(Tokens[0]+1,"SURFACE") )
      { int Pol = ( toupper(Tokens[0][0])=='P' ? 0 : toupper(Tokens[0][0])=='M' ? 1 : 2);
        if (Pol==2) ErrExit("%s:%i: unknown keyword %s",PortFileName,LineNum,Tokens[0]); 
        if (!CurrentPort) ErrExit("%s:%i: %s outside of PORT...ENDPORT",PortFileName,LineNum,Tokens[0]);
        if ( NumTokens!=2 ) ErrExit("%s:%i: syntax error",PortFileName,LineNum);
        CurrentSurface[Pol]=G->GetSurfaceByLabel(Tokens[1]);
        if( !CurrentSurface[Pol] )
         ErrExit("%s:%i: no object/surface %s in geometry %s",PortFileName,LineNum,Tokens[1],G->GeoFileName);
      }
     else if ( !StrCaseCmp(Tokens[0]+1,"EDGES") )
      { int Pol = ( toupper(Tokens[0][0])=='P' ? 0 : toupper(Tokens[0][0])=='M' ? 1 : 2);
        if (Pol==2) ErrExit("%s:%i: unknown keyword %s",PortFileName,LineNum,Tokens[0]); 
        if (!CurrentPort) ErrExit("%s:%i: %s outside of PORT...ENDPORT",PortFileName,LineNum,Tokens[0]);
        RWGSurface *S=CurrentSurface[Pol];
        if (!S) 
         ErrExit("%s:%i: %s without preceding %cOBJECT/%cSURFACE",PortFileName,LineNum,Tokens[0],PolStr[Pol][0]);
        for(int nt=1; nt<NumTokens; nt++)
         { int ne;
           sscanf(Tokens[nt],"%i",&ne);
           if (ne<0 || ne>=S->NumExteriorEdges)
            ErrExit("%s:%i: surface %s has no exterior edge %i",PortFileName,LineNum,S->Label,ne);
           RWGPortEdge *PE = new RWGPortEdge(S->Index, -1-ne, PortList->Ports.size(), Pol);
           CurrentPort->PortEdges[Pol].push_back(PE);
           PortList->PortEdges.push_back(PE);
         }
      }
     else if ( !StrCaseCmp(Tokens[0]+1,"POLYGON") )
      { int Pol = ( toupper(Tokens[0][0])=='P' ? 0 : toupper(Tokens[0][0])=='M' ? 1 : 2);
        if (Pol==2) ErrExit("%s:%i: unknown keyword %s",PortFileName,LineNum,Tokens[0]);
        if (!CurrentPort) ErrExit("%s:%i: %s outside of PORT...ENDPORT",PortFileName,LineNum,Tokens[0]);
        if ( NumTokens%3 != 1 )
         ErrExit("%s:%i: number of arguments to %s must be a multiple of 3",PortFileName,LineNum,Tokens[0]);
        if ( NumTokens < 7 )
         ErrExit("%s:%i: %s requires at least 2 vertices",PortFileName,LineNum,Tokens[0]);
        int NumPolygonVertices = (NumTokens-1)/3;
        dVec PolygonVertices(NumTokens-1);
        for(int nt=1; nt<NumTokens; nt++) sscanf(Tokens[nt],"%le",&(PolygonVertices[nt-1]));
        for(int ns=0; ns<G->NumSurfaces; ns++)
         { iVec neList = FindEdgesInPolygon(G->Surfaces[ns], PolygonVertices);
           for(int nne=0; nne<neList.size(); nne++)
            { RWGPortEdge *PE = new RWGPortEdge(ns, -1-neList[nne], PortList->Ports.size(), Pol);
              CurrentPort->PortEdges[Pol].push_back(PE);
              PortList->PortEdges.push_back(PE);
            }
           Log("Port %i, %s edge: found %i exterior edges of %s in polygon: ",
                PortList->Ports.size()+1,PolStr[Pol], neList.size(), G->Surfaces[ns]->Label);
           for(int nne=0; nne<neList.size(); nne++)
            LogC("%i ",neList[nne]);
         }
      }
     else if ( !StrCaseCmp(Tokens[0]+1,"REFPOINT") )
      { int Pol = ( toupper(Tokens[0][0])=='P' ? 0 : toupper(Tokens[0][0])=='M' ? 1 : 2);
        if (Pol==2) ErrExit("%s:%i: unknown keyword %s",PortFileName,LineNum,Tokens[0]);
        if (!CurrentPort) ErrExit("%s:%i: %s outside of PORT...ENDPORT",PortFileName,LineNum,Tokens[0]);
        if (NumTokens!=4)
         ErrExit("%s:%i: syntax error",PortFileName,LineNum);
        for(int n=0; n<3; n++)
         sscanf(Tokens[1+n],"%le",CurrentPort->RefPoint[Pol] + n);
      }
     else
      ErrExit("%s:%i: unknown keyword %s",PortFileName,LineNum,Tokens[0]);

   }; // while( fgets(buffer, 1000, f) )
  fclose(f);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int np=0; np<PortList->Ports.size(); np++)
   { RWGPort *Port = PortList->Ports[np];
     if( Port->RefPoint[0][0]==HUGE_VAL || Port->RefPoint[1][0]==HUGE_VAL )
      AutoSelectRefPoints(G, Port);
   }

  /***************************************************************/
  /* compute perimeters ******************************************/
  /***************************************************************/
  for(int np=0; np<PortList->Ports.size(); np++)
   { RWGPort *Port = PortList->Ports[np];
     for(int Pol=0; Pol<2; Pol++)
      { Port->Perimeter[Pol]=0.0;
        for(int ie=0; ie<Port->PortEdges[Pol].size(); ie++)
         { RWGPortEdge *PE = Port->PortEdges[Pol][ie];
           Port->Perimeter[Pol] += G->Surfaces[PE->ns]->ExteriorEdges[-(1+PE->ne)]->Length;
         }
      }
   }

  /***************************************************************/
  /* write summary of port list to log file **********************/
  /***************************************************************/
  for(int np=0; np<PortList->Ports.size(); np++)
   for(int Pol=0; Pol<2; Pol++)
    { RWGPort *Port = PortList->Ports[np];
      int NPE = Port->PortEdges[Pol].size();
      double Perimeter=Port->Perimeter[Pol];
      double *X = Port->RefPoint[Pol];
      Log("Port %2i (%s): perimeter %e, X0={%g,%g,%g}, %i edges=[",np,PolStr[Pol],Perimeter,X[0],X[1],X[2],NPE);
      for(int npe=0; npe<NPE; npe++)
       { RWGPortEdge *PE = Port->PortEdges[Pol][npe];
         LogC("%s(%i)%c",G->Surfaces[PE->ns]->Label,-(1+PE->ne),(npe==NPE-1) ? ']' : ' ');
       }
    }

  return PortList;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddPortContributionsToRHS(RWGGeometry *G, RWGPortList *PortList,
                               cdouble *PortCurrents, cdouble Omega, HVector *KN)
{
  int NumPortEdges   = PortList->PortEdges.size();
  int NumEdgeBFPairs = NumPortEdges * G->TotalEdges;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double PortRMinMax[6]={1.0e89,1.0e89,1.0e89,-1.0e89,-1.0e89,-1.0e89};
  double *PortRMin = PortRMinMax+0, *PortRMax=PortRMinMax+3;
  for(int nPE=0; nPE<PortList->PortEdges.size(); nPE++)
   { RWGPortEdge *PE = PortList->PortEdges[nPE];
     RWGSurface *S   = G->Surfaces[PE->ns];
     RWGEdge *E      = S->GetEdgeByIndex(PE->ne);
     double *V[3];
     V[0] = S->Vertices + 3*E->iQP;
     V[1] = S->Vertices + 3*E->iV1;
     V[2] = S->Vertices + 3*E->iV2;
     for(int nv=0; nv<3; nv++)
      for(int Mu=0; Mu<3; Mu++)
       { PortRMin[Mu] = fmin(PortRMin[Mu], V[nv][Mu]);
         PortRMax[Mu] = fmax(PortRMin[Mu], V[nv][Mu]);
       }
   }
  HMatrix XMatrix(3,2,PortRMinMax);
  double RhoMinMax[2], zMinMax[2];
  GetRzMinMax(G, &XMatrix, RhoMinMax, zMinMax);
  Log("Initializing ScalarGF interpolator for Rho range (%e,%e)",RhoMinMax[0],RhoMinMax[1]);
  bool PPIsOnly=true;
  bool Subtract=true;
  bool RetainSingularTerms=false;
  G->Substrate->InitScalarGFInterpolator(Omega, RhoMinMax[0], RhoMinMax[1], 0.0, 0.0, 
                                         PPIsOnly, Subtract, RetainSingularTerms);

  int Order=9;
  char *s=getenv("SCUFF_PORTRHS_ORDER");
  if (s && 1==sscanf(s,"%i",&Order))
   Log("Setting panel-cubature order=%i for RF port RHS entries",Order);
  
#ifndef USE_OPENMP
  Log("Adding port contributions to RHS"); 
  Tic();
#else
  int NumThreads=GetNumThreads();
  Log("Adding port contributions to RHS (%i threads)",NumThreads);
  Tic();
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int nEBFP=0; nEBFP<NumEdgeBFPairs; nEBFP++)
   { 
     LogPercent(nEBFP, NumEdgeBFPairs);

     int nPortEdge       = nEBFP/G->TotalEdges;
     int nBFEdge         = nEBFP%G->TotalEdges;

     RWGPortEdge *PE     = PortList->PortEdges[nPortEdge];
     int nPort           = PE->nPort;
     cdouble PortCurrent = PortCurrents[nPort];
     if (PortCurrent==0.0) continue;
     int nsPort          = PE->ns;
     int nePort          = PE->ne;
     int Pol             = PE->Pol;
     double Sign         = Signs[Pol];
     double Perimeter    = PortList->Ports[nPort]->Perimeter[Pol];

     int nsBF, neBF;
     if (!G->ResolveEdge(nBFEdge, &nsBF, &neBF, 0))
      ErrExit("%s:%i: internal error(%i,%i,%i)\n",nBFEdge,nsBF,neBF);
     
     cdouble ME;
#if 1
     GetMOIMatrixElement(G, G->Substrate, nsPort, nePort, nsBF, neBF, Omega, &ME, Order);
#endif
#if 0
     GetGCMEArgStruct GCMEArgs, *Args=&GCMEArgs;
     InitGetGCMEArgs(Args);
     Args->nsa = nsPort;
     Args->nsb = nsBF;
     Args->NumRegions = 1;
     Args->k[0] = Omega;
     Args->NeedGC = true;
     Args->FIBBICache = (nsPort==nsBF) ? G->FIBBICaches[nsPort] : 0;
     cdouble GabArray[2][NUMGCMES];
     cdouble ikCabArray[2][NUMGCMES];
     GetGCMatrixElements(G, Args, nePort, neBF, GabArray, ikCabArray);
     ME=GabArray[0][0];
#endif
#if 0
  GetEEIArgStruct MyGetEEIArgs, *GetEEIArgs=&MyGetEEIArgs;
  InitGetEEIArgs(GetEEIArgs);
  GetEEIArgs->Sa  = G->Surfaces[nsPort];
  GetEEIArgs->nea = nePort;
  GetEEIArgs->Sb  = G->Surfaces[nsBF];
  GetEEIArgs->neb = neBF;
  GetEEIArgs->k   = Omega;
  GetEdgeEdgeInteractions(GetEEIArgs);
  ME=GetEEIArgs->GC[0];
#endif
     cdouble Weight = -1.0*II*Omega*Sign*PortCurrent/Perimeter;
     KN->AddEntry(nBFEdge,Weight*ME);

   }; // for(int nEBFP=0; nEBFP<NumEdgeBFPairs; nEBFP++)
  double Time=Toc();
  Log("...RHS port contributions done in %e s\n",Time);
  printf("Time = %e s\n",Time);

}

/***************************************************************/
/* contribution of port currents to panel source densities     */
/***************************************************************/
void AddPortContributionsToPSD(RWGGeometry *G, RWGPortList *PortList,
                               cdouble *PortCurrents, cdouble Omega, HMatrix *PSD)
{ 
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
  for(int nPE=0; nPE=PortList->PortEdges.size(); nPE++)
   { 
     RWGPortEdge *PE = PortList->PortEdges[nPE];
     cdouble Weight  = PortCurrents[PE->nPort] / (PortList->Ports[PE->nPort]->Perimeter[PE->Pol]);
     if (Weight==0.0) continue;
     RWGSurface *S   = G->Surfaces[PE->ns];
     RWGEdge *E      = S->GetEdgeByIndex(PE->ne);
     double Length   = E->Length; 
     double *QP      = S->Vertices + 3*(E->iQP);
     int PanelIndex  = G->PanelIndexOffset[S->Index] + E->iPPanel;
     RWGPanel *Pan   = S->Panels[E->iPPanel];

     // prefactor is *negative* for positive port edges 
     cdouble PreFac  = Signs[PE->Pol]*Weight*Length/(2.0*Panel->Area);
     double XmQ[3];
     VecSub(Panel->Centroid, QP, XmQ);
     PSD->AddEntry( PanelIndex,  4, 2.0*PreFac/(II*Omega));  // rho
     PSD->AddEntry( PanelIndex,  5, PreFac*XmQ[0] );  // K_x
     PSD->AddEntry( PanelIndex,  6, PreFac*XmQ[1] );  // K_y
     PSD->AddEntry( PanelIndex,  7, PreFac*XmQ[2] );  // K_z

   }; // for(int nPE=0...)
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
#if 0
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
#endif

/***************************************************************/
/***************************************************************/
/***************************************************************/
void PlotPortsInGMSH(RWGGeometry *G, RWGPortList *PortList, const char *format, ...)
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

  int NumPorts = (int)PortList->Ports.size();
  fprintf(f,"View.LineWidth = 5;\n");
  fprintf(f,"View.LineType  = 1;\n");
  fprintf(f,"View.CustomMax = %i;\n",+(NumPorts+1));
  fprintf(f,"View.CustomMin = %i;\n",-(NumPorts+1));
  fprintf(f,"View.RangeType = 2;\n");
  fprintf(f,"View.ShowScale = 0;\n");

  /***************************************************************/
  /* loop over all ports on all surfaces *************************/
  /***************************************************************/
  for(int nPort=0; nPort<PortList->Ports.size(); nPort++)
   for(int Pol=0; Pol<2; Pol++)
    { 
      fprintf(f,"View \"Port %i %s terminal\" {\n",nPort+1, PolStr[Pol]);

      RWGPort *Port     = PortList->Ports[nPort];
      double *RefPoint  = Port->RefPoint[Pol];
      int Value         = (Pol ? 1 : -1 ) * (nPort+1);

      /*--------------------------------------------------------------*/
      /*- scalar points for ref points                               -*/
      /*--------------------------------------------------------------*/
      fprintf(f,"SP(%e,%e,%e) {%i};\n", RefPoint[0],RefPoint[1],RefPoint[2],Value);

      /*--------------------------------------------------------------*/ 
      /*- scalar lines for port edges                                 */
      /*--------------------------------------------------------------*/
      for(unsigned nPE=0; nPE<Port->PortEdges[Pol].size(); nPE++)
       { RWGPortEdge *PE = Port->PortEdges[Pol][nPE];
         RWGSurface *S   = G->Surfaces[PE->ns];
         RWGEdge *E      = S->GetEdgeByIndex(PE->ne);
         double *V1      = S->Vertices + 3*E->iV1;
         double *V2      = S->Vertices + 3*E->iV2;
         fprintf(f,"SL(%e,%e,%e,%e,%e,%e) {%i,%i};\n",
                    V1[0],V1[1],V1[2],V2[0],V2[1],V2[2],Value,Value);
       }
      fprintf(f,"};\n");
   }
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
