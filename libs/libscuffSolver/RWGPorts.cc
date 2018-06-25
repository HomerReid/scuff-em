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
#include <unistd.h>

#include <libhrutil.h>
#include <libhmat.h>
#include <libscuff.h>
#include <libscuffInternals.h>
#include <libSGJC.h>

#include <vector>
#include <map>
using namespace std;

#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif
#ifdef HAVE_LIBGDSII
  #include <libGDSII.h>
  using namespace libGDSII;
#endif

#include "scuffSolver.h"

using namespace scuff;

#define II cdouble(0.0,1.0)
static const double PolSigns[2]={1.0,-1.0};
static const char *PolNames[2]={"POSITIVE", "NEGATIVE"};
 
namespace scuff {

void GetRzMinMax(RWGGeometry *G, HMatrix *XMatrix, double RhoMinMax[2], double zMinMax[2]);

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
  for(int ne=0; ne<S->NumExteriorEdges; ne++)
   { 
     RWGEdge *E = S->HalfRWGEdges[ne];
     double *V1 = S->Vertices + 3*E->iV1;
     double *V2 = S->Vertices + 3*E->iV2;

     if ( PointInPolygon(V1, PolygonVertices) && PointInPolygon(V2, PolygonVertices) )
      NewEdges.push_back(ne);
   }
  return NewEdges;
}

/***************************************************************/
// general-purpose error checking for keywords in port file    */
/***************************************************************/
int CheckKeywordSyntax(const char *PortFileName, int LineNum, RWGPort *CurrentPort,
                       const char *Keyword, int NumTokens=0, int ReqNumTokens=0)
{
  int Pol = ( toupper(Keyword[0])=='P' ? _PLUS : toupper(Keyword[0])=='M' ? _MINUS : -1);
  if (Pol==-1) ErrExit("%s:%i: unknown keyword %s",PortFileName,LineNum,Keyword);
  if (!CurrentPort) ErrExit("%s:%i: %s outside of PORT...ENDPORT",PortFileName,LineNum,Keyword);
  if ( NumTokens!=ReqNumTokens ) ErrExit("%s:%i: syntax error",PortFileName,LineNum);
  return Pol;
}

void UpdateRMinMax(double RMin[3], double RMax[3], double *V)
{ for(int i=0; i<3; i++)
   { RMin[i] = fmin(RMin[i],V[i]);
     RMax[i] = fmax(RMax[i],V[i]);
   }
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
RWGPortEdgeList AddPointPort(RWGGeometry *G, double *X)
{
  RWGPortEdgeList PortEdges;

  // get distance to closest vertex on any surface
  double MinDistance=HUGE_VAL;
  for(int ns=0; ns<G->NumSurfaces; ns++)
   for(int np=0; np<G->Surfaces[ns]->NumPanels; np++)
    for(int iv=0; iv<3; iv++)
     { int nv = G->Surfaces[ns]->Panels[np]->VI[iv];
       double ThisDistance=VecDistance(X, G->Surfaces[ns]->Vertices + 3*nv);
       MinDistance = fmin(MinDistance, ThisDistance);
     }
 
  for(int ns=0; ns<G->NumSurfaces; ns++)
   for(int np=0; np<G->Surfaces[ns]->NumPanels; np++)
    for(int iv=0; iv<3; iv++)
     {  
       RWGSurface *S = G->Surfaces[ns];
       int nv = S->Panels[np]->VI[iv];
       double ThisDistance=VecDistance(X, S->Vertices + 3*nv);
       if ( fabs(ThisDistance-MinDistance) > 1.0e-6*ThisDistance )
        continue;

       int ne = S->Panels[np]->EI[iv];
       if (ne>=0)
        { 
          RWGEdge *NewEdge = (RWGEdge *)mallocEC(sizeof(RWGEdge));
          NewEdge->iQP = S->Panels[np]->VI[iv];
          NewEdge->iV1 = S->Panels[np]->VI[(iv+1)%3];
          NewEdge->iV2 = S->Panels[np]->VI[(iv+2)%3];
          NewEdge->iPPanel = np;
          NewEdge->PIndex = iv;
          NewEdge->iQM = NewEdge->iMPanel = NewEdge->MIndex = -1;
          
          double *QP = S->Vertices + 3*NewEdge->iQP;
          double *V1 = S->Vertices + 3*NewEdge->iV1;
          double *V2 = S->Vertices + 3*NewEdge->iV2;
          NewEdge->Length = VecDistance(V1, V2);
          NewEdge->Centroid[0] = 0.5*(V1[0]+V2[0]);
          NewEdge->Centroid[1] = 0.5*(V1[1]+V2[1]);
          NewEdge->Centroid[2] = 0.5*(V1[2]+V2[2]);
          NewEdge->Radius = VecDistance(QP, NewEdge->Centroid);
          ne = NewEdge->Index = -1 * (++(S->NumHalfRWGEdges));
          S->HalfRWGEdges=(RWGEdge **)reallocEC(S->HalfRWGEdges, (S->NumHalfRWGEdges)*sizeof(RWGEdge *));
          S->HalfRWGEdges[S->NumHalfRWGEdges-1] = NewEdge;
        }
       PortEdges.push_back(new RWGPortEdge(ns, ne, 0, 0, +1.0));
     }

  return PortEdges;
}

/***************************************************************/
/* port file syntax example                                    */
/*  PORT                                                       */
/*   POSITIVE  A1 A2 A3 B1 B2 B3 ... Z1 Z2 Z3                  */
/*   NEGATIVE  A1 A2 A3 B1 B2 B3 ... Z1 Z2 Z3                  */
/*  ENDPORT                                                    */
/***************************************************************/
RWGPortList *ParsePortFile(RWGGeometry *G, const char *PortFileName)
{
  /***************************************************************/
  /* try to open the file ****************************************/
  /***************************************************************/
  FILE *f=fopenPath(getenv("SCUFF_DATA_PATH"),PortFileName,"r");
  if (f==0)
   ErrExit("could not open file %s",PortFileName);
  
  RWGPortList *PortList = new RWGPortList;

  /***************************************************************/
  /* read through the file and parse lines one-at-a-time *********/
  /***************************************************************/
  RWGPort *CurrentPort=0;
  int CurrentPortIndex=0;
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
      { 
        if (CurrentPort!=0) ErrExit("%s:%i: syntax error (missing ENDPORT?)",PortFileName,LineNum);
        CurrentPort = new RWGPort;
        CurrentPort->Perimeter[_PLUS]=CurrentPort->Perimeter[_MINUS]=0.0;
      }
     /*--------------------------------------------------------------*/  
     else if ( !StrCaseCmp(Tokens[0],"ENDPORT") )
      { 
        PortList->Ports.push_back(CurrentPort);
        CurrentPort=0;
        CurrentPortIndex++;
      }
     /*--------------------------------------------------------------*/
     else if (   !StrCaseCmp(Tokens[0],"POSITIVE") || !StrCaseCmp(Tokens[0],"NEGATIVE")
              || !StrCaseCmp(Tokens[0],"PPOLYGON") || !StrCaseCmp(Tokens[0],"MPOLYGON")
             )
      { 
        if ( (NumTokens<4) || (NumTokens%3 != 1) )
         ErrExit("%s:%i: number of arguments to %s must be a multiple of 3 and >=3",PortFileName,LineNum,Tokens[0]);

        dVec VertexCoordinates(NumTokens-1);
        for(int nt=1; nt<NumTokens; nt++) sscanf(Tokens[nt],"%le",&(VertexCoordinates[nt-1]));
        int NumVertices = VertexCoordinates.size() / 3;

        int Pol = (toupper(Tokens[0][0])=='P') ? _PLUS : _MINUS;
        if (NumVertices==1)
         { 
           // point-like port 
           RWGPortEdgeList NewPEs=AddPointPort(G, &(VertexCoordinates[0]));
           for(unsigned npe=0; npe<NewPEs.size(); npe++)
            { RWGPortEdge *PE   = NewPEs[npe];
              PE->nPort = CurrentPortIndex;
              PE->Pol   = Pol;
              CurrentPort->PortEdges[Pol].push_back(PE);
              PortList->PortEdges.push_back(PE);
            }
         }
        else
         { // exterior edge port
           for(int ns=0; ns<G->NumSurfaces; ns++)
            { iVec neList = FindEdgesInPolygon(G->Surfaces[ns], VertexCoordinates);
              for(unsigned nne=0; nne<neList.size(); nne++)
               { RWGPortEdge *PE = new RWGPortEdge(ns, -1-neList[nne], CurrentPortIndex, Pol, -1.0);
                 CurrentPort->PortEdges[Pol].push_back(PE);
                 PortList->PortEdges.push_back(PE);
               }
            }
         }
      }
     /*--------------------------------------------------------------*/
     else
      ErrExit("%s:%i: unknown keyword %s",PortFileName,LineNum,Tokens[0]);

   } // while( fgets(buffer, 1000, f) )
  fclose(f);

  /***************************************************************/
  /* compute bounding box and perimeters  ************************/
  /***************************************************************/
  double *RMin = PortList->RMinMax + 0;
  double *RMax = PortList->RMinMax + 3;
  RMin[0]=RMin[1]=RMin[2] = +1.23e+45;
  RMax[0]=RMax[1]=RMax[2] = -1.23e+45;
  for(unsigned npe=0; npe<PortList->PortEdges.size(); npe++)
   { RWGPortEdge *PE = PortList->PortEdges[npe];
     RWGSurface *S   = G->Surfaces[PE->ns];
     RWGEdge *E      = S->GetEdgeByIndex(PE->ne);
     PortList->Ports[PE->nPort]->Perimeter[PE->Pol] += E->Length;
     UpdateRMinMax(RMin, RMax, S->Vertices + 3*E->iQP);
     UpdateRMinMax(RMin, RMax, S->Vertices + 3*E->iV1);
     UpdateRMinMax(RMin, RMax, S->Vertices + 3*E->iV2);
   }

  /***************************************************************/
  /* write summary of port list to log file **********************/
  /***************************************************************/
  for(unsigned np=0; np<PortList->Ports.size(); np++)
   for(int Pol=_PLUS; Pol<=_MINUS; Pol++)
    { RWGPort *Port = PortList->Ports[np];
      int NPE = Port->PortEdges[Pol].size();
      double Perimeter=Port->Perimeter[Pol];
      Log("Port %2i (%s): perimeter %e, %i edges%s",np,PolNames[Pol],Perimeter,NPE,NPE==0 ? "" : "=[");
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
RWGPortList *ReadGDSIIPorts(RWGGeometry *G, const char *GDSIIFileName, int Layer)
{
#ifndef HAVE_LIBGDSII
  (void) G; (void) GDSIIFileName; (void) Layer;
  ErrExit("SCUFF-EM must be compiled with libGDSII to process GDSII files");
  return 0;
#else

  // Loop over text strings on the given GDSII layer (or all layers if Layer==-1).
  TextStringList TextStrings = GetTextStrings(GDSIIFileName, Layer);
  int MaxPort=0; 
  multimap<int,dVec> PMPolygons[2];
  for(size_t nt=0; nt<TextStrings.size(); nt++)
   {
     // loop for text strings of the form PORT 3+ or PORT 5-
     char *Text    = TextStrings[nt].Text;
     int nPort;
     char Sign;
     if ( 2 != sscanf(Text,"PORT %i%c",&nPort,&Sign) || !strchr("+pP-mM",Sign) ) continue;
     int PM = (strchr("-mM",Sign) ? 1 : 0);
     if (nPort>MaxPort) MaxPort=nPort;

     // find a polygon on the same layer containing the reference point;
     // if none found, take the reference point itself to be the polygon
     int TextLayer = TextStrings[nt].Layer;
     dVec TextXY   = TextStrings[nt].XY;
     PolygonList Polygons=GetPolygons(GDSIIFileName, TextLayer);
     bool FoundPolygon=false;
     for(size_t np=0; !FoundPolygon && np<Polygons.size(); np++)
      if (libGDSII::PointInPolygon( Polygons[np], TextXY[0], TextXY[1]))
       { PMPolygons[PM].insert( pair<int,dVec>(nPort,Polygons[np]) );
         FoundPolygon=true;
       }
     if (!FoundPolygon)
      PMPolygons[PM].insert( pair<int,dVec>(nPort,TextXY) );
   }
  Log("%s: found %i port definitions (total %lu port-terminal polygons)",MaxPort,PMPolygons[0].size()+PMPolygons[1].size());

  // Write what we found to a `.ports` file
  char PortFileName[100];
  snprintf(PortFileName,100,"/tmp/%s.ports.XXXXXX",GetFileBase(GDSIIFileName));
  FILE *f = (mkstemp(PortFileName)==-1) ? 0 : fopen(PortFileName,"w");
  if (!f) ErrExit("%s:%i: could not create temporary file %s",__FILE__,__LINE__,PortFileName);
  for(int nPort=1; nPort<=MaxPort; nPort++)
   for(int PM=0; PM<2; PM++)
    { if (PM==0) fprintf(f,"PORT \n");
      pair< multimap<int,dVec>::iterator , multimap<int,dVec>::iterator > Range=PMPolygons[PM].equal_range(nPort);
      int NumPolygons=0;
      for(multimap<int,dVec>::iterator p=Range.first; p!=Range.second; p++, NumPolygons++)
       { fprintf(f," %cPOLYGON ",PM==0 ? 'P' : 'M');
         for(size_t n=0; n<p->second.size(); n++) fprintf(f,"%g ",p->second[n]);
         fprintf(f,"\n");
       }
      if (PM==0 && NumPolygons==0)
       ErrExit("%s: no positive terminal defined for port %i/%i",GDSIIFileName,nPort,MaxPort);
      if (PM==1) fprintf(f,"ENDPORT \n");
    }
  fclose(f);
  printf("Wrote port list to %s.\n",PortFileName);

  RWGPortList *PortList=ParsePortFile(G, PortFileName);
  if (!CheckEnv("SCUFF_RETAIN_GDSII_PORTFILE"))
   unlink(PortFileName);
  return PortList;

#endif
}

/***************************************************************/
/* PortBFInteractionMatrix[nbf, np]                            */
/*  = interaction of basis function #nbf with port #np         */
/*                                                             */
/*         ( overlap of basis function #nbf with     )         */
/*  = -1 x ( the electric field produced by port #np )         */
/*         ( driven by a unit-strength current       )         */
/*                                                             */
/* so column #np of this matrix can be used as the RHS vector  */
/* in a BEM scattering calculation to get the response of the  */
/* the system to a unit-strength current at port #np. (Thus,   */
/* the columns of this matrix are the vectors called           */
/* \mathbf{r}_p (p=1,\cdots,NumPorts) in the memo).            */
/***************************************************************/
void scuffSolver::AssemblePortBFInteractionMatrix(cdouble Omega)
{ 
  // don't need to recompute the matrix if has been assembled at this frequency already;
  if (Omega==OmegaPBFI) return;
  OmegaPBFI=Omega;

  Log("Computing port<-->BF interaction matrix");

  int NBF = G->TotalBFs;
  if (PBFIMatrix==0) PBFIMatrix = new HMatrix(NBF,      NumPorts, LHM_COMPLEX);
  if (PPIMatrix==0)  PPIMatrix  = new HMatrix(NumPorts, NumPorts, LHM_COMPLEX);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (G->Substrate)
   {  
     HMatrix XMatrix(3,2,PortList->RMinMax);
     double RhoMinMax[2], zMinMax[2];
     GetRzMinMax(G, &XMatrix, RhoMinMax, zMinMax);
     Log("  initializing ScalarGF interpolator for Rho range (%e,%e)",RhoMinMax[0],RhoMinMax[1]);
     bool PPIsOnly=true;
     bool Subtract=true;
     bool RetainSingularTerms=false;
     G->Substrate->InitScalarGFInterpolator(Omega, RhoMinMax[0], RhoMinMax[1], 0.0, 0.0,
                                            PPIsOnly, Subtract, RetainSingularTerms);
   }

  int Order=-1; CheckEnv("SCUFF_PORTBFI_ORDER",&Order);

  int NumPortEdges = PortList->PortEdges.size();
  int NBFPE        = NBF + (PPIMatrix ? NumPortEdges : 0);
  int NumPairs     = NumPortEdges * NBFPE;
  
  int NumThreads             = GetNumThreads();
  bool AllocateThreadBuffers = (NumThreads>1 && !CheckEnv("SCUFF_NO_THREAD_BUFFERS"));

  HMatrix **PBFPIByThread = 0;
  if (AllocateThreadBuffers)
   { PBFPIByThread = new HMatrix *[NumThreads];
     for(int nt=0; nt<NumThreads; nt++)
      PBFPIByThread[nt] = new HMatrix(NBFPE, NumPorts, LHM_COMPLEX);
   }

  Log("  getting %i port-BF interactions (%i threads)",NumPairs, NumThreads);
  Tic();
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int nPair=0; nPair<NumPairs; nPair++)
   { 
     if (G->LogLevel>=SCUFF_VERBOSE2) LogPercent(nPair, NumPairs);

     int nPortEdge       = nPair/NBFPE;
     int nBFPE           = nPair%NBFPE;

     RWGPortEdge *PE     = PortList->PortEdges[nPortEdge];
     int nPort           = PE->nPort;
     int nsLeft          = PE->ns;
     int neLeft          = PE->ne;
     double Weight       = -1.0*PE->Sign*PolSigns[PE->Pol] / PortList->Ports[nPort]->Perimeter[PE->Pol];

     int nsRight, neRight, nPortPrime=-1;
     if (nBFPE < NBF)
      { if (!G->ResolveEdge(nBFPE, &nsRight, &neRight, 0))
         ErrExit("%s:%i: internal error(%i,%i,%i)\n",nBFPE,nsRight,neRight);
      }
     else
      { RWGPortEdge *PEP  = PortList->PortEdges[nBFPE - NBF];
        nPortPrime        = PEP->nPort;
        nsRight           = PEP->ns;
        neRight           = PEP->ne;
        Weight           *= -1.0*PEP->Sign*PolSigns[PEP->Pol] / PortList->Ports[nPortPrime]->Perimeter[PEP->Pol];
      }

     cdouble ME;
     GetMOIMatrixElement(G, nsLeft, neLeft, nsRight, neRight, Omega, &ME, Order);

     if (AllocateThreadBuffers)
      PBFPIByThread[GetThreadNum()]->AddEntry(nBFPE, nPort, Weight*ME);
     else if (nBFPE < NBF)
      {
#ifdef USE_OPENMP
#pragma omp critical
#endif
        PBFIMatrix->AddEntry(nBFPE, nPort, Weight*ME);
      }
     else // (nBFPE > NBF)
      {
#ifdef USE_OPENMP
#pragma omp critical
#endif
        PPIMatrix->AddEntry(nPortPrime, nPort, Weight*ME);
      }
   } // for(int nPEBF=0...
  double Time=Toc();
  Log("...RHS port contributions done in %e s",Time);

  if (AllocateThreadBuffers)
   { PBFIMatrix->Zero();
     if (PPIMatrix) PPIMatrix->Zero();
     for(int nt=0; nt<NumThreads; nt++)
      { for(int SourcePort=0; SourcePort<NumPorts; SourcePort++)
         { for(int nbf=0; nbf<NBF; nbf++)
            PBFIMatrix->AddEntry(nbf, SourcePort, PBFPIByThread[nt]->GetEntry(nbf,SourcePort));
           for(int npe=0; npe<(NBFPE-NBF); npe++)
            { int DestPort= PortList->PortEdges[npe]->nPort;
               PPIMatrix->AddEntry(DestPort, SourcePort, PBFPIByThread[nt]->GetEntry(NBF+npe,SourcePort));
            }
         }
        delete PBFPIByThread[nt];
      }
     delete[] PBFPIByThread;
   }
  
}

/***************************************************************/
/* if X lies on a meshed surface, compute the charge and       */
/* current densities at X.                                     */
/* iwSigmaK[0] = iw*Sigma                                      */
/* iwSigmaK[1] = Kx                                            */
/* iwSigmaK[2] = Ky                                            */
/* iwSigmaK[3] = Kz                                            */
/***************************************************************/
void scuffSolver::EvalSourceDistribution(const double X[3], cdouble iwSigmaK[4])
{
  iwSigmaK[0]=iwSigmaK[1]=iwSigmaK[2]=iwSigmaK[3]=0.0;

  // get BF contribution
  int ns, np;
  cdouble KNX[6], iwSigmaTauX[2];
  G->EvalSourceDistribution(X, KN, KNX, iwSigmaTauX, &ns, &np);
  if (ns==-1) return; // X does not lie on a meshed surface

  iwSigmaK[0] = iwSigmaTauX[0];
  iwSigmaK[1] = KNX[0];
  iwSigmaK[2] = KNX[1];
  iwSigmaK[3] = KNX[2];

  if ( CheckEnv("SCUFF_RFPSD_PORTONLY") )
   memset(iwSigmaK, 0, 4*sizeof(cdouble));
  if ( CheckEnv("SCUFF_RFPSD_BFONLY") )
   return;

  // add port contribution, if any, by looking at the three edges of the
  // panel containing X and checking if any of them belong to a port
  if (PortCurrents==0) return;
  RWGSurface *S    = G->Surfaces[ns];
  RWGPanel *P      = S->Panels[np];
  int NumPortEdges = PortList->PortEdges.size();
  for(int nei=0; nei<3; nei++)
   {
     int ne = P->EI[nei];
     if (ne>=0) continue; // not a half-RWG edge, can't be a port edge

     int npe;
     RWGPortEdge *PE=0;
     for(npe=0; npe<NumPortEdges; npe++)
      { PE = PortList->PortEdges[npe];
        if ( (PE->ns == ns) && (PE->ne == ne ) )
         break;
      }
     if (npe==NumPortEdges) continue; // is a half-RWG edge, but not a port edge

     int nPort      = PE->nPort;
     int Pol        = PE->Pol;
     double L       = PortList->Ports[nPort]->Perimeter[Pol];
     cdouble Weight = PE->Sign*PolSigns[Pol]*PortCurrents[nPort] / L;
     if(Weight==0.0) continue;

     RWGEdge *E       = S->GetEdgeByIndex(ne);
     double *Q        = S->Vertices + 3*E->iQP;
     double DivbAlpha = E->Length / P->Area;
     double bAlpha[3];
     VecSub(X, Q, bAlpha);
     VecScale(bAlpha, 0.5*DivbAlpha);
     iwSigmaK[0]    += Weight*DivbAlpha;
     iwSigmaK[1]    += Weight*bAlpha[0];
     iwSigmaK[2]    += Weight*bAlpha[1];
     iwSigmaK[3]    += Weight*bAlpha[2];
   }
}

/***************************************************************/
/* contribution of port currents to panel source densities     */
/***************************************************************/
HMatrix *scuffSolver::GetPanelSourceDensities(HMatrix *PSDMatrix)
{ 
  /***************************************************************/
  /* get basis-function contributions     ************************/
  /***************************************************************/
  cdouble Omega = G->StoredOmega;
  PSDMatrix=G->GetPanelSourceDensities(Omega, KN, PSDMatrix);

  if ( CheckEnv("SCUFF_RFPSD_PORTONLY") )
   PSDMatrix->Zero();
  if ( CheckEnv("SCUFF_RFPSD_BFONLY") )
   return PSDMatrix;

  /***************************************************************/
  /* add port contributions **************************************/
  /***************************************************************/
  for(unsigned nPE=0; nPE<PortList->PortEdges.size(); nPE++)
   { 
     RWGPortEdge *PE = PortList->PortEdges[nPE];
     cdouble Weight  = PortCurrents[PE->nPort] / (PortList->Ports[PE->nPort]->Perimeter[PE->Pol]);
     if (Weight==0.0) continue;
     RWGSurface *S   = G->Surfaces[PE->ns];
     RWGEdge *E      = S->GetEdgeByIndex(PE->ne);
     double Length   = E->Length; 
     double *QP      = S->Vertices + 3*(E->iQP);
     int PanelIndex  = G->PanelIndexOffset[S->Index] + E->iPPanel;
     RWGPanel *Panel = S->Panels[E->iPPanel];
     cdouble PreFac  = PE->Sign*PolSigns[PE->Pol]*Weight*Length/(2.0*Panel->Area);
     double XmQ[3];
     VecSub(Panel->Centroid, QP, XmQ);
     PSDMatrix->AddEntry( PanelIndex,  4, 2.0*PreFac/(II*Omega));  // rho
     PSDMatrix->AddEntry( PanelIndex,  5, PreFac*XmQ[0] );  // K_x
     PSDMatrix->AddEntry( PanelIndex,  6, PreFac*XmQ[1] );  // K_y
     PSDMatrix->AddEntry( PanelIndex,  7, PreFac*XmQ[2] );  // K_z
   } // for(int nPE=0...)
 
  return PSDMatrix;
}

} // namespace scuff
