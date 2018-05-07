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
 * RFSolver.cc  -- main code file for SCUFF-EM RF solver
 *
 * homer reid   -- 9/2011 - 4/2018
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
#include <libIncField.h>

#include "RFSolver.h"

#define II cdouble(0.0,1.0)

using namespace scuff;

namespace scuff {

/***************************************************************/
/* RF solver class constructor #1: construct from a .scuffgeo  */
/* file plus a .ports file                                     */
/***************************************************************/
RFSolver::RFSolver(const char *scuffgeoFileName, const char *portFileName)
{ 
  InitSolver();
  scuffgeoFile = scuffgeoFileName ? strdup(scuffgeoFileName) : 0;
  portFile     = portFileName     ? strdup(portFileName)     : 0;
}

/***************************************************************/
/* constructor helper routine that just handles minor stuff;   */
/* most actual initialization is handled by InitGeometry below,*/
/* which is called lazily just-in-time as needed               */
/***************************************************************/
void RFSolver::InitSolver()
{
  /*--------------------------------------------------------------------------*/
  /* initialize internal variables needed to perform RF calculations (which   */
  /* are allocated lazily)                                                    */
  /*--------------------------------------------------------------------------*/
  G            = 0;
  PortList     = 0;
  EEPTable     = 0;
  NumPorts     = 0;
  FileBase     = 0;

  Omega        = -1.0;
  M            = 0;
  MClean       = false;
  PBFIMatrix   = 0; 
  PPIMatrix    = 0;
  PBFIClean    = false;
  KN           = 0;
  PortCurrents = 0;

  DisableSystemBlockCache = false;
  OmegaCache = HUGE_VAL;
  TBlocks = 0;
  UBlocks = 0;

  RetainContributions = CONTRIBUTION_ALL;
  SubstrateFile = 0;
  scuffgeoFile  = 0;
  portFile      = 0;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
RFSolver::~RFSolver()
{
  if (TBlocks)
   { for(int ns=0; ns<G->NumSurfaces; ns++)
      delete TBlocks[ns];
     delete TBlocks;
   }
  if (UBlocks)
   { int NS = G->NumSurfaces, NU=(NS-1)*NS/2;
     for(int nu=0; nu<NU; nu++)
      delete UBlocks[nu];
     delete UBlocks;
   }

  if (PortCurrents) delete PortCurrents;
  if (KN)           delete KN;
  if (PBFIMatrix)   delete PBFIMatrix;
  if (PPIMatrix)    delete PPIMatrix;
  if (M)            delete M;
  if (FileBase)     free(FileBase);
  if (EEPTable)     delete EEPTable;
  if (PortList)     delete PortList;
  if (G)            delete G;
}

/********************************************************************/
/* routines for building up geometries line-by-line from python     */
/* scripts; these just make a note of whatever feature the user     */
/* added, with the actual initialization done later by InitGeometry()*/
/********************************************************************/
void RFSolver::AddSubstrateLayer(double zInterface, cdouble Epsilon, cdouble Mu)
{ 
  if (G) ErrExit("can't modify substrate after geometry has been initialized");

  char Line[100];
  if (isinf(real(Epsilon)))
   snprintf(Line,100,"%e GROUNDPLANE\n",zInterface);
  else if (Mu==1.0)
   snprintf(Line,100,"%e CONST_EPS_%s\n",zInterface,z2s(Epsilon));
  else
   snprintf(Line,100,"%e CONST_EPS_%s_MU_%s\n",zInterface,z2s(Epsilon),z2s(Mu));
  SubstrateLayers.insert(std::pair<double, char *>(-zInterface, strdup(Line)));
}

void RFSolver::AddGroundPlane(double zGP)
{ AddSubstrateLayer(zGP, HUGE_VAL); }

void RFSolver::SetSubstratePermittivity(cdouble Epsilon)
{ AddSubstrateLayer(0.0, Epsilon); }

void RFSolver::SetSubstrateThickness(double h)
{ AddGroundPlane(-h); }

void RFSolver::SetSubstrateFile(const char *_SubstrateFile)
{ 
  if (G) ErrExit("can't modify substrate after geometry has been initialized");
  SubstrateFile = strdup(_SubstrateFile);
}

void RFSolver::SetGeometryFile(const char *_scuffgeoFile)
{
  if (MeshFiles.size() > 0)
   { Warn("can't add geometry file after adding metal traces");
     return;
   }
  if (_scuffgeoFile==0) ErrExit("%s:%i: internal error",__FILE__,__LINE__);
  if (scuffgeoFile)
   { Warn("overwriting existing geometry file %s with new file %s",scuffgeoFile,_scuffgeoFile);
     free(scuffgeoFile);
   }
  scuffgeoFile = vstrdup(_scuffgeoFile);
}

void RFSolver::AddMetalTraceMesh(const char *MeshFile, const char *Transformation)
{
  char *ErrMsg=0;

  if (G) 
   ErrMsg=vstrdup("can't add metal traces after geometry has been initialized");
  if (scuffgeoFile)
   ErrMsg=vstrdup("can't add metal traces to existing geometry file %s",scuffgeoFile);

  // check that mesh is valid
  if (!ErrMsg) 
   { RWGSurface *S = new RWGSurface(MeshFile);
     if (S->ErrMsg) ErrMsg = strdup(S->ErrMsg);
     //delete S;
   }

  // check that transformation is valid
  if (!ErrMsg && Transformation)
   { GTransformation *GT = new GTransformation(Transformation, &ErrMsg);
     //delete GT;
   }

  if (ErrMsg)
   { Warn("AddMetalTraceMesh: %s (ignoring)",ErrMsg);
     free(ErrMsg);
     return;
   }

  MeshFiles.push_back(strdup(MeshFile));
  MeshTransforms.push_back( Transformation ? strdup(Transformation) : 0);
}

void RFSolver::SetPortFile(const char *_portFile)
{ char *ErrMsg=0;
  if (G)
   ErrMsg=vstrdup("can't add port file after geometry has been initialized");
  else if (PortTerminalVertices.size() > 0)
   ErrMsg=vstrdup("can't add port file after calling AddPort() or AddPortTerminal()");
  if (ErrMsg)
   { Warn("SetPortFile: %s (ignoring)",ErrMsg);
     free(ErrMsg);
     return;
   }
  if (portFile) free(portFile);
  portFile = strdup(_portFile);
}

void RFSolver::AddPort(const dVec PVertices, const dVec MVertices)
{ char *ErrMsg=0;
  if (G)
   ErrMsg=vstrdup("can't add port after geometry has been initialized");
  else if (portFile)
   ErrMsg=vstrdup("can't add port after calling SetPortFile()");
  if (ErrMsg)
   { Warn("AddPort: %s (ignoring)",ErrMsg);
     free(ErrMsg);
     return;
   }

  // the following could happen if the user calls AddPortTerminal('P',...) followed by AddPort()
  if ( (PortTerminalVertices.size() % 2) ==1 )
   PortTerminalVertices.push_back(dVec());

  PortTerminalVertices.push_back(PVertices);
  PortTerminalVertices.push_back(MVertices);
}

void RFSolver::AddPort(const dVec PVertices)
{ AddPort(PVertices, dVec()); }

void RFSolver::AddPortTerminal(char PM, const dVec Vertices)
{ 
  bool NewTerminalIsNegative     = ( (toupper(PM)=='M') || (PM=='-') );
  bool CurrentTerminalIsNegative = ((PortTerminalVertices.size() % 2) != 1);
  bool SignChange = (NewTerminalIsNegative != CurrentTerminalIsNegative);
  bool AddNewTerminal = (PortTerminalVertices.size()==0 || SignChange);

  if (AddNewTerminal)
   PortTerminalVertices.push_back(Vertices);
  else
   { PortTerminalVertices.back().push_back(HUGE_VAL);
     for(size_t n=0; n<Vertices.size(); n++) 
      PortTerminalVertices.back().push_back(Vertices[n]);
   }
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RFSolver::InitGeometry()
{ 
  if (G) return;

  FileBase = strdup( scuffgeoFile ? GetFileBase(scuffgeoFile) : "pyscuff");

  /*--------------------------------------------------------------*/
  /*- write .scuffgeo file from user's specifications, then try  -*/
  /*- to create the geometry                                     -*/
  /*--------------------------------------------------------------*/
  bool ownsGeoFile = false;
  if (!scuffgeoFile)
   { 
     ownsGeoFile = true;
     char buffer[100];
     snprintf(buffer,100,".pyscuffgeo.XXXXXX");
     scuffgeoFile=vstrdup(buffer);
     FILE *f=fdopen( mkstemp(scuffgeoFile), "w");
     if (!f) ErrExit("could not write file %s",scuffgeoFile);

     if (MeshFiles.size()==0) ErrExit("no metal traces specified");
     if (MeshFiles.size() != MeshTransforms.size()) ErrExit("%s:%i: internal error");
     for(size_t n=0; n<MeshFiles.size(); n++)
      { fprintf(f,"OBJECT %s_%lu \n MESHFILE %s \n",GetFileBase(MeshFiles[n]),n,MeshFiles[n]);
        free(MeshFiles[n]);
        if (MeshTransforms[n]) { fprintf(f," %s\n",MeshTransforms[n]); free(MeshTransforms[n]); }
        fprintf(f,"ENDOBJECT\n");
      }
     MeshFiles.clear();
     MeshTransforms.clear();
     fclose(f);
   }
  RWGGeometry::UseHRWGFunctions=false;
  G = new RWGGeometry(scuffgeoFile);

  /*--------------------------------------------------------------*/
  /*- write .ports file from user's specifications ---------------*/
  /*--------------------------------------------------------------*/
  bool ownsPortFile = false;
  if (!portFile)
   { 
     ownsPortFile = true;
     char buffer[100];
     snprintf(buffer,100,".pyscuffports.XXXXXX");
     portFile=vstrdup(buffer);
     FILE *f=fdopen( mkstemp(portFile), "w");
     if (!f) ErrExit("could not write file %s",portFile);

     for(size_t nPort=0; nPort<PortTerminalVertices.size()/2; nPort++)
      { fprintf(f,"PORT");
        for(int Pol=_PLUS; Pol<=_MINUS; Pol++)
         for(size_t nv=0; nv<PortTerminalVertices[2*nPort+Pol].size(); nv++)
          { double V = PortTerminalVertices[2*nPort+Pol][nv];
            if (nv==0 || isinf(V)) fprintf(f,"\n   %s ",(Pol==_PLUS ? "POSITIVE" : "NEGATIVE"));
            if (isinf(V)) continue;
            fprintf(f,"%e ",V);
          }
        fprintf(f,"\nENDPORT\n\n");
      }
     fclose(f);
     PortTerminalVertices.clear();
   }
  if (!portFile) ErrExit("no ports specified");
  PortList = ParsePortFile(G, portFile);
  NumPorts = PortList->Ports.size();
  if (NumPorts==0) Warn("no ports found");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (SubstrateLayers.size()>0)
   { char *SubstrateDescription=0;
     for(std::map<double, char *>::iterator it = SubstrateLayers.begin(); it!=SubstrateLayers.end(); it++)
      SubstrateDescription=vstrappend(SubstrateDescription,"%s\n",it->second);
     SubstrateLayers.clear();
     G->Substrate = CreateLayeredSubstrate(SubstrateDescription);
     free(SubstrateDescription);
   }

  EEPTable = CheckEnv("SCUFF_IGNORE_EQUIVALENT_EDGES") ? 0 : new EquivalentEdgePairTable(G);
 
  //if(ownsGeoFile)
  // unlink(scuffgeoFile);
  //if(ownsPortFile)
  // unlink(portFile);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RFSolver::UpdateSystemMatrix()
{ 
  if (!G) InitGeometry();

  // lazy allocation of accelerator
  if (TBlocks==0)
   { int NS = G->NumSurfaces, NU=(NS-1)*NS/2;
     TBlocks = new HMatrix *[NS];
     UBlocks = new HMatrix *[NU];
     for(int nsa=0; nsa<NS; nsa++)
      { int NBFA = G->Surfaces[nsa]->NumBFs;
        int Mate = G->Mate[nsa];
        TBlocks[nsa] = (Mate!=-1) ? TBlocks[Mate] : new HMatrix(NBFA, NBFA, LHM_COMPLEX);
        for(int nsb=nsa+1, nb=0; nsb<NS; nsb++, nb++)
         { int NBFB = G->Surfaces[nsb]->NumBFs;
           UBlocks[nb] = new HMatrix(NBFA, NBFB, LHM_COMPLEX);
         }
      }
   }

  // TBlocks recomputed only if frequency has changed
  if (Omega!=OmegaCache)
   for(int ns=0; ns<G->NumSurfaces; ns++)
    if (G->Mate[ns]==-1)
     AssembleMOIMatrixBlock(G, ns, ns, Omega, TBlocks[ns], 0, 0, EEPTable);

  // UBlocks recomputed if frequency has changed or surfaces were moved
  for(int nsa=0, nb=0; nsa<G->NumSurfaces; nsa++)
   for(int nsb=nsa+1; nsb<G->NumSurfaces; nsb++, nb++)
    if (Omega!=OmegaCache || G->SurfaceMoved[nsa] || G->SurfaceMoved[nsb])
     AssembleMOIMatrixBlock(G, nsa, nsb, Omega, UBlocks[nb], 0, 0, EEPTable);

  // stamp blocks into M matrix
  for(int nsa=0, nb=0; nsa<G->NumSurfaces; nsa++)
   { int OffsetA = G->BFIndexOffset[nsa];
     M->InsertBlock(TBlocks[nsa], OffsetA, OffsetA);
     for(int nsb=nsa+1; nsb<G->NumSurfaces; nsb++, nb++)
      {  int OffsetB = G->BFIndexOffset[nsb];
         M->InsertBlock(UBlocks[nb],OffsetA, OffsetB);
         M->InsertBlockTranspose(UBlocks[nb],OffsetB, OffsetA);
      }
   }

  OmegaCache=Omega;
}

void RFSolver::AssembleSystemMatrix(double Freq)
{ 
  if (!G) InitGeometry();

  if (M==0) M=G->AllocateBEMMatrix();
  Omega = Freq * FREQ2OMEGA;

  Log("Assembling BEM matrix at f=%g GHz...",Freq);
  if (DisableSystemBlockCache)
   AssembleMOIMatrix(G, Omega, M, EEPTable);
  else
   UpdateSystemMatrix();
  Log("Factorizing...");
  M->LUFactorize();
  PBFIClean=false;
}

void RFSolver::Solve(cdouble *CallerPortCurrents)
{
  if (M==0)
   ErrExit("RFSolver: AssembleSystemMatrix() must be called before Solve()");

  /*--------------------------------------------------------------*/
  /*- lazy allocation of internal variables ----------------------*/
  /*--------------------------------------------------------------*/
  int NBF = G->TotalBFs;
  if (PBFIMatrix==0)   PBFIMatrix   = new HMatrix(NBF,      NumPorts, LHM_COMPLEX);
  if (PPIMatrix==0)    PPIMatrix    = new HMatrix(NumPorts, NumPorts, LHM_COMPLEX);
  if (KN==0)           KN           = new HVector(NBF,                LHM_COMPLEX);
  if (PortCurrents==0) PortCurrents = new cdouble[NumPorts];
  
  /*--------------------------------------------------------------*/
  /*- (re)compute port-BF interaction matrix as necessary --------*/
  /*--------------------------------------------------------------*/
  AssemblePortBFInteractionMatrix();

  if (CallerPortCurrents) // if zero, assumes PortCurrents has already been initialized
   memcpy(PortCurrents, CallerPortCurrents, NumPorts*sizeof(cdouble));

  // form RHS vector as weighted linear combination of contributions from each port
  KN->Zero();
  for(int np=0; np<NumPorts; np++)
   for(int nbf=0; nbf<NBF; nbf++)
    KN->AddEntry(nbf, PortCurrents[np]*PBFIMatrix->GetEntry(nbf,np));

  // solve the system
  Log("LU-solving...");
  M->LUSolve(KN);
}

void RFSolver::Solve(cdouble PortCurrent, int WhichPort)
{
  if (PortCurrents==0) PortCurrents = new cdouble[NumPorts];
  memset(PortCurrents, 0, NumPorts*sizeof(cdouble));
  PortCurrents[WhichPort]=PortCurrent;
  Solve(0);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RFSolver::PlotGeometry(const char *PPFormat, ...)
{
  if (!G) InitGeometry();

  /***************************************************************/
  /* open the file ***********************************************/
  /***************************************************************/
  char PPFileName[1000];
  if (!PPFormat)
   snprintf(PPFileName,1000,"%s.pp",FileBase);
  else
   { va_list ap;
     va_start(ap,PPFormat);
     vsnprintfEC(PPFileName,997,PPFormat,ap);
     va_end(ap);
   }

  FILE *f=fopen(PPFileName,"w");
  if (!f)
   { Warn("could not open file %s (skipping geometry plot)",PPFileName);
     return;
   }
  fclose(f);

  /***************************************************************/
  /* plot geometry ***********************************************/
  /***************************************************************/
  G->WritePPMesh(PPFileName,FileBase);
  
  /***************************************************************/
  /* plot ports **************************************************/
  /***************************************************************/
  f=fopen(PPFileName,"a");
  for(unsigned nPort=0; nPort<PortList->Ports.size(); nPort++)
   for(int Pol=_PLUS; Pol<=_MINUS; Pol++)
    { 
      RWGPort *Port = PortList->Ports[nPort];
      if (Port->PortEdges[Pol].size() == 0) continue;

      fprintf(f,"View \"Port %i %s terminal\" {\n",nPort+1, (Pol==_PLUS ? "positive" : "negative") );
      /*--------------------------------------------------------------*/
      /*- scalar lines and arrows for port edges                      */
      /*--------------------------------------------------------------*/
      for(unsigned nPE=0; nPE<Port->PortEdges[Pol].size(); nPE++)
       { RWGPortEdge *PE = Port->PortEdges[Pol][nPE];
         RWGSurface *S   = G->Surfaces[PE->ns];
         RWGEdge *E      = S->GetEdgeByIndex(PE->ne);
         double *V1      = S->Vertices + 3*E->iV1;
         double *V2      = S->Vertices + 3*E->iV2;
         int Value       = (Pol==_PLUS ? 1 : -1 ) * (nPort+1);
         fprintf(f,"SL(%e,%e,%e,%e,%e,%e) {%i,%i};\n",
                    V1[0],V1[1],V1[2],V2[0],V2[1],V2[2],Value,Value);

         // arrow to indicate direction of current flow
         double *X0   = E->Centroid;
         double *ZHat = S->Panels[E->iPPanel]->ZHat;
         double *P0   = S->Panels[E->iPPanel]->Centroid;
         double V1mV2[3]; VecSub(V1, V2, V1mV2);
         double Dir[3];   VecCross(ZHat, V1mV2, Dir);
         double X0P[3];   VecScaleAdd(X0, 0.1, Dir, X0P);
         bool DirPointsIntoPanel = (VecDistance(X0P,P0) < VecDistance(X0,P0));
         bool DirShouldPointIntoPanel = (Pol==_MINUS);
         if (PE->Sign==-1.0) DirShouldPointIntoPanel = !DirShouldPointIntoPanel;
         if (DirPointsIntoPanel!=DirShouldPointIntoPanel) VecScale(Dir, -1.0);
         fprintf(f,"VP(%e,%e,%e) {%e,%e,%e};\n",X0[0],X0[1],X0[2],Dir[0],Dir[1],Dir[2]);
       }
      fprintf(f,"};\n");
      fprintf(f,"View[PostProcessing.NbViews-1].CenterGlyphs=1;\n");
      fprintf(f,"View[PostProcessing.NbViews-1].LineWidth = 5;\n");
      fprintf(f,"View[PostProcessing.NbViews-1].LineType  = 1;\n");
      fprintf(f,"View[PostProcessing.NbViews-1].CustomMax = %i;\n",+(NumPorts+1));
      fprintf(f,"View[PostProcessing.NbViews-1].CustomMin = %i;\n",-(NumPorts+1));
      //fprintf(f,"View[PostProcessing.NbViews-1].RangeType = 2;\n");
      fprintf(f,"View[PostProcessing.NbViews-1].ShowScale = 0;\n");
      fprintf(f,"View.Light = 0;\n");
   }

  /***************************************************************/
  /* plot substrate layers and/or ground plane *******************/
  /***************************************************************/
  if (G->Substrate)
   {  
     LayeredSubstrate *S = G->Substrate;
     double RMax[3] = {-HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
     double RMin[3] = { HUGE_VAL,  HUGE_VAL,  HUGE_VAL};
     for(int ns=0; ns<G->NumSurfaces; ns++)
      for(int i=0; i<3; i++)
       { RMax[i] = fmax(RMax[i], G->Surfaces[ns]->RMax[i]);
         RMin[i] = fmin(RMin[i], G->Surfaces[ns]->RMin[i]);
       }
     int NI = S->NumInterfaces;
     if (!isinf(S->zGP)) NI++;
     for(int ni=0; ni<NI; ni++)
      { double z;
        if (ni==S->NumInterfaces)
         { z=S->zGP;
           fprintf(f,"View \"Ground plane\" {\n");
         }
        else
         { z=S->zInterface[ni];
           fprintf(f,"View \"%s upper surfaces\" {\n",S->MPLayer[ni+1]->Name);
         }
        fprintf(f,"SQ(%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e) {%i,%i,%i,%i};\n};\n",
                   RMin[0], RMin[1], z, RMax[0], RMin[1], z, RMax[0], RMax[1], z, RMin[0], RMax[1], z, ni,ni,ni,ni);
      }
   }
  fclose(f);

  fprintf(stdout,"RF geometry plotted to file %s.\n",PPFileName);
}

void RFSolver::PlotGeometry()
 { PlotGeometry(0); }

} // namespace scuff
