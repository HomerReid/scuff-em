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
 * scuffSolver.cc  -- main code file for SCUFF-EM high-level interface
 *
 * homer reid      -- 9/2011 - 4/2018
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
#include <libscuffInternals.h>

#include "scuffSolver.h"

#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif

using namespace scuff;

namespace scuff {

/***************************************************************/
/* RF solver class constructor #1: construct from a .scuffgeo  */
/* file plus a .ports file                                     */
/***************************************************************/
scuffSolver::scuffSolver(const char *Name)
{ 
  SolverName = strdup( Name ? Name : "pyscuff");
  InitSolver();
}

/***************************************************************/
/* constructor helper routine that just handles minor stuff;   */
/* most actual initialization is handled by InitGeometry below,*/
/* which is called lazily just-in-time as needed               */
/***************************************************************/
void scuffSolver::InitSolver()
{
  /*--------------------------------------------------------------------------*/
  /* most internal variables are initialized to zero and allocated lazily later*/
  /*--------------------------------------------------------------------------*/
  G              = 0;
  PortList       = 0;
  NumPorts       = 0;

  M              = 0;
  PBFIMatrix     = 0; 
  PPIMatrix      = 0;
  KN             = 0;
  CachedPortCurrents   = 0;
  CachedIF             = 0;
  OmegaSIE       = CACHE_DIRTY;
  OmegaPBFI      = CACHE_DIRTY;

  TBlocks        = 0;
  UBlocks        = 0;

  Medium         = 0;
  SubstrateFile  = 0;
  SubstrateInitialized = false;
  scuffgeoFile   = 0;
  portFile       = 0;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
scuffSolver::~scuffSolver()
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

  if (CachedPortCurrents) delete CachedPortCurrents;
  if (CachedIF)     delete CachedIF;
  if (KN)           delete KN;
  if (PBFIMatrix)   delete PBFIMatrix;
  if (PPIMatrix)    delete PPIMatrix;
  if (M)            delete M;
  if (PortList)     delete PortList;
  if (G)            delete G;
  if (SolverName)   free(SolverName);
}

/********************************************************************/
/* routines for building up geometries line-by-line from python     */
/* scripts; these just make a note of whatever feature the user     */
/* added, with the actual initialization done later by              */
/* InitGeometry() and/or InitSubstrate()                            */
/********************************************************************/
void scuffSolver::AddSubstrateLayer(double zInterface, cdouble Epsilon, cdouble Mu)
{ 
  SubstrateInitialized=false;

  char Line[100];
  if (std::isinf(real(Epsilon)))
   snprintf(Line,100,"%e GROUNDPLANE\n",zInterface);
  else if (Mu==1.0)
   snprintf(Line,100,"%e CONST_EPS_%s\n",zInterface,z2s(Epsilon));
  else
   snprintf(Line,100,"%e CONST_EPS_%s_MU_%s\n",zInterface,z2s(Epsilon),z2s(Mu));
  SubstrateLayers.insert(std::pair<double, char *>(-zInterface, strdup(Line)));
}

void scuffSolver::AddGroundPlane(double zGP)
{ AddSubstrateLayer(zGP, HUGE_VAL); }

void scuffSolver::SetSubstratePermittivity(cdouble Epsilon)
{ AddSubstrateLayer(0.0, Epsilon); }

void scuffSolver::SetSubstrateThickness(double h)
{ AddGroundPlane(-h); }

void scuffSolver::SetSubstrateFile(const char *_SubstrateFile)
{ 
  SubstrateInitialized=false;
  SubstrateFile = strdup(_SubstrateFile);
}

void scuffSolver::InitSubstrate()
{
  if (SubstrateInitialized) return;
  SubstrateInitialized=true;
  OmegaSIE=OmegaPBFI=CACHE_DIRTY;

  if (!G) ErrExit("%s:%i: internal error");
  if (G->Substrate) free(G->Substrate);

  if (SubstrateFile)
   { G->Substrate = new LayeredSubstrate(SubstrateFile);
     free(SubstrateFile);
     SubstrateFile=0;
   }
  else if (SubstrateLayers.size()>0)
   { char *SubstrateDescription=0;
     for(std::map<double, char *>::iterator it = SubstrateLayers.begin(); it!=SubstrateLayers.end(); it++)
      { SubstrateDescription=vstrappend(SubstrateDescription,"%s\n",it->second);
        free(it->second);
      }
     G->Substrate = CreateLayeredSubstrate(SubstrateDescription);
     free(SubstrateDescription);
     SubstrateLayers.clear();
   }
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void scuffSolver::SetGeometryFile(const char *_scuffgeoFile)
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

void scuffSolver::AddLatticeVector(dVec L)
{ if (G) 
   { Warn("can't add lattice vectors after geometry has been initialized");
     return;
   }
  if ( L.size()<1 || L.size()>3) ErrExit("invalid lattice vector");
  if ( L.size()==3 && L[2]!=0.0 ) ErrExit("lattice vectors must have zero z-component");
  while( L.size()!=3 ) L.push_back(0.0); // pad with zeros to dimension 3
  LatticeVectors.push_back(L);
}

void scuffSolver::SetMediumMaterial(char *MediumMaterial)
{ Medium=strdup(MediumMaterial); }

void scuffSolver::SetMediumPermittivity(double Epsilon)
{ Medium=vstrdup("CONST_EPS_%g",Epsilon); }

void scuffSolver::AddObject(const char *MeshFile, const char *Material, const char *Label, const char *Transformation)
{
  char *ErrMsg=0;

  if (G) 
   ErrMsg=vstrdup("can't add metal traces after geometry has been initialized");
  else if (scuffgeoFile)
   ErrMsg=vstrdup("can't add metal traces to existing geometry file %s",scuffgeoFile);

  // check mesh
  if (!ErrMsg) 
   { RWGSurface *S = new RWGSurface(MeshFile);
     if (S->ErrMsg) ErrMsg = strdup(S->ErrMsg);
     delete S;
   }

  // check material 
  if (!ErrMsg && Material)
   { MatProp *MP = new MatProp(Material);
     if (MP->ErrMsg)
      ErrMsg=strdup(MP->ErrMsg);
     delete MP;
   }

  // check transformation
  if (!ErrMsg && Transformation)
   { GTransformation *GT = new GTransformation(Transformation, &ErrMsg);
     delete GT;
   }

  if (ErrMsg)
   { Warn("AddObject: %s (ignoring)",ErrMsg);
     free(ErrMsg);
     return;
   }

  MeshFiles.push_back(strdup(MeshFile));
  Materials.push_back(Material ? strdup(Material) : 0);
  MeshLabels.push_back(Label ? strdup(Label) : 0);
  MeshTransforms.push_back( Transformation ? strdup(Transformation) : 0);
}

void scuffSolver::AddMetalTraceMesh(const char *MeshFile, const char *Label, const char *Transformation)
{ AddObject(MeshFile, 0, Label, Transformation); }

void scuffSolver::SetGDSIIPortFile(const char *GDSIIPortFile, iVec Layers)
{
  if (G) 
   { Warn("can't add port file after geometry has been initialized");
     return;
   }
  if (portFile)
   { Warn("Overriding existing port file %s",portFile);
     free(portFile);
   }
  portFile=strdup(GDSIIPortFile);
  GDSIIPortLayers=iVec(Layers);
}

void scuffSolver::SetGDSIIPortFile(const char *GDSIIPortFile, int Layer)
{ SetGDSIIPortFile(GDSIIPortFile, iVec(1,Layer) ); }

void scuffSolver::SetPortFile(const char *_portFile)
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

void scuffSolver::AddPort(const dVec PVertices, const dVec MVertices)
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

void scuffSolver::AddPort(const dVec PVertices)
{ AddPort(PVertices, dVec()); }

void scuffSolver::AddPortTerminal(char PM, const dVec Vertices)
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
void scuffSolver::InitGeometry()
{ 
  if (G) return;

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
     
     if (Medium)
      { fprintf(f,"MEDIUM\n MATERIAL %s\nENDMEDIUM\n\n",Medium);
        free(Medium);
      }

     if (LatticeVectors.size()>0)
      { fprintf(f,"LATTICE\n");
        for(size_t n=0; n<LatticeVectors.size(); n++)
         fprintf(f," %e %e %e \n",LatticeVectors[n][0],LatticeVectors[n][1],LatticeVectors[n][2]);
        fprintf(f,"ENDLATTICE\n");
      }
     LatticeVectors.clear();

     if (MeshFiles.size()==0) ErrExit("no metal traces specified");
     if (MeshFiles.size() != Materials.size()) ErrExit("%s:%i: internal error");
     if (MeshFiles.size() != MeshTransforms.size()) ErrExit("%s:%i: internal error");
     if (MeshFiles.size() != MeshLabels.size()) ErrExit("%s:%i: internal error");
     for(size_t n=0; n<MeshFiles.size(); n++)
      { char Buffer[100], *Label = MeshLabels[n];
        if (!Label)
         { Label=Buffer; snprintf(Label,100,"%s_%lu",GetFileBase(MeshFiles[n]),n); }
        fprintf(f,"OBJECT %s\n MESHFILE %s\n",Label,MeshFiles[n]);
        if (Materials[n]) fprintf(f," MATERIAL %s\n",Materials[n]);
        if (MeshTransforms[n]) fprintf(f," %s\n",MeshTransforms[n]);
        fprintf(f,"ENDOBJECT\n");
        free(MeshFiles[n]);
        if (Materials[n]) free(Materials[n]);
        if (MeshLabels[n]) free(MeshLabels[n]);
        if (MeshTransforms[n]) free(MeshTransforms[n]);
      }
     MeshFiles.clear();
     Materials.clear();
     MeshLabels.clear();
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
            if (nv==0 || std::isinf(V)) fprintf(f,"\n   %s ",(Pol==_PLUS ? "POSITIVE" : "NEGATIVE"));
            if (std::isinf(V)) continue;
            fprintf(f,"%g ",V);
          }
        fprintf(f,"\nENDPORT\n\n");
      }
     fclose(f);
     PortTerminalVertices.clear();
   }
  if (!portFile) ErrExit("no ports specified");
  PortList = (GDSIIPortLayers.size()>0 ? ReadGDSIIPorts(G, portFile, GDSIIPortLayers)
                                       : ParsePortFile(G, portFile)
             );
  NumPorts = PortList->Ports.size();
  //if (NumPorts==0) Warn("no ports found");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  InitSubstrate();
 
  if(ownsGeoFile && !CheckEnv("PYSCUFF_RETAIN_GEOFILE")) unlink(scuffgeoFile);
  if(ownsPortFile) unlink(portFile);
}

/***************************************************************/
/* Identify equivalent surface-surface pairs to allow reuse of */
/*  off-diagonal BEM matrix blocks.                            */
/*                                                             */
/* (SAlpha, SBeta) is equivalent to (SA, SB) if the following  */
/*  three conditions are all satified.                         */
/*                                                             */
/*  1. SAlpha and SA are both the result of geometrical        */
/*     transformations applied to a common surface S1:         */
/*       SAlpha = TAlpha(S1)                                   */
/*       SA     = TA    (S1)                                   */
/*    where TAlpha, TA are geometrical transforms (which may   */
/*    be the identity transform)                               */
/*                                                             */
/*  2. SBeta and SB are both the result of geometrical         */
/*     transformations applied to a common surface S2:         */
/*       SBeta  = TBeta(S2)                                    */
/*       SB     = TB    S2)                                    */
/*                                                             */
/*  3. We have TAlpha^{-1} * TBeta = TA^{-1} TB                */
/*      or                                                     */
/*             TAlpha^{-1} * TBeta * TB^{-1} * TA = identity   */
/*                                                             */
/* If (SAlpha, SBeta) <==> (SA, SB) then the return value is   */
/* the index of the parent off-diagonal block (SA,B).          */
/*                                                             */
/*  M[SAlpha, SBeta] = index of off-diagonal block (SA,SB)     */
/*                                                             */
/* (with NS==number of surfaces)                               */
/***************************************************************/
int OffDiagonalBlockIndex(int NS, int ns, int nsp)
{ return ns*NS - ns*(ns+1)/2 + nsp - ns - 1; }

int FindEquivalentSurfacePair(RWGGeometry *G, int nsAlpha, int nsBeta,
                              bool *pFlipped=0, int *pnsA=0, int *pnsB=0)
{
  if (CheckEnv("SCUFF_IGNORE_EQUIVALENT_SURFACES")) return -1;

  int *Mate = G->Mate;

  int nsAlphaMate = Mate[nsAlpha], nsBetaMate  = Mate[nsBeta];
  if (nsAlphaMate==-1 && nsBetaMate==-1)
   return -1;

  if (nsAlphaMate==-1) nsAlphaMate = nsAlpha;
  if (nsBetaMate==-1)  nsBetaMate  = nsBeta;

  for(int nsA=0; nsA<=nsAlpha; nsA++)
   for(int nsB=1; nsB<=nsBeta; nsB++)
    {
      if (nsA==nsAlpha && nsB==nsBeta) continue;

      int nsAMate = (Mate[nsA]==-1) ? nsA: Mate[nsA];
      int nsBMate = (Mate[nsB]==-1) ? nsB: Mate[nsB];

      if (nsAlphaMate!=nsAMate || nsBetaMate!=nsBMate)
       continue;

      RWGSurface *SAlpha = G->Surfaces[nsAlpha];
      RWGSurface *SBeta  = G->Surfaces[nsBeta];
      RWGSurface *SA     = G->Surfaces[nsA];
      RWGSurface *SB     = G->Surfaces[nsB];
          
      GTransformation T; // identity transformation
    
      if (SAlpha->GT)   T = T - *(SAlpha->GT);
      if (SAlpha->OTGT) T = T - *(SAlpha->OTGT);
    
      if (SBeta->OTGT)  T = T + *(SBeta->OTGT);
      if (SBeta->GT)    T = T + *(SBeta->GT);
    
      if (SB->GT)       T = T - *(SB->GT);
      if (SB->OTGT)     T = T - *(SB->OTGT);
    
      if (SA->OTGT)     T = T + *(SA->OTGT);
      if (SA->GT)       T = T + *(SA->GT);
    
      double Lengthscale = fmin( VecDistance(SA->RMax, SA->RMin),
                                 VecDistance(SB->RMax, SB->RMin)
                              );
      if ( T.IsIdentity(Lengthscale) )
       { bool Flipped=false;
         if (nsB<nsA)
          { int temp = nsA; nsA=nsB; nsB=temp;
             SA=G->Surfaces[nsA];
             SB=G->Surfaces[nsB];
             Flipped=true;
          }
         if (G->LogLevel>=SCUFF_VERBOSE2)
          Log(" %10s<-->%-10s === %10s<-->%-10s",SAlpha->Label,SBeta->Label,SA->Label,SB->Label);
         if (pFlipped) *pFlipped=Flipped;
         if (pnsA)  *pnsA=nsA;
         if (pnsB)  *pnsB=nsB;
         return OffDiagonalBlockIndex(G->NumSurfaces, nsA, nsB);
       }
    }
  return -1;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void scuffSolver::EnableSystemBlockCache()
{ 
  if (G==0)
   { Warn("EnableSysteBlockCache() must be called after InitGeometry (skipping)");
     return;
   }

  if (TBlocks!=0) return; // already enabled

  int NS = G->NumSurfaces, NU=(NS-1)*NS/2;
  if (NS==1)
   { Log("Block caching not meaningful for single-surface geometries (skipping)");
     return;
   }
  else 
   Log("Initializing system block cache (%i/%i diagonal/off-diagonal blocks)",NS,NU);

  TBlocks = new HMatrix *[NS];
  for(int ns=0; ns<NS; ns++)
   { int NBF = G->Surfaces[ns]->NumBFs, Mate=G->Mate[ns];
     TBlocks[ns] = (Mate!=-1) ? TBlocks[Mate] : new HMatrix(NBF, NBF, LHM_COMPLEX);
   }

  UBlocks = new HMatrix *[NU];
  memset(UBlocks, 0, NU*sizeof(HMatrix *)); // allocated lazily

  OmegaSIE = CACHE_DIRTY;
}

/***************************************************************/
/* note: don't call this routine directly; call AssembleSystemMatrix() instead */
/***************************************************************/
void scuffSolver::UpdateSystemMatrix(cdouble Omega)
{ 
  if (TBlocks==0)
   { Warn("EnableSystemBlockCache() must be called before UpdateSystemMatrix (skipping)");
     return;
   }

  // TBlocks recomputed only if frequency has changed
  if (Omega!=OmegaSIE)
   for(int ns=0; ns<G->NumSurfaces; ns++)
    if (G->Mate[ns]==-1)
     { if(NumPorts==0)
        G->AssembleBEMMatrixBlock(ns, ns, Omega, 0, TBlocks[ns]);
       else 
        AssembleMOIMatrixBlock(G, ns, ns, Omega, TBlocks[ns]);
     }

  // UBlocks recomputed if frequency has changed or surfaces were moved
  for(int nsa=0, nu=0; nsa<G->NumSurfaces; nsa++)
   for(int nsb=nsa+1; nsb<G->NumSurfaces; nsb++, nu++)
    { if (FindEquivalentSurfacePair(G, nsa, nsb)!=-1)
       continue;
      if (Omega!=OmegaSIE || G->SurfaceMoved[nsa] || G->SurfaceMoved[nsb])
       { if (UBlocks[nu]==0) 
          UBlocks[nu]=new HMatrix(G->Surfaces[nsa]->NumBFs, G->Surfaces[nsb]->NumBFs, LHM_COMPLEX);
         if(NumPorts==0)
          G->AssembleBEMMatrixBlock(nsa, nsb, Omega, 0, UBlocks[nu]);
         else
          AssembleMOIMatrixBlock(G, nsa, nsb, Omega, UBlocks[nu]);
       }
    }

  // stamp blocks into M matrix
  for(int nsa=0, nu=0; nsa<G->NumSurfaces; nsa++)
   { int OffsetA = G->BFIndexOffset[nsa];
     M->InsertBlock(TBlocks[nsa], OffsetA, OffsetA);
     for(int nsb=nsa+1; nsb<G->NumSurfaces; nsb++, nu++)
      { 
        int OffsetB = G->BFIndexOffset[nsb];
        bool Flipped=false;
        int nuEquivalent=FindEquivalentSurfacePair(G, nsa, nsb, &Flipped);
        HMatrix *UBlock = UBlocks[(nuEquivalent==-1 ? nu : nuEquivalent)];
        int RowOffset = Flipped ? OffsetB : OffsetA;
        int ColOffset = Flipped ? OffsetA : OffsetB;
        M->InsertBlock(UBlock,RowOffset,ColOffset);
        M->InsertBlockTranspose(UBlock,ColOffset,RowOffset);
      }
   }
  OmegaSIE=Omega;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void scuffSolver::AssembleSystemMatrix(double Omega)
{ 
  if (!G) InitGeometry();
  if (!SubstrateInitialized) InitSubstrate();
  if (M==0) M=G->AllocateBEMMatrix();
  G->UpdateCachedEpsMuValues(Omega);

  Log("Assembling BEM matrix at Omega=%g (f=%g GHz)...",Omega,Omega/FREQ2OMEGA);

  // update interpolator
  if (G->Substrate)
   { double RhoMinMax[2];
     GetRhoMinMax(G, RhoMinMax);
     bool PPIsOnly = true;
     bool Subtract = true;
     bool RetainSingularTerms = false;
     bool Verbose = (G->LogLevel == SCUFF_VERBOSE2);
     G->Substrate->InitScalarGFInterpolator(Omega, RhoMinMax[0], RhoMinMax[1], 0.0, 0.0,
                                            PPIsOnly, Subtract, RetainSingularTerms, Verbose);
   }

  if (TBlocks!=0)
   UpdateSystemMatrix(Omega);
  else if (NumPorts!=0)
   AssembleMOIMatrix(G, Omega, M);
  else 
   G->AssembleBEMMatrix(Omega, M);

  Log("Factorizing...");
  M->LUFactorize();

  OmegaSIE=Omega;
}

void scuffSolver::DoSolve(IncField *IF, cdouble *PortCurrents)
{
  if (M==0 || G->StoredOmega!=OmegaSIE)
   { Warn("scuffSolver: AssembleSystemMatrix() must be called before Solve()");
     return;
   }
  cdouble Omega = G->StoredOmega;

  /*--------------------------------------------------------------*/
  /*- lazy allocation of internal variables ----------------------*/
  /*--------------------------------------------------------------*/
  int NBF = G->TotalBFs;
  if (KN==0) KN = new HVector(NBF,                LHM_COMPLEX);

  // stuff needed for RF calculations
  if (PortCurrents)
   { if (PBFIMatrix==0)   PBFIMatrix   = new HMatrix(NBF,      NumPorts, LHM_COMPLEX);
     if (PPIMatrix==0)    PPIMatrix    = new HMatrix(NumPorts, NumPorts, LHM_COMPLEX);
     AssemblePortBFInteractionMatrix(Omega);
    // form RHS vector as weighted linear combination of contributions from each port
     KN->Zero();
     for(int np=0; np<NumPorts; np++)
      for(int nbf=0; nbf<NBF; nbf++)
       KN->AddEntry(nbf, PortCurrents[np]*PBFIMatrix->GetEntry(nbf,np));
     if (CachedPortCurrents==0) CachedPortCurrents = new cdouble[NumPorts];
     memcpy(CachedPortCurrents, PortCurrents, NumPorts*sizeof(cdouble));
   }
  else if (IF)
   { G->AssembleRHSVector(Omega, IF, KN);
     if (CachedIF){
       delete CachedIF;
     }
     if (dynamic_cast<PlaneWave*>(IF) != nullptr)
     {
      PlaneWave *PW = (PlaneWave *) IF;
      CachedIF=new PlaneWave(*PW);
     } 
     if (dynamic_cast<PointSource*>(IF) != nullptr)
     {
      PointSource *PS = (PointSource *) IF;
      CachedIF=new PointSource(*PS);
     }
   }

  // solve the system
  Log("LU-solving...");
  M->LUSolve(KN);
}

void scuffSolver::Solve(IncField *IF)
{ DoSolve(IF, 0); }

void scuffSolver::Solve(cdouble *PortCurrents)
{ DoSolve(0,PortCurrents); }

void scuffSolver::Solve(int WhichPort, cdouble PortCurrent)
{ if (NumPorts==0) ErrExit("no ports defined");
  cdouble *PortCurrents = new cdouble[NumPorts];
  PortCurrents[WhichPort]=PortCurrent;
  DoSolve(0,PortCurrents);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void scuffSolver::PlotGeometry(const char *PPFormat, ...)
{
  if (!G) InitGeometry();

  /***************************************************************/
  /* open the file ***********************************************/
  /***************************************************************/
  char PPFileName[1000];
  if (!PPFormat)
   snprintf(PPFileName,1000,"%s.pp",SolverName);
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
  G->WritePPMesh(PPFileName,SolverName);
  
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

     //
     for(int nl=0; nl<S->NumLayers; nl++)
      { fprintf(f,"View \"Layer %i (%s)\" {\n",nl,S->MPLayer[nl]->Name);
        double zAbove = S->zInterface[0];
        double zBelow = (nl==S->NumLayers-1) ? S->zGP : S->zInterface[nl+1];
        // top surface
        fprintf(f,"SQ(%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e) {%i,%i,%i,%i};\n",
                   RMin[0], RMin[1], zAbove, RMax[0], RMin[1], zAbove, RMax[0], RMax[1], zAbove, RMin[0], RMax[1], zAbove, nl, nl, nl, nl);
        // side surfaces
        fprintf(f,"SQ(%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e) {%i,%i,%i,%i};\n",
                   RMin[0], RMin[1], zAbove, RMin[0], RMax[1], zAbove, RMin[0], RMax[1], zBelow, RMin[0], RMin[1], zBelow, nl, nl, nl, nl);
        fprintf(f,"SQ(%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e) {%i,%i,%i,%i};\n",
                   RMax[0], RMin[1], zAbove, RMax[0], RMax[1], zAbove, RMax[0], RMax[1], zBelow, RMax[0], RMin[1], zBelow, nl, nl, nl, nl);
        fprintf(f,"SQ(%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e) {%i,%i,%i,%i};\n",
                   RMin[0], RMin[1], zAbove, RMax[0], RMin[1], zAbove, RMax[0], RMin[1], zBelow, RMin[0], RMin[1], zBelow, nl, nl, nl, nl);
        fprintf(f,"SQ(%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e) {%i,%i,%i,%i};\n",
                   RMin[0], RMax[1], zAbove, RMax[0], RMax[1], zAbove, RMax[0], RMax[1], zBelow, RMin[0], RMax[1], zBelow, nl, nl, nl, nl);
        // bottom surface
        if ( zBelow!=S->zGP )
        fprintf(f,"SQ(%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e) {%i,%i,%i,%i};\n",
                   RMin[0], RMin[1], zBelow, RMax[0], RMin[1], zBelow, RMax[0], RMax[1], zBelow, RMin[0], RMax[1], zBelow, nl, nl, nl, nl);
        fprintf(f,"};\n");
      }

     if ( !std::isinf(S->zGP) )
      { 
        int nl = S->NumLayers;
        fprintf(f,"View \"Ground plane\" {\n");
        fprintf(f,"SQ(%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e) {%i,%i,%i,%i};\n",
                   RMin[0], RMin[1], S->zGP, RMax[0], RMin[1], S->zGP, RMax[0], RMax[1], S->zGP, RMin[0], RMax[1], S->zGP, nl, nl, nl, nl);
        fprintf(f,"};\n");
      }
   }
  fclose(f);

  fprintf(stdout,"RF geometry visualization written to GMSH file %s.\n",PPFileName);
}

void scuffSolver::PlotGeometry()
 { PlotGeometry(0); }

void scuffSolver::Transform(const char *SurfaceLabel, const char *Transformation)
{ 
  if (!G) InitGeometry();
  RWGSurface *S=G->GetSurfaceByLabel(SurfaceLabel);
  if (!S) 
   { Warn("Transform: no surface %s found in geometry (skipping)",SurfaceLabel);
     return;
   }
  S->Transform(Transformation);
}

void scuffSolver::UnTransform()
{ if (!G) return;
  G->UnTransform();
}

} // namespace scuff
