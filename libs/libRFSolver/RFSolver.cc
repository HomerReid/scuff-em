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
  /*--------------------------------------------------------------*/
  /* try to read in geometry and port list -----------------------*/
  /*--------------------------------------------------------------*/
  RWGGeometry::UseHRWGFunctions=false;
  G = new RWGGeometry(scuffgeoFileName);
  PortList = ParsePortFile(G, portFileName);
  NumPorts = PortList->Ports.size();
  FileBase = strdup(GetFileBase(scuffgeoFileName));

  InitSolver();

}

/***************************************************************/
/* RF solver class constructor #2: construct from a .GDSII     */
/* file                                                        */
/***************************************************************/
RFSolver::RFSolver(const char *GDSIIFileName)
{ 
  FileBase = strdup(GetFileBase(GDSIIFileName));
  NumPorts = 0; // PortList->Ports.size();

  InitSolver();
}

void RFSolver::InitSolver()
{
  /*--------------------------------------------------------------------------*/
  /* initialize internal variables needed to perform RF calculations (which   */
  /* are allocated lazily)                                                    */
  /*--------------------------------------------------------------------------*/
  Omega        = -1.0;
  M            = 0;
  MClean       = false;
  PBFIMatrix   = 0; 
  PPIMatrix    = 0;
  PBFIClean    = false;
  PortCurrents = 0;
  KN           = 0;

  DisableSystemBlockCache = false;
  OmegaCache = HUGE_VAL;
  TBlocks = 0;
  UBlocks = 0;

  RetainContributions = CONTRIBUTION_ALL;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
RFSolver::~RFSolver()
{
  if (G)            delete G;
  if (PortList)     delete PortList;
  if (FileBase)     free(FileBase);
  if (M)            delete M;
  if (PBFIMatrix)   delete PBFIMatrix;
  if (PortCurrents) delete PortCurrents;
  if (KN)           delete KN;

  int NS = G->NumSurfaces, NU=(NS-1)*NS/2;
  if (TBlocks)
   { for(int ns=0; ns<NS; ns++)
      delete TBlocks[ns];
     delete TBlocks;
   }
  if (UBlocks)
   { for(int nu=0; nu<NU; nu++)
      delete UBlocks[nu];
     delete UBlocks;
   }
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RFSolver::SetSubstrate(const char *SubstrateFile)
{ G->Substrate = new LayeredSubstrate(SubstrateFile); }

void RFSolver::SetSubstrate(const char *EpsStr, double h)
{
  char SubstrateDefinition[1000];
  if (EpsStr==0 && h==0.0) // no substrate
    return;
  else if (EpsStr==0) // ground plane onl-
   snprintf(SubstrateDefinition,1000,"0.0 GROUNDPLANE %e\n",-h);
  else if (h==0.0) // dielectric only (no ground plane)
   snprintf(SubstrateDefinition,1000,"0.0 CONST_EPS_%s\n",EpsStr);
  else // grounded dielectric 
   snprintf(SubstrateDefinition,1000,"0.0 CONST_EPS_%s\n%e GROUNDPLANE\n",EpsStr,-h);
  G->Substrate=CreateLayeredSubstrate(SubstrateDefinition);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RFSolver::UpdateSystemMatrix()
{ 
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
     AssembleMOIMatrixBlock(G, ns, ns, Omega, TBlocks[ns]);

  // UBlocks recomputed if frequency has changed or surfaces were moved
  for(int nsa=0, nb=0; nsa<G->NumSurfaces; nsa++)
   for(int nsb=nsa+1; nsb<G->NumSurfaces; nsb++, nb++)
    if (Omega!=OmegaCache || G->SurfaceMoved[nsa] || G->SurfaceMoved[nsb])
     AssembleMOIMatrixBlock(G, nsa, nsb, Omega, UBlocks[nb]);

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
  if (M==0) M=G->AllocateBEMMatrix();
  Omega = Freq * FREQ2OMEGA;
  Log("Assembling BEM matrix at f=%g GHz...",Freq);
  if (DisableSystemBlockCache)
   AssembleMOIMatrix(G, Omega, M);
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

} // namespace scuff
