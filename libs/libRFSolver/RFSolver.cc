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

  RetainContributions = CONTRIBUTION_ALL;

}

/***************************************************************/
/* RF solver class constructor #2: construct from a .GDSII     */
/* file                                                        */
/***************************************************************/
RFSolver::RFSolver(const char *GDSIIFileName)
{ 
  FileBase = strdup(GetFileBase(GDSIIFileName));
  NumPorts = 0; // PortList->Ports.size();

  Omega        = -1.0;
  M            = 0;
  MClean       = 0;
  PBFIMatrix   = 0;
  PPIMatrix    = 0;
  PBFIClean    = false;
  PortCurrents = 0;
  KN           = 0;

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
  if (PortCurrents) free(PortCurrents);
  if (KN)           delete KN;
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

void RFSolver::AssembleSystemMatrix(double Freq)
{ 
  Omega = Freq * FREQ2OMEGA;
  if (M==0) M=G->AllocateBEMMatrix();
  Log("Assembling BEM matrix at f=%g GHz...",Freq);
  AssembleMOIMatrix(G, Omega, M); 
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
  if (PortCurrents==0) PortCurrents = (cdouble *)mallocEC(NumPorts*sizeof(cdouble));
  if (PBFIMatrix==0)   PBFIMatrix   = new HMatrix(NBF,      NumPorts, LHM_COMPLEX);
  if (PPIMatrix==0)    PPIMatrix    = new HMatrix(NumPorts, NumPorts, LHM_COMPLEX);
  if (KN==0)           KN           = new HVector(NBF,                LHM_COMPLEX);
  
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
  if (PortCurrents==0)
   PortCurrents = (cdouble *)mallocEC(NumPorts*sizeof(cdouble));
  memset(PortCurrents, 0, NumPorts*sizeof(cdouble));
  PortCurrents[WhichPort]=PortCurrent;
  Solve(0);
}

} // namespace scuff
