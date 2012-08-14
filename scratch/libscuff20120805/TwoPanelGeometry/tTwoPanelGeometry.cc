#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <libhrutil.h>
#include "libscuff.h"
#include "libhmat.h"

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[]) 
{  
  EnableAllCPUs();
  InstallHRSignalHandler();
  
  cdouble E0[3]  = {1.0, 0.0, 0.0};
  double nHat[3] = {0.0, 0.0, 1.0};
  PlaneWave PW(E0, nHat);
  cdouble Omega=0.1;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  SetLogFileName("tTwoPanelGeometry.log");

  RWGGeometry::AssignBasisFunctionsToExteriorEdges=false;
  RWGGeometry *G = new RWGGeometry("Full.scuffgeo", SCUFF_VERBOSELOGGING);
  HMatrix *M = G->AllocateBEMMatrix();
  HVector *KN = G->AllocateRHSVector();

  G->AssembleBEMMatrix(Omega, M);
  M->ExportToText("MFull.txt");
  M->LUFactorize();

  G->AssembleRHSVector(Omega, &PW, KN);
  KN->ExportToText("RHSFull.txt");
  M->LUSolve(KN);
  KN->ExportToText("KNFull.txt");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGGeometry::AssignBasisFunctionsToExteriorEdges=true;
  RWGGeometry::IncludeLineCharges=false;
  G = new RWGGeometry("Half.scuffgeo", SCUFF_VERBOSELOGGING);
G->Surfaces[0]->Edges[0]=G->Surfaces[0]->Edges[2];
G->Surfaces[0]->NumEdges=1;
G->Surfaces[1]->NumEdges=1;
G->TotalBFs=2;
G->BFIndexOffset[1]=1;
  M = G->AllocateBEMMatrix();
  KN = G->AllocateRHSVector();

  G->AssembleBEMMatrix(Omega, M);
  M->ExportToText("MHalfWO.txt");
  M->LUFactorize();

  G->AssembleRHSVector(Omega, &PW, KN);
  KN->ExportToText("RHSHalfWO.txt");
  M->LUSolve(KN);
  KN->ExportToText("KNHalfWO.txt");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGGeometry::AssignBasisFunctionsToExteriorEdges=true;
  RWGGeometry::IncludeLineCharges=true;
  G = new RWGGeometry("Half.scuffgeo", SCUFF_VERBOSELOGGING);
G->Surfaces[0]->Edges[0]=G->Surfaces[0]->Edges[2];
G->Surfaces[0]->NumEdges=1;
G->Surfaces[1]->NumEdges=1;
G->TotalBFs=2;
G->BFIndexOffset[1]=1;
  M = G->AllocateBEMMatrix();
  KN = G->AllocateRHSVector();

  G->AssembleBEMMatrix(Omega, M);
  M->ExportToText("MHalfW.txt");
  M->LUFactorize();

  G->AssembleRHSVector(Omega, &PW, KN);
  KN->ExportToText("RHSHalfW.txt");
  M->LUSolve(KN);
  KN->ExportToText("KNHalfW.txt");
}
