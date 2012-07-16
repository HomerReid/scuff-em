/*
 *
 */
#include <stdio.h>
#include <stdlib.h>

#include <libhrutil.h>
#include <libscuff.h>
#include <libIncField.h>
#include "PBCGeometry.h"

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  EnableAllCPUs();
  InstallHRSignalHandler();
  
  double LBV1={1.0, 0.0};
  double LBV2={0.0, 1.0};
  double *LBV[2]={LBV1, LBV2};

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  char *GeoFileName="PECSurface_52.scuffgeo";
  RWGGeometry *G=new RWGGeometry(GeoFileName);
  G->SetLogLevel(SCUFF_VERBOSELOGGING);
  PBCGeometry *PG=new PBCGeometry(G, LBV);

  if (Cache) PreloadCache(Cache);
  HMatrix *M = G->AllocateBEMMatrix();
  PG->AssembleBEMMatrix(Omega, BlochP, M);

  void *pCC=HMatrix::CreateC2MLContext("Hello");
  M->ExportToMATLAB(pCC,"Periodic");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGGeometry *G=new RWGGeometry("PECSurface_52_Lattice.scuffgeo");
  PBCGeometry *PG=new PBCGeometry(G, LBV);

  G->Objects[1]->Transform("DISPLACED 1 0 0");

  G->AssembleBEMMatrix(Omega, M);
  M->ExportToMATLAB(pCC,"Direct");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  HMatrix::CloseC2MLContext(pCC);


}
