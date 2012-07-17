/*
 *
 */
#include <stdio.h>
#include <stdlib.h>

#include <libhrutil.h>
#include <libscuff.h>
#include <libIncField.h>
#include "../PBCGeometry.h"

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  EnableAllCPUs();
  InstallHRSignalHandler();
  
  double LBV1[2]={1.0, 0.0};
  double LBV2[2]={0.0, 1.0};
  double *LBV[2]={LBV1, LBV2};

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  const char *GeoFileName="PECSurface_52.scuffgeo";
  RWGGeometry *G=new RWGGeometry(GeoFileName);
  G->SetLogLevel(SCUFF_VERBOSELOGGING);
  PBCGeometry *PG=new PBCGeometry(G, LBV);

  cdouble Omega=0.1;
  double BlochP[2]={0.0, 0.0};

  PreloadCache("../PECSurface_52.cache");
  HMatrix *M = G->AllocateBEMMatrix();
  PG->AssembleBEMMatrix(Omega, BlochP, M);

  void *pCC=HMatrix::OpenMATLABContext("Hello");
  M->ExportToMATLAB(pCC,"Periodic");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  G=new RWGGeometry("PECSurface_52_Lattice.scuffgeo");
  PG=new PBCGeometry(G, LBV);
  M=G->AllocateBEMMatrix();
  G->Objects[1]->Transform("DISPLACED  1 0 0");
  G->Objects[2]->Transform("DISPLACED -1 0 0");

  G->AssembleBEMMatrix(Omega, M);
  M->ExportToMATLAB(pCC,"Direct");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  HMatrix::CloseMATLABContext(pCC);


}
