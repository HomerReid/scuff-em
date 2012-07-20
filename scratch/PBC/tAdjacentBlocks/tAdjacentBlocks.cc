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
  const char *GeoFileName="E10Sphere_48.scuffgeo";
  SetLogFileName("tAdjacentBlocks.log");
  RWGGeometry *G=new RWGGeometry(GeoFileName);
  G->SetLogLevel(SCUFF_VERBOSELOGGING);
  PBCGeometry *PG=new PBCGeometry(G, LBV);

  cdouble Omega=0.1;
  double BlochP[2]={0.0, 0.0};

  PreloadCache("../Sphere_48.cache");
  HMatrix *M = G->AllocateBEMMatrix();
  PG->AssembleBEMMatrix(Omega, BlochP, M);
  StoreCache("../Sphere_48.cache");

  void *pCC=HMatrix::OpenMATLABContext("Hello");
  M->ExportToMATLAB(pCC,"Periodic");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  G=new RWGGeometry("E10Sphere_48_Lattice9.scuffgeo");
  PG=new PBCGeometry(G, LBV);
  M=G->AllocateBEMMatrix();
  G->Objects[1]->Transform("DISPLACED +1 +1 0");
  G->Objects[2]->Transform("DISPLACED -1 -1 0");
  G->Objects[3]->Transform("DISPLACED +1 -1 0");
  G->Objects[4]->Transform("DISPLACED -1 +1 0");
  G->Objects[5]->Transform("DISPLACED +1  0 0");
  G->Objects[6]->Transform("DISPLACED -1  0 0");
  G->Objects[7]->Transform("DISPLACED  0 +1 0");
  G->Objects[8]->Transform("DISPLACED  0 -1 0");

  G->AssembleBEMMatrix(Omega, M);
  M->ExportToMATLAB(pCC,"Direct");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  HMatrix::CloseMATLABContext(pCC);


}
