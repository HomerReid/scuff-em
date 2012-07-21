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
  G=new RWGGeometry("E10Sphere_48_Lattice17.scuffgeo");
#if 0
  G->Objects[0]->ContainingObject=0;
  G->Objects[1]->ContainingObject=0;
  G->Objects[2]->ContainingObject=0;
  G->Objects[3]->ContainingObject=0;
  G->Objects[4]->ContainingObject=0;
  G->Objects[5]->ContainingObject=0;
  G->Objects[6]->ContainingObject=0;
  G->Objects[7]->ContainingObject=0;
  G->Objects[8]->ContainingObject=0;
  G->Objects[9]->ContainingObject=0;
  G->Objects[10]->ContainingObject=0;
  G->Objects[11]->ContainingObject=0;
  G->Objects[12]->ContainingObject=0;
  G->Objects[13]->ContainingObject=0;
  G->Objects[14]->ContainingObject=0;
  G->Objects[15]->ContainingObject=0;
  G->Objects[16]->ContainingObject=0;
#endif
 // PG=new PBCGeometry(G, LBV);
  M=G->AllocateBEMMatrix();
#if 0
  G->Objects[ 1]->Transform("DISPLACED -2 +2 0");
  G->Objects[ 2]->Transform("DISPLACED -1 +2 0");
  G->Objects[ 3]->Transform("DISPLACED +0 +2 0");
  G->Objects[ 4]->Transform("DISPLACED +1 +2 0");
  G->Objects[ 5]->Transform("DISPLACED +2 +2 0");
  G->Objects[ 6]->Transform("DISPLACED -2 -2 0");
  G->Objects[ 7]->Transform("DISPLACED -1 -2 0");
  G->Objects[ 8]->Transform("DISPLACED +0 -2 0");
  G->Objects[ 9]->Transform("DISPLACED +1 -2 0");
  G->Objects[10]->Transform("DISPLACED +2 -2 0");
  G->Objects[11]->Transform("DISPLACED -2 -1 0");
  G->Objects[12]->Transform("DISPLACED -2  0 0");
  G->Objects[13]->Transform("DISPLACED -2 +1 0");
  G->Objects[14]->Transform("DISPLACED -2 -1 0");
  G->Objects[15]->Transform("DISPLACED -2  0 0");
  G->Objects[16]->Transform("DISPLACED -2 +1 0");
#endif
  G->AssembleBEMMatrix(Omega, M);
  M->ExportToMATLAB(pCC,"Direct");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  HMatrix::CloseMATLABContext(pCC);


}
