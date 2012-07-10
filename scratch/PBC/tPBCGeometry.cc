/*
 *
 */
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
  /*--------------------------------------------------------------*/
  /*- process options  -------------------------------------------*/
  /*--------------------------------------------------------------*/
  char *Geometry;
  double LBV[2][2]={ {1.0, 0.0}, {0.0, 1.0} };
  double P[2]={0.0, 0.0};
  cdouble Omega;
  /* name, type, #args, max_instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"Geometry",  PA_STRING,  1, 1, (void *)&GeoFileName, 0,        ".scuffgeo file"},
     {"LBV1",      PA_DOUBLE,  2, 1, (void *)LBV[0],       0,        "lattice basis vector 1"},
     {"LBV2",      PA_DOUBLE,  2, 1, (void *)LBV[1],       0,        "lattice basis vector 2"},
     {"P",         PA_DOUBLE,  2, 1, (void *)P,            0,        "bloch wavevector"},
     {"Omega",     PA_CDOUBLE, 1, 1, (void *)&Omega,       0,        "angular frequency"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (Geometry==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");
  SetLogFileName("tPBC.log");
  RWGGeometry::SetLogLevel(SCUFF_VERBOSELOGGING);
  
  /*--------------------------------------------------------------*/
  /*- add straddlers to all objects ------------------------------*/
  /*--------------------------------------------------------------*/

  /*--------------------------------------------------------------*/
  /*- initialize the PBC accelerator -----------------------------*/
  /*--------------------------------------------------------------*/
  PBCAccelerator *PBCA=CreatePBCAccelerator(G);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  HMatrix *M = G->AllocateBEMMatrix();
  AssembleBEMMatrix_PBC(G, Omega, LBV, P, PBCA, M);
  M->LUFactorize();

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble E0[3]  = {1.0, 0.0, 0.0};
  double nHat[3] = {0.0, 0.0, 1.0};
  PlaneWave PW(E0, nHat);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  HMatrix *KN = G->AllocateRHSVector();
  AssembleRHSVector(Omega, &PW, KN);

  M->LUSolve(KN);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  GetFields_PBC(G, Omega, LBV, P, PBCA, M);

}
