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
   { {"LBV1",      PA_DOUBLE,  2, 1, (void *)LBV[0],       0,        "lattice basis vector 1"},
     {"LBV2",      PA_DOUBLE,  2, 1, (void *)LBV[1],       0,        "lattice basis vector 2"},
     {"P",         PA_DOUBLE,  2, 1, (void *)P,            0,        "bloch wavevector"},
     {"Omega",     PA_CDOUBLE, 1, 1, (void *)&Omega,       0,        "angular frequency"},
     {"Geometry",  PA_STRING,  1, 1, (void *)&GeoFileName, 0,        ".scuffgeo file"},
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
  G->TotalBFs=G->TotalPanels=0;
  for(int no=0; no<G->NumObjects; no++)
   { RWGObject *O = G->Objects[no];
     AddStraddlers(O, LBV);
     G->TotalBFs    += O->NumBFs;
     G->TotalPanels += O->NumPanels;
     if ( no+1 < G->NumObjects )
      { G->BFIndexOffset[no+1]=G->BFIndexOffset[no] + O->NumBFs;
        G->PanelIndexOffset[no+1]=G->PanelIndexOffset[no] + O->NumPanels;
      };
   };

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

}
