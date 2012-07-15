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

  /*--------------------------------------------------------------*/
  /*- process options  -------------------------------------------*/
  /*--------------------------------------------------------------*/
  char *GeoFileName=0;
  double LBV1[2]={1.0, 0.0}; 
  double LBV2[2]={0.0, 1.0}; 
  double *LBV[2]={ LBV1, LBV2 };
  double BlochP[2]={0.0, 0.0};
  cdouble Omega;
  char *Cache=0;
  char *EPFile=0;
  double ZValue=1.0;
  /* name, type, #args, max_instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"Geometry",  PA_STRING,  1, 1, (void *)&GeoFileName, 0,        ".scuffgeo file"},
     {"LBV1",      PA_DOUBLE,  2, 1, (void *)LBV[0],       0,        "lattice basis vector 1"},
     {"LBV2",      PA_DOUBLE,  2, 1, (void *)LBV[1],       0,        "lattice basis vector 2"},
     {"BlochP",    PA_DOUBLE,  2, 1, (void *)BlochP,       0,        "bloch wavevector"},
     {"Omega",     PA_CDOUBLE, 1, 1, (void *)&Omega,       0,        "angular frequency"},
     {"ZValue",    PA_DOUBLE,  1, 1, (void *)&ZValue,      0,        "Z value"},
     {"Cache",     PA_STRING,  1, 1, (void *)&Cache,       0,        ".scuffcache file"},
     {"EPFile",    PA_STRING,  1, 1, (void *)&EPFile,      0,        "list of eval points"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFileName==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");
  SetLogFileName("tPBC.log");

  RWGGeometry *G=new RWGGeometry(GeoFileName);
  G->SetLogLevel(SCUFF_VERBOSELOGGING);
  PBCGeometry *PG=new PBCGeometry(G, LBV);

  if (Cache) PreloadCache(Cache);
  HMatrix *M = G->AllocateBEMMatrix();
  PG->AssembleBEMMatrix(Omega, BlochP, M);
  //PG->G->AssembleBEMMatrix(Omega, M);
  if (Cache) StoreCache(Cache);
  M->LUFactorize();

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble E0[3]  = {1.0, 0.0,  0.0};
  double nHat[3] = {0.0, 0.0, -1.0};
  PlaneWave PW(E0, nHat);
  HVector *KN = G->AllocateRHSVector();
  PG->AssembleRHSVector(Omega, &PW, KN);

PG->G->PlotSurfaceCurrents(KN, Omega, "%s.RHS.pp", GetFileBase(GeoFileName));

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  M->LUSolve(KN);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
#if 0
  int NX=61;
  int NY=9;
  HVector *XCoords=LinSpace(-LBV[0][0], 2.0*LBV[0][0], NX); 
  HVector *YCoords=LinSpace(-LBV[1][1], 2.0*LBV[1][1], NY);
  HMatrix *XMatrix=new HMatrix(NX*NY, 3);
  int nr, nx, ny; 
  for(nr=nx=0; nx<NX; nx++)
   for(ny=0; ny<NY; ny++, nr++)
    { XMatrix->SetEntry(nr, 0, XCoords->GetEntryD(nx));
      XMatrix->SetEntry(nr, 1, XCoords->GetEntryD(ny));
      XMatrix->SetEntry(nr, 2, ZValue);
    };
#endif
  HMatrix *XMatrix=new HMatrix(EPFile); 
  if (XMatrix->ErrMsg)
   ErrExit(XMatrix->ErrMsg);
 
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  PG->G->PlotSurfaceCurrents(KN, Omega, "%s.KN.pp", GetFileBase(GeoFileName));

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  HMatrix *FMatrix=PG->GetFields(&PW, KN, Omega, BlochP, XMatrix);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  char Buffer[100];
  snprintf(Buffer,100,"%s.out",GetFileBase(GeoFileName));
  FILE *f=CreateUniqueFile(Buffer, 1);
  for(int nr=0; nr<XMatrix->NR; nr++)
   fprintf(f,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n",
              XMatrix->GetEntryD(nr,0), XMatrix->GetEntryD(nr,1), XMatrix->GetEntryD(nr,2),
              real(FMatrix->GetEntry(nr,0)), imag(FMatrix->GetEntry(nr,0)),
              real(FMatrix->GetEntry(nr,1)), imag(FMatrix->GetEntry(nr,1)),
              real(FMatrix->GetEntry(nr,2)), imag(FMatrix->GetEntry(nr,2)),
              real(FMatrix->GetEntry(nr,3)), imag(FMatrix->GetEntry(nr,3)),
              real(FMatrix->GetEntry(nr,4)), imag(FMatrix->GetEntry(nr,4)),
              real(FMatrix->GetEntry(nr,5)), imag(FMatrix->GetEntry(nr,5)));
  fclose(f);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  printf("Thank you for your support.\n");

}
