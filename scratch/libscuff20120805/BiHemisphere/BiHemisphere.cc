#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <libhrutil.h>
#include "libscuff.h"

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[]) 
{  
  EnableAllCPUs();
  InstallHRSignalHandler();
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  char *GeoFileName=0;
  char *Cache=0;
  char *EPFile=0;
  cdouble Omega=0.01;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"geometry", PA_STRING,  1, 1, (void *)&GeoFileName,  0,  ".scuffgeo file"},
     {"cache",    PA_STRING,  1, 1, (void *)&Cache,        0,  ".scuffcache file"},
     {"Omega",    PA_CDOUBLE, 1, 1, (void *)&Omega,        0,  "angular frequency"},
     {"EPFile",   PA_STRING,  1, 1, (void *)&EPFile,       0,  "list of evaluation points"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFileName==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  SetLogFileName("BiHemisphere.log");
  RWGGeometry *G = new RWGGeometry(GeoFileName, SCUFF_VERBOSELOGGING);
  HMatrix *M = G->AllocateBEMMatrix();
  HVector *KN = G->AllocateRHSVector();
  
  cdouble E0[3]  = {0.0, 0.0, 1.0};
  double nHat[3] = {1.0, 0.0, 0.0};
  PlaneWave PW(E0, nHat);

  PreloadCache(Cache);
  G->AssembleBEMMatrix(Omega, M);
  StoreCache(Cache);
  M->LUFactorize();

  G->AssembleRHSVector(Omega, &PW, KN);
  M->LUSolve(KN);

  HMatrix *XMatrix = new HMatrix(EPFile);
  HMatrix *FMatrix= G->GetFields(&PW, KN, Omega, XMatrix);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  char buffer[1000];
  sprintf(buffer,"%s.",GetFileBase(GeoFileName));
  strcat(buffer,GetFileBase(EPFile));
  strcat(buffer,".out");
  FILE *f=fopen(buffer,"w");
  for(int nr=0; nr<XMatrix->NR; nr++)
   { fprintf(f,"%e %e %e ", XMatrix->GetEntryD(nr,0), XMatrix->GetEntryD(nr,1), XMatrix->GetEntryD(nr,2));
     fprintf(f,"%e %e ", real(FMatrix->GetEntry(nr,0)), imag(FMatrix->GetEntry(nr,0)));
     fprintf(f,"%e %e ", real(FMatrix->GetEntry(nr,1)), imag(FMatrix->GetEntry(nr,1)));
     fprintf(f,"%e %e ", real(FMatrix->GetEntry(nr,2)), imag(FMatrix->GetEntry(nr,2)));
     fprintf(f,"%e %e ", real(FMatrix->GetEntry(nr,3)), imag(FMatrix->GetEntry(nr,3)));
     fprintf(f,"%e %e ", real(FMatrix->GetEntry(nr,4)), imag(FMatrix->GetEntry(nr,4)));
     fprintf(f,"%e %e ", real(FMatrix->GetEntry(nr,5)), imag(FMatrix->GetEntry(nr,5)));
     fprintf(f,"\n");
   };
  fclose(f);

}
