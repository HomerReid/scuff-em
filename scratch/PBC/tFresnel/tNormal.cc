/*
 * tNormal.cc
 */
#include <stdio.h>
#include <stdlib.h>

#include <libhrutil.h>
#include <libscuff.h>
#include <libIncField.h>
#include "../PBCGeometry.h"

using namespace scuff;

extern int GetFieldObject;

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
  cdouble Omega=0.1;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"geometry", PA_STRING,  1, 1, (void *)&GeoFileName,  0,  ".scuffgeo file"},
     {"omega",    PA_CDOUBLE, 1, 1, (void *)&Omega,        0,  ".scuffgeo file"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFileName==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");
  
  /*--------------------------------------------------------------*/
  /*- create the RWGGeometry and PBCGeometry structures          -*/
  /*--------------------------------------------------------------*/
  SetLogFileName("%s.log",GeoFileName);
  RWGGeometry *G=new RWGGeometry(GeoFileName);

  double LBV1[2]={1.0, 0.0};
  double LBV2[2]={0.0, 1.0};
  double *LBV[2]={LBV1, LBV2};
  PBCGeometry *PG=new PBCGeometry(G, LBV);

  HMatrix *M  = PG->AllocateBEMMatrix();
  HVector *KN = PG->AllocateRHSVector();

  char CacheFileName[100];
  sprintf(CacheFileName,"%s.cache",G->Objects[0]->MeshFileName);
  PreloadCache(CacheFileName);

  double kVacuum = real(Omega); 

  cdouble EpsMedium, MuMedium;
  G->Objects[0]->MP->GetEpsMu(Omega, &EpsMedium, &MuMedium );
  double nMedium = real( csqrt2(EpsMedium*MuMedium) );
  double kMedium = real( nMedium * Omega );

  /*--------------------------------------------------------------*/
  /*- assemble the periodic BEM matrix                           -*/
  /*--------------------------------------------------------------*/
  double BlochP[2]={0.0, 0.0};
  PG->AssembleBEMMatrix(Omega, BlochP, M);
  M->LUFactorize();

  /*--------------------------------------------------------------*/
  /*- create the incident field and solve the scattering problem -*/
  /*--------------------------------------------------------------*/
  cdouble E0[3]={0.0, 1.0, 0.0};
  double nHat[3]={0.0, 0.0, -1.0};
  PlaneWave PW(E0, nHat);
  PG->AssembleRHSVector(Omega, &PW, KN);
  M->LUSolve(KN);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int n;
  HVector *ZVector=LinSpace(1.0, 5.0, 25); 
  HMatrix *XAbove=new HMatrix(ZVector->N, 3);
  HMatrix *XBelow=new HMatrix(ZVector->N, 3);
  for(n=0; n<ZVector->N; n++)
   { 
     XAbove->SetEntry(n, 2,  ZVector->GetEntry(n));
     XBelow->SetEntry(n, 2, -ZVector->GetEntry(n));
   };

GetFieldObject=-1;
   HMatrix *FAbove = PG->GetFields(0, KN, Omega, BlochP, XAbove);
GetFieldObject=0;
   HMatrix *FBelow = PG->GetFields(0, KN, Omega, BlochP, XBelow);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  SetDefaultCD2SFormat("%15.8e %15.8e");
  FILE *f = fopen("tNormal.out","w");
  for(n=0; n<ZVector->N; n++)
   { fprintf(f,"%e %e %e ", XAbove->GetEntryD(n, 0), XAbove->GetEntryD(n, 1), XAbove->GetEntryD(n, 2));
     fprintf(f,"%s ", CD2S(FAbove->GetEntry(n,0)));
     fprintf(f,"%s ", CD2S(FAbove->GetEntry(n,1)));
     fprintf(f,"%s ", CD2S(FAbove->GetEntry(n,2)));
     fprintf(f,"%s ", CD2S(FAbove->GetEntry(n,3)));
     fprintf(f,"%s ", CD2S(FAbove->GetEntry(n,4)));
     fprintf(f,"%s ", CD2S(FAbove->GetEntry(n,5)));
     fprintf(f,"\n");
   };
  fprintf(f,"\n\n");

  for(n=0; n<ZVector->N; n++)
   { fprintf(f,"%e %e %e ", XBelow->GetEntryD(n, 0), XBelow->GetEntryD(n, 1), XBelow->GetEntryD(n, 2));
     fprintf(f,"%s ", CD2S(FBelow->GetEntry(n,0)));
     fprintf(f,"%s ", CD2S(FBelow->GetEntry(n,1)));
     fprintf(f,"%s ", CD2S(FBelow->GetEntry(n,2)));
     fprintf(f,"%s ", CD2S(FBelow->GetEntry(n,3)));
     fprintf(f,"%s ", CD2S(FBelow->GetEntry(n,4)));
     fprintf(f,"%s ", CD2S(FBelow->GetEntry(n,5)));
     fprintf(f,"\n");
   };
  fclose(f);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  printf("Thank you for your support.\n");
  StoreCache(CacheFileName);

}
