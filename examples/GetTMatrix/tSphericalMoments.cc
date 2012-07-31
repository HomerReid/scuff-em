/*
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include <libscuff.h>
#include "GetTMatrix.h"

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  EnableAllCPUs();

  /*--------------------------------------------------------------*/
  /*- process command-line options -------------------------------*/
  /*--------------------------------------------------------------*/
  char *GeoFileName=0;  // .scuffgeo file
  int lMax=5;           // maximum L-value of spherical wave computed
  cdouble Omega=0;      // angular frequency at which to run the computation
  char *Cache=0;        // scuff cache file 
  int l=1, m=0;
  int Type=0;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"geometry",  PA_STRING,  1, 1, (void *)&GeoFileName,  0,  ".scuffgeo file"},
     {"l",         PA_INT,     1, 1, (void *)&l,            0,  "l value"},
     {"m",         PA_INT,     1, 1, (void *)&m,            0,  "m value"},
     {"electric",  PA_BOOL,    0, 1, (void *)&Type,         0,  "use electric monopole (default is magnetic)"},
     {"lMax",      PA_INT,     1, 1, (void *)&lMax,         0,  "maximum l-value"},
     {"Omega",     PA_CDOUBLE, 1, 1, (void *)&Omega,        0,  "angular frequency"},
     {"Cache",     PA_STRING,  1, 1, (void *)&Cache,        0,  "scuff cache file"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFileName==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  SetLogFileName("%s.log",GetFileBase(GeoFileName));
  RWGGeometry *G=new RWGGeometry(GeoFileName);
  HVector *KN = G->AllocateRHSVector();

  SphericalWave SW(l, m, Type);
  SW.Omega=Omega;
  SW.Eps=SW.Mu=1.0;
  G->ExpandCurrentDistribution(&SW, KN);

  int nm, NumMoments = 2*(lMax+1)*(lMax+1);
  HVector *AVector=GetSphericalMoments(G->Objects[0], Omega, lMax, KN, 0);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  SetDefaultCD2SFormat("%+15.8e %+15.8e");
  char FileName[100];
  sprintf(FileName,"%s.%i.%i.SphericalMoments",GetFileBase(GeoFileName),l,m);
  FILE *f=vfopen("%s.%i.%i.SphericalMoments","w",GetFileBase(GeoFileName),l,m);
  const char *TypeChar="ME";
  for(nm=l=0; l<=lMax; l++)
   for(m=-l; m<=l; m++)
    for(Type=0; Type<=1; Type++, nm++)
     fprintf(f,"%i %3i  %3i  %s\n",Type,l,m,CD2S(AVector->GetEntry(nm)));
  fclose(f);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  char buffer[200];
  sprintf(buffer,"sort -g -k1 -k3 -k2 %s > doom",FileName);
  system(buffer);
  sprintf(buffer,"/bin/mv doom %s",FileName);
  system(buffer);

}
