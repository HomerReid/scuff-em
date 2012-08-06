/*
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <libhrutil.h>
#include <libhmat.h>
#include <libMatProp.h>
#include "libscuff.h"
#include "RWGComposite.h"

using namespace scuff;
void GetCPUAffinity();
extern int PSOnly;

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  char *GeoFileName=0;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"geometry",  PA_STRING,  1, 1, (void *)&GeoFileName, 0, ".scuffgeo file"},
     {"psonly",    PA_INT,     1, 1, (void *)&PSOnly,      0, "psonly"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFileName==0)
   OSUsage(argv[0],OSArray,"--Geometry option is mandatory");

  SetLogFileName("tComposite.log");
  Log("tComposite running on %s",GetHostName());
  GetCPUAffinity();

  /*--------------------------------------------------------------*/
  /* create the RWGGeometry from the .scuffgeo file               */
  /*--------------------------------------------------------------*/
  //RWGGeometry *G=new RWGGeometry(f,"BiHemisphere.scuffgeo");
  FILE *f=fopen(GeoFileName,"r");
  char Line[100];
  int LineNum=1;
  fgets(Line,100,f);
EnableAllCPUs(); 
GetCPUAffinity();
  RWGComposite *C=new RWGComposite(f,"BiHemisphere",&LineNum);
  C->SubRegionMPs[0]=new MatProp(MP_VACUUM);
  if (C->ErrMsg)
   ErrExit("%s: %i: %s\n","BiHemisphere.scuffgeo",LineNum,C->ErrMsg);

  //SetLogLevel(SCUFF_VERBOSELOGGING);
  PreloadCache("BiHemisphere.scuffcache");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int NTBF=C->TotalBFs;
  HMatrix *M  = new HMatrix(NTBF, NTBF, LHM_COMPLEX);
  HVector *KN = new HVector(NTBF, LHM_COMPLEX);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble E0[3]  = {0.0, 0.0, 1.0};
  double nHat[3] = {1.0, 0.0, 0.0};
  PlaneWave PW(E0, nHat);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double Omega=0.001;

  ACCMBArgStruct MyACCMBArgs, *ACCMBArgs=&MyACCMBArgs;
  ACCMBArgs->CA=C;
  ACCMBArgs->CB=C;
  ACCMBArgs->Omega=Omega;
  ACCMBArgs->B=M;
  ACCMBArgs->RowOffset=0;
  ACCMBArgs->ColOffset=0;

  AssembleCCMatrixBlock(ACCMBArgs);
  //AddEdgePanelContributions(ACCMBArgs);

  M->LUFactorize();
  //StoreCache("BiHemisphere.scuffcache");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  PW.SetFrequency(Omega);
  AssembleRHSVector_Composite(C, Omega, &PW, KN, 1);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  M->LUSolve(KN);

  PlotSurfaceCurrents(C, KN, Omega, "%s.pp",GetFileBase(GeoFileName));
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int NXPoints=6, NZPoints=100;
  double X, XMin=0.0,  XMax=2.0; 
  double Z, ZMin=-3.0, ZMax=3.0; 
  double R;
  HMatrix *EvalPoints0=new HMatrix(NXPoints*NZPoints,3);
  HMatrix *EvalPoints1=new HMatrix(NXPoints*NZPoints,3);
  HMatrix *EvalPoints2=new HMatrix(NXPoints*NZPoints,3);
  int nr0=0, nr1=0, nr2=0;
  for(X=XMin; X<=1.01*XMax; X+= (XMax-XMin)/((double)(NXPoints-1)))
   for(Z=ZMin; Z<=1.01*ZMax; Z+= (ZMax-ZMin)/((double)(NZPoints-1)))
    { 
      R=sqrt(X*X+Z*Z);
      if ( fabs( R - 1.0 ) < 0.05 )
       continue;
      if ( fabs( Z - 0.0 ) < 0.05 )
       continue;
 
      if( R>1.0 )
       { EvalPoints0->SetEntry(nr0,0,X);
         EvalPoints0->SetEntry(nr0,1,0.0);
         EvalPoints0->SetEntry(nr0,2,Z);
         nr0++;
       }
      else if( R<1.0 && Z<0.0)
       { EvalPoints1->SetEntry(nr1,0,X);
         EvalPoints1->SetEntry(nr1,1,0.0);
         EvalPoints1->SetEntry(nr1,2,Z);
         nr1++;
       }
      else if( R<1.0 && Z>0.0)
       { EvalPoints2->SetEntry(nr2,0,X);
         EvalPoints2->SetEntry(nr2,1,0.0);
         EvalPoints2->SetEntry(nr2,2,Z);
         nr2++;
       }
      else
       ErrExit("bawonkatage");
    };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  HMatrix *EP0=new HMatrix(nr0,3);
  HMatrix *EP1=new HMatrix(nr1,3);
  HMatrix *EP2=new HMatrix(nr2,3);
  int nr;

  f=fopen("EvalPoints.dat","w");
  for(nr=0; nr<nr0; nr++)
   { EP0->SetEntry(nr, 0, EvalPoints0->GetEntry(nr,0));
     EP0->SetEntry(nr, 1, EvalPoints0->GetEntry(nr,1));
     EP0->SetEntry(nr, 2, EvalPoints0->GetEntryD(nr,2));
     fprintf(f,"%e %e %e \n", EP0->GetEntryD(nr,0), EP0->GetEntryD(nr,1), EP0->GetEntryD(nr,2));
   };
  fprintf(f,"\n\n");

  for(nr=0; nr<nr1; nr++)
   { EP1->SetEntry(nr, 0, EvalPoints1->GetEntryD(nr,0));
     EP1->SetEntry(nr, 1, EvalPoints1->GetEntryD(nr,1));
     EP1->SetEntry(nr, 2, EvalPoints1->GetEntryD(nr,2));
     fprintf(f,"%e %e %e \n", EP1->GetEntryD(nr,0), EP1->GetEntryD(nr,1), EP1->GetEntryD(nr,2));
   };

  for(nr=0; nr<nr1; nr++)
   { EP2->SetEntry(nr, 0, EvalPoints2->GetEntryD(nr,0));
     EP2->SetEntry(nr, 1, EvalPoints2->GetEntryD(nr,1));
     EP2->SetEntry(nr, 2, EvalPoints2->GetEntryD(nr,2));
     fprintf(f,"%e %e %e \n", EP2->GetEntryD(nr,0), EP2->GetEntryD(nr,1), EP2->GetEntryD(nr,2));
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  HMatrix *FMatrix0=GetFields(C, 0, &PW, KN, Omega, EP0);
  HMatrix *FMatrix1=GetFields(C, 1, &PW, KN, Omega, EP1);
  HMatrix *FMatrix2=GetFields(C, 2, &PW, KN, Omega, EP2);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  f=vfopen("%s.out","w",GetFileBase(GeoFileName)); 
  SetDefaultCD2SFormat("%+10.3e %+10.3e");
  for(nr=0; nr<nr0; nr++)
   { fprintf(f,"%e %e %e ", EP0->GetEntryD(nr,0), EP0->GetEntryD(nr,1), EP0->GetEntryD(nr,2));
     fprintf(f,"%s %s %s %s %s %s \n", CD2S(FMatrix0->GetEntry(nr,0)), CD2S(FMatrix0->GetEntry(nr,1)), 
                                       CD2S(FMatrix0->GetEntry(nr,2)), CD2S(FMatrix0->GetEntry(nr,3)), 
                                       CD2S(FMatrix0->GetEntry(nr,4)), CD2S(FMatrix0->GetEntry(nr,5)));
   };
  fprintf(f,"\n\n");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  for(nr=0; nr<nr1; nr++)
   { fprintf(f,"%e %e %e ", EP1->GetEntryD(nr,0), EP1->GetEntryD(nr,1), EP1->GetEntryD(nr,2));
     fprintf(f,"%s %s %s %s %s %s \n", CD2S(FMatrix1->GetEntry(nr,0)), CD2S(FMatrix1->GetEntry(nr,1)), 
                                       CD2S(FMatrix1->GetEntry(nr,2)), CD2S(FMatrix1->GetEntry(nr,3)), 
                                       CD2S(FMatrix1->GetEntry(nr,4)), CD2S(FMatrix1->GetEntry(nr,5)));
   };
  fprintf(f,"\n\n");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  for(nr=0; nr<nr2; nr++)
   { fprintf(f,"%e %e %e ", EP2->GetEntryD(nr,0), EP1->GetEntryD(nr,1), EP2->GetEntryD(nr,2));
     fprintf(f,"%s %s %s %s %s %s \n", CD2S(FMatrix2->GetEntry(nr,0)), CD2S(FMatrix2->GetEntry(nr,1)), 
                                       CD2S(FMatrix2->GetEntry(nr,2)), CD2S(FMatrix2->GetEntry(nr,3)), 
                                       CD2S(FMatrix2->GetEntry(nr,4)), CD2S(FMatrix2->GetEntry(nr,5)));
   };
  fprintf(f,"\n\n");

}
