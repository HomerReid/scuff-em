/*
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <libhrutil.h>
#include <libhmat.h>
#include <libMatProp.h>
#include <libscuff.h>

using namespace scuff;



/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  /*--------------------------------------------------------------*/
  /* create the RWGGeometry from the .scuffgeo file               */
  /*--------------------------------------------------------------*/
  //RWGGeometry *G=new RWGGeometry(f,"BiHemisphere.scuffgeo");
  FILE *f=fopen("BiHemisphere.scuffgeo");
  char Line[100];
  int LineNum=1;
  fgets(Line,100,f);
  RWGComposite *C=new RWGComposite(f,"BiHemisphere",&LineNum);
  C->EpsMu[0]=new MatProp(MP_VACUUM);
  if (C->ErrMsg)
   ErrExit("%s: %i: %s\n","BiHemisphere.scuffgeo",LineNum,C->ErrMsg);

  SetLogFileName("BiHemisphere.log");
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
  cdouble E0[3]  = {0.0, 0.0, 1.0};  // point source location 
  double nHat[3] = {1.0, 0.0, 0.0};  // point source location 
  PlaneWave PW(E0, nHat);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double Omega=0.001;

  ACCMBArgStruct MyACCMBArgs, *ACCMBArgs=&MyACCMBArgs;
  ACCMBArgs->CA=C;
  ACCMBArgs->CB=C;
  ACCMBArgs->Omega=Omega;
  ACCMBArgs->M=M;
  ACCMBArgs->RowOffset=0;
  ACCMBArgs->ColOffset=0;

  AssembleCCMatrixBlock(ACCMBArgs);

  M->LUFactorize();
  StoreCache("BiHemisphere.scuffcache");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  AssembleRHSVector_Composite(C, Omega, PW, KN, 0);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  M->LUSolve(KN);

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
  int nr0, nr1, nr2
  for(X=XMin; X<=1.01*XMax; X+= (XMax-XMin)/((double)(NXPoints-1)))
   for(Z=ZMin; Z<=1.01*ZMax; Z+= (ZMax-ZMin)/((double)(NZPoints-1)))
    { 
      R=sqrt(X*X+Z*Z);
      if ( fabs( R - 1.0 ) < 0.05 )
       continue;
      if ( fabs( z - 0.0 ) < 0.05 )
       continue;
 
      if( R>1.0 )
       { EvalPoints0->SetEntry(nr0,0,X);
         EvalPoints0->SetEntry(nr0,1,0.0);
         EvalPoints0->SetEntry(nr0,2,Z);
         nr0++;
       }
      else if( R<1.0 && z<0.0)
       { EvalPoints1->SetEntry(nr1,0,X);
         EvalPoints1->SetEntry(nr1,1,0.0);
         EvalPoints1->SetEntry(nr1,2,Z);
         nr1++;
       }
      else if( R<1.0 && z>0.0)
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

  FILE *f=fopen("EvalPoints.dat","w");
  for(nr=0; nr<nr0; nr++)
   { EP0->SetEntry(nr, 0, EvalPoints0->GetEntry(nr,0));
     EP0->SetEntry(nr, 1, EvalPoints0->GetEntry(nr,1));
     EP0->SetEntry(nr, 2, EvalPoints0->GetEntry(nr,2));
     fprintf(f,"%e %e %e \n", EP0->GetEntry(nr,0), EP0->GetEntry(nr,1), EP0->GetEntry(nr,2));
   };
  fprintf(f,"\n\n");

  for(nr=0; nr<nr1; nr++)
   { EP1->SetEntry(nr, 0, EvalPoints1->GetEntry(nr,0));
     EP1->SetEntry(nr, 1, EvalPoints1->GetEntry(nr,1));
     EP1->SetEntry(nr, 2, EvalPoints1->GetEntry(nr,2));
     fprintf(f,"%e %e %e \n", EP1->GetEntry(nr,0), EP1->GetEntry(nr,1), EP1->GetEntry(nr,2));
   };

  for(nr=0; nr<nr1; nr++)
   { EP2->SetEntry(nr, 0, EvalPoints2->GetEntry(nr,0));
     EP2->SetEntry(nr, 1, EvalPoints2->GetEntry(nr,1));
     EP2->SetEntry(nr, 2, EvalPoints2->GetEntry(nr,2));
     fprintf(f,"%e %e %e \n", EP2->GetEntry(nr,0), EP2->GetEntry(nr,1), EP2->GetEntry(nr,2));
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  HMatrix *FMatrix0=GetFields(C, 0, IF, KN, Omega, EP0);
  HMatrix *FMatrix1=GetFields(C, 1, IF, KN, Omega, EP0);
  HMatrix *FMatrix2=GetFields(C, 2, IF, KN, Omega, EP2);

}
