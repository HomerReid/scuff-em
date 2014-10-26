#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "libhrutil.h"
#include "libhmat.h"

#if defined(_WIN32)
#  define srand48 srand
#  define drand48 my_drand48
static double my_drand48(void) {
  return rand() * 1.0 / RAND_MAX;
}
#endif

int main(int argc, char *argv[])
{ 
   
  /*--------------------------------------------------------------*/
  /*- process command-line arguments -----------------------------*/
  /*--------------------------------------------------------------*/
  int Dim=10;
  int Symmetric=0;
  ArgStruct ASArray[]=
   { {"Dim",       PA_INT,  (void *)&Dim,       "10", "dimension"},
     {"Symmetric", PA_BOOL, (void *)&Symmetric,    0, "do symmetric cases too"},
     {0,0,0,0,0}
   };
  ProcessArguments(argc, argv, ASArray);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  HMatrix *DMNorm, *DMSym, *ZMNorm, *ZMHerm, *ZMSym;
  HMatrix *DMCopy, *ZMCopy;

  DMNorm=new HMatrix(Dim,Dim,LHM_REAL,LHM_NORMAL);
  ZMNorm=new HMatrix(Dim,Dim,LHM_COMPLEX,LHM_NORMAL);
  if (Symmetric)
   { DMSym=new HMatrix(Dim,Dim,LHM_REAL,LHM_SYMMETRIC);
     ZMHerm=new HMatrix(Dim,Dim,LHM_COMPLEX,LHM_HERMITIAN);
     ZMSym=new HMatrix(Dim,Dim,LHM_COMPLEX,LHM_SYMMETRIC);
   };

  DMCopy=new HMatrix(Dim,Dim,LHM_REAL);
  ZMCopy=new HMatrix(Dim,Dim,LHM_COMPLEX);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int nr, nc;
  srand48(time(0));

  for(nr=0; nr<Dim; nr++)
   for(nc=0; nc<Dim; nc++)
    DMNorm->SetEntry(nr, nc, drand48());

  for(nr=0; nr<Dim; nr++)
   for(nc=0; nc<Dim; nc++)
    ZMNorm->SetEntry(nr, nc, cdouble(drand48(),drand48()));

  if (Symmetric) 
   { 
     for(nr=0; nr<Dim; nr++)
      for(nc=nr; nc<Dim; nc++)
       DMSym->SetEntry(nr, nc, drand48());
   
     for(nr=0; nr<Dim; nr++)
      for(nc=nr; nc<Dim; nc++)
       ZMHerm->SetEntry(nr, nc, cdouble(drand48(), nc>nr ? drand48() : 0.0 ));
   
     for(nr=0; nr<Dim; nr++)
      for(nc=nr; nc<Dim; nc++)
       ZMSym->SetEntry(nr, nc,  cdouble(drand48(), drand48() ));
   };
   
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
#if 0
  void *pCC=HMatrix::OpenHDF5Context("tInvert.hdf5");

  DMNorm->ExportToHDF5(pCC,"DMNorm");

  if(Symmetric)
   { DMCopy->Copy(DMSym); 
     DMCopy->ExportToHDF5(pCC,"DMSym");
   };

  for(nr=0; nr<Dim; nr++)
   for(nc=0; nc<Dim; nc++)
    DMCopy->SetEntry(nr, nc, real(ZMNorm->GetEntry(nr,nc)));
  DMCopy->ExportToHDF5(pCC,"ReZMNorm");
  for(nr=0; nr<Dim; nr++)
   for(nc=0; nc<Dim; nc++)
    DMCopy->SetEntry(nr, nc, imag(ZMNorm->GetEntry(nr,nc)));
  DMCopy->ExportToHDF5(pCC,"ImZMNorm");

  for(nr=0; nr<Dim; nr++)
   for(nc=0; nc<Dim; nc++)
    DMCopy->SetEntry(nr, nc, real(ZMHerm->GetEntry(nr,nc)));
  DMCopy->ExportToHDF5(pCC,"ReZMHerm");
  for(nr=0; nr<Dim; nr++)
   for(nc=0; nc<Dim; nc++)
    DMCopy->SetEntry(nr, nc, imag(ZMHerm->GetEntry(nr,nc)));
  DMCopy->ExportToHDF5(pCC,"ImZMHerm");

  for(nr=0; nr<Dim; nr++)
   for(nc=0; nc<Dim; nc++)
    DMCopy->SetEntry(nr, nc, real(ZMSym->GetEntry(nr,nc)));
  DMCopy->ExportToHDF5(pCC,"ReZMSym");
  for(nr=0; nr<Dim; nr++)
   for(nc=0; nc<Dim; nc++)
    DMCopy->SetEntry(nr, nc, imag(ZMSym->GetEntry(nr,nc)));
  DMCopy->ExportToHDF5(pCC,"ImZMSym");
#endif

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double Elapsed;
  printf("\n** LU Factorize: \n");

  Tic(); DMNorm->LUFactorize(); Elapsed=Toc(); 
  printf("double,  normal:    %.3g s\n",Elapsed); 

  Tic(); ZMNorm->LUFactorize(); Elapsed=Toc(); 
  printf("complex, normal:    %.3g s\n",Elapsed);

  if (Symmetric)
   { 
     Tic(); DMSym->LUFactorize(); Elapsed=Toc(); 
     printf("double,  symmetric: %.3g s\n",Elapsed);

     Tic(); ZMHerm->LUFactorize(); Elapsed=Toc(); 
     printf("complex, hermitian: %.3g s\n",Elapsed);

     Tic(); ZMSym->LUFactorize(); Elapsed=Toc(); 
     printf("complex, symmetric: %.3g s\n",Elapsed);
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  printf("\n** Invert:   \n");

  Tic(); DMNorm->LUInvert(); Elapsed=Toc(); 
  printf("double,  normal:    %.3g s\n",Elapsed); 

  Tic(); ZMNorm->LUInvert(); Elapsed=Toc(); 
  printf("complex, normal:    %.3g s\n",Elapsed);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  printf("\n** Invert 2:   \n");

  HMatrix *DMUnit = new HMatrix(Dim, Dim); 
  DMUnit->Zero();
  for(int n=0; n<Dim; n++)
   DMUnit->SetEntry(n,n,1.0);
  Tic(); DMNorm->LUSolve(DMUnit); Elapsed=Toc(); 
  printf("double,  normal:    %.3g s\n",Elapsed); 

  HMatrix *ZMUnit = new HMatrix(Dim, Dim, LHM_COMPLEX); 
  ZMUnit->Zero();
  for(int n=0; n<Dim; n++)
   ZMUnit->SetEntry(n,n,1.0);
  Tic(); ZMNorm->LUSolve(ZMUnit); Elapsed=Toc(); 
  printf("complex, normal:    %.3g s\n",Elapsed);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
#if 0
  DMNorm->ExportToHDF5(pCC,"InvDMNorm");

  DMCopy->Copy(DMSym); DMCopy->ExportToHDF5(pCC,"InvDMSym");

  for(nr=0; nr<Dim; nr++)
   for(nc=0; nc<Dim; nc++)
    DMCopy->SetEntry(nr, nc, real(ZMNorm->GetEntry(nr,nc)));
  DMCopy->ExportToHDF5(pCC,"ReInvZMNorm");
  for(nr=0; nr<Dim; nr++)
   for(nc=0; nc<Dim; nc++)
    DMCopy->SetEntry(nr, nc, imag(ZMNorm->GetEntry(nr,nc)));
  DMCopy->ExportToHDF5(pCC,"ImInvZMNorm");

  for(nr=0; nr<Dim; nr++)
   for(nc=0; nc<Dim; nc++)
    DMCopy->SetEntry(nr, nc, real(ZMHerm->GetEntry(nr,nc)));
  DMCopy->ExportToHDF5(pCC,"ReInvZMHerm");
  for(nr=0; nr<Dim; nr++)
   for(nc=0; nc<Dim; nc++)
    DMCopy->SetEntry(nr, nc, imag(ZMHerm->GetEntry(nr,nc)));
  DMCopy->ExportToHDF5(pCC,"ImInvZMHerm");

  for(nr=0; nr<Dim; nr++)
   for(nc=0; nc<Dim; nc++)
    DMCopy->SetEntry(nr, nc, real(ZMSym->GetEntry(nr,nc)));
  DMCopy->ExportToHDF5(pCC,"ReInvZMSym");
  for(nr=0; nr<Dim; nr++)
   for(nc=0; nc<Dim; nc++)
    DMCopy->SetEntry(nr, nc, imag(ZMSym->GetEntry(nr,nc)));
  DMCopy->ExportToHDF5(pCC,"ImInvZMSym");

  HMatrix::CloseHDF5Context(pCC);
#endif

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  printf("Thank you for your support.\n");
  
}
