#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/timeb.h>

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#ifdef USE_OPENMP
#  include <omp.h>
#endif

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
  HMatrix *DMNorm=0, *DMSym=0, *ZMNorm=0, *ZMHerm=0, *ZMSym=0;
  HMatrix *DMCopy=0, *ZMCopy=0;

  DMNorm=new HMatrix(Dim,Dim,LHM_REAL,LHM_NORMAL);
  ZMNorm=new HMatrix(Dim,Dim,LHM_COMPLEX,LHM_NORMAL);
  if (Symmetric)
   { DMSym=new HMatrix(Dim,Dim,LHM_REAL,LHM_SYMMETRIC);
     ZMHerm=new HMatrix(Dim,Dim,LHM_COMPLEX,LHM_HERMITIAN);
     ZMSym=new HMatrix(Dim,Dim,LHM_COMPLEX,LHM_SYMMETRIC);
   };

  DMCopy=new HMatrix(Dim,Dim,LHM_REAL);
  ZMCopy=new HMatrix(Dim,Dim,LHM_COMPLEX);

  SetConsoleLogging();

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int nr, nc;
  srand48(time(0));

  int NumThreads=GetNumThreads();
printf("%i threads.\n",NumThreads);
  struct random_data **RDs=(struct random_data **)mallocEC(NumThreads*sizeof(struct random_data *));
  int StateSize=16;
  for(int nt=0; nt<NumThreads; nt++)
   { struct timeb tp;
     ftime(&tp);
     unsigned long Seed=1000*tp.time + tp.millitm;
     char *State=(char *)mallocEC(StateSize*sizeof(char));
     RDs[nt]=(struct random_data *)mallocEC(sizeof(struct random_data));
     initstate_r(Seed,State,StateSize,RDs[nt]);
   };

#ifdef USE_OPENMP
  printf("Using %i threads for matrix assembly...\n",NumThreads);
#pragma omp parallel for num_threads(NumThreads)
#endif
  for(nr=0; nr<Dim; nr++)
   { 
     int nt=0;
#ifdef USE_OPENMP
     nt = omp_get_thread_num();
#endif
     LogPercent(nr,Dim,10);
     for(nc=0; nc<Dim; nc++)
      { int ri, ii;
        random_r(RDs[nt], &ri);
        random_r(RDs[nt], &ii);
        int rd=((double)ri)/((double)RAND_MAX);
        int id=((double)ii)/((double)RAND_MAX);
        cdouble z=cdouble( 2*(rd-0.5), 2*(id-0.5) );
        DMNorm->SetEntry(nr, nc, real(z));
        ZMNorm->SetEntry(nr, nc, z);
      };
   };

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
