#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include <libhrutil.h>

typedef struct ThreadData
{ 
  int nt, nThread;
  unsigned long NMax;
  double PartialSum;

} ThreadData;


/***************************************************************/
/***************************************************************/
/***************************************************************/
void *BigExpensiveRoutine(void *v)
{ 
  ThreadData *TD=(ThreadData *)v;

  unsigned long n;
  int nt=0;

  TD->PartialSum=0.0;
  for(n=1; n<TD->NMax; n++)
   { 
     nt++;
     if ( nt==TD->nThread ) nt=0;
     if ( nt!=TD->nt ) continue;

     TD->PartialSum += 1.0/((double)(n*n));

   };

  return 0;
  
} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{ 
  /*--------------------------------------------------------------*/
  /*- process command-line arguments -----------------------------*/
  /*--------------------------------------------------------------*/
  double dNMax=0.0;
  int nThread=0;
  ArgStruct ASArray[]=
   { 
     {"NMax",      PA_DOUBLE, (void *)&dNMax,       0, "NMax parameter"},
     {"nThread",   PA_INT,    (void *)&nThread,     0, "number of threads"},
     {0,0,0,0,0}
   };
  ProcessArguments(argc, argv, ASArray);

  unsigned long NMax = (unsigned long)dNMax;
  if (NMax==0)
   ASUsage(argv[0],ASArray,"--nmax option is mandatory");

  if (nThread==0)
   nThread=GetNumProcs();

  /*--------------------------------------------------------------*/
  /*- fire off threads  ------------------------------------------*/
  /*--------------------------------------------------------------*/
  Tic();
  ThreadData MyTD;
  double TotalSum=0.0;
  int nt;

//#pragma omp parallel for private(MyTD), schedule(static,1), num_threads(nThread)
//#pragma omp parallel for private(MyTD), schedule(auto), num_threads(nThread)
#pragma omp parallel for private(MyTD), schedule(dynamic), num_threads(nThread)
  for(nt=0; nt<nThread; nt++)
   { 
     MyTD.nt=nt;
     MyTD.nThread=nThread;
     MyTD.NMax=NMax;

     BigExpensiveRoutine( (void *)&MyTD);
     TotalSum+=MyTD.PartialSum;
     
   };


  printf("Time elapsed: %e\n",Toc());

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double ExactAnswer=M_PI*M_PI/6.0;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  printf("%18s   %18s   %5s\n","     Summation    ","       Exact      ","  RD  ");
  printf("%18s---%18s---%5s\n","------------------","------------------","------");
  printf("%18.12e   %18.12e   %5.2e\n",
          TotalSum,ExactAnswer, fabs(TotalSum-ExactAnswer)/(ExactAnswer));

}
