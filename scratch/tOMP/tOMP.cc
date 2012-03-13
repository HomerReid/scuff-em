#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <pthread.h>

#include <libhrutil.h>

typedef struct ThreadData
{ 
  int nt, nThread;
  double NMax;
  double PartialSum;

} ThreadData;


/***************************************************************/
/* this routine sums the inverse squares of the integers from  */
/* n=1 to NMax, retaining only the contributions of integers   */
/* satisfying the congruence n==nt modulo nThread.             */
/***************************************************************/
void *BigExpensiveRoutine(void *v)
{ 
  ThreadData *TD=(ThreadData *)v;

  double n;
  int nt=0;

  TD->PartialSum=0.0;
  for(n=1.0; n<TD->NMax; n+=1.0)
   { 
     nt++;
     if ( nt==TD->nThread ) nt=0;
     if ( nt!=TD->nt ) continue;

     TD->PartialSum += 1.0/(n*n);

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
  double NMax=0.0;
  int nThread=0;
  int UsePThreads=0;
  ArgStruct ASArray[]=
   { 
     {"NMax",      PA_DOUBLE, (void *)&NMax,        0, "NMax parameter"},
     {"nThread",   PA_INT,    (void *)&nThread,     0, "number of threads"},
     {"pthreads",  PA_BOOL,   (void *)&UsePThreads, 0, "use pthreads instead of openMP"},
     {0,0,0,0,0}
   };
  ProcessArguments(argc, argv, ASArray);

  if (nThread==0)
   nThread=GetNumProcs();

  if (NMax==0.0)
   ASUsage(argv[0],ASArray,"--nmax option is mandatory");

  /*--------------------------------------------------------------*/
  /*- initialize thread data structures --------------------------*/
  /*--------------------------------------------------------------*/
  int nt;
  ThreadData TDs[nThread];

  for(nt=0; nt<nThread; nt++)
   { TDs[nt].nt=nt;
     TDs[nt].nThread=nThread;
     TDs[nt].NMax=NMax;
   };

  /*--------------------------------------------------------------*/
  /*- use either pthreads or openmp to parallellize --------------*/
  /*--------------------------------------------------------------*/
  Tic();
  if ( UsePThreads )
   { 
     pthread_t Threads[nThread];

     for(nt=0; nt<nThread; nt++)
      pthread_create( &(Threads[nt]), 0, BigExpensiveRoutine, (void *)(&(TDs[nt])));

     for(nt=0; nt<nThread; nt++)
      pthread_join(Threads[nt],0);

   }
  else
   {
     #pragma omp parallel for schedule(static,1), num_threads(nThread)
     for(nt=0; nt<nThread; nt++)
      { 
        BigExpensiveRoutine( (void *)(TDs+nt) );
      };
   };
  Toc();

  double TotalSum=0.0;
  for(nt=0; nt<nThread; nt++)
   TotalSum+=TDs[nt].PartialSum;

  /*--------------------------------------------------------------*/
  /*- print results ----------------------------------------------*/
  /*--------------------------------------------------------------*/
  double ExactAnswer=M_PI*M_PI/6.0;

  printf("\n**\n** Time elapsed: %e s\n**\n",Toc());

  printf("%18s   %18s   %5s\n","     Summation    ","       Exact      ","  RD  ");
  printf("%18s---%18s---%5s\n","------------------","------------------","------");
  printf("%18.12e   %18.12e   %5.2e\n",
          TotalSum,ExactAnswer, fabs(TotalSum-ExactAnswer)/(ExactAnswer));

}
