#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <libhrutil.h>
#include <libhmat.h>

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{ 
  srand48(time(0));

  /*--------------------------------------------------------------*/
  /*- process command-line arguments -----------------------------*/
  /*--------------------------------------------------------------*/
  int NStart, NStop, DeltaN, Real;
  ArgStruct ASArray[]=
   { {"NStart",    PA_INT,  (void *)&NStart,    "1000",  "NStart"},
     {"NStop",     PA_INT,  (void *)&NStop,     "10000", "NStop"},
     {"DeltaN",    PA_INT,  (void *)&DeltaN,    "100",   "DeltaN"}, 
     {"Real",      PA_BOOL, (void *)&Real,        "0",   "use real-valued matrices"}, 
     {0,0,0,0,0}
   };
  ProcessArguments(argc, argv, ASArray);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int nr, NR, nc, NC;
  HMatrix *A, *B, *AB;
  HVector *DAB;
  double Time1, Time2, Error;
  printf("(%6s,%6s)   %10s   %10s   %6s\n","  NR  ","  NC  "," Slow way "," Fast way ","Error ");
  printf("(%6s,%6s)---%10s---%10s---%6s\n","------","------","----------","----------","------");
  for(NR=NStart; NR<NStop; NR+=DeltaN)
   { 
     // we will do the test with non-square matrices just for fun
     NC = 3*NR / 4;

     if (Real)
      { A   = new HMatrix(NR, NC, LHM_COMPLEX);
        B   = new HMatrix(NC, NR, LHM_COMPLEX);
        AB  = new HMatrix(NR, NR, LHM_COMPLEX); 
        DAB = new HVector(NR, LHM_COMPLEX);

        for(nr=0; nr<NR; nr++)
         for(nc=0; nc<NC; nc++)
          { A->SetEntry( nr, nc, cdouble(drand48(), drand48()) );
            B->SetEntry( nc, nr, cdouble(drand48(), drand48()) );
          };
      }
     else
      { A   = new HMatrix(NR, NC);
        B   = new HMatrix(NC, NR);
        AB  = new HMatrix(NR, NR);
        DAB = new HVector(NR);

        for(nr=0; nr<NR; nr++)
         for(nc=0; nc<NC; nc++)
          { A->SetEntry( nr, nc, drand48());
            B->SetEntry( nc, nr, drand48());
          };
      };
     
     Tic();
     A->Multiply(B, AB);
     Time1=Toc();

     Tic();
     A->GetMatrixProductDiagonal(B, DAB);
     Time2=Toc();

     for(Error=0.0, nr=0; nr<NR; nr++)
      Error += norm( DAB->GetEntry(nr) - AB->GetEntry(nr, nr) ); 
     Error=sqrt(Error) / ((double)NR);

     printf("(%6i,%6i)   %+10.3f   %+10.3f   %6.3e\n",NR,NC,Time1, Time2, Error);

     delete A;
     delete B;
     delete AB;
     delete DAB;
   };

}    
