#include "TaylorDuffy.h"
using namespace scuff;

int main(int argc, char *argv[])
{
  /* panel vertices */
  int WhichCase = 2; // 2 common vertices
  double V1[3]  =  { 0.0,  0.0,  0.0 };
  double V2[3]  =  { 0.1,  0.0,  0.0 };
  double V3[3]  =  { 0.05, 0.1,  0.0 };
  double V3P[3] =  { 0.07, -0.08, 0.03 };

  double *Q     = V3;  // source/sink vertex, triangle 1
  double *QP    = V3P; // source/sink vertex, triangle 2

  /* specification of which integrals we want */
  int NumPKs        = 2;
  int PIndex[2]     = {TD_UNITY, TD_PMCHWC};
  int KIndex[2]     = {TD_RP,    TD_RP};
  cdouble KParam[2] = {-1.0,    -3.0};
 
  /* output buffers */
  cdouble Result[2], Error[2];
  
  /* fill in argument structure with problem description */
  TaylorDuffyArgStruct MyTDArgs, *TDArgs=&MyTDArgs;
  InitTaylorDuffyArgs(TDArgs);
  TDArgs->WhichCase = WhichCase;
  TDArgs->V1        = V1;
  TDArgs->V2        = V2;
  TDArgs->V3        = V3;
  TDArgs->V3P       = V3P;
  TDArgs->Q         = Q;
  TDArgs->QP        = QP;
  TDArgs->NumPKs    = NumPKs;
  TDArgs->PIndex    = PIndex;
  TDArgs->KIndex    = KIndex;
  TDArgs->KParam    = KParam;
  TDArgs->Result    = Result;
  TDArgs->Error     = Error;  

  /* specify desired error tolerance */
  TDArgs->RelTol    = 1.0e-10;   // request 10-digit accuracy
  TDArgs->MaxEval   = 25;        // upper limit on integrand samples

  /* calculate the integral */
  TaylorDuffy( TDArgs );

  /* print the results */
  printf("Integrand sampled at %i points.\n",TDArgs->nCalls);
  printf("Integral 1: {%+.8e, %+.8e} (estimated error {%.1e.,%.1e}) \n",
          real(Result[0]),imag(Result[0]),real(Error[0]),imag(Error[0]));
  printf("Integral 2: {%+.8e, %+.8e} (estimated error {%.1e.,%.1e}) \n",
          real(Result[1]),imag(Result[1]),real(Error[1]),imag(Error[1]));

  /* uncomment the following line to report computation time */
//#define MEASURE_RUNTIME
#ifdef MEASURE_RUNTIME
#define REPETITIONS 100
  Tic();
  for(int n=0; n<REPETITIONS; n++)
   TaylorDuffy( TDArgs );
  double TimeElapsed = Toc() / REPETITIONS;
  printf("Computation time: %e microseconds\n",1.0e6*TimeElapsed);
#endif

}
