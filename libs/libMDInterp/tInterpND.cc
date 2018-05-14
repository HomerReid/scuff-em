/*
 * tInterpND
 *
 * homer reid      -- 3/14/2011
 */
#include <stdio.h>
#include <stdlib.h>

#include <libhrutil.h>
#include <libhmat.h>

#include "libMDInterp.h"

#include <readline/readline.h>
#include <readline/history.h>

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define NFUN 3
void Phi(double *X, void *UserData, double *PhiVD)
{
  int D = *((int *)UserData);
  static double Alpha[4]={0.2,0.3,0.4,0.5};
  static double X0[4]={-0.1,0.2,0.35,0.5};
  double XBar[4];

  double ExpArg[2] = {0.0, 0.0};
  for(int d=0; d<D; d++)
   { XBar[d]=X[d]-X0[d];
     ExpArg[0] += Alpha[d]*XBar[d];
     ExpArg[1] += Alpha[d]*XBar[d]*XBar[d];
   };

  int NVD=(1<<D);

  double ExpFac[2];
  ExpFac[0] = exp(-ExpArg[0]);
  ExpFac[1] = exp(-ExpArg[1]);

  iVec Twos(D,2);
  LOOP_OVER_IVECS(nVD, sigmaVec, Twos)
   { double Factor[2]={1.0, 1.0};
     for(int d=0; d<D; d++)
      if (sigmaVec[d]==1)
       { Factor[0]*=-1.0*Alpha[d];
         Factor[1]*=-2.0*Alpha[d]*XBar[d];
       };
     PhiVD[0*NVD + nVD] = Factor[0]*ExpFac[0];
     PhiVD[1*NVD + nVD] = Factor[1]*ExpFac[1];
   }

  static double C[4][4];
  static bool CInitialized=false;
  if (!CInitialized)
   { CInitialized=true;
     srandom(0);
     for(int n0=0; n0<4; n0++)
      for(int n1=0; n1<4; n1++)
       C[n0][n1]=randU(-2.0,2.0);
   }
  int nf=2;
  memset(PhiVD + nf*NVD, 0, NVD*sizeof(double));
  for(int n0=0; n0<4; n0++)
   for(int n1=0; n1<4; n1++)
    { PhiVD[nf*NVD + 0] += C[n0][n1] * pow(X[0],n0) * pow(X[1],n1);
      PhiVD[nf*NVD + 1] += (n0==0 ? 0.0 : n0*C[n0][n1] * pow(X[0],n0-1) * pow(X[1],n1));
      PhiVD[nf*NVD + 2] += (n1==0 ? 0.0 : n1*C[n0][n1] * pow(X[0],n0) * pow(X[1],n1-1));
      PhiVD[nf*NVD + 3] += (n0==0 || n1==0 ? 0.0 : n0*n1*C[n0][n1] * pow(X[0],n0-1) * pow(X[1],n1-1));
    }

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  InstallHRSignalHandler();
  InitializeLog(argv[0]);

  int D=2;
  dVec XMin(2);  XMin[0]=-1.0;  XMin[1]=-2.0;
  dVec XMax(2);  XMax[0]=+2.0;  XMax[1]=+3.0;
  dVec X0(2);    X0[0]=X0[1]=HUGE_VAL;
  double RelTol=1.0e-4;
  bool Verbose=false;
  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { {"XMin",           PA_DOUBLE,  2, 1,       (void *)&(XMin[0]),  0,             ""},
     {"XMax",           PA_DOUBLE,  2, 1,       (void *)&(XMax[0]),  0,             ""},
     {"X0",             PA_DOUBLE,  2, 1,       (void *)&(X0[0]),    0,             ""},
     {"RelTol",         PA_DOUBLE,  1, 1,       (void *)&RelTol,     0,             ""},
     {"Verbose",        PA_BOOL,    0, 1,       (void *)&Verbose,    0,             ""},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  InterpND *Interp = new InterpND(Phi, (void *)&D, NFUN, XMin, XMax, RelTol, Verbose);
  printf("{%lu,%lu} points\n",Interp->XGrids[0].size(),Interp->XGrids[1].size());

  if (!isinf(X0[0]))
   { 
     int NFVD=NFUN*Interp->NVD;
     double *PhiVDExact = new double[NFVD], *PhiVDInterp = new double[NFVD];
     Phi(&(X0[0]), (void *)&D, PhiVDExact);
     Interp->EvaluateVD(&(X0[0]),PhiVDInterp);
     Compare(PhiVDExact, PhiVDInterp, NFVD, "Exact", "Interp");
   }
  double Error=Interp->PlotInterpolationError(Phi, (void *)&D, "/tmp/tInterpND.out");
  printf("MaxError = %e\n",Error);

}
