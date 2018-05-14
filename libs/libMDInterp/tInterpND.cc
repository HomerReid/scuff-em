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
#define NFUN 2
#define X10 0.2
#define X20 0.3
#define X30 0.4
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
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  InstallHRSignalHandler();
  InitializeLog(argv[0]);

  int D=2;
  int NF=2;
  dVec XMin(2);  XMin[0]=-1.0;  XMin[1]=-2.0;
  dVec XMax(2);  XMax[1]=+2.0;  XMax[2]=+3.0;
  double RelTol=1.0e-4;
  bool Verbose=false;
  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { {"XMin",           PA_DOUBLE,  2, 1,       (void *)&(XMin[0]),  0,             ""},
     {"XMax",           PA_DOUBLE,  2, 1,       (void *)&(XMax[0]),  0,             ""},
     {"RelTol",         PA_DOUBLE,  1, 1,       (void *)&RelTol,     0,             ""},
     {"Verbose",        PA_BOOL,    0, 1,       (void *)&Verbose,    0,             ""},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  InterpND *Interp = new InterpND(Phi, (void *)&D, NF, XMin, XMax, RelTol, Verbose);
  double Error=Interp->PlotInterpolationError(Phi, (void *)&D, "/tmp/tInterpND.out");
  printf("{%lu,%lu} points, maxerr=%e\n",Interp->XGrids[0].size(),Interp->XGrids[1].size(),Error);
}
