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
  int D=2;
  int NF=2;
  dVec X0Min(2), X0Max(2);
  X0Min[0] = -1.0;  X0Max[0] = 2.0;
  X0Min[1] = -5.0;  X0Max[1] = 5.0;
  double Error;
  Error=GetInterpolationError(Phi, (void *) &D, NF, 0, 0.1, 0.1, X0Min, X0Max, 1.0e-8, "/tmp/tInterpND.log");
  printf("d=0: Error=%e\n",Error);
  Error=GetInterpolationError(Phi, (void *)&D, NF,  1, 1.1, 0.1, X0Min, X0Max, 1.0e-8, "/tmp/tInterpND.log");
  printf("d=1: Error=%e\n",Error);
}
