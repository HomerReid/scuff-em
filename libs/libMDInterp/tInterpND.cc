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
void printiVec(iVec v, char *s)
{ printf("%s: {",s);
  for(size_t d=0; d<v.size(); d++) printf("%i ",v[d]);
  printf("} ");
}
void printdVec(dVec v, char *s)
{ printf("%s: {",s);
  for(size_t d=0; d<v.size(); d++) printf("%g ",v[d]);
  printf("} ");
}
/***************************************************************/
/***************************************************************/
/***************************************************************/
#define NFUN 3
void Phi(dVec X, void *UserData, double *PhiVD, iVec dXMax)
{
  (void)UserData;
  int D = X.size();
  static double Alpha[4]={0.2,0.3,0.4,0.5};
  static double X0[4]={-0.1,0.2,0.35,0.5};
  double XBar[4];

  double ExpArg[2] = {0.0, 0.0};
  for(int d=0; d<D; d++)
   { XBar[d]=X[d]-X0[d];
     ExpArg[0] += Alpha[d]*XBar[d];
     ExpArg[1] += Alpha[d]*XBar[d]*XBar[d];
   }
  
  double ExpFac[2];
  ExpFac[0] = exp(-ExpArg[0]);
  ExpFac[1] = exp(-ExpArg[1]);

  int NumVDs=1;
  for(int d=0; d<D; d++)
   NumVDs*=dXMax[d];
  //printf("NVD=%i\n",NVD);

 {
  LOOP_OVER_IVECS(nVD, dX, dXMax)
   { double Factor[2]={1.0, 1.0};
     for(int d=0; d<D; d++)
      if (dX[d])
       { Factor[0]*=-1.0*Alpha[d];
         Factor[1]*=-2.0*Alpha[d]*XBar[d];
       }
     PhiVD[0*NumVDs + nVD] = Factor[0]*ExpFac[0];
     PhiVD[1*NumVDs + nVD] = Factor[1]*ExpFac[1];
   }
 }

  int nf=2;
  static double C[4][4][4];
  static bool CInitialized=false;
  if (!CInitialized)
   { CInitialized=true;
     srandom(0);
     for(int n2=0; n2<4; n2++)
      for(int n1=0; n1<4; n1++)
       for(int n0=0; n0<4; n0++)
        C[n2][n1][n0]=randU(-2.0,2.0);
   }
  memset(PhiVD + nf*NumVDs, 0, NumVDs*sizeof(double));
  double Phi2VD[8]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  int n2Max = (D>=3 ? 4 : 1);
  int n1Max = (D>=2 ? 4 : 1);
  double XF[3]={1.0, 1.0, 1.0};
  XF[0]=X[0];
  if (D>1) XF[1]=X[1];
  if (D>2) XF[2]=X[2];
  for(int n2=0; n2<n2Max; n2++)
   for(int n1=0; n1<n1Max; n1++)
    for(int n0=0; n0<4; n0++)
     { Phi2VD[0] += C[n2][n1][n0] * pow(XF[2],n2) * pow(XF[1],n1) * pow(XF[0],n0);
       Phi2VD[1] += n0==0 ? 0.0 : (n0*C[n2][n1][n0] * pow(XF[2],n2) * pow(XF[1],n1) * pow(XF[0],n0-1));
       Phi2VD[2] += n1==0 ? 0.0 : (n1*C[n2][n1][n0] * pow(XF[2],n2) * pow(XF[1],n1-1) * pow(XF[0],n0));
       Phi2VD[3] += n0*n1==0 ? 0.0 : n0*n1*C[n2][n1][n0] * pow(XF[2],n2) * pow(XF[1],n1-1) * pow(XF[0],n0-1);
       Phi2VD[4] += n2==0 ? 0.0 : n2*C[n2][n1][n0] * pow(XF[2],n2-1) * pow(XF[1],n1) * pow(XF[0],n0);
       Phi2VD[5] += n0*n2==0 ? 0.0 : n2*n0*C[n2][n1][n0] * pow(XF[2],n2-1) * pow(XF[1],n1) * pow(XF[0],n0-1);
       Phi2VD[6] += n1*n2==0 ? 0.0 : n2*n1*C[n2][n1][n0] * pow(XF[2],n2-1) * pow(XF[1],n1-1) * pow(XF[0],n0);
       Phi2VD[7] += n0*n1*n2==0 ? 0.0 : n2*n0*n1*C[n2][n1][n0] * pow(XF[2],n2-1) * pow(XF[1],n1-1) * pow(XF[0],n0-1);
     }

  {
    LOOP_OVER_IVECS(nVD, dX, dXMax)
     { int Index = dX[0];
       if (dX.size() > 1 ) Index+=dX[1]*2;
       if (dX.size() > 2 ) Index+=dX[2]*4;
       PhiVD[2*NumVDs + nVD] = Phi2VD[Index];
     }
  }

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  InstallHRSignalHandler();
  InitializeLog(argv[0]);

  int D=3;
  dVec XMin(3);  XMin[0]=-1.0;  XMin[1]=-2.0; XMin[2]=-3.0;
  dVec XMax(3);  XMax[0]=+2.0;  XMax[1]=+3.0; XMax[2]=+4.0;
  dVec X0(3);    X0[0]=X0[1]=X0[2]=HUGE_VAL;
  double RelTol=1.0e-4;
  bool Verbose=false;
  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { {"XMin",           PA_DOUBLE,  3, 1,       (void *)&(XMin[0]),  0,             ""},
     {"XMax",           PA_DOUBLE,  3, 1,       (void *)&(XMax[0]),  0,             ""},
     {"X0",             PA_DOUBLE,  3, 1,       (void *)&(X0[0]),    0,             ""},
     {"RelTol",         PA_DOUBLE,  1, 1,       (void *)&RelTol,     0,             ""},
     {"Verbose",        PA_BOOL,    0, 1,       (void *)&Verbose,    0,             ""},
     {"D",              PA_INT,     1, 1,       (void *)&D,          0,             ""},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  InterpND *Interp=0;
  if (D==1)
   { dVec xMin(1, XMin[0]), xMax(1,XMax[0]);
     Interp = new InterpND(Phi, (void *)&D, NFUN, xMin, xMax, RelTol, Verbose);
   }
  else if (D==2)
   { dVec xMin(2, XMin[0]), xMax(2,XMax[0]);
     xMin[1]=XMin[1]; xMax[1]=XMax[1];
     Interp = new InterpND(Phi, (void *)&D, NFUN, xMin, xMax, RelTol, Verbose);
   }
  else
   Interp = new InterpND(Phi, (void *)&D, NFUN, XMin, XMax, RelTol, Verbose);

  printf("{%lu ",Interp->XGrids[0].size());
  if (D>=2) printf(" %lu ",Interp->XGrids[1].size());
  if (D>=3) printf(" %lu ",Interp->XGrids[2].size());
  printf(" } points\n");
for(int dd=0; dd<Interp->XGrids.size(); dd++)
 { printdVec(Interp->XGrids[dd],"Interp->XGrids[dd]");
   printf("\n");
}

  if (!isinf(X0[0]))
   { 
     int NFVD=NFUN*Interp->NumVDs;
     double *PhiVDExact = new double[NFVD], *PhiVDInterp = new double[NFVD];
     iVec dXMax(D,1);
     dVec x0(D,X0[0]);
     if (D>1) x0[1]=X0[1];
     if (D>2) x0[2]=X0[2];
     Phi(x0, (void *)&D, PhiVDExact, dXMax);
     //Interp->EvaluateVD(&(X0[0]),PhiVDInterp);
     bool Status=Interp->Evaluate(&(X0[0]),PhiVDInterp);
     printf("point %s grid\n",Status ? "inside" : "outside");
     Compare(PhiVDExact, PhiVDInterp, NFUN, "Exact", "Interp");
   }

 // double Error=Interp->PlotInterpolationError(Phi, (void *)&D, "/tmp/tInterpND.out");
 //printf("MaxError = %e\n",Error);

}
