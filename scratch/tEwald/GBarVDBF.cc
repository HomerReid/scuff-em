/*
 *
 */
#include <stdlib.h>
#include <math.h>

#include "libhrutil.h"

#define II cdouble(0,1)
 
#define NSUM 8
#define NFIRSTROUND 1
#define NMAX 10000

/***************************************************************/
/***************************************************************/
/***************************************************************/
namespace scuff{
void AddGBFContribution(double R[3], cdouble k, double P[2],
                        double Lx, double Ly, cdouble *Sum);
}
using namespace scuff;

/***************************************************************/
/* compute \sum_L e^{iP\dot L} G(r+L)                          */
/* where G(r)=exp(i*k*|r|) / (4*pi*|r|)                     */
/***************************************************************/
void GBarVDBF(double *R, cdouble k, double *P, double **LBV, int RetainFirst9,
              double AbsTol, double RelTol, int *pnCells, cdouble *Sum)
{ 
  int nx, ny;
  double RelDelta, AbsDelta;
  cdouble LastSum[NSUM];
  double MaxRelDelta, MaxAbsDelta;
  double Delta, AbsSum;
  int i, NN, ConvergedIters;
  int nCells=0;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  memset(Sum,0,NSUM*sizeof(cdouble));
  for (nx=-NFIRSTROUND; nx<=NFIRSTROUND; nx++)
   for (ny=-NFIRSTROUND; ny<=NFIRSTROUND; ny++, nCells++)
    { 
      if ( RetainFirst9==0 && (abs(nx)<=1) && (abs(ny)<=1) )
       continue; // skip the innermost 9 grid cells unless RetainFirst9==1 

      AddGBFContribution(R, k, P, 
                         nx*LBV[0][0] + ny*LBV[1][0],
                         nx*LBV[0][1] + ny*LBV[1][1],
                         Sum);
    };
        
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  memcpy(LastSum,Sum,NSUM*sizeof(cdouble));
  ConvergedIters=0;
  for(NN=NFIRSTROUND+1; ConvergedIters<3 && NN<=NMAX; NN++)
   {  
     for(nx=-NN; nx<=NN; nx++)
      for(ny=-NN; ny<=NN; ny++)
       { 
         if ( (abs(nx)<NN) && (abs(ny)<NN) )
          continue;

         AddGBFContribution(R, k, P, 
                            nx*LBV[0][0] + ny*LBV[1][0],
                            nx*LBV[0][1] + ny*LBV[1][1],
                            Sum);

         nCells++;
        };

     /*--------------------------------------------------------------*/
     /* convergence analysis ----------------------------------------*/
     /*--------------------------------------------------------------*/
     MaxAbsDelta=MaxRelDelta=0.0;
     for(i=0; i<NSUM; i++)
      { Delta=abs(Sum[i]-LastSum[i]);
        if ( Delta>MaxAbsDelta )
         MaxAbsDelta=Delta;
        AbsSum=abs(Sum[i]);
        if ( AbsSum>0.0 && (Delta > MaxRelDelta*AbsSum) )
         MaxRelDelta=Delta/AbsSum;
      };
     if ( MaxAbsDelta<AbsTol || MaxRelDelta<RelTol )
      ConvergedIters++;
     else
      ConvergedIters=0;

     memcpy(LastSum,Sum,NSUM*sizeof(cdouble));

   };

  if (pnCells) *pnCells=nCells; 

}
