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
void AddGBFContribution(cdouble k, double P[2], double R[3],
                        double Lx, double Ly, cdouble *Sum);

/***************************************************************/
/* compute \sum_L e^{iP\dot L} G(r+L)                          */
/* where G(r)=exp(i*k*|r|) / (4*pi*|r|)                     */
/***************************************************************/
void GBarVDBF(cdouble k, double *P, double *L1, double *L2, double *R, 
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
#if 1
      if ( (abs(nx)<=1) && (abs(ny)<=1) )
       continue; // skip the innermost 9 grid cells 
#endif

      AddGBFContribution(k, P, R, 
                         nx*L1[0] + ny*L2[0],
                         nx*L1[1] + ny*L2[1],
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

         AddGBFContribution(k, P, R, 
                            nx*L1[0] + ny*L2[0],
                            nx*L1[1] + ny*L2[1], Sum);

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
