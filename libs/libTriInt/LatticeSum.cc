/* 
 * LatticeSum.cc -- evaluate a sum of the form
 *
 *  \sum_U f(U)
 *
 * where f(U) is a user-supplied vector-valued function and where
 * U ranges over all points in a 1D or 2D lattice.
 *
 * Homer Reid 2011--2015
 */
#include <stdlib.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "libhmat.h"
#define II cdouble(0,1)

/***************************************************************/
/* type definition for user-supplied summand function.         */
/* note that SummandFunction should ACCUMULATE the             */
/* contribution of lattice point U to the sum, i.e. it should  */
/* implement something like                                    */
/*  Sum[ns] += f_{ns}[U]                                       */
/* instead of                                                  */
/*  Sum[ns] = f_{ns}[U].                                       */
/***************************************************************/
typedef void (*SummandFunction)(double *U, void *UserData, double *Sum);

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddContribution(SummandFunction AddSummand, void *UserData,
                     int n1, int n2, HMatrix *LBasis, double *Sum)
{ 
  double U[3];
  U[0] = n1*LBasis->GetEntryD(0,0);
  U[1] = n1*LBasis->GetEntryD(1,0);
  U[2] = n1*LBasis->GetEntryD(2,0);
  if (LBasis->NC >=2 )
   { U[0] += n2*LBasis->GetEntryD(0,1);
     U[1] += n2*LBasis->GetEntryD(1,1);
     U[2] += n2*LBasis->GetEntryD(2,1);
   };
  AddSummand(U, UserData, Sum);
}

/***************************************************************/
/* NSum is the dimension of the vector-valued summand function.*/
/* LBasis is a 3xD dimensional vector whose columns are the    */
/* lattice basis vectors.                                      */
/* Return value is number of lattice points summed.            */
/***************************************************************/
int GetLatticeSum(SummandFunction Summand, void *UserData, int nSum,
                  HMatrix *LBasis, double *Sum,
                  double AbsTol, double RelTol, int MaxCells)
{ 
  memset(Sum,0,nSum*sizeof(double));

  /***************************************************************/
  /* start by summing the contributions of a first round         */
  /* of cells near the origin                                    */
  /***************************************************************/
  int LDim = LBasis->NC;
  int nCells = 0;
#define NFIRSTROUND 5
  if (LDim==1)
   for (int n=-NFIRSTROUND; n<=NFIRSTROUND; n++, nCells++)
    AddContribution(Summand, UserData, n, 0, LBasis, Sum);
  else // LDim==2
   for (int n1=-NFIRSTROUND; n1<=NFIRSTROUND; n1++)
    for (int n2=-NFIRSTROUND; n2<=NFIRSTROUND; n2++, nCells++)
     AddContribution(Summand, UserData, n1, n2, LBasis, Sum);

  /***************************************************************/
  /* continue to add contributions of outer cells until converged*/
  /***************************************************************/
  double *LastSum = new double[nSum];
  memcpy(LastSum,Sum,nSum*sizeof(double));
  int ConvergedIters=0;
  for(int NN=NFIRSTROUND+1; ConvergedIters<3 && nCells<MaxCells; NN++)
   {  
     if (LDim==1)
      { AddContribution(Summand, UserData,  NN, 0, LBasis, Sum);
        AddContribution(Summand, UserData, -NN, 0, LBasis, Sum);
        nCells+=2;
      }
     else // LDim==2
      { 
        /*--------------------------------------------------------------*/
        /* sum the contributions of the outer perimeter of the innermost*/
        /* NNxNN square of grid cells                                   */
        /*--------------------------------------------------------------*/
        for(int n=-NN; n<NN; n++)
         { AddContribution(Summand, UserData,   n,  NN, LBasis, Sum);
           AddContribution(Summand, UserData,  NN,  -n, LBasis, Sum);
           AddContribution(Summand, UserData,  -n, -NN, LBasis, Sum);
           AddContribution(Summand, UserData, -NN,   n, LBasis, Sum);
           nCells+=4;
         };
      };

     /*--------------------------------------------------------------*/
     /* convergence analysis ----------------------------------------*/
     /*--------------------------------------------------------------*/
     double MaxRelDelta=0.0, MaxAbsDelta=0.0;
     for(int ns=0; ns<nSum; ns++)
      { double Delta=fabs(Sum[ns]-LastSum[ns]);
        MaxAbsDelta=fmax(Delta, MaxAbsDelta);
        double AbsSum=fabs(Sum[ns]);
        if ( AbsSum>0.0 && (Delta > MaxRelDelta*AbsSum) )
         MaxRelDelta=Delta/AbsSum;
      };
     if ( MaxAbsDelta<AbsTol || MaxRelDelta<RelTol )
      ConvergedIters++;
     else
      ConvergedIters=0;
     memcpy(LastSum,Sum,nSum*sizeof(double));

   };

  delete[] LastSum;

  return nCells;

}
