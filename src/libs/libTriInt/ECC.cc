/*
 * ECC.cc -- embedded clenshaw-curtis cubature in one and two dimensions
 *
 * homer reid  -- 6/2014
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

#include <libhrutil.h>
#include <libSGJC.h>
#include <libTriInt.h>

/***************************************************************/
/* embedded clenshaw-curtis cubature in two dimensions.        */
/* p is an integer in the range {2,3,4,5,6}.                   */
/* the code estimates the integral of Integrand, over the      */
/* rectangle with lower-left corner xMin and upper-right       */
/* corner xMax, using nested Clenshaw-Curtis cubature with     */
/* NF points per dimension and NC points per dimension (where  */
/* NF=2^p+1 and NC=2^{p-1}+1) and returns the difference as    */
/* the error.                                                  */
/*                                                             */
/* AllValues points to a user-allocated buffer with space for  */
/* nFun*NF*NF doubles.                                         */
/*                                                             */
/* If OuterValues is nonzero, it points to a buffer containing */
/* function values at the ``outer'' grid points (that is, the  */
/* points common to both the coarse and fine grids). The       */
/* buffer passed as AllValues for a run at order p may later   */
/* be passed as OuterValues for a run at order p+1             */
/*                                                             */
/* If xySymmetric is true, we assume f(x,y) = f(y,x).          */
/***************************************************************/
void ECC2D(int p, double xMin[2], double xMax[2],
           integrand Integrand, void *UserData, int nFun,
           bool xySymmetric,
           double *AllValues, double *OuterValues,
           double *Result, double *Error)
{ 
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int N=0;    // number of points in full cubature rule
  int NSub=0; // number of points in embedded subrule
  switch (p)
   { case 2: N=5;   NSub=3;  break;
     case 3: N=9;   NSub=5;  break;
     case 4: N=17;  NSub=9;  break;
     case 5: N=33;  NSub=17; break;
     case 6: N=65;  NSub=33; break;
     default: ErrExit("unsupported cubature order in ECC2D");
   };
  double *FullQR=GetCCRule(N);
  double *OuterQR=GetCCRule(NSub);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double xAvg[2], xDelta[2];
  xAvg[0] = 0.5*(xMax[0]+xMin[0]);   xDelta[0] = 0.5*(xMax[0]-xMin[0]);
  xAvg[1] = 0.5*(xMax[1]+xMin[1]);   xDelta[1] = 0.5*(xMax[1]-xMin[1]);
  double Jacobian = xDelta[0]*xDelta[1];

  /*--------------------------------------------------------------*/
  /*- populate the table of AllValues. ---------------------------*/
  /*- 'nav' and 'nov' are indices that run sequentially through  -*/
  /*- the sets of all values and outer values, respectively.     -*/
  /*--------------------------------------------------------------*/
  for(int nav=0, nov=0, n1=0; n1<N; n1++)
   for(int n2=0; n2<N; n2++, nav++)
    { 
      if (OuterValues && (n1%2)==0 && (n2%2)==0 )
       { memcpy( AllValues + nFun*nav, OuterValues + nFun*nov, nFun*sizeof(double) );
         nov++;
       }
      else if (xySymmetric && n1>n2)
       { int navPrime = n2*N + n1;
         memcpy( AllValues + nFun*nav, AllValues + nFun*navPrime, nFun*sizeof(double) );
       }
      else
       { double x[2];
         x[0] = xAvg[0] + xDelta[0]*FullQR[2*n1 + 0];
         x[1] = xAvg[1] + xDelta[1]*FullQR[2*n2 + 0];
         Integrand(2, x, UserData, nFun, AllValues + nFun*nav );
       };
    };

  /*--------------------------------------------------------------*/
  /*- evaluate the two cubature rules.----------------------------*/
  /*--------------------------------------------------------------*/
  memset(Result, 0, nFun*sizeof(double));
  memset(Error,  0, nFun*sizeof(double));
  for(int nav=0, n1=0; n1<N; n1++)
   for(int n2=0; n2<N; n2++, nav++)
    { 
      if ( (n1%2)==0 && (n2%2)==0 )
       { double Weight = OuterQR[n1 + 1]*OuterQR[n2 + 1]*Jacobian;
         for(int nf=0; nf<nFun; nf++) 
          Error[nf] += Weight*AllValues[nFun*nav + nf];
       };

      double Weight = FullQR[2*n1 + 1]*FullQR[2*n2 + 1]*Jacobian;
      for(int nf=0; nf<nFun; nf++) 
       Result[nf] += Weight*AllValues[nFun*nav + nf];
    };
 
  for(int nf=0; nf<nFun; nf++)
   Error[nf] = fabs(Error[nf] - Result[nf]);

}

/***************************************************************/
/* embedded clenshaw-curtis cubature in one dimension.         */
/* p is an integer in the range {2,3,4,5,6}.                   */
/* the code estimates the integral of Integrand, over the      */
/* interval [xMin, Max], using nested Clenshaw-Curtis cubature */
/* with NF points and NC points (where NF=2^p+1 and            */
/* NC=2^{p-1}+1) and returns the difference as the error.      */
/*                                                             */
/* AllValues points to a user-allocated buffer with space for  */
/* nFun*NF doubles.                                            */
/*                                                             */
/* If OuterValues is nonzero, it points to a buffer containing */
/* function values at the ``outer'' grid points (that is, the  */
/* points common to both the coarse and fine grids). The       */
/* buffer passed as AllValues for a run at order p may later   */
/* be passed as OuterValues for a run at order p+1.            */
/***************************************************************/
void ECC(int p, double xMin, double xMax,
         integrand Integrand, void *UserData, int nFun,
         double *AllValues, double *OuterValues,
         double *Result, double *Error)
{ 
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int N;    // number of points in full cubature rule
  int NSub; // number of points in embedded subrule
  switch (p)
   { case 2: N=5;   NSub=3;  break;
     case 3: N=9;   NSub=5;  break;
     case 4: N=17;  NSub=9;  break;
     case 5: N=33;  NSub=17; break;
     case 6: N=65;  NSub=33; break;
     default: ErrExit("unsupported cubature order in ECC2D");
   };
  double *FullQR=GetCCRule(N);
  double *OuterQR=GetCCRule(NSub);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double xAvg   = 0.5*(xMax+xMin);
  double xDelta = 0.5*(xMax-xMin);
  double Jacobian = xDelta;

  /*--------------------------------------------------------------*/
  /*- populate the table of AllValues. ---------------------------*/
  /*- 'nav' and 'nov' are indices that run sequentially through  -*/
  /*- the sets of all values and outer values, respectively.     -*/
  /*--------------------------------------------------------------*/
  for(int nov=0, n=0; n<N; n++)
   { 
     if (OuterValues && (n%2)==0 )
      { memcpy( AllValues + nFun*n, OuterValues + nFun*nov, nFun*sizeof(double) );
        nov++;
      }
     else
      { double x = xAvg + xDelta*FullQR[2*n + 0];
        Integrand(1, &x, UserData, nFun, AllValues + nFun*n );
      };
   };

  /*--------------------------------------------------------------*/
  /*- evaluate the two cubature rules.----------------------------*/
  /*--------------------------------------------------------------*/
  memset(Result, 0, nFun*sizeof(double));
  memset(Error,  0, nFun*sizeof(double));
  for(int n=0; n<N; n++)
   { 
     if ( (n%2)==0 )
      { double Weight = OuterQR[n + 1]*Jacobian;
        for(int nf=0; nf<nFun; nf++) 
         Error[nf] += Weight*AllValues[nFun*n + nf];
      };

     double Weight = FullQR[2*n + 1]*Jacobian;

     for(int nf=0; nf<nFun; nf++) 
      Result[nf] += Weight*AllValues[nFun*n + nf];
   };
 
  for(int nf=0; nf<nFun; nf++)
   Error[nf] = fabs(Error[nf] - Result[nf]);
}
