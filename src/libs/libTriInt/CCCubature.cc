/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This file is part of SCUFF-EM.
 *
 * SCUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * SCUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * CCCubature.cc -- clenshaw-curtis quadrature routines
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

#define MAXDIM 10

/***************************************************************/
/***************************************************************/
/***************************************************************/
int CCCubature(int Order, unsigned fdim, integrand f, void *fdata,
	       unsigned dim, const double *xmin, const double *xmax, 
	       size_t maxEval, double reqAbsError, double reqRelError,
               error_norm norm, double *Integral, double *Error)
{
  if (Order==0)
   return pcubature(fdim, f, fdata, dim, xmin, xmax, maxEval,
                    reqAbsError, reqRelError, norm, Integral, Error);
 
  if (Order<0)
   return RRCubature(-Order, 0, fdim, f, fdata, dim, xmin, xmax, Integral, Error);

  double *CCQR = GetCCRule(Order);
  if (!CCQR) 
   ErrExit("invalid CCRule order (%i) in CCCubature",Order);

  if (dim>MAXDIM) 
   ErrExit("dimension too high in CCCubature");

  double uAvg[MAXDIM], uDelta[MAXDIM];
  for(unsigned d=0; d<dim; d++)
   { uAvg[d]   = 0.5*(xmax[d] + xmin[d]);
     uDelta[d] = 0.5*(xmax[d] - xmin[d]);
   };

  int ncp[MAXDIM];
  memset(ncp, 0, dim*sizeof(int));

  double *Integrand=Error;
  memset(Integral, 0, fdim*sizeof(double));
  bool Done=false;
  int nCalls=0;
  while(!Done)
   { 
     // get d-dimensional cubature point and weight
     double u[MAXDIM], w=1.0;
     for(unsigned nd=0; nd<dim; nd++)
      { u[nd]  = uAvg[nd] - uDelta[nd]*CCQR[2*ncp[nd] + 0];
           w  *=            uDelta[nd]*CCQR[2*ncp[nd] + 1];
      }; 

     // advance to next d-dimensional cubature point
     for(unsigned nd=0; nd<dim; nd++)
      { ncp[nd] = (ncp[nd]+1)%Order;
        if(ncp[nd]) break;
        if(nd==(dim-1)) Done=true;
      };

     f(dim, u, fdata, fdim, Integrand);
     VecPlusEquals(Integral, w, Integrand, fdim);
     nCalls++;
   };

  return nCalls;
}


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
   { case 1: N=3;   NSub=1;  break;
     case 2: N=5;   NSub=3;  break;
     case 3: N=9;   NSub=5;  break;
     case 4: N=17;  NSub=9;  break;
     case 5: N=33;  NSub=17; break;
     case 6: N=65;  NSub=33; break;
     case 7: N=129; NSub=65; break;
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
  int N=0;    // number of points in full cubature rule
  int NSub=0; // number of points in embedded subrule
  switch (p)
   { case 1: N=3;   NSub=1;  break;
     case 2: N=5;   NSub=3;  break;
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

/***************************************************************/
/* FDim = number of doubles in integrand vector                */
/* IDim = number of integration variables                      */
/* If Orders is non-NULL, then Orders[0,1,..] are the number   */
/* of sample points in each dimension.                         */
/* If Orders is NULL, then Order is used for all dimensions.   */
/***************************************************************/
int RRCubature(int Order, int *Orders, 
               int FDim, integrand f, void *UserData,
	       int IDim, const double *Lower, const double *Upper,
               double *Integral, double *Error)
{
  int N[MAXDIM];
  if (Orders)
   memcpy(N, Orders, IDim*sizeof(int));
  else
   for(int d=0; d<IDim; d++)
    N[d]=Order;

  for(int d=0; d<IDim; d++)
   if (N[d]<=0)
    ErrExit("invalid order %i in RRCubature",N[d]);

  double Delta[MAXDIM];
  double Weight=1.0;
  for(int d=0; d<IDim; d++)
   { Delta[d] = (Upper[d] - Lower[d]) / ((double)N[d]);
     Weight*=Delta[d];
   };

  int ncp[MAXDIM];
  memset(ncp, 0, IDim*sizeof(int));

  double *Integrand=Error;
  memset(Integral, 0, FDim*sizeof(double));
  bool Done=false;
  int nCalls=0;
  while(!Done)
   { 
     // get d-dimensional cubature point
     double x[MAXDIM];
     for(int d=0; d<IDim; d++)
      x[d] = Lower[d] + ncp[d]*Delta[d];

     // advance to next d-dimensional cubature point
     for(int d=0; d<IDim; d++)
      { ncp[d] = (ncp[d]+1)%N[d];
        if(ncp[d]) break;
        if(d==(IDim-1)) Done=true;
      };

     f(IDim, x, UserData, FDim, Integrand);
     VecPlusEquals(Integral, Weight, Integrand, FDim);
     nCalls++;
   };

  return nCalls;
}
