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
 * GetEPIs_BruteForce.cc -- compute edge-panel integrals by brute force
 * 
 * homer reid -- 6/2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>
#include <libSGJC.h>
#include <libPolyFit.h>

#include "libscuff.h"

using namespace scuff;

#define ABSTOL 1.0e-8    // absolute tolerance
#define RELTOL 1.0e-4    // relative tolerance

#define II cdouble(0.0,1.0)

/***************************************************************/
/* data structure used to pass data to integrand routines      */
/***************************************************************/
typedef struct EPIBFData
 { 
   double V0[3], A[3], B[3];
   double *L1, *L2;
   cdouble k;
 } EPIBFData;

/***************************************************************/
/***************************************************************/
/***************************************************************/
static void EPIBFIntegrand(unsigned ndim, const double *x, void *parms,
                           unsigned nfun, double *fval)
{
  EPIBFData *EPIBFD=(EPIBFData *)parms;
  
  double u=x[0];
  double v=u*x[1];
  double up=x[2];

  double X[3], XP[3], R[3];
  double r;

  X[0] = EPIBFD->V0[0] + u*EPIBFD->A[0] + v*EPIBFD->B[0];
  X[1] = EPIBFD->V0[1] + u*EPIBFD->A[1] + v*EPIBFD->B[1];
  X[2] = EPIBFD->V0[2] + u*EPIBFD->A[2] + v*EPIBFD->B[2];

  XP[0] = EPIBFD->L1[0] + up*(EPIBFD->L2[0] - EPIBFD->L1[0]);
  XP[1] = EPIBFD->L1[1] + up*(EPIBFD->L2[1] - EPIBFD->L1[1]);
  XP[2] = EPIBFD->L1[2] + up*(EPIBFD->L2[2] - EPIBFD->L1[2]);

  VecSub(X, XP, R);
  r=VecNorm(R);

  cdouble ik=II*EPIBFD->k;
  cdouble Phi=exp(ik*r)/(4.0*M_PI*r);

  cdouble *zf=(cdouble *)fval;
  zf[0] = u*Phi;

} 

/***************************************************************/
/* compute edge-panel integrals using brute-force technique    */
/* (adaptive cubature over the panel and the edge).            */
/***************************************************************/
cdouble GetEPI_BruteForce(double **V, double **L, cdouble k, int PlotFits)
{ 
  /***************************************************************/
  /* detect common vertices **************************************/
  /***************************************************************/
  int ncv=0;
  if ( VecEqualFloat(L[0], V[0]) || VecEqualFloat(L[0], V[1]) || VecEqualFloat(L[0],V[2]))
   ncv++;
  if ( VecEqualFloat(L[1], V[0]) || VecEqualFloat(L[1], V[1]) || VecEqualFloat(L[1],V[2]))
   ncv++;

  /***************************************************************/
  /* setup for call to cubature routine    ***********************/
  /***************************************************************/
  EPIBFData MyEPIBFData, *EPIBFD=&MyEPIBFData;
 
  memcpy(EPIBFD->V0, V[0], 3*sizeof(double));
  VecSub(V[1], V[0], EPIBFD->A);
  VecSub(V[2], V[1], EPIBFD->B);
  EPIBFD->L1=L[0];
  EPIBFD->L2=L[1];

  EPIBFD->k = k;

  double Lower[3]={0.0, 0.0, 0.0};
  double Upper[3]={1.0, 1.0, 1.0};

  int fDim=2;
  cdouble Result, Error;

  /***************************************************************/
  /* switch off based on whether or not there are any common     */
  /* vertices                                                    */
  /***************************************************************/
  if (ncv==0)
   {
     adapt_integrate(fDim, EPIBFIntegrand, (void *)EPIBFD, 3, Lower, Upper,
                     0, ABSTOL, RELTOL, (double *)&Result, (double *)&Error);
     return Result;
   }
  else
   {
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     double ZHat[3];
     VecCross(EPIBFD->A, EPIBFD->B, ZHat);
     VecNormalize(ZHat);

     /*--------------------------------------------------------------*/
     /* if there are common vertices then we estimate the edge-panel */
     /* integrals using a limiting process in which we displace the  */
     /* panel through a distance Z in the direction of its surface   */
     /* normal and try to extrapolate to Z=0                         */
     /*--------------------------------------------------------------*/
     int nz, NZ=10;
     double Z[NZ], GR[NZ], GI[NZ];
     double DeltaZ=VecDistance(L[1],L[0])/100.0;

     for(nz=0; nz<NZ; nz++)
      { 
        Z[nz]=((double)(nz+1))*DeltaZ;
        VecScaleAdd(V[0], Z[nz], ZHat, EPIBFD->V0);

       adapt_integrate(fDim, EPIBFIntegrand, (void *)EPIBFD, 3, Lower, Upper,
                       100000, ABSTOL, RELTOL, (double *)&Result, (double *)&Error);
 
        GR[nz]=real(Result);
        GI[nz]=imag(Result);
      };
 
     PolyFit *PF=new PolyFit(Z, GR, NZ, 4);
     if (PlotFits) 
      PF->PlotFit(Z, GR, NZ, 0.0, Z[NZ-1], "real(<fa|G|fb>)");
     real(Result) = PF->f(0.0);
     delete PF;

     PF=new PolyFit(Z, GI, NZ, 4);
     if (PlotFits) 
      PF->PlotFit(Z, GI, NZ, 0.0, Z[NZ-1], "imag(<fa|G|fb>)");
     imag(Result) = PF->f(0.0);
     delete PF;

     return Result;

   }; // if (ncv==0 ... else)

}
