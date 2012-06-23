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
#include "libscuffInternals.h"

using namespace scuff;

#define ABSTOL 1.0e-8    // absolute tolerance
#define RELTOL 1.0e-4    // relative tolerance

#define II cdouble(0.0,1.0)

bool VecEqualFloat(const double *a, const double *b) 
{
   return (     float(a[0]) == float(b[0])
            &&  float(a[1]) == float(b[1])
            &&  float(a[2]) == float(b[2])
          );
}


/***************************************************************/
/* data structure used to pass data to integrand routines      */
/***************************************************************/
typedef struct EPIBFData
 { 
   double *V0, A[3], B[3];
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

  X[0] = EPIBFD->V[0] + u*EPIBFD->A[0] + v*EPIBFD->B[0];
  X[1] = EPIBFD->V[1] + u*EPIBFD->A[1] + v*EPIBFD->B[1];
  X[2] = EPIBFD->V[2] + u*EPIBFD->A[2] + v*EPIBFD->B[2];

  XP[0] = EPIBFD->L0[0] + up*(EPIBFD->L1[0] - EPIBFD->L0[0]);
  XP[1] = EPIBFD->L0[1] + up*(EPIBFD->L1[1] - EPIBFD->L0[1]);
  XP[2] = EPIBFD->L0[2] + up*(EPIBFD->L1[2] - EPIBFD->L0[2]);

  VecSub(X, XP, R);
  r=VecNorm(R);

  cdouble ik=II*EPIBFD->k;
  cdouble Phi=exp(ikr)/(4.0*M_PI*r);

  cdouble *zf=(cdouble *)fval;
  zf[0] = u*Phi;
} 

/***************************************************************/
/* compute edge-panel integrals using brute-force technique    */
/* (adaptive cubature over the panel and the edge).            */
/***************************************************************/
void GetEPIs_BruteForce(GetEPIArgStruct *Args, int PlotFits)
{ 
  /***************************************************************/
  /* extract fields from argument structure **********************/
  /***************************************************************/
  RWGObject *Oa = Args->Oa;
  int npa = Args->npa;
  RWGPanel *Pa = Oa->Panels[npa];

  RWGObject *Ob = Args->Ob;
  int iL1 = Args->iL1;
  int iL2 = Args->iL2;

  /***************************************************************/
  /* detect common vertices **************************************/
  /***************************************************************/
  double VA = Oa->Vertices + 3*Pa->
  double VL1=
  int ncv=0;
  if ( Oa==Ob )
   { if ( VeciL1==P->VI[0] || iL1==P->VI[1] || iL1==P->VI[2]) 
      ncv++;
     if (iL2==P->VI[0] || iL2==P->VI[1] || iL2==P->VI[2]) 
      ncv++;
   };

  /***************************************************************/
  /* setup for call to cubature routine    ***********************/
  /***************************************************************/
  EPIBFData MyEPIBFData, *EPIBFD=&MyEPIBFData;
 
  memcpy(EPIBFD->V0, Vam, 
  VecSub(Va[1], Va[0], EPIBFD->A);
  VecSub(Va[2], Va[1], EPIBFD->B);
  EPIBFD->Q  = Qa;

  // note that for Va[0] and Qb we make copies of the 
  // entries (not just the pointers) because we may need
  // to displace them, see below.
  memcpy(EPIBFD->V0P,Vb[0],3*sizeof(double));
  VecSub(Vb[1], Vb[0], EPIBFD->AP);
  VecSub(Vb[2], Vb[1], EPIBFD->BP);
  memcpy(EPIBFD->QP,Qb,3*sizeof(double));
   
  EPIBFD->k = k;
  EPIBFD->NeedGradient = NeedGradient;

  double Lower[4]={0.0, 0.0, 0.0, 0.0};
  double Upper[4]={1.0, 1.0, 1.0, 1.0};

  int fDim=(EPIBFD->NeedGradient) ? 16 : 4;
  double Result[fDim], Error[fDim];

  /***************************************************************/
  /* switch off based on whether or not there are any common     */
  /* vertices                                                    */
  /***************************************************************/
  if (ncv==0)
   {
     adapt_integrate(fDim, EPIBFIntegrand, (void *)EPIBFD, 4, Lower, Upper,
                     0, ABSTOL, RELTOL, Result, Error);
     Args->EPI = cdouble(Result[0], Result[1]);
   }
  else
   {
     /*--------------------------------------------------------------*/
     /* if there are common vertices then we estimate the edge-panel */
     /* integrals using a limiting process in which we displace the  */
     /* panel through a distance Z in the direction of its surface   */
     /* normal and try to extrapolate to Z=0                         */
     /*--------------------------------------------------------------*/
     int nz, NZ=10;
     double Z[NZ], GR[NZ], GI[NZ];
     double DeltaZ=Oa->Panels[npa]->Radius/100.0;
     double *ZHat=Oa->Panels[npa]->ZHat;

     for(nz=0; nz<NZ; nz++)
      { 
        Z[nz]=((double)(nz+1))*DeltaZ;
        VecScaleAdd(Vb[0], Z[nz], ZHat, EPIBFD->V0P);
        VecScaleAdd(Qb,    Z[nz], ZHat, EPIBFD->QP);

 //       adapt_integrate(fDim, EPIBFIntegrand, (void *)EPIBFD, 4, Lower, Upper,
 //                       0, ABSTOL, RELTOL, Result, Error);
       adapt_integrate_log(fDim, EPIBFIntegrand, (void *)EPIBFD, 4, Lower, Upper,
                           100000, ABSTOL, RELTOL, Result, Error, "SGJC.log", 15);

        GR[nz]=Result[0];
        GI[nz]=Result[1];
        CR[nz]=Result[2];
        CI[nz]=Result[3];
      };
 
     PolyFit *PF=new PolyFit(Z, GR, NZ, 4);
     if (PlotFits) 
      PF->PlotFit(Z, GR, NZ, 0.0, Z[NZ-1], "real(<fa|G|fb>)");
     real(Args->H[0])=PF->f(0.0);
     delete PF;

     PF=new PolyFit(Z, GI, NZ, 4);
     if (PlotFits) 
      PF->PlotFit(Z, GI, NZ, 0.0, Z[NZ-1], "imag(<fa|G|fb>)");
     imag(Args->H[0])=PF->f(0.0);
     delete PF;

     PF=new PolyFit(Z, CR, NZ, 4);
     if (PlotFits) 
      PF->PlotFit(Z, CR, NZ, 0.0, Z[NZ-1], "imag(<fa|C|fb>)");
     real(Args->H[1])=PF->f(0.0);
     delete PF;

     PF=new PolyFit(Z, CI, NZ, 4);
     if (PlotFits) 
      PF->PlotFit(Z, CI, NZ, 0.0, Z[NZ-1], "imag(<fa|C|fb>)");
     imag(Args->H[1])=PF->f(0.0);
     delete PF;
     
   }; // if (ncv==0 ... else)

}
