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
 * qIntegral.cc -- evaluate integrals over Fourier wavevectors
 *
 * homer reid   -- 3/2017-10/2017
 *
 */
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>

#include "libhrutil.h"
#include "libMDInterp.h"
#include "libSGJC.h"
#include "libSpherical.h"
#include "libSubstrate.h"

#define II cdouble(0.0,1.0)
#ifndef ZVAC
#define ZVAC 376.73031346177
#endif

#define SQRT2 1.41421356237309504880

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct qIntegrandData 
 {
   LayeredSubstrate *Substrate;
   cdouble Omega;
   double rkLayer2;
   double q0;
   bool Propagating;
   qFunction UserFunction;
   void *UserData;
   int nCalls;

 } qIntegrandData;

int qIntegrand(unsigned ndim, const double *uVector,
               void *pqIData, unsigned fdim, double *fval)
{
  qIntegrandData *Data        = (qIntegrandData *)pqIData;
  LayeredSubstrate *Substrate = Data->Substrate;
  cdouble Omega               = Data->Omega;
  double rkLayer2             = Data->rkLayer2;
  double q0                   = Data->q0;
  bool Propagating            = Data->Propagating;
  qFunction UserFunction      = Data->UserFunction;
  void *UserData              = Data->UserData;
  Data->nCalls++;

  #define UMIN 1.0e-3
  double u=uVector[0], Jacobian = 1.0, qMag;
  if (Propagating)
   { 
     if ( u<=UMIN )
      u=UMIN;
     if ( u>=(1.0-UMIN) )
      u=1.0-UMIN;
     double qz = u*q0;
     Jacobian  = q0*qz/(2.0*M_PI);
     qMag      = sqrt(rkLayer2 - qz*qz);
   }
  else // evanescent
   { 
     if (u>1.0-UMIN)
      return 0;    // integrand vanishes at qMag=iqz=infinity
     if (u<=UMIN)
      u=UMIN;
     double Denom = 1.0 / (1.0-u);
     qMag         = q0 + u*Denom;
     Jacobian     = qMag*Denom*Denom/(2.0*M_PI);
   };

  double q2D[2];
  if (ndim==2)
   { double Theta = uVector[1];
     q2D[0] = qMag*cos(Theta);
     q2D[1] = qMag*sin(Theta);
     Jacobian/=(2.0*M_PI);
   }
  else
   { q2D[0] = qMag;
     q2D[1] = 0.0;
   };
 
  cdouble *Integrand = (cdouble *)fval;
  memset(Integrand, 0, fdim*sizeof(double));
  UserFunction(Substrate, q2D, Omega, UserData, Integrand);
  VecScale(Integrand, Jacobian, fdim/2);

  return 0;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void LayeredSubstrate::qIntegrate(cdouble Omega,
                                  qFunction UserFunction,
                                  void *UserData,
                                  cdouble *Integral, int FDim,
                                  bool ThetaIndependent)
{
  double Lower[2]={0.0,0.0}, Upper[2]={1.0,2.0*M_PI};
  int IDim = ThetaIndependent ? 1 : 2;
  qIntegrandData MyData, *Data=&MyData;
  Data->Substrate    = this;
  Data->Omega        = Omega;
  Data->UserFunction = UserFunction;
  Data->UserData     = UserData;

  // compute wavevectors in all layers and sort by real part
  UpdateCachedEpsMu(Omega);
  HMatrix *kLayer   = new HMatrix(NumLayers, 2, LHM_REAL);
  for(int n=0; n<NumLayers; n++)
   { cdouble k=sqrt(EpsLayer[n]*MuLayer[n])*Omega;
     kLayer->SetEntry(n, 0, real(k));
     kLayer->SetEntry(n, 1, imag(k));
   };
  kLayer->Sort(0);

  cdouble *IPartial = new cdouble[FDim];
  cdouble *Error    = new cdouble[FDim];

  /***************************************************************/
  /* evaluate 'propagating' integrals                            */
  /***************************************************************/
  memset(Integral, 0, FDim*sizeof(cdouble));
  for(int n=0; n<NumLayers; n++)
   {
     if (LayerOnly!=-1 && LayerOnly!=n)
      continue;

     double knm1       = (n==0) ? 0.0 : kLayer->GetEntryD(n-1,0);
     double kn         = kLayer->GetEntryD(n,0);
     if (kn==knm1)
      continue;
     Data->rkLayer2    = kn*kn;
     Data->q0          = sqrt(kn*kn - knm1*knm1);
     Data->nCalls      = 0;
     Data->Propagating = true;

     memset(Times,0,NUMTIMES*sizeof(double));
     hcubature(2*FDim, qIntegrand, (void *)Data, IDim,
               Lower, Upper, qMaxEval, qAbsTol, qRelTol, 
               ERROR_PAIRED, (double *)IPartial, (double *)Error);
     VecPlusEquals(Integral, 1.0, IPartial, FDim);

     if (LogLevel>=LIBSUBSTRATE_VERBOSE)
      Log("Integral %i [%g,%g]: %i calls {%e,%e} {%e,%e}", n+1,knm1,kn,Data->nCalls,
           real(IPartial[0]),imag(IPartial[0]),real(IPartial[FDim-1]),imag(IPartial[FDim-1]));

     if (LogLevel>=LIBSUBSTRATE_VERBOSE2)
      for(int nt=0; nt<NUMTIMES; nt++)
       printf("%10s  %f \n",TimeNames[nt],Times[nt]);
   };

  /***************************************************************/
  /* evaluate 'evanescent' integral                              */
  /***************************************************************/
  if (LayerOnly==-1 || LayerOnly==NumLayers)
   { 
     Data->nCalls      = 0;
     Data->Propagating = false;
     Data->q0          = kLayer->GetEntryD(NumLayers-1,0);
     memset(Times,0,NUMTIMES*sizeof(double));

     hcubature(2*FDim, qIntegrand, (void *)Data, IDim,
               Lower, Upper, qMaxEval, qAbsTol, qRelTol,
               ERROR_PAIRED, (double *)IPartial, (double *)Error);
     VecPlusEquals(Integral, 1.0, IPartial, FDim);

     if (LogLevel>=LIBSUBSTRATE_VERBOSE)
      Log("Integral %i [%g,inf]: %i calls {%e,%e} {%e,%e}",NumLayers,Data->q0,Data->nCalls,
           real(IPartial[0]),imag(IPartial[0]), real(IPartial[FDim-1]),imag(IPartial[FDim-1]));

     if (LogLevel>=LIBSUBSTRATE_VERBOSE2)
      for(int nt=0; nt<NUMTIMES; nt++)
       printf("%10s  %f \n",TimeNames[nt],Times[nt]);
   };
  
  delete[] Error;
  delete[] IPartial;
  delete kLayer;
  
}
