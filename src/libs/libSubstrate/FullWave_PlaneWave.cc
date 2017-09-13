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
 * FullWave_PWDecomposition.cc -- compute the full-wave green's function
 *                             -- above a multi-layer dielectric substrate
 *                             -- using the plane-wave-decomposition method
 *
 * homer reid                  -- 3/2017-9/2017
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

/***************************************************************/
/* r[Mu][Nu] = amplitude of Mu-polarized upward-traveling wave */
/*             due to unit-strength Nu-polarized downward-traveling*/
/*             wave impinging on uppermost substrate surface   */
/***************************************************************/
void Substrate::GetReflectionCoefficients(double Omega, double *q2D,
                                          cdouble r[2][2])
{
  UpdateCachedEpsMu(Omega);

  // assume substrate does not mix polarizations
  r[POL_TE][POL_TM] = r[POL_TM][POL_TE] = 0.0; 

  /***************************************************************/
  /* handle the simple case of Fresnel scattering from a single  */
  /* semi-infinite dielectric half-space                         */
  /***************************************************************/
  if (NumLayers==0 && zGP!=HUGE_VAL)
   { r[POL_TE][POL_TE] = -1.0;
     r[POL_TM][POL_TM] = +1.0;
     return;
   }
  else if (NumLayers==1 && zGP==HUGE_VAL)
   { 
     cdouble k2Above = EpsMedium*Omega*Omega;
     cdouble k2Below = EpsLayer[0]*MuLayer[0]*Omega;
     double qMag2 = q2D[0]*q2D[0] + q2D[1]*q2D[1];
     cdouble qzAbove = sqrt( k2Above - qMag2 );
     cdouble qzBelow = sqrt( k2Below - qMag2 );
     cdouble Eps = EpsLayer[0] / EpsMedium;
     cdouble Mu  = MuLayer[0];
     r[POL_TE][POL_TE] = (Mu*qzAbove - qzBelow) / (Mu*qzAbove + qzBelow);
     r[POL_TM][POL_TM] = (Eps*qzAbove - qzBelow) / (Eps*qzAbove + qzBelow);
     return;
   };

  ErrExit("substrate configuration not yet supported");

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct PWDData
 {
   LayeredSubstrate *Substrate;
   HMatrix *XMatrix;
   double Omega;
   bool Polar;
   bool Accumulate;
   bool Propagating;
   int nCalls;

 } PWDData;

/***************************************************************/
/* Integrand[ 36*nx + 0*9 + 3*Mu + Nu ] = G^{EE}_{Mu,Nu}       */
/* Integrand[ 36*nx + 1*9 + 3*Mu + Nu ] = G^{EM}_{Mu,Nu}       */
/* Integrand[ 36*nx + 2*9 + 3*Mu + Nu ] = G^{ME}_{Mu,Nu}       */
/* Integrand[ 36*nx + 3*9 + 3*Mu + Nu ] = G^{MM}_{Mu,Nu}       */
/***************************************************************/
int SubstrateDGFIntegrand_PlaneWave(unsigned ndim, const double *q,
                                    void *UserData, unsigned fdim,
                                    double *fval)
{

  PWDData *Data               = (PWDData *)UserData;
  Data->nCalls++;

  LayeredSubstrate *Substrate = Data->Substrate;
  HMatrix *XMatrix            = Data->XMatrix;
  cdouble Omega               = Data->Omega;
  bool Polar                  = Data->Polar;

  int IDim = 36*XMatrix->NR;
  if (Data->Accumulate == false)
   memset(Integrand, 0, IDim*sizeof(cdouble));

  // Polar = true --> we have already integrated out 
  //                  q_Theta to yield Bessel functions,
  //                  and what we are evaluating here 
  //                  is just the integrand of the 1-dimensional
  //                  q_r integral 
  //
  // Polar = false--> we are evaluating the 2-dimensional
  //                  (qx,qy) integral
  //
  double q2, qMag, Jacobian=1.0;
  cdouble One, Cos, Sin, Cos2, Sin2, CosSin;
  if (Polar)
   { 
     if (Data->Propagating)
      qMag     = q[0]*real(k0);
     else 
      { double Denom = 1.0/(1.0-q[0]);
        qMag = real(k0) * (1.0 + q[0]*Denom);
        Jacobian = Denom*Denom;
      };
   }
  else
   { q2       = q[0]*q[0] + q[1]*q[1];
     qMag     = sqrt(q2);
     double qxHat = (qMag==0.0) ? 1.0 : q[0] / qMag;
     double qyHat = (qMag==0.0) ? 0.0 : q[1] / qMag;
     One      = 1.0;
     Cos      = qxHat;
     Sin      = qyHat;
     Cos2     = qxHat*qxHat;
     CosSin   = qxHat*qyHat;
     Sin2     = qyHat*qyHat;
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble k02 = EpsMedium*Omega*Omega;
  cdouble qz2 = k02 - q2;
  if (qz2==0.0)
   return;
  cdouble qz = sqrt(qz2);
  if ( imag(qz)<0.0 )
   qz*=-1.0;

  cdouble rr[4];
  GetReflectionCoefficients(Omega, q, rr);
  cdouble rTE = rr[POL_TE][POL_TE];
  cdouble rTM = rr[POL_TM][POL_TM];

  cdouble MTE[3][3], MTM[3][3];
  MTE[0][2] = MTE[1][2] = MTE[2][0] = MTE[2][1] = MTE[2][2] = 0.0;

  bool TwoPointDGF = (XMatrix->NC>=6);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int nx=0; nx<XMatrix->NR; nx++)
   { 
     double XDest[3], XSourceBuffer[3];
     double *XSource = (TwoPointDGF) ? XSourceBuffer : XDest;
     XMatrix->GetEntriesD(nx,"0:2",XDest);
     if (TwoPointDGF)
      XMatrix->GetEntriesD(nx,"3:5",XSource);

     if ( abs(imag(qz*(XSource[2] + XDest[2])) > 40.0 ) )
      continue;

     double R[3];
     VecSub(XDest, XSource, R);
     double Rho=sqrt( R[0]*R[0] + R[1]*R[1] );
     double xHat = (Rho==0.0) ? 1.0 : R[0] / Rho;
     double yHat = (Rho==0.0) ? 0.0 : R[1] / Rho;

     double qDotRho=0.0;
     if (Polar)
      { cdouble J[3];
        double qRho = qMag*Rho;
        double TPQ = 2.0*M_PI*qMag;
        AmosBessel('J', qRho, 0.0, 3, false, J, Data->Workspace);
        cdouble J1oqRho = (qRho==0.0 ? 0.0 : J[1]/qRho);
        cdouble Bracket = (J[0] - 2.0*J1oqRho - J[2]);
        One     = TPQ*J[0];
        Cos     = II*TPQ*J[1]*xHat;
        Sin     = II*TPQ*J[1]*yHat;
        Cos2    = TPQ*(0.5*Bracket*xHat*xHat + J1oqRho);
        CosSin  = TPQ*0.5*Bracket*xHat*yHat;
        Sin2    = TPQ*(0.5*Bracket*yHat*yHat + J1oqRho);
      }
     else
      qDotRho = q[0]*R[0] + q[1]*R[1];
     
     MTE[0][0] = Sin2;
     MTE[1][1] = Cos2;
     MTE[0][1] = MTE[1][0] = -1.0*CosSin;

     MTM[0][0] = -qz2*Cos2 / k02;
     MTM[1][1] = -qz2*Sin2 / k02;
     MTM[2][2] = q2*One  / k02;
     MTM[0][1] = MTM[1][0] = -qz2*CosSin / k02;
     MTM[2][0] =      qMag*qz*Cos / k02;
     MTM[0][2] = -1.0*MTM[2][0];
     MTM[2][1] =      qMag*qz*Sin / k02;
     MTM[1][2] = -1.0*MTM[2][1];

     cdouble ExpArg = II*( qDotRho + qz*(XSource[2]+XDest[2]) );
     cdouble Factor = II*exp(ExpArg) / (8.0*M_PI*M_PI*qz);
  
     for(int Mu=0; Mu<3; Mu++)
      for(int Nu=0; Nu<3; Nu++)
       { Integrand[36*nx + 0*9 + 3*Mu + Nu]
          += Jacobian * Factor * (rTE*MTE[Mu][Nu] + rTM*MTM[Mu][Nu]);
         Integrand[36*nx + 1*9 + 3*Mu + Nu]
          += Jacobian * Factor * 0.0; // (rTE*MTE[Mu][Nu] + rTM*MTM[Mu][Nu]);
         Integrand[36*nx + 2*9 + 3*Mu + Nu]
          += Jacobian * Factor * 0.0; // (rTE*MTE[Mu][Nu] + rTM*MTM[Mu][Nu]);
         Integrand[18*nx + 3*9 + 3*Mu + Nu]
          += Jacobian * Factor * (rTM*MTE[Mu][Nu] + rTE*MTM[Mu][Nu]);
       };

     if (LogFile)
      { fprintf(LogFile,"%e %e ",qMag/real(k0),XSource[2]);
        fprintf(LogFile,"%e %e %e %e ",real(rTE),imag(rTE),real(rTM),imag(rTM));
        fprintf(LogFile,"%e ",imag(Factor*rTE*MTE[0][0]));
        fprintf(LogFile,"%e ",imag(Factor*rTE*MTE[1][1]));
        fprintf(LogFile,"%e ",imag(Factor*rTE*MTE[2][2]));
        fprintf(LogFile,"%e ",imag(Factor*rTM*MTM[0][0]));
        fprintf(LogFile,"%e ",imag(Factor*rTM*MTM[1][1]));
        fprintf(LogFile,"%e ",imag(Factor*rTM*MTM[2][2]));
        fprintf(LogFile,"\n");
      };


   };
}

/***************************************************************/
/* get half-space DGFs via plane-wave decomposition approach   */
/***************************************************************/
void LayeredSubstrate::GetSubstrateDGF_PlaneWave(cdouble Omega,
                                                 HMatrix *XMatrix,
                                                 HMatrix *GMatrix)
{ 
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int NX = XMatrix->NR;
  int IDim = 36*NX;
  static int IDimSave=0;
  static cdouble *Integral1=0, *Integral2=0, *Error=0;
  if (IDimSave!=IDim)
   { IDimSave = IDim;
     Integral1 = (cdouble *)reallocEC(Integral1, 3 * IDim * sizeof(cdouble));
     Integral2 = Integral1 + IDim;
     Error     = Integral2 + IDim;
   };
  UpdateCachedEpsMu(Omega);

  HalfSpaceData MyData, *Data=&MyData;
  Data->XMatrix    = XMatrix;
  Data->Omega      = Omega;
  Data->Epsilon    = EpsLayer[0];
  Data->Mu         = MuLayer[0];
  Data->Polar      = true;
  Data->Accumulate = false;

  double Lower, Upper;
  Lower=0.0;
  Upper=1.0;
  Data->nCalls      = 0;
  Data->Propagating= true;
  hcubature(2*IDim, SubstrateDGFIntegrand_PlaneWave, (void *)Data, 1,
            &Lower, &Upper, MaxEvals, AbsTol, RelTol, 
            ERROR_INDIVIDUAL, (double *)Integral1, (double *)Error);
  Log(" small-q integral: %i calls",Data->nCalls);

  Lower=0.0;
  Upper=1.0;
  Data->nCalls      = 0;
  Data->Propagating = false;
  hcubature(2*IDim, SubstrateDGFIntegrand_PlaneWave, (void *)Data, 1,
            &Lower, &Upper, MaxEvals, AbsTol, RelTol,
            ERROR_INDIVIDUAL, (double *)Integral2, (double *)Error);
  Log(" large-q integral: %i calls",Data->nCalls);

 
  for(int nx=0; nx<NX; nx++)
   for(int ng=0; ng<18; ng++)
    GMatrix->SetEntry(nx, ng, Integral1[18*nx + ng] + Integral2[18*nx + ng]);
  
}
