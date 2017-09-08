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
 * FullWave_SurfaceCurrents.cc -- compute the full-wave green's function
 *                             -- above a multi-layer dielectric substrate
 *                             -- using the method of induced surface currents
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

#define II cdouble(0.0,1.0)
#ifndef ZVAC
#define ZVAC 376.73031346177
#endif

#if 0

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct SCData 
 {
   LayeredSubstrate *Substrate;
   HMatrix *XMatrix;
   double Omega;
   bool Polar;
   bool Accumulate;
   bool Propagating;
   int nCalls;

 } SCData;

/***************************************************************/
/* Integrand[ 18*nx + 0*9 + 3*Mu + Nu ] = G^{E}_{Mu,Nu}        */
/* Integrand[ 18*nx + 1*9 + 3*Mu + Nu ] = G^{M}_{Mu,Nu}        */
/***************************************************************/
int SubstrateDGFIntegrand_SC(unsigned ndim, const double *q,
                             void *UserData, unsigned fdim, double *fval)
{

  SCData *Data                = (SCData *)UserData;
  Data->nCalls++;

  LayeredSubstrate *Substrate = Data->Substrate;
  HMatrix *XMatrix            = Data->XMatrix;
  cdouble Omega               = Data->Omega;
  bool Polar                  = Data->Polar;

  cdouble EpsMedium           = Data->Substrate->EpsMedium

  int IDim = 18*XMatrix->NR;
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
       { Integrand[18*nx + 0*9 + 3*Mu + Nu] 
          += Jacobian * Factor * (rTE*MTE[Mu][Nu] + rTM*MTM[Mu][Nu]);
         Integrand[18*nx + 1*9 + 3*Mu + Nu] 
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
void LayeredSubstrate::GetHalfSpaceDGFs_SC(cdouble Omega,
                                           HMatrix *XMatrix,
                                           HMatrix *GMatrix)
{ 
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int NX = XMatrix->NR;
  int IDim = 18*NX;
  static int IDimSave=0;
  static cdouble *Integral1=0, *Integral2=0, *Error=0;
  if (IDimSave!=IDim)
   { IDimSave = IDim;
     Integral1 = (cdouble *)reallocEC(Integral1, IDim * sizeof(cdouble));
     Integral2 = (cdouble *)reallocEC(Integral2, IDim * sizeof(cdouble));
     Error     = (cdouble *)reallocEC(Error, IDim * sizeof(cdouble));
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Log("Evaluating qr integral for DGFs at %i points...",NX);

  double Lower, Upper;
  HalfSpaceData MyData, *Data=&MyData;
  Data->XMatrix    = XMatrix;
  Data->Omega      = Omega; 
  Data->Epsilon    = Epsilon;
  Data->Mu         = Mu;     
  Data->Polar      = true;
  Data->Accumulate = false;

  Lower=0.0;
  Upper=1.0;
  Data->nCalls      = 0;
  Data->Propagating= true;
  hcubature(2*IDim, HalfSpaceDGFIntegrand_Polar, (void *)Data, 1,
            &Lower, &Upper, MaxEvals, AbsTol, RelTol, 
            ERROR_INDIVIDUAL, (double *)Integral1, (double *)Error);
  Log(" small-q integral: %i calls",Data->nCalls);

  Lower=0.0;
  Upper=1.0;
  Data->nCalls      = 0;
  Data->Propagating = false;
  hcubature(2*IDim, HalfSpaceDGFIntegrand_Polar, (void *)Data, 1, 
            &Lower, &Upper, MaxEvals, AbsTol, RelTol, 
            ERROR_INDIVIDUAL, (double *)Integral2, (double *)Error);
  Log(" large-q integral: %i calls",Data->nCalls);

 
  for(int nx=0; nx<NX; nx++)
   for(int ng=0; ng<18; ng++)
    GMatrix->SetEntry(nx, ng, Integral1[18*nx + ng] + Integral2[18*nx + ng]);
  
}

#endif

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetGC0Twiddle(cdouble k2, double q2D[3], cdouble qz, double Sign,
                   cdouble GT[3][3], cdouble CT[3][3])
{
  double qx  = q2D[0];
  double qy  = q2D[1];

  GT[0][0] = 1.0 - qx*qx/k2;
  GT[1][1] = 1.0 - qy*qy/k2;
  GT[2][2] = 1.0 - qz*qz/k2;

  GT[0][1] = GT[1][0] =      -1.0*qx*qy/k2;
  GT[0][2] = GT[2][0] = -1.0*Sign*qx*qz/k2;
  GT[1][2] = GT[2][1] = -1.0*Sign*qy*qz/k2;

  CT[0][0] = CT[1][1] = CT[2][2] = 0.0;
  CT[0][1] = -1.0*Sign;   CT[1][0] = +1.0*Sign;

  if (qz==0.0)
   CT[0][2]=CT[2][0]=CT[1][2]=CT[2][1]=0.0;
  else
   { CT[0][2] =      qy/qz;  CT[2][0]=-1.0*CT[0][2];
     CT[1][2] = -1.0*qx/qz;  CT[2][1]=-1.0*CT[1][2];
   };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddScriptG0Twiddle(cdouble Omega, cdouble Eps, cdouble Mu, double q2D[2],
                        double zDmS, cdouble ScriptG0Twiddle[6][6])
{
  cdouble GT[3][3], CT[3][3];
  cdouble k2=Eps*Mu*Omega*Omega;
  cdouble qz = sqrt(k2 - q2D[0]*q2D[0] - q2D[1]*q2D[1]);
  double Sign = (zDmS > 0.0) ? 1.0 : -1.0;
  GetGC0Twiddle(k2, q2D, qz, Sign, GT, CT);
  cdouble ExpFac = exp(II*qz*fabs(zDmS));
  cdouble EEPrefac = -0.5*Omega*Mu*ZVAC*ExpFac/qz;
  cdouble EMPrefac = +0.5*ExpFac;
  cdouble MEPrefac = -0.5*ExpFac;
  cdouble MMPrefac = -0.5*Omega*Eps*ExpFac/(qz*ZVAC);
  for(int i=0; i<3; i++)
   for(int j=0; j<3; j++)
    { ScriptG0Twiddle[0+i][0+j] += EEPrefac * GT[i][j];
      ScriptG0Twiddle[0+i][3+j] += EMPrefac * CT[i][j];
      ScriptG0Twiddle[3+i][0+j] += MEPrefac * CT[i][j];
      ScriptG0Twiddle[3+i][3+j] += MMPrefac * GT[i][j];
    };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetScriptG0Twiddle(cdouble Omega, cdouble Eps, cdouble Mu, double q2D[2],
                        double zDmS, cdouble ScriptG0Twiddle[6][6])
{ memset(ScriptG0Twiddle, 0, 36*sizeof(cdouble));
  AddScriptG0Twiddle(Omega, Eps, Mu, q2D, zDmS, ScriptG0Twiddle);
}

/**********************************************************************/
/* assemble the matrix that operates on the vector of external-field  */
/* Fourier coefficients in all regions to yield the vector of         */
/* surface-current Fourier coefficients on all material interfaces    */
/**********************************************************************/
void LayeredSubstrate::AssembleMF2SMatrix(cdouble Omega, double q2D[2],
                                          HMatrix *MF2S)
{
  UpdateCachedEpsMu(Omega);
  MF2S->Zero();
  for(int a=0; a<NumInterfaces; a++)
   {
     int RowOffset = 4*a;
     for(int b=a-1; b<=a+1; b++)
      { 
        if (b<0 || b>=NumInterfaces) continue;
        int ColOffset = 4*b;

        // contributions of surface currents on interface z_b
        // to tangential-field matching equations at interface z_a
        cdouble ScriptG0Twiddle[6][6];
        double Sign=-1.0;
        if ( b==(a-1) )
         GetScriptG0Twiddle(Omega, EpsLayer[a], MuLayer[a], q2D, zInterface[a]-zInterface[b], ScriptG0Twiddle);
        else if ( b==(a+1) )
         GetScriptG0Twiddle(Omega, EpsLayer[a+1], MuLayer[a+1], q2D, zInterface[a]-zInterface[b], ScriptG0Twiddle);
        else // (b==a)
         { GetScriptG0Twiddle(Omega, EpsLayer[a],   MuLayer[a],   q2D, +1.0e-12, ScriptG0Twiddle);
           AddScriptG0Twiddle(Omega, EpsLayer[a+1], MuLayer[a+1], q2D, -1.0e-12, ScriptG0Twiddle);
           Sign=1.0;
         };

        for(int EH=0; EH<2; EH++)
         for(int KN=0; KN<2; KN++)
          for(int i=0; i<2; i++)
           for(int j=0; j<2; j++)
            MF2S->AddEntry(RowOffset+2*EH+0, ColOffset+2*KN+0, Sign*ScriptG0Twiddle[3*EH+i][3*KN+j]);
      };
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void LayeredSubstrate::GetScriptGTwiddle_SC(cdouble Omega, double q2D[2], double zDest, double zSource,
                                            cdouble ScriptGTwiddle[2][6][6])
{
  UpdateCachedEpsMu(Omega);
 
  /**********************************************************************/
  /* assemble F->S matrix ***********************************************/
  /**********************************************************************/
  int NumVars = 4*NumInterfaces;
  HMatrix MF2S(NumVars, NumVars, LHM_COMPLEX);
  AssembleMF2SMatrix(Omega, q2D, &MF2S);
  MF2S.LUFactorize();

  /**********************************************************************/
  /* assemble RHS vector for each (pource point, polarization, orientation) */
  /**********************************************************************/
  int nrSource = GetRegionIndex(zSource);
  cdouble EpsSource=EpsLayer[nrSource];
  cdouble  MuSource=MuLayer[nrSource];
  cdouble ScriptG0TSource[2][6][6];
  if (nrSource>0)
   GetScriptG0Twiddle(Omega, EpsSource, MuSource, q2D, zInterface[nrSource-1]-zSource, ScriptG0TSource[0]);
  if (nrSource<NumInterfaces)
   GetScriptG0Twiddle(Omega, EpsSource, MuSource, q2D, zInterface[nrSource]-zSource, ScriptG0TSource[1]);

  int nrDest = GetRegionIndex(zDest);
  cdouble EpsDest=EpsLayer[nrDest];
  cdouble  MuDest=MuLayer[nrDest];
  cdouble ScriptG0TDest[2][6][6];
  if (nrDest>0)
   GetScriptG0Twiddle(Omega, EpsDest, MuDest, q2D, zDest - zInterface[nrDest-1], ScriptG0TDest[0]);
  if (nrDest<NumInterfaces)
   GetScriptG0Twiddle(Omega, EpsDest, MuDest, q2D, zDest - zInterface[nrDest], ScriptG0TDest[1]);

  HVector KNTwiddle(NumVars, LHM_COMPLEX);
  for(int Nu=0; Nu<6; Nu++)
   {
     // fill in RHS vector describing fields sourced by point in region #nrSource
     KNTwiddle.Zero();
     for (int Which=0; Which<=1; Which++)
      { int nr = nrSource-1+Which;
        if (nr<0 || nr>=NumInterfaces) continue;
        int RowOffset = 4*nr;  
        double Sign   = (Which == 1 ? -1.0 : 1.0);
        for(int EH=0; EH<2; EH++)
         for(int i=0; i<2; i++)
          KNTwiddle.SetEntry(RowOffset + 2*EH + i, Sign*ScriptG0TSource[Which][3*EH+i][Nu]);
      };

     // solve for surface currents on all layers
     MF2S.LUSolve(&KNTwiddle); 

     // convolve surface currents with GF to yield fields at eval point
     for (int Which=0; Which<=1; Which++)
      { int nr = nrDest-1+Which;
        if (nr<0 || nr>=NumInterfaces) continue;
        int RowOffset = 4*nr;
        double Sign   = (Which == 1 ? +1.0 : -1.0);
        for(int Mu=0; Mu<6; Mu++)
         for(int EH=0; EH<2; EH++)
          for(int i=0; i<2; i++)
           ScriptGTwiddle[Which][Mu][Nu] += Sign*ScriptG0TDest[Which][Mu][3*EH+i]*KNTwiddle.GetEntry(RowOffset+2*EH+i);
      }; // for (Which=0 ...)

   }; // for(int Nu=0; Nu<6; Nu++)

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
//void LayeredSubstrate::GetHalfSpaceDGFs_SC(cdouble Omega, HMatrix *XMatrix, HMatrix *GMatrix)
//{}
