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
 * SommerfeldIntegrand.cc -- compute the integrand of the integral
 *                        -- over q (in-plane wavevector magnitude)
 *                        -- that we evaluate using Sommerfeld     
 *                        -- integration techniques to get the
 *                        -- real-space DGF
 *
 * homer reid  -- 3/2017-9/2017
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

/********************************************************************/
/* If environment variable SCUFF_GTWIDDLE_ONLY is set, zero out all */
/* components of the gTwiddle vector except those listed in         */
/* SCUFF_GTWIDDLE_ONLY                                              */
/* example: export SCUFF_GTWIDDLE_ONLY="4 7 12"                     */
/********************************************************************/
void ZerogTwiddle(cdouble *gTwiddle)
{
  static bool Retain[NUMGTWIDDLE], Init=false;
  
  if (!Init)
   { 
     Init=true;
     for(int ng=0; ng<NUMGTWIDDLE; ng++) Retain[ng]=true;
     char *s=getenv("SCUFF_GTWIDDLE_ONLY");
     if (!s) return;
     for(int ng=0; ng<NUMGTWIDDLE; ng++) Retain[ng]=false;
     char *Tokens[25];
     int NumTokens=Tokenize(s, Tokens, 25);
     for(int nt=0; nt<NumTokens; nt++)
      { int ng;
        sscanf(Tokens[nt],"%i",&ng);
        Retain[ng]=true;
        Log("Retaining #%i\n",ng);
      };
   };

  for(int ng=0; ng<NUMGTWIDDLE; ng++)
   if (!Retain[ng])
    gTwiddle[ng]=0.0;
  
}


/********************************************************************/
/* Call GetGTwiddle to compute the full 6x6 Fourier-space DGF, then */
/* extract from it the 18 distinct scalar functions gTwiddle(q)     */
/* that get paired with Bessel-function factors to form the         */
/* integrand of the Sommerfeld integral.                            */
/*                                                                  */
/* gTwiddleVD[0][0] = gTwiddle factors                              */
/* gTwiddleVD[0][1] = dg/dzDest   (if dzDest=true)                  */
/* gTwiddleVD[1][0] = dg/dzSource (if dzSource=true)                */
/* gTwiddleVD[1][1] = d2g/dzDestdzSource (if dzDest=dzSource=true)  */
/********************************************************************/
void LayeredSubstrate::gTwiddleFromGTwiddle(cdouble Omega, cdouble q,
                                            double zDest, double zSource,
                                            HMatrix *RTwiddle, HMatrix *WMatrix,
                                            HMatrix *STwiddle, cdouble *gTwiddleVD[2][2],
                                            bool dzDest, bool dzSource)
{
  cdouble GBuffer[4*36];
  HMatrix   GTwiddle(6,6,LHM_COMPLEX,GBuffer+0*36);
  HMatrix     dGTdzD(6,6,LHM_COMPLEX,GBuffer+1*36);
  HMatrix     dGTdzS(6,6,LHM_COMPLEX,GBuffer+2*36);
  HMatrix d2GTdzDdzS(6,6,LHM_COMPLEX,GBuffer+3*36);

  cdouble q2D[2];
  q2D[0]=q;
  q2D[1]=0.0;

  for(int ndzDest=0; ndzDest <= (dzDest ? 1 : 0); ndzDest++)
   for(int ndzSource=0; ndzSource <= (dzSource ? 1 : 0); ndzSource++)
    {   
     GetScriptGTwiddle(Omega, q2D, zDest, zSource,
                       RTwiddle, WMatrix, STwiddle, &GTwiddle, 
                       (ndzDest==1), (ndzSource==1));

     cdouble *gTwiddle = gTwiddleVD[ndzDest][ndzSource];

     gTwiddle[_EE0P] = GTwiddle.GetEntry(0+1,0+1);
     gTwiddle[_EE0Z] = GTwiddle.GetEntry(0+2,0+2);
     gTwiddle[_EE1A] = GTwiddle.GetEntry(0+0,0+2);
     gTwiddle[_EE2A] = GTwiddle.GetEntry(0+0,0+0) - gTwiddle[_EE0P];

     gTwiddle[_EM0P] = GTwiddle.GetEntry(0+0,3+1);
     gTwiddle[_EM1A] = GTwiddle.GetEntry(0+1,3+2);
     gTwiddle[_EM1B] = GTwiddle.GetEntry(0+2,3+1);
     gTwiddle[_EM2A] = -1.0*(GTwiddle.GetEntry(0+1,3+0) + gTwiddle[_EM0P]);

     gTwiddle[_MM0P] = GTwiddle.GetEntry(3+1,3+1);
     gTwiddle[_MM0Z] = GTwiddle.GetEntry(3+2,3+2);
     gTwiddle[_MM1A] = GTwiddle.GetEntry(3+0,3+2);
     gTwiddle[_MM2A] = GTwiddle.GetEntry(3+0,3+0) - gTwiddle[_MM0P];

     // the following are only independent of the preceding
     // if zDest, zSource lie in different substrate layers
     // so in principle their calculation could be omitted
     // for the same-layer case
     gTwiddle[_EE1B] = GTwiddle.GetEntry(0+2,0+0);
     gTwiddle[_ME0P] = GTwiddle.GetEntry(3+1,0+0);
     gTwiddle[_ME1A] = GTwiddle.GetEntry(3+2,0+1);
     gTwiddle[_ME1B] = GTwiddle.GetEntry(3+1,0+2);
     gTwiddle[_ME2A] = -1.0*(GTwiddle.GetEntry(3+0,0+1) + gTwiddle[_ME0P]);
     gTwiddle[_MM1B] = GTwiddle.GetEntry(3+2,3+0);

     /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
     ZerogTwiddle(gTwiddle);
     /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
   };

}

/********************************************************************/
/********************************************************************/
/********************************************************************/
void LayeredSubstrate::gTwiddleHardCoded(cdouble Omega, cdouble q,
                                         double zDest, double zSource,
                                         cdouble *gTwiddleVD[2][2],
                                         bool dzDest, bool dzSource)
{
  UpdateCachedEpsMu(Omega);

  cdouble *gTwiddle=gTwiddleVD[0][0];
  if (ForceFreeSpace)
   { 
     // vacuum DGF
     cdouble k2  = Omega*Omega, q2=q*q, qz2 = k2-q2;
     if ( abs(qz2) < 1.0e-12 )
      { q*=(1.0+1.0e-5); q2=q*q; qz2=k2-q2; };
     cdouble qz=csqrt2(qz2);
     cdouble ExpFac = exp(II*qz*fabs(zDest-zSource));
     double Sign = (zDest>=zSource) ? 1.0 : -1.0;
     cdouble EEFactor = -0.5*Omega*ZVAC*ExpFac/qz;
     cdouble EMFactor = +0.5*ExpFac;
     cdouble MEFactor = -0.5*ExpFac;
     cdouble MMFactor = -0.5*Omega*ExpFac/(ZVAC*qz);
     gTwiddle[_EE0P] = EEFactor;
     gTwiddle[_EE0Z] = EEFactor*(1.0-qz2/k2);
     gTwiddle[_EE1A] = -1.0*EEFactor*Sign*q*qz/k2;
     gTwiddle[_EE1B] = -1.0*EEFactor*Sign*q*qz/k2;
     gTwiddle[_EE2A] = EEFactor*-1.0*q2/k2;
     gTwiddle[_EM0P] = EMFactor*-1.0*Sign;
     gTwiddle[_EM1A] = -1.0*EMFactor*q/qz;
     gTwiddle[_EM1B] = +1.0*EMFactor*q/qz;
     gTwiddle[_EM2A] = 0.0;
     gTwiddle[_ME0P] = MEFactor*-1.0*Sign;
     gTwiddle[_ME1A] = -1.0*MEFactor*q/qz;
     gTwiddle[_ME1B] = +1.0*MEFactor*q/qz;
     gTwiddle[_ME2A] = 0.0;
     gTwiddle[_MM0P] = MMFactor;
     gTwiddle[_MM0Z] = MMFactor*(1.0-qz2/k2);
     gTwiddle[_MM1A] = -1.0*MMFactor*Sign*q*qz/k2;
     gTwiddle[_MM1B] = -1.0*MMFactor*Sign*q*qz/k2;
     gTwiddle[_MM2A] = MMFactor*-1.0*q2/k2;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
     ZerogTwiddle(gTwiddle);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
   }
  else if (NumInterfaces==1)
   {
     // dielectric half-space
#if 0
     cdouble Eps = EpsLayer[GetLayerIndex(zDest)];
     cdouble k2  = Eps*Eps*Omega*Omega, q2=q*q, qz2 = k2-q2;
     if ( abs(qz2) < 1.0e-12 )
      { q*=(1.0+1.0e-5); q2=q*q; qz2=k2-q2; };
     cdouble qz=csqrt2(qz2);
     cdouble ExpFac = exp(II*qz*fabs(zDest-zSource));
     double Sign = (zDest>=zSource) ? 1.0 : -1.0;
     cdouble EEFactor = -0.5*Omega*ZVAC*ExpFac/qz;
     cdouble EMFactor = +0.5*ExpFac;
     cdouble MEFactor = -0.5*ExpFac;
     cdouble MMFactor = -0.5*Omega*ExpFac/(ZVAC*qz);
     gTwiddle[_EE0P] = EEFactor;
     gTwiddle[_EE0Z] = EEFactor*(1.0-qz2/k2);
     gTwiddle[_EE1A] = -1.0*EEFactor*Sign*q*qz/k2;
     gTwiddle[_EE1B] = -1.0*EEFactor*Sign*q*qz/k2;
     gTwiddle[_EE2A] = EEFactor*-1.0*q2/k2;
     gTwiddle[_EM0P] = EMFactor*-1.0*Sign;
     gTwiddle[_EM1A] = -1.0*EMFactor*q/qz;
     gTwiddle[_EM1B] = +1.0*EMFactor*q/qz;
     gTwiddle[_EM2A] = 0.0;
     gTwiddle[_ME0P] = MEFactor*-1.0*Sign;
     gTwiddle[_ME1A] = -1.0*MEFactor*q/qz;
     gTwiddle[_ME1B] = +1.0*MEFactor*q/qz;
     gTwiddle[_ME2A] = 0.0;
     gTwiddle[_MM0P] = MMFactor;
     gTwiddle[_MM0Z] = MMFactor*(1.0-qz2/k2);
     gTwiddle[_MM1A] = -1.0*MMFactor*Sign*q*qz/k2;
     gTwiddle[_MM1B] = -1.0*MMFactor*Sign*q*qz/k2;
     gTwiddle[_MM2A] = MMFactor*-1.0*q2/k2;
#endif
   };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int SommerfeldIntegrand(unsigned ndim, const double *x,
                        void *UserData, unsigned fdim, double *fval)
{

  SommerfeldIntegrandData *SID = (SommerfeldIntegrandData *)UserData;
  LayeredSubstrate *S          = SID->Substrate;
  cdouble Omega                = SID->Omega;
  double q0                    = SID->q0;
  bool uTransform              = SID->uTransform;
  HMatrix *XMatrix             = SID->XMatrix;
  HMatrix *RTwiddle            = SID->RTwiddle;
  HMatrix *WMatrix             = SID->WMatrix; 
  HMatrix *STwiddle            = SID->STwiddle;
  bool dzSource                = SID->dzSource;
  bool dzDest                  = SID->dzDest;
  SID->NumPoints++;
  FILE *byqFile                = SID->byqFile;
 
  memset(fval, 0, fdim*sizeof(double));

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble q=x[0];
  double Jac=1.0;
  if(ndim==2)
   { 
     q += II*x[1];
   }
  else if (uTransform)
   { double u=x[0];
     if (u==1.0)
      return 0;
     q = q0 + u/(1.0-u);
     Jac = 1.0/( (1.0-u)*(1.0-u) );
   };

  /***************************************************************/
  /* get gTwiddle functions **************************************/
  /***************************************************************/
  cdouble *zfval = (cdouble *)fval;
  cdouble gBuffer[4*NUMGTWIDDLE];
  cdouble *gTwiddleVD[2][2]={{gBuffer+0*NUMGTWIDDLE, gBuffer+1*NUMGTWIDDLE},
                             {gBuffer+2*NUMGTWIDDLE, gBuffer+3*NUMGTWIDDLE}};
  cdouble J[4];
  double ScaleFac=1.0;
  double LastzDest=1.234e56, LastzSource=2.345e67, LastRho=3.456e78;
  for(int nx=0; nx<XMatrix->NR; nx++)
   { 
     double Rhox    = XMatrix->GetEntryD(nx, 0) - XMatrix->GetEntryD(nx,3);
     double Rhoy    = XMatrix->GetEntryD(nx, 1) - XMatrix->GetEntryD(nx,4);
     double Rho     = sqrt(Rhox*Rhox + Rhoy*Rhoy);
     double zDest   = XMatrix->GetEntryD(nx,2);
     double zSource = XMatrix->GetEntryD(nx,5);

     /***************************************************************/
     /* get gTwiddle factors ****************************************/
     /***************************************************************/
     if ( nx==0 || !EqualFloat(zDest,LastzDest) || !EqualFloat(zSource,LastzSource) )
      { LastzDest=LastzDest;
        LastzSource=LastzSource;
        if (S->HardCoded)
         S->gTwiddleHardCoded(Omega, q, zDest, zSource,
                              gTwiddleVD, dzDest, dzSource);
        else
         S->gTwiddleFromGTwiddle(Omega, q, zDest, zSource, 
                                 RTwiddle, WMatrix, STwiddle,
                                 gTwiddleVD, dzDest, dzSource);
      };
     
     /***************************************************************/
     /* get bessel-function factors  ********************************/
     /***************************************************************/
     if (nx==0 || !EqualFloat(Rho, LastRho) )
      { LastRho=Rho;
        cdouble qRho=q*Rho;
        double Workspace[4*3];
        bool Scale = (abs(qRho)<1.0e-4 ? false : true);
        ScaleFac = Scale ? exp(fabs(imag(qRho))) : 1.0;
        AmosBessel('J', qRho, 0, 3, Scale, J, Workspace);
        J[3] = ( (abs(qRho) < 1.0e-4) ? 0.5*(1.0 - qRho*qRho/8.0)
                                      : J[1]/qRho
               );
        J[1]*=II;
        J[2]*=-1.0;
      };

     /***************************************************************/
     /* assemble integrand vector ***********************************/
     /***************************************************************/
     cdouble Factor = q*Jac*ScaleFac / (2.0*M_PI);
     cdouble *Integrand = zfval;
     for(int dzd=0; dzd<=(dzDest ? 1 : 0); dzd++)
      for(int dzs=0; dzs<=(dzSource ? 1 : 0); dzs++)
       { 
         cdouble *gTwiddle=gTwiddleVD[dzd][dzs];

         Integrand[_EE0P] = Factor*gTwiddle[_EE0P]*J[0];
         Integrand[_EE0Z] = Factor*gTwiddle[_EE0Z]*J[0];
         Integrand[_EE1A] = Factor*gTwiddle[_EE1A]*J[1];
         Integrand[_EE1B] = Factor*gTwiddle[_EE1B]*J[1];
         Integrand[_EE2A] = Factor*gTwiddle[_EE2A]*J[2];
         Integrand[_EE2B] = Factor*gTwiddle[_EE2A]*J[3];
    
         Integrand[_EM0P] = Factor*gTwiddle[_EM0P]*J[0];
         Integrand[_EM1A] = Factor*gTwiddle[_EM1A]*J[1];
         Integrand[_EM1B] = Factor*gTwiddle[_EM1B]*J[1];
         Integrand[_EM2A] = Factor*gTwiddle[_EM2A]*J[2];
         Integrand[_EM2B] = Factor*gTwiddle[_EM2A]*J[3];
    
         Integrand[_ME0P] = Factor*gTwiddle[_ME0P]*J[0];
         Integrand[_ME1A] = Factor*gTwiddle[_ME1A]*J[1];
         Integrand[_ME1B] = Factor*gTwiddle[_ME1B]*J[1];
         Integrand[_ME2A] = Factor*gTwiddle[_ME2A]*J[2];
         Integrand[_ME2B] = Factor*gTwiddle[_ME2A]*J[3];
    
         Integrand[_MM0P] = Factor*gTwiddle[_MM0P]*J[0];
         Integrand[_MM0Z] = Factor*gTwiddle[_MM0Z]*J[0];
         Integrand[_MM1A] = Factor*gTwiddle[_MM1A]*J[1];
         Integrand[_MM1B] = Factor*gTwiddle[_MM1B]*J[1];
         Integrand[_MM2A] = Factor*gTwiddle[_MM2A]*J[2];
         Integrand[_MM2B] = Factor*gTwiddle[_MM2A]*J[3];

         if (byqFile)
          { fprintf(byqFile,"%e %e %e %e %e %i %i ",real(q),imag(q),Rho,zDest,zSource,dzd,dzs);
            fprintVecCR(byqFile,Integrand,NUMGSCALAR);
          };

         Integrand += NUMGSCALAR;

       };

   }; // for(int nx=0; nx<XMatrix->NR; nx++)

  return 0;

}
