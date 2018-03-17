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
 * gFrak.cc -- compute the quantities gFrak(Rho,zDest,zSource)
 *          -- that we need to compute the substrate DGF
 *
 * homer reid  -- 3/2017-2/2018
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
                                            cdouble *gTwiddleVD[2][2],
                                            cdouble *Workspace,
                                            bool dzDest, bool dzSource)
{
  cdouble GBuffer[36];
  HMatrix GTwiddle(6,6,LHM_COMPLEX,GBuffer);

  cdouble q2D[2];
  q2D[0]=q;
  q2D[1]=0.0;

  for(int ndzDest=0; ndzDest <= (dzDest ? 1 : 0); ndzDest++)
   for(int ndzSource=0; ndzSource <= (dzSource ? 1 : 0); ndzSource++)
    {   
     GetScriptGTwiddle(Omega, q2D, zDest, zSource, &GTwiddle, Workspace,
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

   };

}

/***************************************************************/
/* fetch Bessel-function factors                               */
/*  JdJFactors[0][i] = J[i]                                    */
/*  JdJFactors[1][i] = d/dRho J_i                              */
/*                                                             */
/* J[_J0] = J_0(q*Rho)                                         */
/* J[_J1] = I*J_1(q*Rho)                                       */
/* J[_J2] = -1*J_2(q*Rho)                                      */
/* J[_JJ] = J_1(q*Rho) / (q*Rho)                               */
/*                                                             */
/***************************************************************/
#define _J0 0
#define _J1 1
#define _J2 2
#define _JJ 3
void GetJdJFactors(cdouble q, double Rho, cdouble JdJFactors[2][4],
                   bool NeedRhoDerivatives, int NuMax)
{
  cdouble qRho = q*Rho;
  cdouble J0, J1, J2,  J1oqRho;
  cdouble dJ0, dJ1, dJ2, dJ1oqRho;
  if ( abs(qRho)<1.0e-4 ) // series expansions for small arguments
   { cdouble qRho2=qRho*qRho;
     J0       = 1.0 - qRho2/4.0;
     J1oqRho  = 0.5*(1.0-qRho2/8.0);
     J1       = J1oqRho*qRho;
     J2       = 0.125*qRho2*(1.0-qRho2/12.0);
     dJ0      = -0.5*q*qRho;
     dJ1oqRho = -0.125*q*qRho;
     dJ1      = q*J1oqRho + qRho*dJ1oqRho;
     dJ2      = q*qRho/4.0 - qRho2*qRho*Rho/24.0;
   }
  else if (abs(qRho)>1.0e4) // asymptotic forms for large argument
   { double rqRho = real(qRho);
     double JPreFac = sqrt(2.0/(M_PI*rqRho));
     J0       = JPreFac * cos(rqRho - 0.25*M_PI);
     J1       = JPreFac * cos(rqRho - 0.75*M_PI);
     J2       = JPreFac * cos(rqRho - 1.25*M_PI);
     J1oqRho  = J1/qRho;
     dJ0      = -0.5*J0/Rho - q*JPreFac*sin(rqRho-0.25*M_PI);
     dJ1      = -0.5*J1/Rho - q*JPreFac*sin(rqRho-0.75*M_PI);
     dJ1oqRho =    dJ1/qRho - J1oqRho/Rho;
     dJ2      = -0.5*J2/Rho - q*JPreFac*sin(rqRho-1.25*M_PI);
   }
  else
   { double Workspace[12];
     cdouble JFactors[3]={0.0, 0.0, 0.0};
     int NumOrders = NuMax+1;
     if (NeedRhoDerivatives && NuMax<2)
      NumOrders++;
     if (NumOrders>3) NumOrders=3;
     AmosBessel('J', qRho, 0.0, NumOrders, false, JFactors, Workspace);
     J0=JFactors[0];
     J1=JFactors[1];
     J2=JFactors[2];
     J1oqRho = JFactors[1]/qRho;
     dJ0 = -q*J1;
     dJ1 = 0.5*q*(J0 - J2);
     dJ1oqRho = dJ1/qRho - J1oqRho/Rho;
     dJ2 = q*(J1 - 2.0*J2/qRho);
   };

  JdJFactors[0][_J0] = J0;
  JdJFactors[0][_J1] = II*J1;
  JdJFactors[0][_J2] = -1.0*J2;
  JdJFactors[0][_JJ] = J1oqRho;

  JdJFactors[1][_J0] = dJ0;
  JdJFactors[1][_J1] = II*dJ1;
  JdJFactors[1][_J2] = -1.0*dJ2;
  JdJFactors[1][_JJ] = dJ1oqRho;

}

/***************************************************************/
/* integrand passed to SommerfeldIntegrator to get gFrak() at  */
/* an array of user-specified points.                          */
/***************************************************************/
typedef struct gFrakIntegrandData
 {
   LayeredSubstrate *Substrate;
   cdouble Omega;
   double q0;
   bool uTransform;
   HMatrix *XMatrix;
   cdouble *Workspace;
   bool dRho, dzDest, dzSource;
   bool Accumulate;
   bool SubtractQS;
   bool EEOnly;    
   bool ScalarPPIs;
   FILE *byqFile;
   int NumPoints;
 } gFrakIntegrandData;

int gFrakIntegrand(unsigned ndim, const double *x,
                   void *UserData, unsigned fdim, double *fval)
{
  gFrakIntegrandData *Data = (gFrakIntegrandData *)UserData;
  LayeredSubstrate *S      = Data->Substrate;
  cdouble Omega            = Data->Omega;
  double q0                = Data->q0;
  bool uTransform          = Data->uTransform;
  HMatrix *XMatrix         = Data->XMatrix;
  cdouble *Workspace       = Data->Workspace;
  bool dRho                = Data->dRho;
  bool dzDest              = Data->dzDest;
  bool dzSource            = Data->dzSource;
  bool Accumulate          = Data->Accumulate;
  bool SubtractQS          = Data->SubtractQS;
  bool EEOnly              = Data->EEOnly;
  bool ScalarPPIs          = Data->ScalarPPIs;
  FILE *byqFile            = Data->byqFile;
  Data->NumPoints++;
 
  if (!Accumulate)
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
  if (q==0.0)
   q=1.0e-6;
  cdouble q2=q*q;
  cdouble Factor = q*Jac/(2.0*M_PI);
  int NumgFrak = EEOnly ? 6 : NUMGFRAK;

  /***************************************************************/
  /* assemble gFrakTwiddle vector for all spatial evaluation points */
  /***************************************************************/
  cdouble *Integrand = (cdouble *)fval;
  cdouble gTVDBuffer[4*NUMGTWIDDLE];
  cdouble *gTwiddleVD[2][2]=
   { {gTVDBuffer+0*NUMGTWIDDLE, gTVDBuffer+1*NUMGTWIDDLE},
     {gTVDBuffer+2*NUMGTWIDDLE, gTVDBuffer+3*NUMGTWIDDLE}
   };
  cdouble JdJ[2][4];
  double LastRho=1.234e56, LastzDest=2.345e67, LastzSource=3.456e78;
  bool LastdRho=false, LastdzDest=false,  LastdzSource=false;
  for(int nx=0, nn=0; nx<XMatrix->NR; nx++)
   { 
     double Rhox    = XMatrix->GetEntryD(nx, 0) - XMatrix->GetEntryD(nx,3);
     double Rhoy    = XMatrix->GetEntryD(nx, 1) - XMatrix->GetEntryD(nx,4);
     double Rho     = sqrt(Rhox*Rhox + Rhoy*Rhoy);
     double zDest   = XMatrix->GetEntryD(nx,2);
     double zSource = XMatrix->GetEntryD(nx,5);

     /***************************************************************/
     /* get gTwiddle factors ****************************************/
     /***************************************************************/
     bool NeedgTwiddle = (   !EqualFloat(zDest,LastzDest)
                          || (dzDest && !LastdzDest)
                          || !EqualFloat(zSource,LastzSource)
                          || (dzSource && !LastdzSource)
                         );

     if (NeedgTwiddle)
      { LastzDest=zDest;      LastdzDest=dzDest;
        LastzSource=zSource;  LastdzSource=dzSource;
        S->gTwiddleFromGTwiddle(Omega, q, zDest, zSource, gTwiddleVD, 
                                Workspace, dzDest, dzSource);
      };
     
     /***************************************************************/
     /* get bessel-function factors  ********************************/
     /***************************************************************/
     bool NeedJ = ( !EqualFloat(Rho,LastRho) || (dRho && !LastdRho) );
     if (NeedJ)
      { LastRho  = Rho;
        LastdRho = dRho;
        GetJdJFactors(q, Rho, JdJ, dRho);
      };

     /***************************************************************/
     /* assemble integrand vector ***********************************/
     /***************************************************************/
     for(int dr=0; dr<=(dRho ? 1 : 0); dr++)
      for(int dzd=0; dzd<=(dzDest ? 1 : 0); dzd++)
       for(int dzs=0; dzs<=(dzSource ? 1 : 0); dzs++)
        { 
          cdouble *J=JdJ[dr];
          cdouble *gTwiddle=gTwiddleVD[dzd][dzs];

          cdouble *gFTwiddle = Integrand + NumgFrak*(nn++);

          if (ScalarPPIs)
           { int DestLayer  = S->GetLayerIndex(zDest);
             cdouble EpsRel = S->EpsLayer[DestLayer], MuRel=S->MuLayer[DestLayer];
             cdouble qz = sqrt(EpsRel*MuRel*Omega*Omega - q2);
             gFTwiddle[_EE0P] += Factor*gTwiddle[_EE0P]*J[0];
             gFTwiddle[_EE2A] += Factor*gTwiddle[_EE2A]*J[0]/q2;
             gFTwiddle[_ME1A] += Factor*gTwiddle[_EE0P]*J[1]*q/(II*Omega*MuRel*ZVAC);
             gFTwiddle[_ME0P] += -II*qz*Factor*gTwiddle[_EE0P]*J[0]/(II*Omega*MuRel*ZVAC);
             gFTwiddle[_EM1A] -= Factor*gTwiddle[_MM0P]*J[1]*q/(II*Omega*EpsRel/ZVAC);
             gFTwiddle[_EM0P] += -II*qz*Factor*gTwiddle[_MM0P]*J[0]/(II*Omega*EpsRel/ZVAC);
             gFTwiddle[_MM0P] -= Factor*gTwiddle[_MM0P]*J[0];
             gFTwiddle[_MM2A] += Factor*gTwiddle[_MM2A]*J[0]/q2;
             continue;
           }

          gFTwiddle[_EE0P] += Factor*gTwiddle[_EE0P]*J[0];
          gFTwiddle[_EE0Z] += Factor*gTwiddle[_EE0Z]*J[0];
          gFTwiddle[_EE1A] += Factor*gTwiddle[_EE1A]*J[1];
          gFTwiddle[_EE1B] += Factor*gTwiddle[_EE1B]*J[1];
          gFTwiddle[_EE2A] += Factor*gTwiddle[_EE2A]*J[2];
          gFTwiddle[_EE2B] += Factor*gTwiddle[_EE2A]*J[3];
          if (EEOnly) continue;
     
          gFTwiddle[_EM0P] += Factor*gTwiddle[_EM0P]*J[0];
          gFTwiddle[_EM1A] += Factor*gTwiddle[_EM1A]*J[1];
          gFTwiddle[_EM1B] += Factor*gTwiddle[_EM1B]*J[1];
          gFTwiddle[_EM2A] += Factor*gTwiddle[_EM2A]*J[2];
          gFTwiddle[_EM2B] += Factor*gTwiddle[_EM2A]*J[3];
     
          gFTwiddle[_ME0P] += Factor*gTwiddle[_ME0P]*J[0];
          gFTwiddle[_ME1A] += Factor*gTwiddle[_ME1A]*J[1];
          gFTwiddle[_ME1B] += Factor*gTwiddle[_ME1B]*J[1];
          gFTwiddle[_ME2A] += Factor*gTwiddle[_ME2A]*J[2];
          gFTwiddle[_ME2B] += Factor*gTwiddle[_ME2A]*J[3];
     
          gFTwiddle[_MM0P] += Factor*gTwiddle[_MM0P]*J[0];
          gFTwiddle[_MM0Z] += Factor*gTwiddle[_MM0Z]*J[0];
          gFTwiddle[_MM1A] += Factor*gTwiddle[_MM1A]*J[1];
          gFTwiddle[_MM1B] += Factor*gTwiddle[_MM1B]*J[1];
          gFTwiddle[_MM2A] += Factor*gTwiddle[_MM2A]*J[2];
          gFTwiddle[_MM2B] += Factor*gTwiddle[_MM2A]*J[3];

          if (byqFile)
           { fprintf(byqFile,"%e %e %e %e %e %i %i %i ",real(q),imag(q),Rho,zDest,zSource,dr,dzd,dzs);
             fprintVecCR(byqFile,gFTwiddle,NumgFrak);
           };
        };

   }; // for(int nx=0; nx<XMatrix->NR; nx++)

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (SubtractQS)
   { cdouble *zfval=(cdouble *)fval;
     int zfdim = fdim/2;
     VecScale(zfval, -1.0, zfdim);
     Data->Omega*=-1.0;
     Data->Accumulate=true;
     Data->SubtractQS=false;
     gFrakIntegrand(ndim, x, UserData, fdim, fval);
     Data->Omega*=-1.0;
     Data->Accumulate=false;
     Data->SubtractQS=true;
     VecScale(zfval, -1.0, zfdim);
   }

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  return 0;

} // gFrakIntegrand(...)

void LayeredSubstrate::GetgFrakTwiddle(cdouble Omega, cdouble q, double Rho,
                                       double zDest, double zSource,
                                       cdouble*gFrakTwiddle,
                                       bool SubtractQS, bool EEOnly,
                                       bool ScalarPPIs)
{
  double XDS[6]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  XDS[0]=Rho;
  XDS[2]=zDest;
  XDS[5]=zSource;
  HMatrix XMatrix(1,6,LHM_REAL,XDS);
  
  gFrakIntegrandData MyData, *Data=&MyData;

  Data->Substrate  = this;
  Data->Omega      = Omega;
  Data->q0         = 0.0;
  Data->uTransform = false;
  Data->XMatrix    = &XMatrix;
  Data->Workspace  = 0;
  Data->dRho       = false;
  Data->dzDest     = false;
  Data->dzSource   = false;
  Data->Accumulate = false;
  Data->SubtractQS = SubtractQS;
  Data->EEOnly     = EEOnly;
  Data->EEOnly     = ScalarPPIs;
  Data->byqFile    = 0;
  Data->NumPoints  = 0;

  int zfdim = EEOnly ? 6 : NUMGFRAK;
  int ndim=2;
  gFrakIntegrand(ndim, (double *)&q, (void *)Data, 2*zfdim, (double *)gFrakTwiddle);
}

/***************************************************************/
/* compute a, c parameters for the contour-integral portion of */
/* the Sommerfeld integral via the procedure described in      */
/* Golubovic et al, "Efficient Algorithms for Computing        */
/* Sommerfeld Intergral Tails," IEEE Transactions on Antennas  */
/* and Propagation **60** 2409 (2012).                         */
/* DOI: 10.1109/TAP.2012.2189718                               */
/*--------------------------------------------------------------*/
void GetacSommerfeld(LayeredSubstrate *S, cdouble Omega,
                     double Rho, double zDest, double zSource,
                     double *pa, double *pc)
{
  int NumLayers     = S->NumLayers;
  cdouble *EpsLayer = S->EpsLayer;
  cdouble *MuLayer  = S->MuLayer;

  double MaxRealn2 = 0.0;
  for(int nl=0; nl<NumLayers; nl++)
   MaxRealn2 = fmax(MaxRealn2, real(EpsLayer[nl]*MuLayer[nl]));
  double a = 1.5*sqrt(MaxRealn2)*real(Omega);
  double c = real(Omega);
   
  if ( Rho>fabs(zDest-zSource) && real(Omega*Rho)>1.0 )
   c = 1.0/Rho;

  *pa=a;
  *pc=c;
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int LayeredSubstrate::GetgFrak(cdouble Omega, HMatrix *XMatrix, cdouble *gFrak,
                               cdouble *Workspace,
                               bool SubtractQS, bool EEOnly, bool ScalarPPIs,
                               bool dRho, bool dzDest, bool dzSource)
{
  int NX         = XMatrix->NR;
  int NF         = NUMGFRAK * (dRho?2:1) * (dzDest?2:1) * (dzSource?2:1);
  int zfdim      = NX * NF;
  int fdim       = 2*zfdim;
  cdouble *Error = new cdouble[zfdim];

  bool OwnsWorkspace = (Workspace==0);
  if (OwnsWorkspace)
   Workspace = CreateScriptGTwiddleWorkspace();
 
  char *s=getenv("SOMMERFELD_BYQFILE");
  FILE *byqFile = (s ? fopen(s,"w") : 0);
  s=getenv("SOMMERFELD_INTEGRATOR_BYPASS");
  bool Bypass = (s && s[0]=='1');
  s=getenv("SOMMERFELD_INTEGRATOR_XNU");
  int xNu = ( s ? s[0]-'0' : 0);
  s=getenv("SOMMERFELD_FORCE_SUBTRACT");
  if (s) SubtractQS = (s[0]=='1');

  const char *STStr = (SubtractQS ? "subtracting" : "not subtracting");

  gFrakIntegrandData MyData, *Data=&MyData;
  Data->Substrate   = this;
  Data->Omega       = Omega;
  Data->q0          = 0.0;
  Data->uTransform  = false;
  Data->XMatrix     = XMatrix;
  Data->Workspace   = Workspace;
  Data->dRho        = dRho;
  Data->dzSource    = dzDest;
  Data->dzDest      = dzSource;
  Data->Accumulate  = false;
  Data->SubtractQS  = SubtractQS;
  Data->EEOnly      = EEOnly;
  Data->ScalarPPIs  = ScalarPPIs;
  Data->byqFile     = byqFile;
  Data->NumPoints   = 0;
  
  /***************************************************************/
  /* if requested, bypass SommerfeldIntegrator() and just do     */
  /* an adaptive quadrature over the entire positive q axis      */
  /***************************************************************/
  if (Bypass)
   { Data->uTransform=true;
     double uMin=0.0, uMax=1.0;
     int ndim=1;
     Log("Computing gFrak via simple quadrature (%s) ...",STStr);
     pcubature(fdim, gFrakIntegrand, (void *)Data,
               ndim, &uMin, &uMax, qMaxEval, qAbsTol, qRelTol, ERROR_PAIRED,
               (double *)gFrak, (double *)Error);
   }
  /***************************************************************/
  /* otherwise use the Sommerfeld integrator *********************/
  /***************************************************************/
  else
   { 
     bool Verbose=true;
     double Rhox  = XMatrix->GetEntryD(0,0) - XMatrix->GetEntryD(0,3);
     double Rhoy  = XMatrix->GetEntryD(0,1) - XMatrix->GetEntryD(0,4);
     double Rho   = sqrt(Rhox*Rhox + Rhoy*Rhoy);
     double zDest = XMatrix->GetEntryD(0,2), zSource=XMatrix->GetEntryD(0,5);
     double a, c;
     GetacSommerfeld(this, Omega, Rho, zDest, zSource, &a, &c);

     Log("Computing gFrak via Sommerfeld(a=%g,c=%g,xNu=%i,%s)...",a,c,xNu,STStr);
     SommerfeldIntegrate(gFrakIntegrand, (void *)Data, zfdim,
                         a, c, xNu, Rho, qMaxEvalA, qMaxEvalB,
                         qAbsTol, qRelTol, gFrak, Error, Verbose);
   }
  Log("...%i points",Data->NumPoints);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (Data->byqFile) fclose(Data->byqFile);
  delete[] Error;
  if (OwnsWorkspace)
   DestroyScriptGTwiddleWorkspace(Workspace);
  return Data->NumPoints;
}
