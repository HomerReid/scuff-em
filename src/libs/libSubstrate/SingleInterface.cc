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

/***************************************************************/
/* return value is updated prefactor for Phi terms *************/
/***************************************************************/
cdouble AddVSITerm(cdouble k, double Rho, double z, double h,
                   cdouble Eta, bool PPIsOnly,
                   int n, cdouble EtaFac, cdouble *V)
{
  double zn = z + (n==0 ? 0.0 : 2.0*n*h);
  double r2 = Rho*Rho + zn*zn;
  double r=sqrt(r2);
  cdouble ikr=II*k*r;
  double Sign = (n%2) ? -1.0 : 1.0;
  cdouble ExpFac = Sign*ZVAC*exp(ikr) / (4.0*M_PI*r);

  // contributions to V^A_parallel
  if (n<=1) V[_VSI_APAR] += ExpFac;

  // contributions to V^Phi and derivatives
  ExpFac *= EtaFac;
  V[_VSI_PHI] += ExpFac;
  if (!PPIsOnly)
   { ExpFac *= (ikr-1.0)/r2;
     V[_VSI_DRPHI] += Rho*ExpFac;
     if (fabs(z)>1.0e-12) V[_VSI_DZPHI] += zn*ExpFac;
   };

  return Eta*EtaFac;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int GetVSICorrection(cdouble Omega, double Rho, double z,
                     cdouble Eps, double h, cdouble *V,
                     bool PPIsOnly, bool RetainSingularTerms,
                     int MaxTerms, double RelTol, double AbsTol)
{
  cdouble k      = Omega;
  cdouble Eta    = (Eps-1.0)/(Eps+1.0);

  int NumPotentials = (PPIsOnly ? 2 : NUMVSI);
  memset(V, 0, NumPotentials*sizeof(cdouble));

  // n=0 terms (singular at Rho=0)
  if (RetainSingularTerms)
   AddVSITerm(k, Rho, z, h, Eta, PPIsOnly, 0, 1.0-Eta, V);

  if (isinf(h) || MaxTerms==1) return 1;

  // n=1 (first ground-plane image)
  cdouble EtaFac = AddVSITerm(k, Rho, z, h, Eta, PPIsOnly, 1, 1.0-Eta*Eta, V);
  if (MaxTerms==2) return 1;

  /*--------------------------------------------------------------*/
  /*- second and higher ground-plane image terms -----------------*/
  /*--------------------------------------------------------------*/
  double RefVals[3];
  RefVals[0] = abs(V[_VSI_PHI]);
  RefVals[1] = abs(V[_VSI_DRPHI]);
  RefVals[2] = abs(V[_VSI_DZPHI]);
  int ConvergedIters=0;
  for(int nTermPairs=1; nTermPairs<MaxTerms/2; nTermPairs++)
   { 
     cdouble Vn[NUMVSI]={0.0, 0.0, 0.0, 0.0, 0.0};
     EtaFac = AddVSITerm(k, Rho, z, h, Eta, PPIsOnly, 2*nTermPairs,   EtaFac, Vn);
     EtaFac = AddVSITerm(k, Rho, z, h, Eta, PPIsOnly, 2*nTermPairs+1, EtaFac, Vn);
     VecPlusEquals(V, 1.0, Vn, NumPotentials);
    
     bool AllConverged=true;
     for(int p=0; p<(PPIsOnly ? 1 : 3); p++)
      { int Index = (p==0 ? _VSI_PHI : p==1 ? _VSI_DRPHI : _VSI_DZPHI);
        double absTerm = abs(Vn[Index]);
        if ( (absTerm > AbsTol) && (absTerm>RelTol*RefVals[Index]) )
         AllConverged=false;
      }
     if (AllConverged)
      ConvergedIters++;
     else
      ConvergedIters=0;
     if (ConvergedIters==3)
      return 2*nTermPairs;
   }

  return 0; // never get here
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetVSITwiddle(cdouble q, cdouble u0, cdouble u,
                   double Rho, double z, cdouble Eps, double h,
                   cdouble *VTVector, bool PPIsOnly)
{
  if (q==0.0)
   { memset(VTVector, 0, (PPIsOnly ? 2 : NUMVSI)*sizeof(cdouble));
     return;
   };

  // fetch bessel-function factors
  cdouble JdJFactors[2][4];
  bool NeedRhoDerivatives = false;
  int NuMax = (PPIsOnly ? 0 : 1);
  GetJdJFactors(q, Rho, JdJFactors, NeedRhoDerivatives, NuMax);
  cdouble J0Fac = q*JdJFactors[0][0] / (2.0*M_PI);
  cdouble J1Fac = PPIsOnly ? 0.0 : II*q*q*JdJFactors[0][1] / (2.0*M_PI);
  
  // fetch hyperbolic factors
  cdouble uh = (isinf(h) ? 0.0 : u*h), uh2=uh*uh;
  cdouble uTanhFac, uCothFac, SinhFac, CoshFac, dzSinhFac;
  cdouble euz = (z<0.0) ? exp(u*z) : 1.0, emuz=1.0/euz;
  if ( isinf(h) || real(uh)>20.0 )
   { uTanhFac = uCothFac = u;
     SinhFac = CoshFac = euz;
     dzSinhFac = u*euz;
   }
  else if (abs(uh)<1.0e-3)
   { uTanhFac  = uh2/h;
     uCothFac  = (1.0 + uh2/3.0)/h;
     SinhFac   = 1.0 + z/h + u*u*z*(2.0*h*h + 3.0*h*z + z*z)/(6.0*h);
     CoshFac   = 1.0 + (z*h + 0.5*z*z)*u*u;
     dzSinhFac = 1.0/h + u*u*(2.0*h*h + 6.0*h*z + 3.0*z*z)/(6.0*h);
   }
  else
   { cdouble euh = exp(u*h), emuh = 1.0/euh;
     cdouble TanhFac = (euh-emuh)/(euh+emuh);
     uTanhFac  = u*TanhFac;
     uCothFac  = u/TanhFac;
     SinhFac   = (euh*euz - emuh*emuz)/(euh - emuh);
     CoshFac   = (euh*euz + emuh*emuz)/(euh + emuh);
     dzSinhFac = u*(euh*euz + emuh*emuz)/(euh - emuh);
   }
  cdouble Num = u0+uTanhFac, DTE=u0+uCothFac, DTM=Eps*u0+uTanhFac, DTETM=DTE*DTM;
 
  VTVector[_VSI_APAR ] = ZVAC*J0Fac/DTE;
  VTVector[_VSI_PHI  ] = ZVAC*J0Fac*Num/DTETM;
  VTVector[_VSI_APERP] = -1.0*ZVAC*(Eps-1.0)*J1Fac/DTETM;
  VTVector[_VSI_DRPHI] = ZVAC*J1Fac*Num/DTETM;
  VTVector[_VSI_DZPHI] = VTVector[_VSI_PHI];
  if ( z >= 0.0 )
   { cdouble ExpFac = exp(-u0*z);
     VTVector[_VSI_APAR ]*=ExpFac;
     VTVector[_VSI_PHI  ]*=ExpFac;
     VTVector[_VSI_APERP]*=ExpFac;
     VTVector[_VSI_DRPHI]*=ExpFac;
     VTVector[_VSI_DZPHI]*=-1.0*u0*ExpFac;
   }
  else
   { VTVector[_VSI_APAR]  *= SinhFac;
     VTVector[_VSI_PHI ]  *= SinhFac;
     VTVector[_VSI_APERP] *= CoshFac;
     VTVector[_VSI_DRPHI] *= SinhFac;
     VTVector[_VSI_DZPHI] *= u*CoshFac/SinhFac;
   };

  if (fabs(z)<1.0e-12) VTVector[_VSI_DZPHI]=0.0;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct SingleInterfaceData
 { 
   LayeredSubstrate *S;
   cdouble Omega;
   cdouble Eps;
   double h;
   HMatrix *XMatrix;
   bool RetainTerm[2];
   bool PPIsOnly;
   FILE *byqFile; 
   int NumCalls;

 } SingleInterfaceData;

int qIntegrand_SingleInterface(unsigned ndim, const double *x, void *UserData,
                               unsigned fdim, double *fval)
{
  (void) fdim;
  SingleInterfaceData *SIData = (SingleInterfaceData *)UserData;
  LayeredSubstrate *S  = SIData->S;
  cdouble Omega        = SIData->Omega;
  cdouble Eps          = SIData->Eps;
  double h             = SIData->h;
  HMatrix *XMatrix     = SIData->XMatrix;
  bool *RetainTerm     = SIData->RetainTerm;
  bool PPIsOnly        = SIData->PPIsOnly;
  FILE *byqFile        = SIData->byqFile;
  SIData->NumCalls++;

  S->UpdateCachedEpsMu(Omega);

  cdouble q  = x[0] + ((ndim==2) ? II*x[1] : 0.0);
  cdouble q2=q*q, u0 = sqrt(q2 - Omega*Omega), u = sqrt(q2 - Eps*Omega*Omega);

  int NumPotentials = PPIsOnly ? 2 : NUMVSI;
  int NX            = XMatrix->NR;
  HMatrix VTMatrix(NumPotentials, NX, (cdouble *)fval);
  for(int nx=0; nx<NX; nx++)
   { 
     double Rhox    = XMatrix->GetEntryD(nx,0) - XMatrix->GetEntryD(nx,3);
     double Rhoy    = XMatrix->GetEntryD(nx,1) - XMatrix->GetEntryD(nx,4);
     double Rho     = sqrt(Rhox*Rhox + Rhoy*Rhoy);
     double zDest   = XMatrix->GetEntryD(nx,2);
     double zSource = XMatrix->GetEntryD(nx,5);
     double zDelta  = zDest - zSource;
                       
     if (!EqualFloat(zSource, S->zInterface[0]))
      ErrExit("SingleInterface routines require sources on interface");

     cdouble VTVector[NUMVSI]  = {0.0,0.0,0.0,0.0,0.0};
     cdouble VT0Vector[NUMVSI] = {0.0,0.0,0.0,0.0,0.0};
     if(RetainTerm[0])
      GetVSITwiddle(q, u0, u, Rho, zDelta, Eps, h, VTVector, PPIsOnly);
     if(RetainTerm[1])
      GetVSITwiddle(q, u0, u0, Rho, zDelta, Eps, h, VT0Vector, PPIsOnly);

     VT0Vector[_VSI_APERP]=0.0; // we don't do the subtraction for this one
     for(int i=0; i<NumPotentials; i++)
      VTMatrix.SetEntry(i,nx,VTVector[i]-VT0Vector[i]);
 
     if (byqFile)
      { fprintf(byqFile,"%e %e %e %e ",real(q),imag(q),Rho,zDelta);
        fprintVec(byqFile,VTVector,NumPotentials);
        fprintVecCR(byqFile,VT0Vector,NumPotentials);
      }
   };
  
  return 0;
}

/***************************************************************/
/* if !NeedDerivatives:                                        */
/*  VMatrix[{0,1,2},nx] = VA_InPlane, VA_Z, VPhi for point #nx */
/*                                                             */
/* if NeedDerivatives we additionally have                     */
/*  VMatrix[{3,4},nx]   = GradPhi_InPlane, GradPhi_z           */
/***************************************************************/
int LayeredSubstrate::GetSingleInterfacePotentials(cdouble Omega,
                                                   HMatrix *XMatrix,
                                                   HMatrix *VMatrix,
                                                   bool PPIsOnly,
                                                   bool Subtract,
                                                   bool RetainSingularTerms,
                                                   bool CorrectionOnly)
{
  UpdateCachedEpsMu(Omega);

  double Rhox    = XMatrix->GetEntryD(0,0) - XMatrix->GetEntryD(0,3);
  double Rhoy    = XMatrix->GetEntryD(0,1) - XMatrix->GetEntryD(0,4);
  double Rho     = sqrt(Rhox*Rhox + Rhoy*Rhoy);
  double zDest   = XMatrix->GetEntryD(0,2);
  double zSource = XMatrix->GetEntryD(0,5);
  double q0, qR;
  GetacSommerfeld(this, Omega, Rho, zDest, zSource, &q0, &qR);

  int MaxEvalA = qMaxEvalA, MaxEvalB = qMaxEvalB;
  
  char *s=getenv("SOMMERFELD_BYQFILE");
  FILE *byqFile = (s ? fopen(s,"w") : 0);
  int xNu = 0;
  s=getenv("SOMMERFELD_INTEGRATOR_XNU");
  if (s)
   { sscanf(s,"%i",&xNu);
     printf("setting xNu=%i\n",xNu);
   }
  s=getenv("SOMMERFELD_FORCE_SUBTRACT");
  if (s) Subtract = (s[0]=='1');
  s=getenv("SOMMERFELD_MAXEVALA");
  if (s) sscanf(s,"%i",&MaxEvalA);
  s=getenv("SOMMERFELD_MAXEVALB");
  if (s) sscanf(s,"%i",&MaxEvalB);

  cdouble Eps = EpsLayer[1];
  double h    = zInterface[0] - zGP;

  SingleInterfaceData SIData;
  SIData.S             = this;
  SIData.Omega         = Omega;
  SIData.Eps           = Eps;
  SIData.h             = h;
  SIData.XMatrix       = XMatrix;
  SIData.PPIsOnly      = PPIsOnly;
  SIData.RetainTerm[0] = !CorrectionOnly;
  SIData.RetainTerm[1] = (Subtract || CorrectionOnly);
  SIData.byqFile       = byqFile;

  int NX = XMatrix->NR;
  int NumPotentials = PPIsOnly ? 2 : NUMVSI;
  int zfdim = NumPotentials * NX;
  cdouble *Error = new cdouble[zfdim]; 
  SIData.NumCalls=0;
  SommerfeldIntegrate(qIntegrand_SingleInterface, (void *)&SIData, zfdim,
                      q0, qR, xNu, Rho, MaxEvalA, MaxEvalB,
                      qAbsTol, qRelTol, VMatrix->ZM, Error);
  if (byqFile) fclose(byqFile);
  if (CorrectionOnly) VMatrix->Scale(-1.0);
  delete[] Error;

  bool NeedCorrection = Subtract || (!Subtract && !RetainSingularTerms);
  if (NeedCorrection)
   for(int nx=0; nx<NX; nx++)
    { Rhox     = XMatrix->GetEntryD(nx,0) - XMatrix->GetEntryD(nx,3);
      Rhoy     = XMatrix->GetEntryD(nx,1) - XMatrix->GetEntryD(nx,4);
      Rho      = sqrt(Rhox*Rhox + Rhoy*Rhoy);
      double z = XMatrix->GetEntryD(nx,2) - XMatrix->GetEntryD(nx,5);
      cdouble VCorr[NUMVSI];
      cdouble *V=(cdouble *)VMatrix->GetColumnPointer(nx);
      int NumTerms;
      if (Subtract)
       { NumTerms=GetVSICorrection(Omega, Rho, z, Eps, h, VCorr, PPIsOnly,
                                   RetainSingularTerms, 1000, qRelTol, qAbsTol);
         VecPlusEquals(V,1.0,VCorr,NumPotentials);
       }
      else if (!Subtract && !RetainSingularTerms)
       { NumTerms=GetVSICorrection(Omega, Rho, z, Eps, h, VCorr, PPIsOnly,
                                   true, 1, qRelTol, qAbsTol);
         VecPlusEquals(V,-1.0,VCorr,NumPotentials);
       }
      Log("VSICorrection({w,Rho,z}={%g,%g,%g}: %i terms summed",real(Omega),Rho,z,NumTerms);
    }
    
  return SIData.NumCalls;
}

int LayeredSubstrate::GetSingleInterfacePotentials(cdouble Omega, double Rho,
                                                   double zDest, cdouble *VSI,
                                                   bool PPIsOnly, bool Subtract,
                                                   bool RetainSingularTerms,
                                                   bool CorrectionOnly)
{
  double XDS[6]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  XDS[0]=Rho;
  XDS[2]=zDest;
  HMatrix XMatrix(1,6,XDS);

  int NumPotentials = (PPIsOnly ? 2 : NUMVSI);
  HMatrix VMatrix(NumPotentials, 1, VSI);

  return GetSingleInterfacePotentials(Omega, &XMatrix, &VMatrix, PPIsOnly, Subtract,
                                      RetainSingularTerms, CorrectionOnly);
}
