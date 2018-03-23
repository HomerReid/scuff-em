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
 * MOI.cc -- specialized routines in libSubstrate for handling the
 *        -- case of metal directly atop an insulating substrate
 *        -- (either infinite-thickness or grounded)
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
cdouble AddSGFTerm_MOI(cdouble k, double Rho, double z, double h,
                       cdouble Eta, bool PPIsOnly,
                       int n, cdouble EtaFac, cdouble *V)
{
  double zn = z + (n==0 ? 0.0 : 2.0*n*h);
  double r2 = Rho*Rho + zn*zn;
  if (r2==0.0)
   return 0.0;
  double r=sqrt(r2);
  cdouble ikr=II*k*r;
  double Sign = (n%2) ? -1.0 : 1.0;
  cdouble ExpFac = Sign*exp(ikr) / (4.0*M_PI*r);

  // contributions to V^A_parallel
  if (n<=1) V[_SGF_APAR] += ExpFac;

  // contributions to V^Phi and derivatives
  ExpFac *= EtaFac;
  V[_SGF_PHI] += ExpFac;
  if (!PPIsOnly)
   { ExpFac *= (ikr-1.0)/r2;
     V[_SGF_DRPHI] += Rho*ExpFac;
     V[_SGF_DZPHI] += zn*ExpFac;
   }

  return Eta*EtaFac;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int GetSGFCorrection_MOI(cdouble Omega, double Rho, double z,
                         cdouble Eps, double h, cdouble *V,
                         bool PPIsOnly, bool RetainSingularTerms,
                         int MaxTerms, double RelTol, double AbsTol)
{
  cdouble k      = Omega;
  cdouble Eta    = (Eps-1.0)/(Eps+1.0);

  int NumSGFs = (PPIsOnly ? 2 : NUMSGFS_MOI);
  memset(V, 0, NumSGFs*sizeof(cdouble));

  // n=0 terms (singular at Rho=0)
  if (RetainSingularTerms)
   AddSGFTerm_MOI(k, Rho, z, h, Eta, PPIsOnly, 0, 1.0-Eta, V);

  if (isinf(h) || MaxTerms==1) return 1;

  // n=1 (first ground-plane image)
  cdouble EtaFac = AddSGFTerm_MOI(k, Rho, z, h, Eta, PPIsOnly, 1, 1.0-Eta*Eta, V);
  if (MaxTerms==2) return 1;

  /*--------------------------------------------------------------*/
  /*- second and higher ground-plane image terms -----------------*/
  /*--------------------------------------------------------------*/
  double RefVals[3];
  RefVals[0] = abs(V[_SGF_PHI]);
  RefVals[1] = abs(V[_SGF_DRPHI]);
  RefVals[2] = abs(V[_SGF_DZPHI]);
  int ConvergedIters=0;
  for(int nTermPairs=1; nTermPairs<MaxTerms/2; nTermPairs++)
   { 
     cdouble Vn[NUMSGFS_MOI]={0.0, 0.0, 0.0, 0.0, 0.0};
     EtaFac = AddSGFTerm_MOI(k, Rho, z, h, Eta, PPIsOnly, 2*nTermPairs,   EtaFac, Vn);
     EtaFac = AddSGFTerm_MOI(k, Rho, z, h, Eta, PPIsOnly, 2*nTermPairs+1, EtaFac, Vn);
     VecPlusEquals(V, 1.0, Vn, NumSGFs);
    
     bool AllConverged=true;
     for(int p=0; p<(PPIsOnly ? 1 : 3); p++)
      { int Index = (p==0 ? _SGF_PHI : p==1 ? _SGF_DRPHI : _SGF_DZPHI);
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
void GetVTwiddle_MOI(cdouble q, cdouble u0, cdouble u,
                     double Rho, double z, cdouble Eps, double h,
                     cdouble *VTVector, bool PPIsOnly)
{
  if (q==0.0)
   { memset(VTVector, 0, (PPIsOnly ? 2 : NUMSGFS_MOI)*sizeof(cdouble));
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
     CoshFac   = 1.0 + (h + 0.5*z)*z*u*u;
     dzSinhFac = 1.0/h + u*u*(2.0*h*h + 6.0*h*z + 3.0*z*z)/(6.0*h);
   }
  else
   { cdouble euh = exp(uh), emuh = 1.0/euh;
     cdouble TanhFac = (euh-emuh)/(euh+emuh);
     uTanhFac  = u*TanhFac;
     uCothFac  = u/TanhFac;
     SinhFac   = (euh*euz - emuh*emuz)/(euh - emuh);
     CoshFac   = (euh*euz + emuh*emuz)/(euh + emuh);
     dzSinhFac = u*(euh*euz + emuh*emuz)/(euh - emuh);
   }
  cdouble Num   =     u0 + uTanhFac;
  cdouble DTE   =     u0 + uCothFac;
  cdouble DTM   = Eps*u0 + uTanhFac;
  cdouble DTETM = DTE*DTM;
 
  VTVector[_SGF_APAR ] = J0Fac/DTE;
  VTVector[_SGF_PHI  ] = J0Fac*Num/DTETM;
  VTVector[_SGF_AZ   ] = -1.0*(Eps-1.0)*J1Fac/DTETM;
  VTVector[_SGF_DRPHI] = J1Fac*Num/DTETM;
  VTVector[_SGF_DZPHI] = VTVector[_SGF_PHI];
  if ( z >= 0.0 )
   { cdouble ExpFac = exp(-u0*z);
     VTVector[_SGF_APAR ] *= ExpFac;
     VTVector[_SGF_PHI  ] *= ExpFac;
     VTVector[_SGF_AZ   ] *= ExpFac;
     VTVector[_SGF_DRPHI] *= ExpFac;
     VTVector[_SGF_DZPHI] *= -1.0*u0*ExpFac;
   }
  else
   { VTVector[_SGF_APAR ] *= SinhFac;
     VTVector[_SGF_PHI  ] *= SinhFac;
     VTVector[_SGF_AZ   ] *= CoshFac;
     VTVector[_SGF_DRPHI] *= SinhFac;
     VTVector[_SGF_DZPHI] *= dzSinhFac;
   };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct MOIData
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

 } MOIData;

int qIntegrand_MOISGFs(unsigned ndim, const double *x, void *UserData,
                       unsigned fdim, double *fval)
{
  (void) fdim;
  MOIData *Data        = (MOIData *)UserData;
  LayeredSubstrate *S  = Data->S;
  cdouble Omega        = Data->Omega;
  cdouble Eps          = Data->Eps;
  double h             = Data->h;
  HMatrix *XMatrix     = Data->XMatrix;
  bool *RetainTerm     = Data->RetainTerm;
  bool PPIsOnly        = Data->PPIsOnly;
  FILE *byqFile        = Data->byqFile;
  Data->NumCalls++;

  S->UpdateCachedEpsMu(Omega);

  cdouble q  = x[0] + ((ndim==2) ? II*x[1] : 0.0);
  cdouble q2=q*q, u0 = sqrt(q2 - Omega*Omega), u = sqrt(q2 - Eps*Omega*Omega);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
if ( real(u0)<0.0 || real(u)<0.0 )
 printf("q2,w2,u0,u={%e,%e} {%e} {%e,%e} {%e,%e}\n",real(q2),imag(q2),real(Omega),real(u0),imag(u0),real(u),imag(u));
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  int NumSGFs = PPIsOnly ? 2 : NUMSGFS_MOI;
  int NX      = XMatrix->NR;
  HMatrix VTMatrix(NumSGFs, NX, (cdouble *)fval);
  for(int nx=0; nx<NX; nx++)
   { 
     double Rhox    = XMatrix->GetEntryD(nx,0) - XMatrix->GetEntryD(nx,3);
     double Rhoy    = XMatrix->GetEntryD(nx,1) - XMatrix->GetEntryD(nx,4);
     double Rho     = sqrt(Rhox*Rhox + Rhoy*Rhoy);
     double zDest   = XMatrix->GetEntryD(nx,2);
     double zSource = XMatrix->GetEntryD(nx,5);
     double zDelta  = zDest - zSource;

     if (!EqualFloat(zSource, S->zInterface[0]))
      ErrExit("MOI routines require sources on interface");

     cdouble VTVector[NUMSGFS_MOI]  = {0.0,0.0,0.0,0.0,0.0};
     cdouble VT0Vector[NUMSGFS_MOI] = {0.0,0.0,0.0,0.0,0.0};
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
#if 0
static int Howdage=0;
Howdage++;
printf("%10i %e %e {%e,%e} ",Howdage,Rho,zDelta,real(q),imag(q));
#endif
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
     if(RetainTerm[0])
      GetVTwiddle_MOI(q, u0, u, Rho, zDelta, Eps, h, VTVector, PPIsOnly);
     if(RetainTerm[1] && (zDelta > -1.0e-10) )
      GetVTwiddle_MOI(q, u0, u0, Rho, zDelta, Eps, h, VT0Vector, PPIsOnly);

     VT0Vector[_SGF_AZ]=0.0; // we don't do the subtraction for this one
     for(int i=0; i<NumSGFs; i++)
      VTMatrix.SetEntry(i,nx,VTVector[i]-VT0Vector[i]);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
#if 0
for(int i=0; i<NumSGFs; i++)
 printf("{%+.2e,%+.2e}-{%+.2e,%+.2e} ",real(VTVector[i]),imag(VTVector[i]),real(VT0Vector[i]),imag(VT0Vector[i]));
printf("\n");
#endif
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
 
     if (byqFile)
      { fprintf(byqFile,"%e %e %e %e ",real(q),imag(q),Rho,zDelta);
        fprintVec(byqFile,VTVector,NumSGFs);
        fprintVecCR(byqFile,VT0Vector,NumSGFs);
      }
   };
  
  return 0;
}

/***************************************************************/
/* VMatrix[{0,1},nx] = VA_Parallel, VPhi for point #nx         */
/*                                                             */
/* if !PPIsOnly we additionally have                           */
/*  VMatrix[{2,3,4},nx] = VA_z, dPhi/dRho, dPhi/dz             */
/***************************************************************/
int LayeredSubstrate::GetScalarGFs_MOI(cdouble Omega,
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
  if (fabs(zSource-zInterface[0])>1.0e-12)
   ErrExit("MOI routines require sources on dielectric interface");
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

  MOIData Data;
  Data.S             = this;
  Data.Omega         = Omega;
  Data.Eps           = Eps;
  Data.h             = h;
  Data.XMatrix       = XMatrix;
  Data.PPIsOnly      = PPIsOnly;
  Data.RetainTerm[0] = !CorrectionOnly;
  Data.RetainTerm[1] = (Subtract || CorrectionOnly);
  Data.byqFile       = byqFile;

  int NX = XMatrix->NR;
  int NumSGFs = PPIsOnly ? 2 : NUMSGFS_MOI;
  int zfdim = NumSGFs * NX;
  cdouble *Error = new cdouble[zfdim];
  Data.NumCalls=0;
  SommerfeldIntegrate(qIntegrand_MOISGFs, (void *)&Data, zfdim,
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
      if (z < -1.0e-10) continue; // don't do subtraction for eval points in dielectric
      cdouble VCorr[NUMSGFS_MOI];
      cdouble *V=(cdouble *)VMatrix->GetColumnPointer(nx);
      int NumTerms=0;
      if (Subtract)
       { NumTerms=GetSGFCorrection_MOI(Omega, Rho, z, Eps, h, VCorr, PPIsOnly,
                                       RetainSingularTerms, 1000, qRelTol, qAbsTol);
         VecPlusEquals(V,1.0,VCorr,NumSGFs);
       }
      else if (!Subtract && !RetainSingularTerms)
       { NumTerms=GetSGFCorrection_MOI(Omega, Rho, z, Eps, h, VCorr, PPIsOnly,
                                       true, 1, qRelTol, qAbsTol);
         VecPlusEquals(V,-1.0,VCorr,NumSGFs);
       }
      Log("VSICorrection({w,Rho,z}={%g,%g,%g}: %i terms summed",real(Omega),Rho,z,NumTerms);
    }
    
  return Data.NumCalls;
}

int LayeredSubstrate::GetScalarGFs_MOI(cdouble Omega, double Rho,
                                       double zDest, cdouble *V,
                                       bool PPIsOnly, bool Subtract,
                                       bool RetainSingularTerms,
                                       bool CorrectionOnly)
{
  double XDS[6]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  XDS[0]=Rho;
  XDS[2]=zDest;
  HMatrix XMatrix(1,6,XDS);

  int NumSGFs = (PPIsOnly ? 2 : NUMSGFS_MOI);
  HMatrix VMatrix(NumSGFs, 1, V);

  return GetScalarGFs_MOI(Omega, &XMatrix, &VMatrix, PPIsOnly, Subtract,
                          RetainSingularTerms, CorrectionOnly);
}
