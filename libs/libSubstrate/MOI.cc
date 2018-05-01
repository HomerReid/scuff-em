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

void InitScalarGFOptions(ScalarGFOptions *Options)
{ 
  Options->PPIsOnly            = false;
  Options->Subtract            = true;
  Options->RetainSingularTerms = false;
  Options->CorrectionOnly      = false;
  Options->UseInterpolator     = true;
  Options->NeedDerivatives     = 0;
  Options->MaxTerms            = 1000;
  Options->RelTol              = 1.0e-6;
  Options->AbsTol              = 1.0e-12;
}

static const ScalarGFOptions DefaultSGFOptions
 = {false, true, false, false, true, 0, 1000, 1.0e-6, 1.0e-12};

/***************************************************************/
/* Add a single term to the series defining the ScalarGF       */
/* corrections.                                                */
/* Return value is updated prefactor for Phi terms.            */
/***************************************************************/
cdouble AddSGFTerm_MOI(const ScalarGFOptions *Options,
                       cdouble k, double Rho, double z, double h,
                       cdouble Eta, int n, cdouble EtaFac, cdouble *VVector)
{ 
  if (Options==0) Options = &DefaultSGFOptions;
  bool PPIsOnly = Options->PPIsOnly;
  int NumSGFs   = (PPIsOnly ? 2 : NUMSGFS_MOI);
  int dRhoMax   = (Options->NeedDerivatives & NEED_DRHO) ? 1 : 0;
  int dzMax     = (Options->NeedDerivatives & NEED_DZ) ? 1 : 0;

  double zn = z + (n==0 ? 0.0 : 2.0*n*h);
  double r2 = Rho*Rho + zn*zn;
  if (r2==0.0)
   return 0.0;
  double r=sqrt(r2);
  cdouble ikr=II*k*r, ikr2=ikr*ikr, ikr3=ikr2*ikr;
  double Sign = (n%2) ? -1.0 : 1.0;
  cdouble ExpFac[4];
  ExpFac[0] = Sign*exp(ikr) / (4.0*M_PI*r);
  ExpFac[1] = ExpFac[0] * (ikr - 1.0)/r2;
  ExpFac[2] = ExpFac[0] * (ikr2 - 3.0*ikr + 3.0)/(r2*r2);
  ExpFac[3] = ExpFac[0] * (ikr3 - 6.0*ikr2 + 15.0*ikr - 15.0)/(r2*r2*r2);

  for(int dz=0, nvd=0; dz<=dzMax; dz++)
   for(int dRho=0; dRho<=dRhoMax; dRho++, nvd++)
    { 
      cdouble *V = VVector + nvd*NumSGFs;
      double zFac = (dz ? zn : 1.0), RhoFac = (dRho ? Rho : 1.0), zRho=zFac*RhoFac;
      int Index = dRho+dz;

      // contributions to V^A_parallel
      if (n<=1) V[_SGF_APAR] += zRho*ExpFac[Index];

      // contributions to V^Phi
      V[_SGF_PHI] += EtaFac*zRho*ExpFac[Index];

      if (PPIsOnly) continue;

      // contributions to V^Phi derivatives
      V[_SGF_DRPHI] += EtaFac*(zRho*Rho*ExpFac[Index+1] + dRho*zFac*ExpFac[Index]);
      V[_SGF_DZPHI] += EtaFac*(zRho*zn *ExpFac[Index+1] + dz*RhoFac*ExpFac[Index]);
    }

  return Eta*EtaFac;
}

/***************************************************************/
/* Compute the real-space correction to the scalar GFs, i.e.   */
/* the result of integrating VTwiddle(q,u0,u0) over all q.     */
/***************************************************************/
int GetSGFCorrection_MOI(cdouble Omega, double Rho, double z,
                         cdouble Eps, double h, cdouble *V,
                         const ScalarGFOptions *Options)
{
  if (!Options) Options    = &DefaultSGFOptions;
  bool PPIsOnly            = Options->PPIsOnly;
  bool RetainSingularTerms = Options->RetainSingularTerms;
  int NeedDerivatives      = Options->NeedDerivatives;
  int MaxTerms             = Options->MaxTerms;
  double RelTol            = Options->RelTol;
  double AbsTol            = Options->AbsTol;

  cdouble k   = Omega;
  cdouble Eta = (Eps-1.0)/(Eps+1.0);

  int NumSGFs = (PPIsOnly ? 2 : NUMSGFS_MOI);
  int NumSGFVDs = NumSGFs;
  if (NeedDerivatives & NEED_DRHO) NumSGFVDs*=2; 
  if (NeedDerivatives & NEED_DZ  ) NumSGFVDs*=2; 
  memset(V, 0, NumSGFs*sizeof(cdouble));

  // n=0 terms (singular at Rho=0)
  if (RetainSingularTerms)
   AddSGFTerm_MOI(Options, k, Rho, z, h, Eta, 0, 1.0-Eta, V);

  if (isinf(h) || MaxTerms==1) return 1;

  // n=1 (first ground-plane image)
  cdouble EtaFac = AddSGFTerm_MOI(Options, k, Rho, z, h, Eta, 1, 1.0-Eta*Eta, V);

  if (Eta==0.0 || MaxTerms==2) return 1;

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
     cdouble Vn[4*NUMSGFS_MOI];
     memset(Vn, 0, NumSGFVDs);
     EtaFac = AddSGFTerm_MOI(Options, k, Rho, z, h, Eta, 2*nTermPairs,   EtaFac, Vn);
     EtaFac = AddSGFTerm_MOI(Options, k, Rho, z, h, Eta, 2*nTermPairs+1, EtaFac, Vn);
     VecPlusEquals(V, 1.0, Vn, NumSGFVDs);
    
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
/* Compute the integrand of the integral over q that defines   */
/* the scalar GFs.                                             */
/***************************************************************/
void GetVTwiddle_MOI(cdouble q, cdouble u0, cdouble u,
                     double Rho, double z, cdouble Eps, double h,
                     cdouble *VTVector, const ScalarGFOptions *Options=0)
{
  if (!Options) Options  = &DefaultSGFOptions;
  bool PPIsOnly          = Options->PPIsOnly;
  bool NeedRhoDerivative = (Options->NeedDerivatives & NEED_DRHO) ? true : false;
  bool NeedzDerivative   = (Options->NeedDerivatives & NEED_DZ) ? true : false;

  int NumSGFs  = (PPIsOnly ? 2 : NUMSGFS_MOI);
  int NumVDs   = (NeedRhoDerivative ? 2 : 1) * (NeedzDerivative ? 2 : 1);
  int NumSGFVDs = NumSGFs * NumVDs;

  if (q==0.0)
   { memset(VTVector, 0, NumSGFVDs*sizeof(cdouble));
     return;
   };

  // fetch bessel-function factors
  cdouble JdJFactors[2][4];
  int NuMax = (PPIsOnly ? 0 : 1);
  GetJdJFactors(q, Rho, JdJFactors, NeedRhoDerivative, NuMax);
  
  // fetch hyperbolic factors
  cdouble uh = (isinf(h) ? 0.0 : u*h), uh2=uh*uh;
  cdouble uTanhFac, uCothFac, SinhFac, CoshFac, dzSinhFac, dzCoshFac, d2zSinhFac;
  cdouble euz = (z<0.0) ? exp(u*z) : 1.0, emuz=1.0/euz;
  if ( isinf(h) || real(uh)>20.0 )
   { uTanhFac = uCothFac = u;
     SinhFac = CoshFac = euz;
     dzSinhFac = dzCoshFac = u*euz;
     d2zSinhFac= u*u*euz;
   }
  else if (abs(uh)<1.0e-3)
   { uTanhFac   = uh2/h;
     uCothFac   = (1.0 + uh2/3.0)/h;
     SinhFac    = 1.0 + z/h + u*u*z*(2.0*h*h + 3.0*h*z + z*z)/(6.0*h);
     CoshFac    = 1.0 + (h + 0.5*z)*z*u*u;
     dzSinhFac  = 1.0/h + u*u*(2.0*h*h + 6.0*h*z + 3.0*z*z)/(6.0*h);
     dzCoshFac  = (h+z)*u*u;
     d2zSinhFac = u*u*(1.0 + z/h);
   }
  else
   { cdouble euh = exp(uh), emuh = 1.0/euh;
     cdouble TanhFac = (euh-emuh)/(euh+emuh);
     uTanhFac   = u*TanhFac;
     uCothFac   = u/TanhFac;
     SinhFac    = (euh*euz - emuh*emuz)/(euh - emuh);
     CoshFac    = (euh*euz + emuh*emuz)/(euh + emuh);
     dzSinhFac  = u*(euh*euz + emuh*emuz)/(euh - emuh);
     dzCoshFac  = u*(euh*euz - emuh*emuz)/(euh + emuh);
     d2zSinhFac = u*u*SinhFac;
   }
  cdouble Num   =     u0 + uTanhFac;
  cdouble DTE   =     u0 + uCothFac;
  cdouble DTM   = Eps*u0 + uTanhFac;
  cdouble DTETM = DTE*DTM;

  cdouble zFactors[2][NUMSGFS_MOI];
  if (z>=0.0)
   { cdouble ExpFac=exp(-u0*z);
     for(int n=0; n<NUMSGFS_MOI; n++)
      { zFactors[0][n] = (n==_SGF_DZPHI ? -u0 : 1.0)*ExpFac; 
        zFactors[1][n] = -u0*zFactors[0][n];
      }
   }
  else
   { zFactors[0][_SGF_APAR ] = zFactors[0][_SGF_PHI] = zFactors[0][_SGF_DRPHI] = SinhFac;
     zFactors[0][_SGF_AZ   ] = CoshFac;
     zFactors[0][_SGF_DZPHI] = dzSinhFac;
     zFactors[1][_SGF_APAR ] = zFactors[1][_SGF_PHI] = zFactors[1][_SGF_DRPHI] = dzSinhFac;
     zFactors[1][_SGF_AZ   ] = dzCoshFac;
     zFactors[1][_SGF_DZPHI] = d2zSinhFac;
   }

  cdouble RhoFactors[2][NUMSGFS_MOI];
  for(int dRho=0; dRho <= (NeedRhoDerivative ? 1 : 0); dRho++)
   { cdouble J0Fac = q*JdJFactors[dRho][0] / (2.0*M_PI);
     cdouble J1Fac = II*q*q*JdJFactors[dRho][1] / (2.0*M_PI);
     RhoFactors[dRho][_SGF_APAR ] = J0Fac / DTE;
     RhoFactors[dRho][_SGF_PHI  ] = J0Fac * Num/DTETM;
     RhoFactors[dRho][_SGF_AZ   ] = J1Fac * (1.0-Eps)/DTETM;
     RhoFactors[dRho][_SGF_DRPHI] = J1Fac * Num/DTETM;
     RhoFactors[dRho][_SGF_DZPHI] = J0Fac * Num/DTETM;
   }
 
  for(int dz=0, nvd=0; dz<=(NeedzDerivative?1:0); dz++)
   for(int dRho=0; dRho<=(NeedRhoDerivative?1:0); dRho++, nvd++)
    for(int nSGF=0; nSGF<NumSGFs; nSGF++)
     VTVector[ nvd*NumSGFs + nSGF ] = zFactors[dz][nSGF] * RhoFactors[dRho][nSGF];
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
   ScalarGFOptions Options;
   FILE *byqFile;
   int NumCalls;

 } MOIData;

int qIntegrand_MOISGFs(unsigned ndim, const double *x, void *UserData,
                       unsigned fdim, double *fval)
{
  (void) fdim;
  MOIData *Data            = (MOIData *)UserData;
  LayeredSubstrate *S      = Data->S;
  cdouble Omega            = Data->Omega;
  cdouble Eps              = Data->Eps;
  double h                 = Data->h;
  HMatrix *XMatrix         = Data->XMatrix;
  ScalarGFOptions *Options = &(Data->Options);
  FILE *byqFile            = Data->byqFile;
  Data->NumCalls++;

  bool PPIsOnly          = Options->PPIsOnly;
  bool Subtract          = Options->Subtract;
  bool CorrectionOnly    = Options->CorrectionOnly;
  bool NeedRhoDerivative = (Options->NeedDerivatives & NEED_DRHO) ? true : false;
  bool NeedzDerivative   = (Options->NeedDerivatives & NEED_DZ) ? true : false;

  // Term 0 is the full integrand,    VTwiddle(q,u0,u).
  // Term 1 is minus the correction, -VTwiddle(q,u0,u0).
  bool RetainTerm[2];
  RetainTerm[0] = !CorrectionOnly;
  RetainTerm[1] = Subtract || CorrectionOnly;

  S->UpdateCachedEpsMu(Omega);

  cdouble q  = x[0] + ((ndim==2) ? II*x[1] : 0.0);
  cdouble q2=q*q, u0 = sqrt(q2 - Omega*Omega), u = sqrt(q2 - Eps*Omega*Omega);

  int NumSGFs   = (PPIsOnly ? 2 : NUMSGFS_MOI);
  int NumVDs    = (NeedRhoDerivative ? 2 : 1) * (NeedzDerivative ? 2 : 1);
  int NumSGFVDs = NumSGFs * NumVDs;
  int NX        = XMatrix->NR;
  HMatrix VTMatrix(NumSGFVDs, NX, (cdouble *)fval);
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

     cdouble VTVector[4*NUMSGFS_MOI];  memset(VTVector, 0, 4*NUMSGFS_MOI*sizeof(cdouble));
     cdouble VT0Vector[4*NUMSGFS_MOI]; memset(VT0Vector, 0, 4*NUMSGFS_MOI*sizeof(cdouble));
     if(RetainTerm[0])
      GetVTwiddle_MOI(q, u0, u, Rho, zDelta, Eps, h, VTVector, Options);
     if(RetainTerm[1] && (zDelta > -1.0e-10) )
      { GetVTwiddle_MOI(q, u0, u0, Rho, zDelta, Eps, h, VT0Vector, Options);
        if (!PPIsOnly)
         for(int nvd=0; nvd<NumVDs; nvd++)
          VT0Vector[nvd*NumSGFs + _SGF_AZ]=0.0; // we don't do the subtraction for this one
      }

     for(int i=0; i<NumSGFVDs; i++)
      VTMatrix.SetEntry(i,nx,VTVector[i]-VT0Vector[i]);
 
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
                                       const ScalarGFOptions *Options)
{
  if (!Options) Options    = &DefaultSGFOptions;
  bool PPIsOnly            = Options->PPIsOnly;
  bool Subtract            = Options->Subtract;
  bool RetainSingularTerms = Options->RetainSingularTerms;
  bool CorrectionOnly      = Options->CorrectionOnly;
  bool NeedRhoDerivative   = (Options->NeedDerivatives & NEED_DRHO) ? true : false;
  bool NeedzDerivative     = (Options->NeedDerivatives & NEED_DZ) ? true : false;

  UpdateCachedEpsMu(Omega);
  
  /*--------------------------------------------------------------*/
  /*- look at the first point in the list to determine contour    */
  /*- integral parameters for SommerfeldIntegrator                */
  /*--------------------------------------------------------------*/
  double Rhox    = XMatrix->GetEntryD(0,0) - XMatrix->GetEntryD(0,3);
  double Rhoy    = XMatrix->GetEntryD(0,1) - XMatrix->GetEntryD(0,4);
  double Rho     = sqrt(Rhox*Rhox + Rhoy*Rhoy);
  double zDest   = XMatrix->GetEntryD(0,2);
  double zSource = XMatrix->GetEntryD(0,5);
  if (fabs(zSource-zInterface[0])>1.0e-12)
   ErrExit("MOI routines require sources on dielectric interface");
  double q0, qR;
  GetacSommerfeld(this, Omega, Rho, zDest, zSource, &q0, &qR);

  /*--------------------------------------------------------------*/
  /*- look at environment variables to configure further options  */
  /*- FIXME                                                       */
  /*--------------------------------------------------------------*/
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

  /***************************************************************/
  /* If Subtract==true and Eps==1.0, i.e. the substrate consists */
  /* of (at most) just a ground plane, we can skip the q integral*/
  /* the subtracted term accounts for the entire GF in this case.*/
  /***************************************************************/
  bool SkipqIntegral = (Eps==1.0 && Subtract);

  MOIData Data;
  Data.S             = this;
  Data.Omega         = Omega;
  Data.Eps           = Eps;
  Data.h             = h;
  Data.XMatrix       = XMatrix;
  Data.byqFile       = byqFile;
  memcpy(&(Data.Options), Options, sizeof(ScalarGFOptions));

  int NumSGFs   = (PPIsOnly ? 2 : NUMSGFS_MOI);
  int NumVDs    = (NeedRhoDerivative ? 2 : 1) * (NeedzDerivative ? 2 : 1);
  int NumSGFVDs = NumSGFs * NumVDs;
  int NX        = XMatrix->NR;
  int zfdim     = NumSGFVDs * NX;
  cdouble *Error = new cdouble[zfdim];
  Data.NumCalls=0;
  if (SkipqIntegral)
   VMatrix->Zero();
  else 
   SommerfeldIntegrate(qIntegrand_MOISGFs, (void *)&Data, zfdim,
                       q0, qR, xNu, Rho, MaxEvalA, MaxEvalB,
                       qAbsTol, qRelTol, VMatrix->ZM, Error);
  if (byqFile) fclose(byqFile);
  if (CorrectionOnly) VMatrix->Scale(-1.0);
  delete[] Error;

  /***************************************************************/
  /* if we subtracted the near-field term from the Fourier-space */
  /* integrand, compute and add that term in real space.         */
  /* Alternatively, if we *didn't* subtract in the Fourier       */
  /* integral, but we were asked to exclude the most singular    */
  /* terms from the potentials, then we must evaluate and        */
  /* subtract those.                                             */
  /***************************************************************/
  bool NeedCorrection = Subtract || (!Subtract && !RetainSingularTerms);
  if (NeedCorrection)
   for(int nx=0; nx<NX; nx++)
    { Rhox     = XMatrix->GetEntryD(nx,0) - XMatrix->GetEntryD(nx,3);
      Rhoy     = XMatrix->GetEntryD(nx,1) - XMatrix->GetEntryD(nx,4);
      Rho      = sqrt(Rhox*Rhox + Rhoy*Rhoy);
      double z = XMatrix->GetEntryD(nx,2) - XMatrix->GetEntryD(nx,5);
      cdouble VCorr[4*NUMSGFS_MOI];
      cdouble *V=(cdouble *)VMatrix->GetColumnPointer(nx);
      if (Subtract)
       { if (z < -1.0e-10) 
          ErrExit("%s:%i: internal error (singularity subtraction not implemented for points in dielectric)",__FILE__,__LINE__);
         GetSGFCorrection_MOI(Omega, Rho, z, Eps, h, VCorr, Options);
         VecPlusEquals(V,1.0,VCorr,NumSGFVDs);
       }
      else if (!Subtract && !RetainSingularTerms)
       { 
         ScalarGFOptions OptionsPrime;
         memcpy(&OptionsPrime, Options, sizeof(ScalarGFOptions));
         OptionsPrime.RetainSingularTerms=true;
         OptionsPrime.MaxTerms=1;
         GetSGFCorrection_MOI(Omega, Rho, z, Eps, h, VCorr, &OptionsPrime);
         VecPlusEquals(V,-1.0,VCorr,NumSGFVDs);
       }
    }
    
  return Data.NumCalls;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
bool LayeredSubstrate::GetScalarGFs_Interp(cdouble Omega, double Rho, double z,
                                           cdouble *V, const ScalarGFOptions *Options)
{ 
  /*--------------------------------------------------------------*/
  /*- checks to determine whether the interpolator is usable     -*/
  /*--------------------------------------------------------------*/
  if (!Options->UseInterpolator) return false;
  if (!ScalarGFInterpolator) return false;
  if (Omega!=OmegaCache) return false;
 
  dVec XMin = ScalarGFInterpolator->XMin;
  dVec DX   = ScalarGFInterpolator->DX;
  iVec NVec = ScalarGFInterpolator->NVec;
  int D     = XMin.size();
  
  double RhoMin = XMin[0];
  double RhoMax = RhoMin + (NVec[0]-1)*DX[0];
  if ( Rho < RhoMin || Rho>RhoMax ) return false;

  if (D==1)
   { 
     if (!EqualFloat(z, zSGFI)) return false;
   }
  else // (D>1)
   { double zMin = XMin[1];
     double zMax = zMin + (NVec[1]-1)*DX[1];
     if ( z < zMin || z>zMax ) return false;
   }

  int NumSGFs      = (Options->PPIsOnly    ? 2 : NUMSGFS_MOI);
  int NumSGFsTable = (SGFIOptions.PPIsOnly ? 2 : NUMSGFS_MOI);
  if( NumSGFs > NumSGFsTable )
   return false;
  if( SGFIOptions.Subtract            != Options->Subtract)
   return false;
  if( SGFIOptions.RetainSingularTerms != Options->RetainSingularTerms )
   return false;
  if( SGFIOptions.CorrectionOnly      != Options->CorrectionOnly )
   return false;
  if( (Options->NeedDerivatives != 0) )
   return false;

  /*--------------------------------------------------------------*/
  /*- checks passed: get data by interpolation -------------------*/
  /*--------------------------------------------------------------*/
  double RhoZ[2];
  RhoZ[0] = Rho;
  RhoZ[1] = z;
  cdouble PhiVD[NUMSGFS_MOI];
  ScalarGFInterpolator->Evaluate(RhoZ,(double *)PhiVD);
  memcpy(V, PhiVD, NumSGFs * sizeof(cdouble));
  return true;
}

int LayeredSubstrate::GetScalarGFs_MOI(cdouble Omega, double Rho,
                                       double zDest, cdouble *V,
                                       const ScalarGFOptions *Options)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if ( GetScalarGFs_Interp(Omega, Rho, zDest, V, Options) )
   return 0;

  double XDS[6]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  XDS[0]=Rho;
  XDS[2]=zDest;
  HMatrix XMatrix(1,6,XDS);

  int NumSGFVDs = (Options->PPIsOnly ? 2 : NUMSGFS_MOI);
  if (Options->NeedDerivatives & NEED_DRHO) NumSGFVDs*=2;
  if (Options->NeedDerivatives & NEED_DZ  ) NumSGFVDs*=2;
  HMatrix VMatrix(NumSGFVDs, 1, V);

  return GetScalarGFs_MOI(Omega, &XMatrix, &VMatrix, Options);
}

/***************************************************************/
/* wrapper around ScalarGF routines for passage to libMDInterp */
/* constructor                                                 */
/***************************************************************/
typedef struct PhiVDFuncData
 { 
   LayeredSubstrate *S;
   cdouble Omega;
   int Dimension;
   double zFixed;
   ScalarGFOptions Options;

 } PhiVDFuncData;

void PhiVDFunc_ScalarGFs(double *RhoZ, void *UserData, double *PhiVD)
{
  PhiVDFuncData *Data      = (PhiVDFuncData *)UserData;
  LayeredSubstrate *S      = Data->S;
  cdouble Omega            = Data->Omega;
  int Dimension            = Data->Dimension;
  ScalarGFOptions *Options = &(Data->Options);
 
  double Rho = RhoZ[0];
  double z   = (Dimension>1 ? RhoZ[1] : Data->zFixed);

  cdouble V[4*NUMSGFS_MOI];
  S->GetScalarGFs_MOI(Omega, Rho, z, V, Options);

  // reorder the data in the order expected by libMDInterp, which is:
  //  PhiVD[ nPhi * NumVD + nVD ] = value/deriv #nVD for function #nPhi
  // where nPhi = 2*nSGF + 0,1 for real,imag part of scalar GF #nSGF
  int NumSGFs = (Options->PPIsOnly ? 2 : NUMSGFS_MOI);
  int NumVDs  = 2*Dimension;
  for(int nSGF=0; nSGF<NumSGFs; nSGF++)
   for(int nVD=0; nVD<NumVDs; nVD++)
    { PhiVD[ (2*nSGF + 0)*NumVDs + nVD] = real( V[nVD*NumSGFs + nSGF] );
      PhiVD[ (2*nSGF + 1)*NumVDs + nVD] = imag( V[nVD*NumSGFs + nSGF] );
    }
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
InterpND *LayeredSubstrate::InitScalarGFInterpolator(cdouble Omega,
                                                     double RhoMin, double RhoMax,
                                                     double zMin, double zMax,
                                                     bool PPIsOnly, bool Subtract,
                                                     bool RetainSingularTerms,
                                                     double DeltaRho, double Deltaz,
                                                     bool Verbose)
{
  UpdateCachedEpsMu(Omega);

  bool TrivialSubstrate = (EpsLayer[1]==1.0 && isinf(zGP));
  if (TrivialSubstrate)
   return 0;

  if (ScalarGFInterpolator)
   { delete ScalarGFInterpolator;
     ScalarGFInterpolator=0;
   }

  int Dimension = ( EqualFloat(zMin,zMax) ? 1 : 2 );

  dVec RzMin(Dimension), RzMax(Dimension);
  RzMin[0] = RhoMin;  RzMax[0] = RhoMax;
  if (Dimension>=2)
   { RzMin[1] = fmin(zMin,zMax);  RzMax[1] = fmax(zMin,zMax); }
 
  char *s=getenv("SCUFF_SUBSTRATE_DELTARHO");
  if (s) sscanf(s,"%le", &DeltaRho);
  if (Dimension==2)
   { s=getenv("SCUFF_SUBSTRATE_DELTAZ");
     if (s) sscanf(s,"%le", &Deltaz);
   }
  if (DeltaRho==0.0) DeltaRho=0.1;
  if (Deltaz==0.0) Deltaz=0.1;

  iVec NRz(Dimension);
  NRz[0]=(int)(ceil((RhoMax - RhoMin)/DeltaRho));
  if (Dimension>=2)
   NRz[1]=(int)(ceil((zMax - zMin)/Deltaz));

  int zNF = (PPIsOnly ? 2 : NUMSGFS_MOI);
  int NF  = 2*zNF;

  PhiVDFuncData Data;
  Data.S         = this;
  Data.Omega     = Omega;
  Data.Dimension = Dimension;
  Data.zFixed    = zMin;
   
  ScalarGFOptions *Options = &(Data.Options);
  InitScalarGFOptions(Options);
  Options->PPIsOnly=PPIsOnly;
  Options->Subtract=Subtract;
  Options->RetainSingularTerms=RetainSingularTerms;
  Options->NeedDerivatives = NEED_DRHO;
  if (Dimension>=2) Options->NeedDerivatives |= NEED_DZ;
 
  if (Dimension==1)
   Log("Initializing ScalarGF interpolator for Rho=[%g,%g](%i points)",RhoMin,RhoMax,NRz[0]);
  else
   Log("Initializing ScalarGF interpolator for Rho=[%g,%g](%i pts), z=[%g,%g](%i pts)",RhoMin,RhoMax,NRz[0],zMin,zMax,NRz[1]);
  ScalarGFInterpolator = new InterpND(RzMin, RzMax, NRz, NF, PhiVDFunc_ScalarGFs, (void *)&Data, Verbose);
  memcpy(&SGFIOptions, Options, sizeof(ScalarGFOptions));
  if (Dimension==1) zSGFI=zMin;
  return ScalarGFInterpolator;
}

void LayeredSubstrate::DestroyScalarGFInterpolator()
{ if (ScalarGFInterpolator)
   delete ScalarGFInterpolator;
  ScalarGFInterpolator=0;
}
