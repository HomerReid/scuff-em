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
/***************************************************************/
/***************************************************************/
int GetVSICorrection(cdouble Omega, double Rho, double z,
                     cdouble Eps, double h, cdouble *V,
                     bool PPIsOnly, bool RetainSingularTerms,
                     int MaxTerms, double RelTol, double AbsTol)
{
  int NumPotentials = (PPIsOnly ? 2 : 5);
  memset(V, 0, NumPotentials*sizeof(cdouble));

  cdouble k      = Omega;
  cdouble Eta    = (Eps-1.0)/(Eps+1.0);
  cdouble ExpFac, ikr;
  double r, r2;

  /*--------------------------------------------------------------*/
  /* terms singular at Rho=0                                      */
  /*--------------------------------------------------------------*/
  if (RetainSingularTerms)
   { r2 = Rho*Rho + z*z;
     r=sqrt(r2);
     ikr=II*k*r;
     ExpFac = ZVAC*exp(ikr) / (4.0*M_PI*r);
     V[_VSI_APAR]   += ExpFac;
     V[_VSI_PHI ]   += (1.0-Eta)*ExpFac;
     if (!PPIsOnly)
      { ExpFac *= (ikr-1.0)/r2;
        V[_VSI_DRPHI] += Rho*(1.0-Eta)*ExpFac;
        V[_VSI_DZPHI] += z*(1.0-Eta)*ExpFac;
      }
   }

  if (isinf(h)) return 0;

  /*--------------------------------------------------------------*/
  /*- first ground-plane image term-------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble EtaFac = (1.0-Eta*Eta);
  z += 2.0*h;
  r2 = Rho*Rho + z*z;
  r=sqrt(r2);
  ikr=II*k*r;
  ExpFac = ZVAC*exp(ikr) / (4.0*M_PI*r);
  V[_VSI_APAR] -= ExpFac;
  V[_VSI_PHI]  -= EtaFac*ExpFac;
  if (!PPIsOnly)
   { ExpFac *= (ikr-1.0)/r2;
     V[_VSI_DRPHI] -= Rho*EtaFac*ExpFac;
     V[_VSI_DZPHI] -= z*EtaFac*ExpFac;
   }
 
  /*--------------------------------------------------------------*/
  /*- second and higher ground-plane image terms -----------------*/
  /*--------------------------------------------------------------*/
  double RefVals[3]; 
  RefVals[0] = abs(V[_VSI_PHI]);
  RefVals[1] = abs(V[_VSI_DRPHI]);
  RefVals[2] = abs(V[_VSI_DZPHI]);
  int ConvergedIters=0;
  for(int nTerms=1; nTerms<MaxTerms; nTerms++)
   { 
     // contribution of 2nth and (2n+1)th images
     cdouble Terms[3]={0.0, 0.0, 0.0};
     for(int n=0; n<2; n++)
      { z += 2.0*h;
        r2 = Rho*Rho + z*z;
        r=sqrt(r2);
        ikr=II*k*r;
        ExpFac = ZVAC*exp(ikr) / (4.0*M_PI*r);
        EtaFac *= -1.0*Eta;
        Terms[0] += EtaFac*ExpFac;
        if (PPIsOnly)
         { ExpFac *= (ikr-1.0)/r2;
           Terms[1] += Rho*EtaFac*ExpFac;
           Terms[2] += z*EtaFac*ExpFac;
         }
      };

     V[_VSI_PHI] -= Terms[0];
     if (!PPIsOnly)
      { V[_VSI_DRPHI] -= Terms[1];
        V[_VSI_DZPHI] -= Terms[2];
      }
     bool AllConverged=true;
     for(int p=0; p<(PPIsOnly ? 1 : 3); p++)
      { double absTerm = abs(Terms[p]);
        if ( (absTerm > AbsTol) && (absTerm>RelTol*RefVals[p]) )
         AllConverged=false;
      }
     if (AllConverged)
      ConvergedIters++;
     else
      ConvergedIters=0;
     if (ConvergedIters==3)
      return nTerms;
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
   { memset(VTVector, 0, (PPIsOnly ? 2 : 5)*sizeof(cdouble));
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
  cdouble uh = u*h, uh2=uh*uh;
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
   };
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
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct SingleInterfaceData
 { 
   LayeredSubstrate *S;
   cdouble Omega;
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
  HMatrix *XMatrix     = SIData->XMatrix;
  bool *RetainTerm     = SIData->RetainTerm;
  bool PPIsOnly        = SIData->PPIsOnly;
  FILE *byqFile        = SIData->byqFile;
  SIData->NumCalls++;

  S->UpdateCachedEpsMu(Omega);
  cdouble Eps = S->EpsLayer[1];
  double h    = S->zInterface[0] - S->zGP; // substrate thickness

  cdouble q  = x[0] + ((ndim==2) ? II*x[1] : 0.0);
  cdouble q2=q*q, u0 = sqrt(q2 - Omega*Omega), u = sqrt(q2 - Eps*Omega*Omega);

  int NumPotentials = PPIsOnly ? 2 : 5;
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

     VT0Vector[_VSI_APERP]=0.0;
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
                                                   bool PPIsOnly, bool Subtract,
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

  SingleInterfaceData SIData;
  SIData.S             = this;
  SIData.Omega         = Omega;
  SIData.XMatrix       = XMatrix;
  SIData.PPIsOnly      = PPIsOnly;
  SIData.RetainTerm[0] = !CorrectionOnly;
  SIData.RetainTerm[1] = (Subtract || CorrectionOnly);
  SIData.byqFile       = byqFile;

  int NumPotentials = PPIsOnly ? 2 : 5;
  int zfdim = NumPotentials * XMatrix->NR;
  cdouble *Error = new cdouble[zfdim]; 
  SIData.NumCalls=0;
  SommerfeldIntegrate(qIntegrand_SingleInterface, (void *)&SIData, zfdim,
                      q0, qR, xNu, Rho, MaxEvalA, MaxEvalB,
                      qAbsTol, qRelTol, VMatrix->ZM, Error);
  if (byqFile) fclose(byqFile);
  if (CorrectionOnly) VMatrix->Scale(-1.0);
  delete[] Error;
  return SIData.NumCalls;
}

int LayeredSubstrate::GetSingleInterfacePotentials(cdouble Omega, double Rho,
                                                   double zDest, cdouble *VSI, 
                                                   bool PPIsOnly, bool Subtract, 
                                                   bool CorrectionOnly)
{
  double XDS[6]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  XDS[0]=Rho;
  XDS[2]=zDest;
  HMatrix XMatrix(1,6,XDS);

  int NumPotentials = (PPIsOnly ? 2 : 5);
  HMatrix VMatrix(NumPotentials, 1, VSI);

  return GetSingleInterfacePotentials(Omega, &XMatrix, &VMatrix, PPIsOnly, Subtract,
                                      CorrectionOnly);
}
