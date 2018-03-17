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
 * tVSI -- test of potential functions needed for single-interface substrates
 */
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>
#include <unistd.h>

#include "libhrutil.h"
#include "libSubstrate.h"
#include "libscuff.h"

const char *GroundedSiSlab=
 "0.0 CONST_EPS_11.7\n"
 "-1.0 GROUNDPLANE\n";

using namespace scuff;

#define II cdouble(0.0,1.0)

int GetVSICorrection(cdouble Omega, double Rho, double zDest,
                     cdouble Eps, double h, cdouble *V,
                     bool PPIsOnly=true, bool RetainSingularTerms=true,
                     int MaxTerms=1000, 
                     double RelTol=1.0e-6, double AbsTol=1.0e-12);

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  InstallHRSignalHandler();
  InitializeLog(argv[0]);

  /*--------------------------------------------------------------*/
  /*- process command-line arguments -----------------------------*/
  /*--------------------------------------------------------------*/
  char *SubstrateFile=0;
  cdouble Eps=12.0;
  double h=1.0;
  double XDS[6]={1.0, 0.0, 1.0, 0.0, 0.0, 0.0};
  char *XDSFile=0;
  cdouble Omega=0.7;
  bool PPIsOnly=false;
  bool CorrectionOnly=false;
  int MaxTerms=1000;
  double RelTol=1.0e-6;
  double AbsTol=1.0e-12;
  /* name, type, #args, max_instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"SubstrateFile",  PA_STRING,  1, 1, (void *)&SubstrateFile,  0, ".substrate file"},
     {"Eps",            PA_CDOUBLE, 1, 1, (void *)&Eps,            0, "substrate permittivity"},
     {"h",              PA_DOUBLE,  1, 1, (void *)&h,              0, "substrate thickness"},
//
     {"Omega",          PA_CDOUBLE, 1, 1, (void *)&Omega,          0, "angular frequency"},
//
     {"XDS",            PA_DOUBLE,  6, 1, (void *)XDS,             0, "XDS file"},
     {"XDSFile",        PA_STRING,  1, 1, (void *)&XDSFile,        0, "XDS file"},
//
//
     {"PPIsOnly",       PA_BOOL,    0, 1, (void *)&PPIsOnly,       0, ""},
     {"CorrectionOnly", PA_BOOL,    0, 1, (void *)&CorrectionOnly, 0, "retain only correction"},
     {"MaxTerms",       PA_INT,     1, 1, (void *)&MaxTerms,       0, ""},
     {"RelTol",         PA_DOUBLE,  1, 1, (void *)&RelTol,         0, ""},
     {"AbsTol",         PA_DOUBLE,  1, 1, (void *)&AbsTol,         0, ""},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  LayeredSubstrate *S; 
  if (SubstrateFile)
   { S = new LayeredSubstrate(SubstrateFile);
     S->UpdateCachedEpsMu(Omega);
   }
  else
   { char EpsStr[100], SubstrateDefinition[1000];
     if (imag(Eps)==0.0)
      snprintf(EpsStr,100,"%g",real(Eps));
     else
      snprintf(EpsStr,100,"%g+%gi",real(Eps),imag(Eps));
     if (h==0.0)
      snprintf(SubstrateDefinition,1000,"0.0 CONST_EPS_%s\n",EpsStr);
     else
      snprintf(SubstrateDefinition,1000,"0.0 CONST_EPS_%s\n%e GROUNDPLANE\n",EpsStr,-h);
     S=CreateLayeredSubstrate(SubstrateDefinition);
   };
  S->Describe();
  S->UpdateCachedEpsMu(Omega);
  Eps=S->EpsLayer[1];
  h = S->zInterface[0]-S->zGP;

  HMatrix *XMatrix = XDSFile ? new HMatrix(XDSFile) 
                             : new HMatrix(1,6,XDS);
  int NX=XMatrix->NR;
  int NumPotentials = (PPIsOnly ? 2 : 5);

  if (NX==1)
   { 
     double Rhox   = XMatrix->GetEntryD(0,0) - XMatrix->GetEntryD(0,3);
     double Rhoy   = XMatrix->GetEntryD(0,1) - XMatrix->GetEntryD(0,4);
     double Rho    = sqrt(Rhox*Rhox + Rhoy*Rhoy);
     double zDest  = XMatrix->GetEntryD(0,2);

     cdouble VFull[5], VSubt[5]; 
     int NFull = S->GetSingleInterfacePotentials(Omega, Rho, zDest, VFull, PPIsOnly, false, false);
     int NSubt = S->GetSingleInterfacePotentials(Omega, Rho, zDest, VSubt, PPIsOnly, true , false);

     printf("Full (%5i calls): %s    %s\n",NFull,CD2S(VFull[0]),CD2S(VFull[1]));
     printf("Subt (%5i calls): %s    %s\n",NSubt,CD2S(VSubt[0]),CD2S(VSubt[1]));

     cdouble VCorrection[5], VFull2[5];
     int nTerms=GetVSICorrection(Omega, Rho, zDest, Eps, h, VCorrection, 
                                 PPIsOnly, true, MaxTerms, RelTol, AbsTol);
     printf("Corr (%5i terms): %s    %s\n",nTerms,CD2S(VCorrection[0]),CD2S(VCorrection[1]));
     VecAdd(VSubt, VCorrection, VFull2, NumPotentials);

     Compare(VFull, VFull2, NumPotentials, "Full", "Subt+Corr");
     exit(1);
   }

}

#if 0
  FILE *f=fopen("tVSI.out","w");
  for(int nx=0; nx<NX; nx++)
   { double Rhox   = XMatrix->GetEntryD(nx,0) - XMatrix->GetEntryD(nx,3);
     double Rhoy   = XMatrix->GetEntryD(nx,1) - XMatrix->GetEntryD(nx,4);
     double Rho    = sqrt(Rhox*Rhox + Rhoy*Rhoy);
     double zDelta = XMatrix->GetEntryD(nx,2) - XMatrix->GetEntryD(nx,5);

     cdouble VTVector[5];
     bool PPIOnly=false;
     GetVTwiddle(q, u0, u, Rho, zDelta, Eps, h, VTVector, PPIOnly);
     fprintf(f,"%e %e ",Rho,zDelta);
     fprintVecCR(f,VTVector,5);

     if (NX==1) // test derivatives by finite-differencing
      { 
        cdouble dPhiHR[2], dPhiBF1[2], dPhiBF2[2];
        dPhiHR[0] = VTVector[_VSI_DRPHI];
        dPhiHR[1] = VTVector[_VSI_DZPHI];

        cdouble VTP[5], VTM[5];
        PPIOnly=true;

        double dRho = 1.0e-4*fabs(Rho);
        if (dRho==0.0) dRho=1.0e-4;
        GetVTwiddle(q, u0, u, Rho+dRho, zDelta, Eps, h, VTP, PPIOnly);
        GetVTwiddle(q, u0, u, Rho-dRho, zDelta, Eps, h, VTM, PPIOnly);
        dPhiBF1[0] = (VTP[_VSI_PHI] - VTVector[_VSI_PHI])/dRho;
        dPhiBF2[0] = (VTP[_VSI_PHI] - VTM[_VSI_PHI])/(2.0*dRho);

        double dz  = 1.0e-4*fabs(zDelta);
        if (dz==0.0) dz=1.0e-4;
        GetVTwiddle(q, u0, u, Rho, zDelta+dz, Eps, h, VTP, PPIOnly);
        GetVTwiddle(q, u0, u, Rho, zDelta-dz, Eps, h, VTM, PPIOnly);
        dPhiBF1[1] = (VTP[_VSI_PHI] - VTVector[_VSI_PHI])/dz;
        dPhiBF2[1] = (VTP[_VSI_PHI] - VTM[_VSI_PHI])/(2.0*dz);
        Compare(dPhiHR, dPhiBF1, dPhiBF2, 2, "HR", "BF1", "BF2");
      }
   }
  fclose(f);
#endif
