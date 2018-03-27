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
 * tInterpolation -- test of interpolation table for scalar GFs
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

double EtaFD=1.0e-4;

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
//
  cdouble Omega=0.7;
//
  double RhoMin=0.0, RhoMax=2.0;
  int RhoPoints=100;
//
  bool PPIsOnly=false;
  bool Subtract=true;
  bool CorrectionOnly=false;
  bool OmitSingularTerms=true;
//
  /* name, type, #args, max_instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"SubstrateFile",  PA_STRING,  1, 1, (void *)&SubstrateFile,  0, ".substrate file"},
     {"Eps",            PA_CDOUBLE, 1, 1, (void *)&Eps,            0, "substrate permittivity"},
     {"h",              PA_DOUBLE,  1, 1, (void *)&h,              0, "substrate thickness"},
//
     {"Omega",          PA_CDOUBLE, 1, 1, (void *)&Omega,          0, "angular frequency"},
//
     {"RhoMin",         PA_DOUBLE,  1, 1, (void *)&RhoMin,         0, "XDS file"},
     {"RhoMax",         PA_DOUBLE,  1, 1, (void *)&RhoMax,         0, "XDS file"},
     {"RhoPoints",      PA_INT,     1, 1, (void *)&RhoPoints,        0, ""},
//
     {"PPIsOnly",       PA_BOOL,    0, 1, (void *)&PPIsOnly,       0, ""},
     {"Subtract",       PA_BOOL,    0, 1, (void *)&Subtract,       0, ""},
     {"CorrectionOnly", PA_BOOL,    0, 1, (void *)&CorrectionOnly, 0, "retain only correction"},
     {"OmitSingularTerms", PA_BOOL, 0, 1, (void *)&OmitSingularTerms, 0, "omit singular contributions"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
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

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  ScalarGFOptions Options;
  InitScalarGFOptions(&Options);
  Options.PPIsOnly            = PPIsOnly;
  Options.Subtract            = Subtract;
  Options.CorrectionOnly      = CorrectionOnly;
  Options.RetainSingularTerms = !OmitSingularTerms;

  Tic();
  printf("Initializing accelerator...");
  char RPStr[10];
  snprintf(RPStr,10,"%i",RhoPoints);
  setenv("SCUFF_SUBSTRATE_RHOPOINTS",RPStr,1);
  S->InitScalarGFInterpolator(Omega, RhoMin, RhoMax, 0.0, 0.0, PPIsOnly, Subtract, !OmitSingularTerms);
  double InitTime = Toc();
  printf("...%e s (%e s / point)\n",InitTime,InitTime/RhoPoints);
   
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double DRho = (RhoMax-RhoMin)/(RhoPoints-1);
  double MaxRD=0.0, MaxRDRho=0.0;
  double MeanRD[NUMSGFS_MOI]={0.0, 0.0, 0.0, 0.0, 0.0};
  int MaxRDComponent=0, Samples=0;
  double InterpTime=0.0, FullTime=0.0;
  int NumSGFs = PPIsOnly ? 2 : NUMSGFS_MOI;
  FILE *f=fopen("/tmp/tInterpolation.out","w");
  double zDest=0.0;
  for(double Rho=RhoMin; Rho<=RhoMax; Rho+=0.25*DRho)
   { 
     Options.UseInterpolator=true;
     cdouble VInterp[NUMSGFS_MOI];
     Tic();
     S->GetScalarGFs_MOI(Omega, Rho, zDest, VInterp, &Options);
     InterpTime+=Toc();

     Options.UseInterpolator=false;
     cdouble VFull[NUMSGFS_MOI];
     Tic();
     S->GetScalarGFs_MOI(Omega, Rho, zDest, VFull, &Options);
     FullTime+=Toc();
   
     fprintf(f,"%e ",Rho);
     for(int n=0; n<NumSGFs; n++)
      { double RDn = RD(VInterp[n], VFull[n]);
        fprintf(f,"%e ",RDn);
        if (RDn > MaxRD)
         { MaxRD = RDn;
           MaxRDRho = Rho;
           MaxRDComponent=n;
         }
        MeanRD[n]+=RDn;
      }
     fprintVec(f,VInterp,NumSGFs);
     fprintVecCR(f,VFull,NumSGFs);
     Samples++;
   }
  fclose(f);

  printf("Full   time: %e ms\n",FullTime);
  printf("Interp time: %e ms\n",InterpTime);
  printf("Max RD = %e (Rho=%e, n=%i)\n",MaxRD,MaxRDRho,MaxRDComponent);
  for(int n=0; n<NumSGFs; n++)
   printf("mean RD %i: %e\n",n,MeanRD[n]/Samples);
  printf("\nThank you for your support.\n");

}

