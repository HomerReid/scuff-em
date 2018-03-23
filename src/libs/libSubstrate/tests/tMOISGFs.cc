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
 * tMOISGFs -- test of scalar GFs for metal-on-insulator geometries
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
  bool OmitSingularTerms=false;
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
     {"OmitSingularTerms", PA_BOOL, 0, 1, (void *)&OmitSingularTerms, 0, "omit singular contributions"},
//
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
  int NumSGFs = (PPIsOnly ? 2 : 5);

  FILE *f=fopen("/tmp/tVSI.out","w");
  for(int nx=0; nx<NX; nx++)
   { 
     double Rhox   = XMatrix->GetEntryD(nx,0) - XMatrix->GetEntryD(nx,3);
     double Rhoy   = XMatrix->GetEntryD(nx,1) - XMatrix->GetEntryD(nx,4);
     double Rho    = sqrt(Rhox*Rhox + Rhoy*Rhoy);
     double z      = XMatrix->GetEntryD(nx,2) - XMatrix->GetEntryD(nx,5);

     bool Subtract;
     cdouble VFull[5], VSubt[5];

     Subtract=false;
     int NFull = S->GetScalarGFs_MOI(Omega, Rho, z, VFull, PPIsOnly,
                                     Subtract, !OmitSingularTerms);
     printf("Full (%5i calls): %s    %s\n",NFull,CD2S(VFull[0]),CD2S(VFull[1]));

     Subtract=true;
     int NSubt = S->GetScalarGFs_MOI(Omega, Rho, z, VSubt, PPIsOnly,
                                     Subtract, !OmitSingularTerms);
     printf("Subt (%5i calls): %s    %s\n",NSubt,CD2S(VSubt[0]),CD2S(VSubt[1]));

     if (NX==1) 
      Compare(VFull, VSubt, NumSGFs, "Full", "Subt");

     fprintf(f,"%e %e %i %i ",Rho,z,NFull,NSubt);
     fprintVec(f,VFull,NumSGFs);
     fprintVecCR(f,VSubt,NumSGFs);
   }
  fclose(f);

}
