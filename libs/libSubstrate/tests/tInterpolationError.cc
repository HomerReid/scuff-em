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

typedef struct PhiVDFuncData
 { 
   LayeredSubstrate *S;
   cdouble Omega;
   int Dimension;
   double zFixed;
   ScalarGFOptions Options;

 } PhiVDFuncData;

void PhiVDFunc_ScalarGFs(double *RhoZ, void *UserData, double *PhiVD);

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
  cdouble Eps=2.25;
  double h=0.79;
//
  cdouble Omega=0.7;
  double Freq=0.0;
//
  double RhoMin=0.0, RhoMax=2.0;
  double ZMin=0.0,     ZMax=0.0;
//
  double RelTol=1.0e-3;
//
  bool PPIsOnly=true;
  bool NoSubtract=false;
//
  /* name, type, #args, max_instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"SubstrateFile",  PA_STRING,  1, 1, (void *)&SubstrateFile,  0, ".substrate file"},
     {"Eps",            PA_CDOUBLE, 1, 1, (void *)&Eps,            0, "substrate permittivity"},
     {"h",              PA_DOUBLE,  1, 1, (void *)&h,              0, "substrate thickness"},
//
     {"Omega",          PA_CDOUBLE, 1, 1, (void *)&Omega,          0, "angular frequency"},
     {"Freq",           PA_DOUBLE,  1, 1, (void *)&Freq,           0, "frequency in GHz"},
//
     {"RhoMin",         PA_DOUBLE,  1, 1, (void *)&RhoMin,         0, ""},
     {"RhoMax",         PA_DOUBLE,  1, 1, (void *)&RhoMax,         0, ""},
     {"ZMin",           PA_DOUBLE,  1, 1, (void *)&ZMin,           0, ""},
     {"ZMax",           PA_DOUBLE,  1, 1, (void *)&ZMax,           0, ""},
//
     {"RelTol",         PA_DOUBLE,  1, 1, (void *)&RelTol,         0, ""},
//
     {"PPIsOnly",       PA_BOOL,    0, 1, (void *)&PPIsOnly,       0, ""},
     {"NoSubtract",     PA_BOOL,    0, 1, (void *)&NoSubtract,       0, ""},
     {0,0,0,0,0,0,0}
   };
  bool Subtract = !NoSubtract;
  ProcessOptions(argc, argv, OSArray);
  if (Freq!=0.0)
   Omega = Freq * (2.0*M_PI/300.0);

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
  bool RetainSingularTerms = !Subtract;
  S->InitScalarGFInterpolator(Omega, RhoMin, RhoMax, ZMin, ZMax,
                              PPIsOnly, Subtract, RetainSingularTerms);
  InterpND *Interp = S->ScalarGFInterpolator;
  printf("{%lu,%lu} points\n",Interp->XGrids[0].size(),Interp->XGrids[1].size());

  PhiVDFuncData Data;
  Data.S         = S;
  Data.Omega     = Omega;
  Data.Dimension = 2; //Interp->xPoints.size();
  Data.zFixed    = S->zSGFI;
  memcpy(&(Data.Options), &(S->SGFIOptions), sizeof(ScalarGFOptions));

  double MaxError = Interp->PlotInterpolationError(PhiVDFunc_ScalarGFs, (void *)&Data,
                                                   const_cast<char *>("/tmp/tInterpolationError.out"));
/*
  printf("Rho points (%lu) = { ",Interp->xPoints[0].size());
  for(size_t n=0; n<Interp->xPoints[0].size(); n++)
   printf("%e ",Interp->xPoints[0][n]);
  printf("}\n");
  if (ZMin!=ZMax)
   { printf("Z   points (%lu) = { ",Interp->xPoints[1].size());
     for(size_t n=0; n<Interp->xPoints[1].size(); n++)
       printf("%e ",Interp->xPoints[1][n]);
   }
*/
  printf("{%lu",Interp->XGrids[0].size());
  if (Interp->XGrids.size()>1)
   printf(",%lu",Interp->XGrids[1].size());
  printf("} points: Max error = %e \n",MaxError);
  printf("Wrote error data to /tmp/tInterpolationError.out.\n");
  printf("Thank you for your support.\n");

}
