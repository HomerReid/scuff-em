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
#include <unistd.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "libhrutil.h"
#include "libSGJC.h"
#include "libSubstrate.h"
#include "libscuff.h"
#include "libscuffInternals.h"

const char GroundedSiSlab[]=
 "0.0  CONST_EPS_12.0\n"
 "-1.0 GROUNDPLANE\n";

using namespace scuff;

#define II cdouble(0.0,1.0)

// FIXME put me somewhere else
#define NPFC 7

void GetRzMinMax(RWGGeometry *G, HMatrix *XMatrix, double RhoMinMax[2], double zMinMax[2]);

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
  SetDefaultCD2SFormat("%+.8e %+.8e");

  /*--------------------------------------------------------------*/
  /*- process command-line arguments -----------------------------*/
  /*--------------------------------------------------------------*/
  char *GeoFile=0;
  char *SubstrateFile=0;
  cdouble Eps  = -1.0;
  double     h = 0.0;
//
  cdouble Omega=1.0;
//
  double DeltaRho=0.1;
  double Deltaz=0.1;
//
  char *EPFile = 0;
  /* name, type, #args, max_instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"Geometry",      PA_STRING,  1, 1, (void *)&GeoFile,       0, ".scuffgeo file"},
     {"SubstrateFile", PA_STRING,  1, 1, (void *)&SubstrateFile, 0, ".substrate file"},
     {"Eps",           PA_CDOUBLE, 1, 1, (void *)&Eps,        0, ""},
     {"h",             PA_DOUBLE,  1, 1, (void *)&h,          0, ""},
//
     {"Omega",         PA_CDOUBLE, 1, 1, (void *)&Omega,         0, "angular frequency"},
//
     {"EPFile",        PA_STRING,  1, 1, (void *)&EPFile,     0, ""},
//
     {"DeltaRho",      PA_DOUBLE,  1, 1, (void *)&DeltaRho,   0, ""},
     {"Deltaz",        PA_DOUBLE,  1, 1, (void *)&Deltaz,     0, ""},
//
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFile==0)
   OSUsage(argv[0], VERSION, OSArray,"--geometry option is mandatory");
  if (EPFile==0)
   OSUsage(argv[0], VERSION, OSArray,"--epfile option is mandatory");

  RWGGeometry *G = new RWGGeometry(GeoFile);

  /***************************************************************/
  /* read substrate                                              */
  /***************************************************************/
  LayeredSubstrate *Substrate; 
  if (Eps!=-1.0)
   { char EpsStr[100], SubstrateDefinition[1000];
     if (imag(Eps)==0.0)
      snprintf(EpsStr,100,"%g",real(Eps));
     else
      snprintf(EpsStr,100,"%g+%gi",real(Eps),imag(Eps));
     if (h==0.0)
      snprintf(SubstrateDefinition,1000,"0.0 CONST_EPS_%s\n",EpsStr);
     else
      snprintf(SubstrateDefinition,1000,"0.0 CONST_EPS_%s\n%e GROUNDPLANE\n",EpsStr,-h);
     Substrate=CreateLayeredSubstrate(SubstrateDefinition);
   }
  else if (SubstrateFile)
   Substrate = new LayeredSubstrate(SubstrateFile);
  else
   Substrate=CreateLayeredSubstrate(GroundedSiSlab);

  Substrate->Describe();
  Substrate->UpdateCachedEpsMu(Omega);
  Eps = Substrate->EpsLayer[1];
  h   = Substrate->zInterface[0] - Substrate->zGP;
  printf("Eps={%g,%g}, h=%g\n",real(Eps),imag(Eps),h);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  HMatrix *XMatrix = new HMatrix(EPFile); 
  XMatrix->Transpose();

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double RhoMinMax[2], zMinMax[2];
  GetRzMinMax(G, XMatrix, RhoMinMax, zMinMax);
  printf("Rho=[%e,%e]  z=[%e,%e]\n",RhoMinMax[0],RhoMinMax[1],zMinMax[0],zMinMax[1]);
  bool PPIsOnly = false;
  bool Subtract = true;
  bool RetainSingularTerms = false;

  PhiVDFuncData Data;
  Data.S         = Substrate;
  Data.Omega     = Omega;
  Data.Dimension = 2;
  Data.zFixed    = zMinMax[0];
   
  ScalarGFOptions *Options = &(Data.Options);
  InitScalarGFOptions(Options);
  Options->PPIsOnly=PPIsOnly;
  Options->Subtract=Subtract;
  Options->RetainSingularTerms=RetainSingularTerms;
  Options->NeedDerivatives = NEED_DRHO;
  Options->NeedDerivatives = NEED_DZ;

  Substrate->InitScalarGFInterpolator(Omega, RhoMinMax[0], RhoMinMax[1], zMinMax[0], zMinMax[1],
                                      PPIsOnly, Subtract, RetainSingularTerms, DeltaRho, Deltaz);
  double MaxErr=
   Substrate->ScalarGFInterpolator->PlotInterpolationError(PhiVDFunc_ScalarGFs, (void *)&Data, 
                                                           "/tmp/InterpError.dat");
  printf("Max error = %e\n",MaxErr);
  
  dVec RzVec(2);
  double MeanRelError[20], MeanAbsError[20];
  FILE *f=fopen("/tmp/tInterpolation.out","w");
  for(int RhoPoints=10; RhoPoints<=50; RhoPoints+=10)
   for(int zPoints=10; zPoints<=10; zPoints+=10)
    { dVec dRzVec(2);
      dRzVec[0] = (RhoMinMax[1] - RhoMinMax[0]) / (RhoPoints-1);
      dRzVec[1] = (  zMinMax[1] -   zMinMax[0]) / (  zPoints-1);
   
      for(int nr=0; nr<=2; nr++)
       for(int nz=0; nz<=2; nz++)
        { RzVec[0] = RhoMinMax[0] + 0.5*nr*(RhoMinMax[1]-RhoMinMax[0]);
          RzVec[1] =   zMinMax[0] + 0.5*nz*(  zMinMax[1]-  zMinMax[0]);
          double MaxErr=GetInterpolationError(PhiVDFunc_ScalarGFs, (void *)&Data, 2*(PPIsOnly ? 2 : 5),
                                              RzVec, dRzVec, MeanRelError, MeanAbsError);
          printf("(%i,%i) (%i%i) %.2e \n",RhoPoints,zPoints,nr,nz,MaxErr);
          fprintf(f,"%e %e %e %e %i %i ",dRzVec[0],dRzVec[1],RzVec[0],RzVec[1],nr,nz);
          fprintVec(f,MeanRelError,20);
          fprintVecCR(f,MeanAbsError,20);
        }
    }
  fclose(f);
}
