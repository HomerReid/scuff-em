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
dVec GetRhoPoints(PhiVDFunc UserFunc, void *UserData, int NF,
                  double RhoMin, double RhoMax, double DesiredMaxRE)
{ 
  double DeltaRhoMin = (RhoMax - RhoMin) / (1000.0);
  double DeltaRhoMax = (RhoMax-RhoMin) / 2.0;
  double DeltaRho = (RhoMax - RhoMin) / (10.0);
  double Rho=RhoMin;
  dVec RhoPoints(1,Rho);
  while( Rho<RhoMax )
   { 
     dVec XVec(1,Rho);
     bool GotDeltaRho = false;
     while(!GotDeltaRho)
      { dVec dXVec(1,DeltaRho);
        double Err=GetInterpolationError(UserFunc, UserData, NF, XVec, dXVec);
        if ( Err > DesiredMaxRE )
         DeltaRho *= 0.9;
        else if (Err<0.5*DesiredMaxRE)
         DeltaRho *= 1.1;
        else 
         GotDeltaRho = true;

        if (DeltaRho<DeltaRhoMin || DeltaRho>DeltaRhoMax)
         { DeltaRho = (DeltaRho<DeltaRhoMin ? DeltaRhoMin : DeltaRhoMax);
           GotDeltaRho=true;
         }
      }
     Rho+=DeltaRho;
     if (Rho>RhoMax) Rho=RhoMax;
     RhoPoints.push_back(Rho);
   }
  return RhoPoints;
}

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
  double DeltaRho=0.1;
  double ZMin=0.0,     ZMax=0.0;
  double DeltaZ=0.1;
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
     {"DeltaRho",       PA_DOUBLE,  1, 1, (void *)&DeltaRho,       0, ""},
     {"ZMin",           PA_DOUBLE,  1, 1, (void *)&ZMin,           0, ""},
     {"ZMax",           PA_DOUBLE,  1, 1, (void *)&ZMax,           0, ""},
     {"DeltaZ",         PA_DOUBLE,  1, 1, (void *)&DeltaZ,         0, ""},
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
  int Dimension = (ZMin==ZMax) ? 1 : 2;
  int zNF = (PPIsOnly ? 2 : NUMSGFS_MOI);
  int NF  = 2*zNF;

  PhiVDFuncData Data;
  Data.S         = S;
  Data.Omega     = Omega;
  Data.Dimension = Dimension;
  Data.zFixed    = ZMin;
   
  ScalarGFOptions *Options = &(Data.Options);
  InitScalarGFOptions(Options);
  Options->PPIsOnly=PPIsOnly;
  Options->Subtract=Subtract;
  Options->RetainSingularTerms=!Subtract;
  Options->NeedDerivatives = NEED_DRHO;
  if (Dimension>=2) Options->NeedDerivatives |= NEED_DZ;

  dVec dRZVec(Dimension);
  dRZVec[0] = DeltaRho;
  if (Dimension==2) dRZVec[1] = DeltaZ;

  double *MeanRelError = new double[Dimension*NF];
  double *MeanAbsError = new double[Dimension*NF];

  dVec RhoPoints=GetRhoPoints(PhiVDFunc_ScalarGFs, (void *)&Data, NF, RhoMin, RhoMax, RelTol);
  printf("%lu points\n",RhoPoints.size());

  vector<dVec> xPoints(1,RhoPoints);
  InterpND Interp(xPoints, NF, PhiVDFunc_ScalarGFs, (void *)&Data);

  int NVD = (1<<Dimension);
  double *PhiVDExact  = new double[NVD*NF];
  double *PhiExact    = new double[NF];
  double *PhiInterp   = new double[NF];
  double *LastPhi     = new double[NF]; 
  memset(LastPhi, 0, NF*sizeof(double));
  double LastRho=0.0;
  FILE *f=fopen("tInterpError.out","w");
  for(size_t nrd=0; nrd<4*(RhoPoints.size()-1); nrd++)
   {
     int nr = nrd/4;
     int nd = nrd%4;
     double Rho0 = RhoPoints[nr];
     double DRho = 0.25*(RhoPoints[nr+1] - RhoPoints[nr]);
     double Rho = Rho0 + nd*DRho;
     dVec RZVec(1,Rho);

     PhiVDFunc_ScalarGFs(&(RZVec[0]), (void *)&Data, PhiVDExact);
     for(int nf=0; nf<NF; nf++) PhiExact[nf] = PhiVDExact[nf*NVD + 0];
     Interp.Evaluate(&Rho, PhiInterp);

     double MaxRE = 0.0;
     for(int nf=0; nf<NF; nf++) 
      MaxRE = fmax(MaxRE, RD(PhiInterp[nf], PhiExact[nf]));

     fprintVec(f,&(RZVec[0]),Dimension);
     fprintf(f,"%e ",MaxRE);
     fprintVec(f,PhiExact,NF),
     fprintVec(f,PhiInterp,NF);
     for(int nf=0; nf<NF; nf++)
      fprintf(f,"%e ",PhiVDExact[ nf*2 + 1]);
 
     if(nrd>0)
      for(int nf=0; nf<NF; nf++)
       fprintf(f,"%e ",(PhiExact[nf] - LastPhi[nf])/(Rho-LastRho));
     LastRho=Rho;
     memcpy(LastPhi, PhiExact, NF*sizeof(double));
     fprintf(f,"\n");
   }
  fclose(f);

/*

  int PointsPerDim = 3;
  int NumPoints = (Dimension==1 ? PointsPerDim : PointsPerDim*PointsPerDim);
  double RhoSpacing = (RhoMax - RhoMin) / ((double)(PointsPerDim-1));
  double ZSpacing   = (ZMax   - ZMin)   / ((double)(PointsPerDim-1));
  for(int np=0; np<NumPoints; np++)
   { 
     int nr = np%NumPoints;
     int nz = np/NumPoints;

     dVec RZVec(Dimension);
     RZVec[0] = RhoMin + nr*RhoSpacing;
     if (Dimension==2) RZVec[1] = ZMin + nz*ZSpacing;

     double MaxRE=GetInterpolationError(PhiVDFunc_ScalarGFs, (void *)&Data, NF, RZVec, dRZVec,
                                        MeanRelError, MeanAbsError);
     printf("\n");
     printf("(Rho,Z)={");   fprintVec(stdout, &(RZVec[0]), Dimension, "%g"); printf("}: ");
     printf("(dRho,dZ)={"); fprintVec(stdout, &(dRZVec[0]), Dimension,"%g"); printf("}: ");
     printf("{MaxRE}={%e}\n",MaxRE);
     printf("    nf   Mean    Max\n");
     for(int nf=0; nf<Dimension*NF; nf++)
      printf("   %i    %.2e    %.2e\n",nf,MeanRelError[nf],MeanAbsError[nf]);
   }
*/

}
