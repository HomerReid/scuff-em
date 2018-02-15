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
 * tDerivatives.cc -- libSubstrate unit test for derivatives of Sommerfeld integrand
 *                 -- with respect to rho, zDest, zSource
 */

#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>
#include <unistd.h>

#include "libhrutil.h"
#include "libSubstrate.h"
#include "libSpherical.h"
#include "libscuff.h"

#define NUMTESTS 5
const char *TestSubstrates[NUMTESTS]=
 { "0.0 CONST_EPS_11.7\n",
   "0.0 GROUNDPLANE\n",
   "0.0 CONST_EPS_11.7\n",
   "0.0 CONST_EPS_11.7\n -1.0 VACUUM\n",
   "0.0 CONST_EPS_11.7\n -1.0 GROUNDPLANE\n"
 };
const char *TestNames[NUMTESTS]=
 { "Free space",
   "Ground plane",
   "Si half space",
   "Si slab",
   "Si slab with ground plane",
 };

using namespace scuff;

#define II cdouble(0.0,1.0)

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetIntegrand(double rzz[3], cdouble q, SommerfeldIntegrandData *SID,
                  cdouble *Integrand)
{
  double Rho = rzz[0], zDest=rzz[1], zSource=rzz[2];

  double XBuffer[6]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  XBuffer[0]=Rho;
  XBuffer[2]=zDest;
  XBuffer[5]=zSource;
  HMatrix MyXMatrix(1,6,LHM_REAL,XBuffer);
  SID->XMatrix=&MyXMatrix;
 
  int ndim=2;
  int zfdim = NUMGSCALAR;
  if (SID->dRho) zfdim*=2;
  if (SID->dzDest) zfdim*=2;
  if (SID->dzSource) zfdim*=2;
  int dfdim = 2*zfdim;
  double *dfval = (double *)Integrand;
  SommerfeldIntegrand(ndim, (double *)&q, (void *)SID, dfdim, dfval);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
bool RunTest(const char *SubstrateFile, int NumTest,
             cdouble Omega, cdouble q,
             double Rho, double zDest, double zSource, 
             double Eta, double FirstOrder=false)
{
  LayeredSubstrate *S = 0;
 
  if (SubstrateFile) 
   { S = new LayeredSubstrate(SubstrateFile);
     Log("using substrate file %s",SubstrateFile);
   }
  else
   { S = CreateLayeredSubstrate(TestSubstrates[NumTest]);
     Log("using test substrate %i (%s)",NumTest,TestNames[NumTest]);
     if (NumTest==0)
      S->ForceFreeSpace=true;
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SommerfeldIntegrandData *SID=CreateSommerfeldIntegrandData(S);
  SID->Omega=Omega;
  SID->dRho = SID->dzDest = SID->dzSource = true;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble Integrand[8*NUMGSCALAR];
  cdouble *d001HR = Integrand + 1*NUMGSCALAR;
  cdouble *d010HR = Integrand + 2*NUMGSCALAR;
  cdouble *d011HR = Integrand + 3*NUMGSCALAR;
  cdouble *d100HR = Integrand + 4*NUMGSCALAR;
  cdouble *d101HR = Integrand + 5*NUMGSCALAR;
  cdouble *d110HR = Integrand + 6*NUMGSCALAR;
  cdouble *d111HR = Integrand + 7*NUMGSCALAR;

  double rzz[3];
  rzz[0]=Rho;
  rzz[1]=zDest;
  rzz[2]=zSource;

  GetIntegrand(rzz, q, SID, Integrand);

  cdouble dI[3][8*NUMGSCALAR];
  for(int nu=0; nu<3; nu++)
   { 
     double Delta = Eta*fabs(rzz[nu]); 
     if (Delta==0.0) Delta=Eta;
     rzz[nu] += Delta;
     GetIntegrand(rzz, q, SID, dI[nu]);
     rzz[nu] -= Delta;
     if (FirstOrder)
      { for(int nf=0; nf<8*NUMGSCALAR; nf++)
         dI[nu][nf] = (dI[nu][nf] - Integrand[nf])/(Delta);
      }
     else
      { cdouble Scratch[8*NUMGSCALAR];
        rzz[nu] -= Delta;
        GetIntegrand(rzz, q, SID, Scratch);
        rzz[nu] += Delta;
        for(int nf=0; nf<8*NUMGSCALAR; nf++)
         dI[nu][nf] = (dI[nu][nf] - Scratch[nf])/(2.0*Delta);
      };
   };

  Compare(d001HR, dI[2], NUMGSCALAR, "dzSource (HR)", "dzSource (BF)");
  Compare(d010HR, dI[1], NUMGSCALAR, "dzDest (HR)", "dzDest (BF)");
  Compare(d100HR, dI[0], NUMGSCALAR, "dRho (HR)", "dRho (BF)");

  DestroySommerfeldIntegrandData(SID);

  return true;
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
  cdouble Omega       = 1.1;
  cdouble q           = 0.4;
  double Rho          = 0.4;
  double zDest        = 1.0;
  double zSource      = 0.6;
  double Eta          = 1.0e-4;
  int NumTest         = -1;
NumTest         = 0;
  char *SubstrateFile = 0;
  bool FirstOrder     = false;
  /* name, type, #args, max_instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"Omega",         PA_CDOUBLE, 1, 1, (void *)&Omega,      0, ""},
     {"q",             PA_CDOUBLE, 1, 1, (void *)&q,          0, ""},
     {"Rho",           PA_DOUBLE,  1, 1, (void *)&Rho,        0, ""},
     {"zDest",         PA_DOUBLE,  1, 1, (void *)&zDest,      0, ""},
     {"zSource",       PA_DOUBLE,  1, 1, (void *)&zSource,    0, ""},
     {"Eta",           PA_DOUBLE,  1, 1, (void *)&Eta,        0, ""},
     {"NumTest",       PA_INT,     1, 1, (void *)&NumTest,    0, ""},
     {"SubstrateFile", PA_STRING,  1, 1, (void *)&SubstrateFile, 0, ""},
     {"FirstsOrder",   PA_BOOL,    0, 1, (void *)&FirstOrder, 0, ""},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  if (SubstrateFile) // user-specified substrate file
   { 
     RunTest(SubstrateFile, 0, Omega, q, Rho, zDest, zSource, Eta, FirstOrder);
   }
  else if (NumTest!=-1) // one particular test substrate
   { 
     RunTest(0, NumTest, Omega, q, Rho, zDest, zSource, Eta, FirstOrder);
   }
  else // unit test: run all test substrates
   {
     int TestsPassed=0;
     for(NumTest=0; NumTest<NUMTESTS; NumTest++)
      if (RunTest(0, NumTest, Omega, q, Rho, zDest, zSource, Eta))
       TestsPassed++;
     return (TestsPassed==NUMTESTS) ? 0 : 1;
   }

}
