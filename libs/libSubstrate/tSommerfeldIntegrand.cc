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
 * tSommerfeldIntegrand -- test the integrand of the Sommerfeld integral, as
 *                      -- computed by libSubstrate for various cases, against
 *                      -- known results
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
void RunUnitTest(int NumTest, cdouble Omega, cdouble q,
                 double Rho, double zDest, double zSource,
                 FILE *DataFile=0)
                 
{
  Log("Running unit test %i (%s)...",NumTest,TestNames[NumTest]);
  LayeredSubstrate *S = CreateLayeredSubstrate(TestSubstrates[NumTest]);
  int NumInterfaces   = S->NumInterfaces;

  if (NumTest==0)
   S->ForceFreeSpace=true;

  double XBuffer[6]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  HMatrix MyXMatrix(1,6,LHM_REAL,XBuffer);
  XBuffer[0]=Rho;
  XBuffer[2]=zDest;
  XBuffer[5]=zSource;

  SommerfeldIntegrandData MySID, *SID=&MySID;
  SID->Substrate  = S;
  SID->Omega      = Omega;
  SID->q0         = 0.0;
  SID->uTransform = false;
  SID->XMatrix    = &MyXMatrix;
  SID->RTwiddle   = new HMatrix(6, 4*NumInterfaces, LHM_COMPLEX);
  SID->WMatrix    = new HMatrix(4*NumInterfaces, 4*NumInterfaces, LHM_COMPLEX);
  SID->STwiddle   = new HMatrix(4*NumInterfaces, 6,               LHM_COMPLEX);
  SID->byqFile    = 0;

  cdouble TestIntegrand[NUMGSCALAR];
  S->HardCoded = false;
  SommerfeldIntegrand(2, (double *)&q, (void *)SID, 2*NUMGSCALAR,
                      (double *)TestIntegrand);

  cdouble ExactIntegrand[NUMGSCALAR];
  S->HardCoded = true;
  SommerfeldIntegrand(2, (double *)&q, (void *)SID, 2*NUMGSCALAR,
                      (double *)ExactIntegrand);
  if (DataFile==0)
   Compare(ExactIntegrand, TestIntegrand, NUMGSCALAR, "Exact", "Test");
  else
   { fprintf(DataFile,"%e %e %e %e %e %e %e ",real(q),imag(q),real(Omega),imag(Omega),Rho,zSource,zDest);
     fprintVec(DataFile,ExactIntegrand,NUMGSCALAR);
     fprintVecCR(DataFile,TestIntegrand,NUMGSCALAR);
   };

  delete SID->RTwiddle;
  delete SID->WMatrix;
  delete SID->STwiddle;
  delete S;
    
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
  cdouble Omega=1.0;
  cdouble q = 0.4;
  double Rho = 0.4;
  double zDest = 1.0;
  double zSource = 0.6;
  char *qFile=0;
  /* name, type, #args, max_instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"Omega",         PA_CDOUBLE, 1, 1, (void *)&Omega,      0, ""},
     {"q",             PA_CDOUBLE, 1, 1, (void *)&q,          0, ""},
     {"Rho",           PA_DOUBLE,  1, 1, (void *)&Rho,        0, ""},
     {"zDest",         PA_DOUBLE,  1, 1, (void *)&zDest,      0, ""},
     {"zSource",       PA_DOUBLE,  1, 1, (void *)&zSource,    0, ""},
     {"qFile",         PA_STRING,  1, 1, (void *)&qFile,      0, ""},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  if (qFile)
   { HVector *qList = new HVector(qFile);

     FILE *DataFile=fopen("/tmp/tSommerfeldIntegrand.out","w");
     fprintf(DataFile,"# 1-7 qr qi wr wi Rho zDest zSource\n");
     fprintf(DataFile,"# 8,9   re, im gTwiddle[_EE0P] (exact) \n");
     fprintf(DataFile,"# ...\n");
     fprintf(DataFile,"# 52,53 re, im gTwiddle[_EE0P] (HR) \n");

     for(int nq=0; nq<qList->N; nq++)
      RunUnitTest(0,Omega,qList->GetEntry(nq),Rho,zDest,zSource,DataFile);
     fclose(DataFile);
     exit(0);
   };

  RunUnitTest(0, Omega, q, Rho, zDest, zSource);

}
