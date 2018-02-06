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
 * tFullWaveSubstrate -- libSubstrate unit test for full-wave substrate functionality
 */

#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>
#include <unistd.h>

#include "libhrutil.h"
#include "libSubstrate.h"
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
namespace scuff{
void CalcGC(double R[3], cdouble Omega,
            cdouble EpsR, cdouble MuR,
            cdouble GMuNu[3][3], cdouble CMuNu[3][3],
            cdouble GMuNuRho[3][3][3], cdouble CMuNuRho[3][3][3]);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetGExact(cdouble Omega, double XDS[6], bool GroundPlane,
               cdouble ScriptG[6][6])
{
  double R[3];
  VecSub(XDS+0, XDS+3, R);
  if (GroundPlane) R[2]=XDS[2]+XDS[5];
  cdouble G[3][3], C[3][3], dG[3][3][3], dC[3][3][3];
  CalcGC(R, Omega, 1.0, 1.0, G, C, dG, dC);
  for(int P=0; P<2; P++)
   for(int Q=0; Q<2; Q++)
    { cdouble PreFac = II*Omega;
      if (P==0 && Q==0) PreFac*=ZVAC; // EE quadrant
      if (P==1 && Q==0) PreFac*=-1.0; // ME quadrant
      if (P==1 && Q==1) PreFac/=ZVAC; // MM quadrant
      for(int Mu=0; Mu<3; Mu++)
       for(int Nu=0; Nu<3; Nu++)
        { double Sign = 1.0;
          if (GroundPlane && ((Q==0 && Nu<2) || (Q==1 && Nu==2)))
           Sign=-1.0;
          ScriptG[3*P+Mu][3*Q+Nu]=PreFac*Sign*(P==Q ? G[Mu][Nu] : C[Mu][Nu]);
        };
    };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
bool RunUnitTest(int NumTest, cdouble Omega, double *XDS, FILE *DataFile=0)
{
  Log("Running unit test %i (%s)...",NumTest,TestNames[NumTest]);
  LayeredSubstrate *S =CreateLayeredSubstrate(TestSubstrates[NumTest]);
  if (NumTest==0)
   S->ForceFreeSpace=true;
   
  cdouble GTest[6][6];
  S->GetSubstrateDGF(Omega, XDS+0, XDS+3, GTest);

  cdouble GExact[6][6];
  if (NumTest==0 || NumTest==2)
   GetGExact(Omega, XDS, (NumTest==2), GExact);

  if (DataFile)
   { fprintVec(DataFile,XDS,6);
     fprintVec(DataFile,(cdouble *)GExact,36);
     fprintVecCR(DataFile,(cdouble *)GTest,36);
   }
  else
   Compare((cdouble *)GExact, (cdouble *)GTest, 36, "Exact", "Test");

  delete S;
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
  int UnitTest=-1;
  char *SubstrateFile=0;
  cdouble Omega=1.0;
  double XDS[6]={1.0, 0.9, 0.8, 0.7, 0.8, 0.4};
  char *XDSFile=0;
  bool FreeSpace=false;
  bool OmitFreeSpace=false;
  bool HardCoded=false;
  int EntryOnly=-1;
  bool EEOnly=false; 
  bool XYOnly=false;
  /* name, type, #args, max_instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"UnitTest",      PA_INT,     1, 1, (void *)&UnitTest,   0, "run a single unit test"},
     {"SubstrateFile", PA_STRING,  1, 1, (void *)&SubstrateFile, 0, ".substrate file"},
     {"Omega",         PA_CDOUBLE, 1, 1, (void *)&Omega,      0, "angular frequency"},
     {"XDS",           PA_DOUBLE,  6, 1, (void *)XDS,         0, ""},
     {"XDSFile",       PA_STRING,  1, 1, (void *)&XDSFile,    0, ""},
     {"FreeSpace",     PA_BOOL,    1, 1, (void *)&FreeSpace,  0, ""},
     {"OmitFreeSpace", PA_BOOL,    0, 1, (void *)&OmitFreeSpace, 0, ""},
     {"HardCoded",     PA_BOOL,    1, 1, (void *)&HardCoded,  0, ""},
     {"EntryOnly",     PA_INT,     1, 1, (void *)&EntryOnly,  0, ""},
     {"EEOnly",        PA_BOOL,    0, 1, (void *)&EEOnly,     0, ""},
     {"XYOnly",        PA_BOOL,    0, 1, (void *)&XYOnly,     0, ""},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  /***************************************************************/
  /* If the user specified an --XDSFile, compute substrate DGF   */
  /* at all points and write results to data file.               */
  /***************************************************************/
  if (XDSFile)
   { HMatrix *XMatrix=new HMatrix(XDSFile);
     if (XMatrix==0) ErrExit("invalid file %s",XDSFile);
     LayeredSubstrate *S = 0;
     if (SubstrateFile)
      S=new LayeredSubstrate(SubstrateFile);
     else
      S=CreateLayeredSubstrate(TestSubstrates[0]);
     if (FreeSpace) S->ForceFreeSpace=true;
     if (HardCoded) S->HardCoded=true;
     FILE *DataFile=fopen("/tmp/tFullWaveSubstrate.out","w");
     for(int nx=0; nx<XMatrix->NR; nx++)
      { XMatrix->GetEntriesD(nx,":",XDS);
        cdouble ScriptG[6][6];
        S->GetSubstrateDGF(Omega, XDS+0, XDS+3, ScriptG, AUTO, true);
        fprintVec(DataFile,XDS,6);
        fprintVecCR(DataFile,(cdouble *)ScriptG,36);
      }
     fclose(DataFile);
     exit(0);
   }

  /***************************************************************/
  /* Otherwise run unit tests. ***********************************/
  /***************************************************************/
  RunUnitTest(0, Omega, XDS);

}
