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

void AddGamma0(double XD[3], double XS[3], cdouble Omega,
               cdouble EpsRel, cdouble MuRel, cdouble *Gamma0,
               double zGP=-1.0*HUGE_VAL, bool Image=false);

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetGamma0(double XDS[6], cdouble Omega, cdouble EpsRel, cdouble MuRel,
               cdouble Gamma0[6][6])
{
  double R[3];
  VecSub(XDS+0, XDS+3, R);
  cdouble G[3][3], C[3][3], dG[3][3][3], dC[3][3][3];
  CalcGC(R, Omega, EpsRel, MuRel, G, C, dG, dC);
  cdouble k=sqrt(EpsRel*MuRel)*Omega;
  cdouble ZAbs=ZVAC*sqrt(MuRel/EpsRel);
  for(int P=0; P<2; P++)
   for(int Q=0; Q<2; Q++)
    { cdouble PreFac = II*k;
      if (P==0 && Q==0) PreFac*=ZAbs; // EE quadrant
      if (P==1 && Q==0) PreFac*=-1.0; // ME quadrant
      if (P==1 && Q==1) PreFac/=ZAbs; // MM quadrant
      for(int Mu=0; Mu<3; Mu++)
       for(int Nu=0; Nu<3; Nu++)
        Gamma0[3*P+Mu][3*Q+Nu]=PreFac*(P==Q ? G[Mu][Nu] : C[Mu][Nu]);
    };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
#if 0
bool RunUnitTest(int NumTest, cdouble Omega, double *XDS, FILE *DataFile=0)
{
  Log("Running unit test %i (%s)...",NumTest,TestNames[NumTest]);
  LayeredSubstrate *S =CreateLayeredSubstrate(TestSubstrates[NumTest]);
  if (NumTest==0)
   S->ForceFreeSpace=true;
   
  cdouble GTest[6][6];
  S->GetSubstrateDGF(Omega, XDS+0, XDS+3, GTest);

  cdouble GVacuum[6][6];
  if (NumTest==0 || NumTest==2)
   GetGVacuum(Omega, XDS, (NumTest==2), GVacuum);

  if (DataFile)
   { fprintVec(DataFile,XDS,6);
     fprintVec(DataFile,(cdouble *)GVacuum,36);
     fprintVecCR(DataFile,(cdouble *)GTest,36);
   }
  else
   Compare((cdouble *)GVacuum, (cdouble *)GTest, 36, "Exact", "Test");

  delete S;
  return true;
}
#endif

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
  int NumTest=-1;
  char *SubstrateFile=0;
  cdouble Omega=1.0;
  double XDS[6]={1.0, 0.9, 0.8, 0.7, 0.8, 0.4};
  char *XDSFile=0;
  bool Full=false;
  bool OmitFreeSpace=false;
  /* name, type, #args, max_instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"NumTest",       PA_INT,     1, 1, (void *)&NumTest,   0, "run a single unit test"},
     {"SubstrateFile", PA_STRING,  1, 1, (void *)&SubstrateFile, 0, ".substrate file"},
     {"Omega",         PA_CDOUBLE, 1, 1, (void *)&Omega,      0, "angular frequency"},
     {"XDS",           PA_DOUBLE,  6, 1, (void *)XDS,         0, ""},
     {"XDSFile",       PA_STRING,  1, 1, (void *)&XDSFile,    0, ""},
     {"Full",          PA_BOOL,    0, 1, (void *)&Full,       0, ""},
     {"OmitFreeSpace", PA_BOOL,    0, 1, (void *)&OmitFreeSpace, 0, ""},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  //if (NumTest==-1 && SubstrateFile==0
  // RunUnitTest(0, Omega, XDS);

  LayeredSubstrate *S = 0;
  if (SubstrateFile)
   S=new LayeredSubstrate(SubstrateFile);
  else
   { S=CreateLayeredSubstrate(TestSubstrates[NumTest]);
     if (NumTest==0)
      S->ForceFreeSpace=true;
   }
  S->UpdateCachedEpsMu(Omega);
  cdouble EpsRel = S->EpsLayer[1];
  cdouble MuRel  = S->MuLayer[1];
  
  HMatrix *XMatrix = XDSFile ? new HMatrix(XDSFile)
                             : new HMatrix(1,6,LHM_REAL,XDS);

  FILE *DataFile=fopen("/tmp/tFullWaveSubstrate.out","w");
  SetDefaultCD2SFormat("%e %e");
  for(int nx=0; nx<XMatrix->NR; nx++)
   { XMatrix->GetEntriesD(nx,":",XDS);

     cdouble ScriptG[6][6];
     DGFMethod Method=AUTO;
     bool AddHomogeneousDGF=true;
     S->GetSubstrateDGF(Omega, XDS+0, XDS+3, ScriptG, Method, true);
     
     HMatrix MyXMatrix(1,6);
     MyXMatrix.SetEntriesD(0,":",XDS);
     HMatrix MyGMatrix(6,6,LHM_COMPLEX);
     AddGamma0(XDS+0, XDS+3, Omega, EpsRel, MuRel, MyGMatrix.ZM, S->zGP);

     fprintVec(DataFile,XDS,6);
     for(int i=0; i<6; i++)
      for(int j=0; j<6; j++)
       fprintf(DataFile,"%s ",CD2S(ScriptG[i][j]));
     for(int i=0; i<6; i++)
      for(int j=0; j<6; j++)
        fprintf(DataFile,"%s ",CD2S(MyGMatrix.GetEntry(i,j)));
     fprintf(DataFile,"\n");

 //  HMatrix MyGMatrix(36,1,LHM_COMPLEX);
 //  S->GetSubstrateDGF(Omega, &MyXMatrix, &MyGMatrix, AUTO, AddHomogeneousDGF);

 //   cdouble Gamma0BF[6][6];
 //    GetGamma0(XDS, Omega, EpsRel, MuRel, Gamma0BF);

   }
  fclose(DataFile);
}
