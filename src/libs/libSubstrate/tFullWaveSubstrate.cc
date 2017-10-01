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
 * tFullWaveSubstrate.cc -- unit test for libSubstrate
 *
 * homer reid       -- 8/2017
 *
 */
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>
#include <unistd.h>

#include "libhrutil.h"
#include "libSubstrate.h"
#include "libscuff.h"

// infinite silicon half-space
const char SISubstrateFile[]=
 "0.0 CONST_EPS_11.7\n";

namespace scuff {
void CalcGC(double R[3], cdouble Omega,
            cdouble EpsR, cdouble MuR,
            cdouble GMuNu[3][3], cdouble CMuNu[3][3],
            cdouble GMuNuRho[3][3][3], cdouble CMuNuRho[3][3][3]);
}

#define II cdouble(0.0,1.0)

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddScriptG0(LayeredSubstrate *S, cdouble Omega, HMatrix *XMatrix, HMatrix *GMatrix)
{
  for(int nx=0; nx<XMatrix->NR; nx++)
   { 
     double XD[6], *XS=XD+3;
     XMatrix->GetEntriesD(nx, "0:5", XD);
     int nrDest   = S->GetRegionIndex(XD[2]);
     int nrSource = S->GetRegionIndex(XS[2]);
     if (nrDest!=nrSource) continue;

     double R[3];
     VecSub(XD, XS, R);
     cdouble EpsRel, MuRel;
     S->MPLayer[nrDest]->GetEpsMu(Omega, &EpsRel, &MuRel);
     cdouble k=sqrt(EpsRel*MuRel)*Omega;
     cdouble G[3][3], C[3][3], dG[3][3][3], dC[3][3][3];
     scuff::CalcGC(R, Omega, EpsRel, MuRel, G, C, dG, dC);
     for(int Mu=0; Mu<3; Mu++)
      for(int Nu=0; Nu<3; Nu++)
       { GMatrix->AddEntry(nx, 6*(0+Mu) + 0+Nu,
                           II*Omega*MuRel*G[Mu][Nu]);
         GMatrix->AddEntry(nx, 6*(0+Mu) + 3+Nu,
                           II*k*C[Mu][Nu]);
         GMatrix->AddEntry(nx, 6*(3+Mu) + 0+Nu,
                           -1.0*II*k*C[Mu][Nu]);
         GMatrix->AddEntry(nx, 6*(3+Mu) + 3+Nu,
                           II*Omega*EpsRel*G[Mu][Nu]);
       };
   };
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
  char *EPFile=0;
  double Omega=0.1;
  bool FreeSpace=false;
  /* name, type, #args, max_instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"SubstrateFile", PA_STRING,  1, 1, (void *)&SubstrateFile, 0, ".substrate file"},
     {"EPFile",        PA_STRING,  1, 1, (void *)&EPFile,     0, "list of evaluation points"},  
     {"Omega",         PA_CDOUBLE, 1, 1, (void *)&Omega,      0, "angular frequency"},
     {"FreeSpace",     PA_BOOL,    1, 1, (void *)&FreeSpace,  0, ""},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  bool OwnsSubstrateFile=false;
  if (SubstrateFile==0)
   { SubstrateFile=strdup("XXXXXX");
     if ( mkstemp(SubstrateFile) == -1 )
      ErrExit("could not create temporary file");
     FILE *f=fopen(SubstrateFile,"w");
     fprintf(f,SISubstrateFile);
     fclose(f);
     OwnsSubstrateFile=true;
   };

  LayeredSubstrate *S=new LayeredSubstrate(SubstrateFile);
  if (S->ErrMsg)
   ErrExit(S->ErrMsg);
  if (OwnsSubstrateFile) unlink(SubstrateFile);

  if (FreeSpace) 
   { S->ForceFreeSpace=true;
     printf("Doing the free-space case foryaf.\n");
   };
  
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  HMatrix *XMatrix;
  if (EPFile)
   XMatrix = new HMatrix(EPFile);
  else
   { 
     #define NUMPTS  150
     #define XMIN   -1.5
     #define XMAX    1.5
     #define DX      (XMAX - XMIN)/(NUMPTS-1)
     XMatrix = new HMatrix(NUMPTS, 6);
     double XXP[6] = {0.0, 0.0, 0.0, 0.25, 0.5, 1.0};
     for(int nx=0; nx<NUMPTS; nx++)
      { XXP[2] = XMIN + ((double)nx)*DX;
        XMatrix->SetEntriesD(nx, "0:5", XXP);
      };
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  FILE *ff=fopen("tFullWaveSubstrate.gp","w");
  fprintf(ff,"XD(x)=(column(1))\n");
  fprintf(ff,"YD(x)=(column(2))\n");
  fprintf(ff,"ZD(x)=(column(3))\n");
  fprintf(ff,"XS(x)=(column(4))\n");
  fprintf(ff,"YS(x)=(column(5))\n");
  fprintf(ff,"ZS(x)=(column(6))\n");
  const char *EM="EM";
  const char *xyz="xyz";
  for(int pMu=0, nc=7; pMu<6; pMu++)
   for(int qNu=0; qNu<6; qNu++, nc+=2)
    { int p  = pMu/3;
      int Mu = pMu%3;
      int q  = qNu/3;
      int Nu = qNu%3;
      fprintf(ff,"rG%c%c%c%cHR(x)=($%i)\n",EM[p],EM[q],xyz[Mu],xyz[Nu],nc+0);
      fprintf(ff,"iG%c%c%c%cHR(x)=($%i)\n",EM[p],EM[q],xyz[Mu],xyz[Nu],nc+1);
      fprintf(ff,"mG%c%c%c%cHR(x)=(D2($%i,$%i))\n",EM[p],EM[q],xyz[Mu],xyz[Nu],nc+0,nc+1);
      fprintf(ff,"rG%c%c%c%cRef(x)=($%i)\n",EM[p],EM[q],xyz[Mu],xyz[Nu],nc+72);
      fprintf(ff,"iG%c%c%c%cRef(x)=($%i)\n",EM[p],EM[q],xyz[Mu],xyz[Nu],nc+73);
      fprintf(ff,"mG%c%c%c%cRef(x)=(D2($%i,$%i))\n",EM[p],EM[q],xyz[Mu],xyz[Nu],nc+72,nc+73);
    };
  fclose(ff);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  HMatrix *GSL= S->GetSubstrateDGF(Omega, XMatrix, STATIC_LIMIT);

  HMatrix *GSC= S->GetSubstrateDGF(Omega, XMatrix, SURFACE_CURRENT);
  if (!FreeSpace)
   AddScriptG0(S, Omega, XMatrix, GSC);

  FILE *f=fopen("tFullWaveSubstrate.out","w");
  for(int nx=0; nx<XMatrix->NR; nx++)
   { 
     double XXP[6];
     XMatrix->GetEntriesD(nx,"0:5",XXP);
     fprintVec(f,XXP,6);

     cdouble GVector[36];
     GSC->GetEntries(nx,":",GVector);
     fprintVec(f,GVector,36);

     if (FreeSpace)
      { cdouble GMuNu[3][3], CMuNu[3][3], GMuNuRho[3][3][3], CMuNuRho[3][3][3];
        double R[3]; 
        VecSub(XXP+0, XXP+3, R);
        cdouble EpsRel=1.0, MuRel=1.0;
        scuff::CalcGC(R, Omega, EpsRel, MuRel, GMuNu, CMuNu, GMuNuRho, CMuNuRho);
        for(int Mu=0; Mu<3; Mu++)
         for(int Nu=0; Nu<3; Nu++)
          { GVector[(0*3 + Mu)*6 + (0*3 + Nu)]=II*Omega*ZVAC*MuRel*GMuNu[Mu][Nu];
            GVector[(0*3 + Mu)*6 + (1*3 + Nu)]=+1.0*II*Omega*CMuNu[Mu][Nu];
            GVector[(1*3 + Mu)*6 + (0*3 + Nu)]=-1.0*II*Omega*CMuNu[Mu][Nu];
            GVector[(1*3 + Mu)*6 + (1*3 + Nu)]=II*Omega*ZVAC*EpsRel*GMuNu[Mu][Nu];
          };
      }
     else 
      GSL->GetEntries(nx,":",GVector);
     fprintVec(f,GVector,36);

     fprintf(f,"\n");
   };
  fclose(f);
  printf("Thank you for your support.\n");

}
