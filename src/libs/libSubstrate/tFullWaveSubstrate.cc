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
     XMatrix->GetEntriesD(nx, ":", XD);
     if (XMatrix->NC==3)
      { XS[0]=XD[0]; XS[1]=XD[1]; XS[2]=XD[2]; };
     double R[3];
     VecSub(XD, XS, R);
     int nr = S->GetRegionIndex(XS[2]);
     cdouble EpsRel, MuRel;
     S->MPLayer[nr]->GetEpsMu(Omega, &EpsRel, &MuRel);
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
  /*--------------------------------------------------------------*/
  /*- process command-line arguments -----------------------------*/
  /*--------------------------------------------------------------*/
  char *SubstrateFile=0;
  char *EPFile=0;
  double Omega=0.01;
  /* name, type, #args, max_instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"SubstrateFile", PA_STRING,  1, 1, (void *)&SubstrateFile, 0, ".substrate file"},
     {"EPFile",        PA_STRING,  1, 1, (void *)&EPFile, 0, "list of evaluation points"},  
     {"Omega",         PA_CDOUBLE, 1, 1, (void *)&Omega,  0, "angular frequency"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

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

  HMatrix *GPW= S->GetSubstrateDGF(Omega, XMatrix, PLANE_WAVE);
  //HMatrix *GSC= S->GetSubstrateDGF(Omega, XMatrix, SURFACE_CURRENT);

  HMatrix *G0 = new HMatrix(XMatrix->NR, 36, LHM_COMPLEX);
  AddScriptG0(S, Omega, XMatrix, G0);

  //HMatrix *GMatrix=new HMatrix(XMatrix->NR, 36, LHM_COMPLEX);
  HMatrix *GStatic  = S->GetSubstrateDGF(Omega, XMatrix);
  FILE *f=fopen("tFullWaveSubstrate.out","w");
  for(int nx=0; nx<XMatrix->NR; nx++)
   { double XXP[6];
     XMatrix->GetEntriesD(nx,"0:5",XXP);
     fprintVec(f,XXP,6);

     cdouble GVector[36];
     GMatrix->GetEntries(nx,":",GVector);
     fprintVec(f,GVector,36);

     G0Matrix->GetEntries(nx,":",GVector);
     fprintVec(f,GVector,36);

     fprintf(f,"\n");
   };
  fclose(f);
  printf("Thank you for your support.\n");

}
