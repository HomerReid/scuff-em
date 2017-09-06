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

// infinite silicon half-space
const char SISubstrateFile[]=
 "0.0 CONST_EPS_11.7\n";

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
  double Omega=0.0;
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

  HMatrix *GMatrix=new HMatrix(XMatrix->NR, 18, LHM_COMPLEX);
  //S->GetHalfSpaceDGFs_SC(Omega, XMatrix, GMatrix);
  FILE *f=fopen("tFullWaveSubstrate.out","w");
  for(int nx=0; nx<XMatrix->NR; nx++)
   { double XXP[6];
     XMatrix->GetEntriesD(nx,"0:5",XXP);
     fprintVec(f,XXP,6);
     cdouble GVector[18];
     GMatrix->GetEntries(nx,":",GVector);
     fprintVecCR(f,GVector,18);
   };
  fclose(f);
  printf("Thank you for your support.\n");

}
