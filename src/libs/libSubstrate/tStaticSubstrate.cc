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
 * tStaticSubstrate.cc -- unit test for libSubstrate
 *
 * homer reid          -- 8/2017
 *
 */
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>
#include <unistd.h>

#include "libhrutil.h"
#include "libSubstrate.h"

// freestanding silicon slab of thickness 1 length unit
const char SISubstrateFile[]=
 "0.0 CONST_EPS_11.7\n"
 "-1.0 VACUUM\n";

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
  /* name, type, #args, max_instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"SubstrateFile", PA_STRING,  1, 1, (void *)&SubstrateFile, 0, ".substrate file"},
     {"EPFile",        PA_STRING,  1, 1, (void *)&EPFile, 0, "list of evaluation points"},
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
     XMatrix = new HMatrix(NUMPTS, 3);
     for(int n=0; n<NUMPTS; n++)
      XMatrix->SetEntry(n, 2, XMIN + ((double)n)*DX);
   };

  double XS[3]={0.5, 0.0, 1.0};
  FILE *f=fopen("tStaticSubstrate.out","w");
  for(int nx=0; nx<XMatrix->NR; nx++)
   { double XD[3];
     XMatrix->GetEntriesD(nx,"0:2",XD);
     double G0Correction=0.0;
     double PhiE[4]={0.0, 0.0, 0.0, 0.0}, DeltaPhiE[4];
     S->GetDeltaPhiE(XD, XS, DeltaPhiE, &G0Correction);
     AddPhiE0(XD, XS[0], XS[1], XS[2], 1.0, PhiE);
     if (G0Correction) PhiE[0]*=G0Correction;
     VecPlusEquals(PhiE, 1.0, DeltaPhiE, 4);
     fprintVec(f,XD,3);
     fprintVecCR(f,PhiE,4);
   };
  fclose(f);
  printf("Thank you for your support.\n");

}
