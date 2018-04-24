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
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <libhrutil.h>
#include <libhmat.h>
#include <libscuff.h>
#include <libIncField.h>
#include <libSubstrate.h>
#include <RFSolver.h>

#define II cdouble(0.0,1.0)

using namespace scuff;

/***************************************************************/
/* main function   *********************************************/
/***************************************************************/  
int main(int argc, char *argv[])
{
  InstallHRSignalHandler();
  InitializeLog(argv[0]);

  /***************************************************************/
  /** process command-line arguments *****************************/
  /***************************************************************/
  char *GeoFile=0;
//
  char *SubstrateFile=0;
  char *EpsStr = const_cast<char *>("1.0");
  double h     = 0.0;
//
  cdouble Omega=1.0;
  int nsa=0;
  int nsb=0;

  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { {"geometry",       PA_STRING,  1, 1,       (void *)&GeoFile,    0,             "geometry file"},
     {"Omega",          PA_CDOUBLE, 1, 1,       (void *)&Omega,      0,             ""},
     {"nsa",            PA_INT,     1, 1,       (void *)&nsa,        0,             ""},
     {"nsb",            PA_INT,     1, 1,       (void *)&nsb,        0,             ""},
//
     {"SubstrateFile",  PA_STRING,  1, 1,       (void *)&SubstrateFile,   0,        "substrate definition file"},
     {"Eps",            PA_STRING,  1, 1,       (void *)&EpsStr,     0,             "substrate permittivity"},
     {"h",              PA_DOUBLE,  1, 1,       (void *)&h,          0,             "substrate thickness"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFile==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");

  /***************************************************************/
  /* create the geometry                                         */
  /***************************************************************/
  RWGGeometry *G=new RWGGeometry(GeoFile);

  char SubstrateDefinition[100];
  if (h==0.0)
   snprintf(SubstrateDefinition,100,"0.0 CONST_EPS_%s\n",EpsStr);
  else
   snprintf(SubstrateDefinition,100,"0.0 CONST_EPS_%s\n%e GROUNDPLANE\n",EpsStr,-h);
  G->Substrate = CreateLayeredSubstrate(SubstrateDefinition);

  int NEA = G->Surfaces[nsa]->NumEdges, NEB=G->Surfaces[nsb]->NumEdges;
  HMatrix *MSlow = new HMatrix(NEA, NEB, LHM_COMPLEX);
  HMatrix *MFast = new HMatrix(NEA, NEB, LHM_COMPLEX);

  Tic();
  AssembleMOIMatrixBlock(G, nsa, nsb, Omega, MSlow, 0, 0, 0);
  printf("Slow time 1: %e s\n",Toc());

  Tic();
  AssembleMOIMatrixBlock(G, nsa, nsb, Omega, MSlow, 0, 0, 0);
  printf("Slow time 2: %e s\n",Toc());

  Tic();
  EquivalentEdgePairTable *EEPTable = new EquivalentEdgePairTable(G);
  printf("EEPTable:    %e s\n",Toc());

  Tic();
  AssembleMOIMatrixBlock(G, nsa, nsb, Omega, MFast, 0, 0, EEPTable);
  printf("Fast time: %e s\n",Toc());

  double MaxRD=0.0, MeanRD=0.0;
  for(int nea=0; nea<NEA; nea++)
   for(int neb=0; neb<NEB; neb++)
    { 
      cdouble MESlow=MSlow->GetEntry(nea,neb), MEFast=MFast->GetEntry(nea,neb);
      double ThisRD = RD(MESlow, MEFast);
      MaxRD=fmax(MaxRD,ThisRD);
      MeanRD+=ThisRD;
      if(ThisRD>0.1) 
       printf("Bawonkatage! M[%i,%i]={%s} (slow), {%s} (fast)\n",nea,neb,CD2S(MESlow),CD2S(MEFast));
    }
  printf("Mean/max RD %e/%e\n",MeanRD/(NEA*NEB),MaxRD);

}
