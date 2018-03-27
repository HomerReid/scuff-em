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
 *
 */
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>
#include <unistd.h>

#include "libhrutil.h"
#include "libSGJC.h"
#include "libSubstrate.h"
#include "libscuff.h"
#include "libscuffInternals.h"

const char GroundedSiSlab[]=
 "0.0  CONST_EPS_12.0\n"
 "-1.0 GROUNDPLANE\n";

using namespace scuff;

#define II cdouble(0.0,1.0)

/***************************************************************/
/***************************************************************/
/***************************************************************/
int GetMOIMatrixElement(RWGGeometry *G, LayeredSubstrate *Substrate,
                        int nsa, int nea, int nsb, int neb,
                        cdouble Omega, cdouble *MEE,
                        int Order=-1, bool Subtract=false, cdouble *Terms=0);

int GetMOIMatrixElement_NC(RWGGeometry *G, LayeredSubstrate *Substrate,
                           int nsa, int nea, int nsb, int neb,
                           cdouble Omega, cdouble *MEE,
                           int InnerOrder=-1, int OuterOrder=-1, 
                           bool Subtract=false, cdouble *Terms=0);


/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  InstallHRSignalHandler();
  InitializeLog(argv[0]);
  SetDefaultCD2SFormat("%+.8e %+.8e");

  /*--------------------------------------------------------------*/
  /*- process command-line arguments -----------------------------*/
  /*--------------------------------------------------------------*/
  char *GeoFile=0;
  char *SubstrateFile=0;
  cdouble Omega=1.0;
//
  int nsa=0;
  int nea=-1;
  int nsb=0;
  int neb=-1;
//
  int Order      = 9;
  int InnerOrder = 9;
  int OuterOrder = 9;
//
  cdouble Eps  = -1.0;
  double     h = 0.0;
  /* name, type, #args, max_instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"Geometry",      PA_STRING,  1, 1, (void *)&GeoFile,       0, ".scuffgeo file"},
     {"SubstrateFile", PA_STRING,  1, 1, (void *)&SubstrateFile, 0, ".substrate file"},
     {"Omega",         PA_CDOUBLE, 1, 1, (void *)&Omega,         0, "angular frequency"},
//
     {"nsa",           PA_INT,     1, 1, (void *)&nsa,        0, ""},
     {"nea",           PA_INT,     1, 1, (void *)&nea,        0, ""},
     {"nsb",           PA_INT,     1, 1, (void *)&nsb,        0, ""},
     {"neb",           PA_INT,     1, 1, (void *)&neb,        0, ""},
//
     {"Order",         PA_INT,     1, 1, (void *)&Order,      0, ""},
     {"InnerOrder",    PA_INT,     1, 1, (void *)&InnerOrder, 0, ""},
     {"OuterOrder",    PA_INT,     1, 1, (void *)&OuterOrder, 0, ""},
//
     {"Eps",           PA_CDOUBLE, 1, 1, (void *)&Eps,        0, ""},
     {"h",             PA_DOUBLE,  1, 1, (void *)&h,          0, ""},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFile==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");

  RWGGeometry *G = new RWGGeometry(GeoFile);
  RWGSurface *Sa = G->Surfaces[nsa];
  RWGSurface *Sb = G->Surfaces[nsb];
  if (nea==-1)
   nea = irand(0, Sa->NumEdges);
  if (neb==-1)
   neb = irand(0, Sb->NumEdges);
  printf("--nea %i --neb %i\n",nea,neb);

  /***************************************************************/
  /* read substrate                                              */
  /***************************************************************/
  LayeredSubstrate *Substrate; 
  if (Eps!=-1.0)
   { char EpsStr[100], SubstrateDefinition[1000];
     if (imag(Eps)==0.0)
      snprintf(EpsStr,100,"%g",real(Eps));
     else
      snprintf(EpsStr,100,"%g+%gi",real(Eps),imag(Eps));
     if (h==0.0)
      snprintf(SubstrateDefinition,1000,"0.0 CONST_EPS_%s\n",EpsStr);
     else
      snprintf(SubstrateDefinition,1000,"0.0 CONST_EPS_%s\n%e GROUNDPLANE\n",EpsStr,-h);
     Substrate=CreateLayeredSubstrate(SubstrateDefinition);
   }
  else if (SubstrateFile)
   Substrate = new LayeredSubstrate(SubstrateFile);
  else
   Substrate=CreateLayeredSubstrate(GroundedSiSlab);

  Substrate->Describe();
  Substrate->UpdateCachedEpsMu(Omega);
  Eps = Substrate->EpsLayer[1];
  h   = Substrate->zInterface[0] - Substrate->zGP;
  printf("Eps={%g,%g}, h=%g\n",real(Eps),imag(Eps),h);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble MS[1], STerms[8], MNS[1], NSTerms[8], MNC[1], NCTerms[8];
  bool Subtract;

  Subtract=true;
  Tic();
  int SCalls=GetMOIMatrixElement(G, Substrate, nsa, nea, nsb, neb, Omega, MS, Order, Subtract, STerms);
  double STime=Toc();
  printf("    Subtracted: %i calls (%e s)\n",SCalls,STime);
  
  Subtract=false;
  Tic();
  int NSCalls=GetMOIMatrixElement(G, Substrate, nsa, nea, nsb, neb, Omega, MNS, Order, Subtract, NSTerms);
  double NSTime=Toc();
  printf("Not subtracted: %i calls (%e s)\n",NSCalls,NSTime);

  Subtract=true;
  Tic();
  memset(NCTerms, 0, 8*sizeof(cdouble));
  int NCCalls=GetMOIMatrixElement_NC(G, Substrate, nsa, nea, nsb, neb, Omega, MNC,
                                     InnerOrder, OuterOrder, Subtract, NCTerms);
  double NCTime=Toc();
  printf("Nested cubature: %i calls (%e s)\n",NCCalls,NCTime);

  printf("\n**\n** Full ME: \n**\n");
  Compare(MS, MNS, MNC, 1, "Subtracted", "Not subtracted", "Nested");
  
  printf("\n**\n** Subtracted q integral:\n**\n");
  cdouble SQTerms[2];
  SQTerms[0] = STerms[0] - STerms[2];
  SQTerms[1] = STerms[1] - STerms[3];
  printf("A:   {%e,%e} (term/full %e)\n",real(SQTerms[0]),imag(SQTerms[0]),abs(SQTerms[0]/MS[0]));
  printf("Phi: {%e,%e} (term/full %e)\n",real(SQTerms[1]),imag(SQTerms[1]),abs(SQTerms[1]/MS[0]));

  printf("\n**\n** Nonsingular correction:\n**\n");
  printf("A:   {%e,%e} (term/full %e)\n",real(STerms[2]),imag(SQTerms[2]),abs(SQTerms[2]/MS[0]));
  printf("Phi: {%e,%e} (term/full %e)\n",real(STerms[3]),imag(SQTerms[3]),abs(SQTerms[3]/MS[0]));

  printf("\n**\n** Singular correction:\n**\n");
  Compare(STerms+6, STerms+4, NCTerms+2, 2, "Exact", "Cubature", "Nested");

  printf("\n**\n** NS q Integral:\n**\n");
  Compare(STerms, NCTerms, 2, "Subtracted", "Nested");
  
  printf("Fractions of total:\n");
  printf("A   (q integral):    %e \n",abs(SQTerms[0]/MS[0]));
  printf("A   (NS correction): %e \n",abs(STerms[2]/MS[0]));
  printf("A   (singular):      %e \n",abs(STerms[6]/MS[0]));
  printf("A   (total):         %e \n",abs((STerms[0]+STerms[6])/MS[0]));
  printf("Phi (q integral):    %e \n",abs(SQTerms[1]/MS[0]));
  printf("Phi (NS correction): %e \n",abs(STerms[3]/MS[0]));
  printf("Phi (singular):      %e \n",abs(STerms[7]/MS[0]));
  printf("Phi (total):         %e \n",abs((STerms[1]+STerms[7])/MS[0]));
}
