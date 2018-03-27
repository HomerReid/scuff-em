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

// FIXME put me somewhere else
#define NPFC 7

/***************************************************************/
/***************************************************************/
/***************************************************************/
int Get1BFMOIFields(RWGGeometry *G, LayeredSubstrate *Substrate,
                    int ns, int ne, cdouble Omega, double *XDest,
                    cdouble *PFVector, int Order=-1, bool Subtract=true, 
                    cdouble *SingularTerms=0);

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
  int ns=0;
  int ne=-1;
//
  double X[3]  = {1.0, 0.0, 1.0};
  char *EPFile = 0;
  int Order    = 9;
  cdouble Eps  = -1.0;
  double     h = 0.0;
  /* name, type, #args, max_instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"Geometry",      PA_STRING,  1, 1, (void *)&GeoFile,       0, ".scuffgeo file"},
     {"SubstrateFile", PA_STRING,  1, 1, (void *)&SubstrateFile, 0, ".substrate file"},
     {"Omega",         PA_CDOUBLE, 1, 1, (void *)&Omega,         0, "angular frequency"},
//
     {"ns",            PA_INT,     1, 1, (void *)&ns,         0, ""},
     {"ne",            PA_INT,     1, 1, (void *)&ne,         0, ""},
//
     {"X",             PA_DOUBLE,  3, 1, (void *)X,           0, ""},
     {"EPFile",        PA_STRING,  1, 1, (void *)&EPFile,     0, ""},
//
     {"Order",         PA_INT,     1, 1, (void *)&Order,      0, ""},
//
     {"Eps",           PA_CDOUBLE, 1, 1, (void *)&Eps,        0, ""},
     {"h",             PA_DOUBLE,  1, 1, (void *)&h,          0, ""},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFile==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");

  RWGGeometry *G = new RWGGeometry(GeoFile);
  RWGSurface *S=G->Surfaces[ns];
  if (ne==-1)
   ne = irand(0, S->NumEdges);

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
  HMatrix *XMatrix = EPFile ? new HMatrix(EPFile) : new HMatrix(1,3);
  if (!EPFile) XMatrix->SetEntriesD(0,":",X);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  FILE *f=fopen("/tmp/tMOIFields.out","w");
  fprintf(f,"#1-3 x y z\n");
  fprintf(f,"#4-9   iwA_{x,y,z}  (unsubtracted)\n");
  fprintf(f,"#10-17 Phi, GradPhi (unsubtracted)\n");
  fprintf(f,"#18-23 iwA_{x,y,z}  (subtracted)\n");
  fprintf(f,"#24-31 Phi, GradPhi (subtracted)\n");
  fprintf(f,"#32-37 iwA_{x,y,z}  (sub-corr)\n");
  fprintf(f,"#38-45 Phi, GradPhi (sub-corr)\n");
  fprintf(f,"#46-51 iwA_{x,y,z}  (corr)\n");
  fprintf(f,"#52-59 Phi, GradPhi (corr)\n");
  for(int nx=0; nx<XMatrix->NR; nx++)
   { double XDest[3];
     XMatrix->GetEntriesD(nx,"0:2",XDest);
     cdouble PFVectorNS[NPFC], PFVectorS[NPFC];
     bool Subtract;

     Subtract=false;
     Get1BFMOIFields(G, Substrate, ns, ne, Omega, XDest, PFVectorNS, Order, Subtract);

     Subtract=true;
     cdouble SingularTerms[NPFC];
     Get1BFMOIFields(G, Substrate, ns, ne, Omega, XDest, PFVectorS, Order, Subtract, SingularTerms);

     const char *PFNames[NPFC]={"Ax","Ay","Az","Phi","dxPhi","dyPhi","dzPhi"};
     if (nx==0)
      { 
        Compare(PFVectorS, PFVectorNS, NPFC, "Subtracted", "Not subtracted");

        cdouble NSTerms[NPFC];
        VecSub(PFVectorS, SingularTerms, NSTerms, NPFC);
        printf("%5s | %-23s | %-23s | qInt/Total\n","Cmpt","Subtracted q Integral","Correction");
        for(int n=0; n<NPFC; n++)
         printf("%5s | (%+.3e,%+.3e) | (%+.3e,%+.3e) | %e\n",PFNames[n],
                 real(NSTerms[n]),imag(NSTerms[n]),real(SingularTerms[n]),imag(SingularTerms[n]),
                 abs(NSTerms[n]/PFVectorS[n]));
      }

     fprintVec(f,   XDest);
     fprintVec(f,   PFVectorNS, NPFC);
     fprintVec(f,   PFVectorS,  NPFC);
     fprintVec(f,   SingularTerms, NPFC);
   }
  fclose(f);

}
