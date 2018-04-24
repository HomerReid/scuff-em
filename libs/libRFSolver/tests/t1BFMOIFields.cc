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
 * tZMatrix.cc
 * homer reid    -- 4/2018
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
  char *PortFile=0;
//
  char *SubstrateFile=0;
  char *EpsStr = 0;
  double h     = 0.0;
//
  double Freq=1.0;
//
  int ns=0;
  int ne=-1;
//
  char *EPFile=0;
  double X[3]; int nX=0;
  bool Subtract = true;
  
  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { {"geometry",       PA_STRING,  1, 1,       (void *)&GeoFile,    0,             "geometry file"},
//
     {"portfile",       PA_STRING,  1, 1,       (void *)&PortFile,   0,             "port file"},
//
     {"SubstrateFile",  PA_STRING,  1, 1,       (void *)&SubstrateFile,   0,        "substrate definition file"},
     {"Eps",            PA_STRING,  1, 1,       (void *)&EpsStr,     0,             "substrate permittivity"},
     {"h",              PA_DOUBLE,  1, 1,       (void *)&h,          0,             "substrate thickness"},
//
     {"Freq",           PA_DOUBLE,  1, 1,       (void *)&Freq,       0,             "frequency in GHz"},
//
     {"ns",             PA_INT,     1, 1,       (void *)&ns,         0,             ""},
     {"ne",             PA_INT,     1, 1,       (void *)&ne,         0,             ""},
     {"X",              PA_DOUBLE,  3, 1,       (void *)X,          &nX,            ""},
//
     {"Subtract",       PA_BOOL,    0, 1,       (void *)&Subtract,   0,             ""},
//
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFile==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");
  if (PortFile==0)
   OSUsage(argv[0],OSArray,"--PortFileName option is mandatory");

  /***************************************************************/
  /* create the geometry                                         */
  /***************************************************************/
  RWGGeometry::UseHRWGFunctions=false;
  RWGGeometry *G=new RWGGeometry(GeoFile);
 
  HMatrix *M=G->AllocateBEMMatrix();
  HVector *KN=G->AllocateRHSVector();

  /***************************************************************/
  /* process substrate-related options                           */
  /***************************************************************/
  if (SubstrateFile)
   G->Substrate = new LayeredSubstrate(SubstrateFile);
  else if (EpsStr!=0)
   { 
     char SubstrateDefinition[1000];
     if (h==0.0) // no ground plane
      snprintf(SubstrateDefinition,1000,"0.0 CONST_EPS_%s\n",EpsStr);
     else
      snprintf(SubstrateDefinition,1000,"0.0 CONST_EPS_%s\n%e GROUNDPLANE\n",EpsStr,-h);
     G->Substrate=CreateLayeredSubstrate(SubstrateDefinition);
   }

  RWGSurface *Surf = G->Surfaces[ns];
  srandom(time(0));
  if (ne==-1)
   ne = irand(0,Surf->NumEdges);
  double Radius = Surf->Edges[ne]->Radius;
  double *X0 = Surf->Edges[ne]->Centroid;

  /***************************************************************/
  /* parse the port list and plot if requested *******************/
  /***************************************************************/
  RWGPortList *PortList=ParsePortFile(G, PortFile);
  int NumPorts = PortList->Ports.size();

  cdouble Omega=FREQ2OMEGA * Freq;
  int Orders[]={0, 4,7,9,13,16,20,25};
  double Times[8];
  cdouble PFVector[8][NPFC];

  double XPoints[5][3];
  int NumPoints;
  if (nX)
   { XPoints[0][0]=X[0];
     XPoints[0][1]=X[1];
     XPoints[0][2]=X[2];
     NumPoints=1;
   }
  else
   { NumPoints=0;
     for(double r=0.1*Radius; r<=10.1*Radius; r*=10.0)
      { double CosTheta=randU(-1.0,1.0), Phi=randU(0.0,2.0*M_PI);
        double SinTheta=sqrt(1.0-CosTheta*CosTheta);
        XPoints[NumPoints][0] = X0[0] + r*SinTheta*cos(Phi);
        XPoints[NumPoints][1] = X0[1] + r*SinTheta*sin(Phi);
        XPoints[NumPoints][2] = X0[2] + r*CosTheta;
        NumPoints++;
      }
   }
   

  SetDefaultCD2SFormat("{%+.4e,%+.6e}");
  for(int nx=0; nx<NumPoints; nx++)
   { 
     printf("\n\n X=(%g,%g,%g) r=%g: \n", XPoints[nx][0],XPoints[nx][1],XPoints[nx][2],VecDistance(X0,XPoints[nx]));
     for(int no=0; no<8; no++)
      { Tic();
        Get1BFMOIFields(G, ns, ne, Omega, XPoints[nx], PFVector[no], Orders[no], Subtract);
        Times[no]=Toc(); 
        printf("Order %2i: (%.2e s)\n",Orders[no],Times[no]);
        for(int npfc=0; npfc<NPFC; npfc++)
         { if (npfc==1 || npfc==5) continue;
           printf("%s (%.1e)  ",CD2S(PFVector[no][npfc]),RD(PFVector[no][npfc],PFVector[0][npfc]));
         };
        printf("\n");
         
      } 
        
   }

}
