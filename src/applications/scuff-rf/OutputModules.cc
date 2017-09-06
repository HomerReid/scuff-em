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
 * ProcessEPFile.cc -- scuff-rf code file for evaluating radiated fields 
 *                  -- at user-specified evaluation points 
 * 
 * homer reid       -- 9/2011
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

#include "RWGPorts.h"

#define FREQ2OMEGA (2.0*M_PI/300.0)
#define II cdouble(0.0,1.0)

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetPanelPotentials2(RWGSurface *O, int np, int iQ,
                         cdouble Omega, double *X, 
                         cdouble *A, cdouble *CurlA, cdouble *GradPhi);

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetFieldOfDrivenPorts(RWGPort **Ports, int NumPorts, 
                           cdouble *PortCurrents, 
                           cdouble Omega, double *X, cdouble *EH)
{
  int nPort;
  RWGPort *Port;
  RWGSurface *S;
  int nPanel, PanelIndex, iQ;
  cdouble PortCurrent, Weight;
  cdouble A[3], CurlA[3], GradPhi[3];
  cdouble IZK=II*ZVAC*Omega, IZoK=II*ZVAC/Omega;

  memset(EH, 0, 6*sizeof(cdouble));
  for(nPort=0; nPort<NumPorts; nPort++)
   {
     PortCurrent=PortCurrents[nPort];
     if (PortCurrent==0.0) continue;
     Port=Ports[nPort];
     
     /***************************************************************/
     /* contributions of panels on the positive side of the port    */
     /***************************************************************/
     S=Port->PSurface;
     Weight=PortCurrent/(Port->PPerimeter);
     for(nPanel=0; nPanel<Port->NumPEdges; nPanel++)
      { 
        PanelIndex  = Port->PPanelIndices[nPanel];
        iQ          = Port->PPaneliQs[nPanel];

        GetPanelPotentials2(S, PanelIndex, iQ, Omega, X, A, CurlA, GradPhi );

        EH[0] -= Weight*(IZK*A[0] - IZoK*GradPhi[0]); 
        EH[1] -= Weight*(IZK*A[1] - IZoK*GradPhi[1]); 
        EH[2] -= Weight*(IZK*A[2] - IZoK*GradPhi[2]); 
        EH[3] -= Weight*CurlA[0]; 
        EH[4] -= Weight*CurlA[1]; 
        EH[5] -= Weight*CurlA[2]; 

      }; // for nPanel=0; nPanel<Port->NumPEdges; nPanel++)
     
     /***************************************************************/
     /* contributions of panels on the negative side of the port    */
     /***************************************************************/
     S=Port->MSurface;
     Weight=PortCurrent/(Port->MPerimeter);
     for(nPanel=0; nPanel<Port->NumMEdges; nPanel++)
      { 
        PanelIndex  = Port->MPanelIndices[nPanel];
        iQ          = Port->MPaneliQs[nPanel];

        GetPanelPotentials2(S, PanelIndex, iQ, Omega, X, A, CurlA, GradPhi );

        EH[0] += Weight*(IZK*A[0] - IZoK*GradPhi[0]); 
        EH[1] += Weight*(IZK*A[1] - IZoK*GradPhi[1]); 
        EH[2] += Weight*(IZK*A[2] - IZoK*GradPhi[2]); 
        EH[3] += Weight*CurlA[0]; 
        EH[4] += Weight*CurlA[1]; 
        EH[5] += Weight*CurlA[2];

      }; // for nPanel=0; nPanel<Port->NumPEdges; nPanel++)

   }; // for(nPort=0; nPort<NumPorts; nPort++))

}
  
/***************************************************************/
/***************************************************************/
/***************************************************************/
void ProcessEPFile(RWGGeometry *G, HVector *KN, cdouble Omega,
                   RWGPort **Ports, int NumPorts, cdouble *PortCurrents,
                   char *EPFile, char *FileBase)
{
  /***************************************************************/
  /* read in the list of evaluation points  **********************/
  /***************************************************************/
  HMatrix *XMatrix=new HMatrix(EPFile,LHM_TEXT,"--ncol 3");
  if (XMatrix->ErrMsg)
   ErrExit("%s: %s",EPFile,XMatrix->ErrMsg);

  /***************************************************************/
  /* create the output file **************************************/
  /***************************************************************/
  char buffer[1000], FieldFileName[1000];
  snprintf(buffer,1000,"%s.fields",FileBase);
  FILE *FieldFile=CreateUniqueFile(buffer,1,FieldFileName);

  /***************************************************************/
  /* call the usual scuff-EM routine to get the contribution of  */
  /* the main structure to the fields                            */
  /***************************************************************/
  HMatrix *FMatrix=G->GetFields(0,KN,Omega,XMatrix);

  /***************************************************************/
  /* now go through the list of evaluation points and add the    */
  /* contributions of driven ports to the fields at each point   */
  /***************************************************************/
  double X[3];
  cdouble EH[6];
  int np, nc;
  for(np=0; np<XMatrix->NR; np++)
   { 
     X[0]=XMatrix->GetEntryD(np, 0);
     X[1]=XMatrix->GetEntryD(np, 1);
     X[2]=XMatrix->GetEntryD(np, 2);
     fprintf(FieldFile,"%e %e %e ",X[0],X[1],X[2]);

     GetFieldOfDrivenPorts(Ports, NumPorts, PortCurrents, Omega, X, EH);

     for(nc=0; nc<6; nc++)
      { 
        EH[nc] += FMatrix->GetEntry(np, nc);
        fprintf(FieldFile,"%.8e %.8e ",real(EH[nc]),imag(EH[nc]));
      };

     fprintf(FieldFile,"\n");
     fflush(FieldFile);
   };

  fclose(FieldFile);
  printf("Components of radiated fields written to file %s.\n",FieldFileName);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct RFMeshData
 { RWGGeometry *G;
   HVector *KN;
   cdouble Omega;
   RWGPort **Ports;
   int NumPorts;
   cdouble *PortCurrents;
 } RFMeshData;

HMatrix *RFFluxMDF(void *UserData, HMatrix *XMatrix, 
                   const char ***pDataNames)
{
  static const char *DataNames[]=
   { "Ex", "Ey", "Ez", "|E|", "|H|", "Radial Poynting flux" };
  int NumData = (sizeof(DataNames) / sizeof(DataNames[0]));
  *pDataNames = DataNames;

  RFMeshData *Data  = (RFMeshData *)UserData;
  RWGGeometry *G    = Data->G;
  HVector *KN       = Data->KN;
  cdouble Omega     = Data->Omega;
  RWGPort **Ports   = Data->Ports;
  int NumPorts      = Data->NumPorts;
  cdouble *PortCurrents = Data->PortCurrents;

  int NX = XMatrix->NR;
  Log("Generating FVMesh data for %i points...",NX);

  /***************************************************************/
  /* call the usual scuff-EM routine to get the contribution of  */
  /* the main structure to the fields                            */
  /***************************************************************/
  HMatrix *FMatrix=G->GetFields(0,KN,Omega,XMatrix);

  /***************************************************************/
  /* now go through the list of evaluation points and add the    */
  /* contributions of driven ports to the fields at each point,  */
  /* then fill in the appropriate colums of the DataMatrix.      */
  /***************************************************************/
  HMatrix *DataMatrix=new HMatrix(NX, NumData, LHM_REAL);
  DataMatrix->Zero();
  for(int nx=0; nx<NX; nx++)
   { 
     double X[3];
     XMatrix->GetEntriesD(nx, "0:2", X);

     cdouble EH[6], *E=EH+0, *H=EH+3;
     GetFieldOfDrivenPorts(Ports, NumPorts, PortCurrents, Omega, X, EH);
     for(int nc=0; nc<6; nc++)
      EH[nc] += FMatrix->GetEntry(nx, nc);

     double ENorm = sqrt( norm(E[0]) + norm(E[1]) + norm(E[2]) );
     double HNorm = sqrt( norm(H[0]) + norm(H[1]) + norm(H[2]) );

     double PV[3];
     PV[0] = 0.5*real( conj(E[1])*H[2] - conj(E[2])*H[1] );
     PV[1] = 0.5*real( conj(E[2])*H[0] - conj(E[0])*H[2] );
     PV[2] = 0.5*real( conj(E[0])*H[1] - conj(E[1])*H[0] );
     double Pr = (X[0]*PV[0] + X[1]*PV[1] + X[2]*PV[2])/VecNorm(X);

     DataMatrix->SetEntry(nx, 0, EH[0]);
     DataMatrix->SetEntry(nx, 1, EH[1]);
     DataMatrix->SetEntry(nx, 2, EH[2]);
     DataMatrix->SetEntry(nx, 3, ENorm);
     DataMatrix->SetEntry(nx, 4, HNorm);
     DataMatrix->SetEntry(nx, 5, Pr );
   };

  delete FMatrix;
  return DataMatrix;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
const char *PPOptions=
 //"General.ScaleX   = 3;\n"
 //"General.ScaleY   = 3;\n"
 //"General.ScaleZ   = 3;\n"
 //"General.Axes     = 0;\n"
 "View.Light   = 0;\n"
 "View.Visible = 0;\n";

void ProcessFVMesh(RWGGeometry *G, HVector *KN, cdouble Omega,
                   RWGPort **Ports, int NumPorts, cdouble *PortCurrents, 
                   char *FVMesh, char *FileBase)
{
  RFMeshData MyData, *Data = &MyData;
  Data->G            = G;
  Data->KN           = KN;
  Data->Omega        = Omega;
  Data->Ports        = Ports;
  Data->NumPorts     = NumPorts;
  Data->PortCurrents = PortCurrents;

  char OutFileName[1000];
  snprintf(OutFileName,1000,"%s.%s",FileBase,GetFileBase(FVMesh));
  MakeMeshPlot(RFFluxMDF, (void *)Data, FVMesh, PPOptions, 
               OutFileName);
}
