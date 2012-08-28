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
                   char *EPFile)
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
  snprintf(buffer,1000,"%s.fields",GetFileBase(EPFile));
  FILE *FieldFile=CreateUniqueFile(buffer,1,FieldFileName);
  setlinebuf(FieldFile);

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
     fprintf(FieldFile,"%e %e %e",X[0],X[1],X[2]);

     GetFieldOfDrivenPorts(Ports, NumPorts, PortCurrents, Omega, X, EH);

     for(nc=0; nc<6; nc++)
      { 
        EH[nc] += FMatrix->GetEntry(np, nc);
        fprintf(FieldFile,"%.8e %.8e ",real(EH[nc]),imag(EH[nc]));
      };

     fprintf(FieldFile,"\n");

   };

  fclose(FieldFile);
  printf("Components of radiated fields written to file %s.\n",FieldFileName);

}
