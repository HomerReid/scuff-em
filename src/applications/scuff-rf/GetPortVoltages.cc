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
 * GetPortVoltages -- a scuff-EM module to compute the 'voltage' gaps between
 *                 -- the positive and negative terminals of a port
 * 
 * homer reid      -- 9/2011                                        
 *                                         
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <libhrutil.h>
#include <libhmat.h>
#include <libscuff.h>
#include <libSGJC.h>

#include "RWGPorts.h"

#define ABSTOL 1.0e-8
#define RELTOL 1.0e-4
#define FREQ2OMEGA (2.0*M_PI/300.0)
#define II cdouble(0.0,1.0)

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetPanelPotentials(RWGSurface *S, int np, int iQ, cdouble IK,
                        double *X, cdouble *PhiA);

cdouble GetPanelPotential(RWGSurface *S, int np, cdouble IK, double *X);
namespace scuff{

cdouble GetEdgePanelInteraction(double **PV, double **EV, cdouble K);

}
void GetEdgeiwA_Multipole(RWGSurface *S, int ne, cdouble IK, double X[3], cdouble iwA[3]);

/***************************************************************/
/* integrand function used to evaluate the line integral       */
/* iwA \cdot dl from X1 to X2                                  */
/* x[0] = Tau where X = X_1 + Tau(X_2-X_1)                     */
/***************************************************************/
typedef struct iwAIData 
 {
   RWGGeometry *G;
   HVector *KN;
   RWGPort **Ports;
   int NumPorts;
   cdouble *PortCurrents;
   cdouble IK;
   double *X1, *X2;
 } iwAIData;

void iwaIntegrand(unsigned ndim, const double *x, void *params, 
                  unsigned fdim, double *fval)
{
  (void) ndim;
  (void) fdim;

  iwAIData *iwAID = (iwAIData *)params;

  RWGGeometry *G        = iwAID->G;
  HVector *KN           = iwAID->KN;
  RWGPort **Ports       = iwAID->Ports;
  int NumPorts          = iwAID->NumPorts;
  cdouble *PortCurrents = iwAID->PortCurrents;
  cdouble IK            = iwAID->IK;
  double *X1            = iwAID->X1;
  double *X2            = iwAID->X2;

  double Tau=x[0];
  double X[3], X2mX1[3];
  VecSub(X2, X1, X2mX1);
  X[0] = X1[0] + Tau*X2mX1[0];
  X[1] = X1[1] + Tau*X2mX1[1];
  X[2] = X1[2] + Tau*X2mX1[2];

  /*--------------------------------------------------------------*/
  /*- contributions of interior edges ----------------------------*/
  /*--------------------------------------------------------------*/
  int BFIndex, ns, ne;
  cdouble PhiAP[4], PhiAM[4], iwA[3], iwAI;
  RWGSurface *S;
  RWGEdge *E;
  iwAI=0.0;
  for(BFIndex=0, ns=0; ns<G->NumSurfaces; ns++)
   for(S=G->Surfaces[ns], ne=0; ne<S->NumEdges; ne++, BFIndex++)
    { 
      E=S->Edges[ne];
  
      if ( VecDistance(X, E->Centroid) > 3.0*E->Radius )
       GetEdgeiwA_Multipole(S, ne, IK, X, iwA);
      else
       { GetPanelPotentials(S, E->iPPanel, E->PIndex, IK, X, PhiAP);
         GetPanelPotentials(S, E->iMPanel, E->MIndex, IK, X, PhiAM);
         iwA[0] = PhiAP[1] - PhiAM[1];
         iwA[1] = PhiAP[2] - PhiAM[2];
         iwA[2] = PhiAP[3] - PhiAM[3];
       };

      iwAI += KN->GetEntry(BFIndex) * (  (PhiAP[1]-PhiAM[1])*X2mX1[0]
                                        +(PhiAP[2]-PhiAM[2])*X2mX1[1]
                                        +(PhiAP[3]-PhiAM[3])*X2mX1[2]
                                      );
    };

  /*--------------------------------------------------------------*/
  /*- contributions of driven ports ------------------------------*/
  /*--------------------------------------------------------------*/
  int nPort, nPanel, PanelIndex, iQ;
  cdouble PortCurrent, Weight, PhiA[4];
  RWGPort *Port;
  for(nPort=0; nPort<NumPorts; nPort++)
   { 
     PortCurrent=PortCurrents[nPort];
     if (PortCurrent==0.0) 
      continue;
     Port=Ports[nPort];
     
     /*--------------------------------------------------------------*/
     /*- contribution of panels on positive edge of port ------------*/
     /*--------------------------------------------------------------*/
     S = Port->PSurface;
     Weight = PortCurrent / Port->PPerimeter;
     for(nPanel=0; nPanel<Port->NumPEdges; nPanel++)
      { 
        PanelIndex   = Port->PPanelIndices[nPanel];
        iQ           = Port->PPaneliQs[nPanel];
        
        GetPanelPotentials(S, PanelIndex, iQ, IK, X, PhiA);
        iwAI -= Weight * ( PhiA[1]*X2mX1[0] + PhiA[2]*X2mX1[1] + PhiA[3]*X2mX1[2] );
      };
     
     /*--------------------------------------------------------------*/
     /*- contribution of panels on negative edge of port ------------*/
     /*--------------------------------------------------------------*/
     S = Port->MSurface;
     Weight = PortCurrent / Port->MPerimeter;
     for(nPanel=0; nPanel<Port->NumMEdges; nPanel++)
      { 
        PanelIndex   = Port->MPanelIndices[nPanel];
        iQ           = Port->MPaneliQs[nPanel];
        
        GetPanelPotentials(S, PanelIndex, iQ, IK, X, PhiA);
        iwAI += Weight * ( PhiA[1]*X2mX1[0] + PhiA[2]*X2mX1[1] + PhiA[3]*X2mX1[2] );
      };

   };
  
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  *(cdouble *)fval=iwAI;

}

/***************************************************************/
/* compute the line integral \int iwA dl over the straight     */
/* line connecting X1 to X2                                    */
/***************************************************************/
cdouble iwAIntegral(RWGGeometry *G, HVector *KN, 
                    RWGPort **Ports, int NumPorts, cdouble *PortCurrents,
                    cdouble IK, double *X1, double *X2, double RefVal)
{
  iwAIData MyiwAIData, *iwAID=&MyiwAIData;

  iwAID->G=G;
  iwAID->KN=KN;
  iwAID->Ports=Ports;
  iwAID->NumPorts=NumPorts;
  iwAID->PortCurrents=PortCurrents;
  iwAID->IK=IK;
  iwAID->X1=X1;
  iwAID->X2=X2;

  cdouble iwAI, Error;

  double Lower=0.0;
  double Upper=1.0;
  adapt_integrate(2, iwaIntegrand, (void *)iwAID, 1, 
                  &Lower, &Upper,
		  0, RELTOL*fabs(RefVal), RELTOL,
		  (double *)&iwAI, (double *)&Error);

  return iwAI;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ComputePanelCharges(RWGGeometry *G, HVector *KN,
                         RWGPort **Ports, int NumPorts, cdouble *PortCurrents,
                         cdouble Omega, cdouble *PanelCharges)
{
  int NumPanels = G->TotalPanels;
  memset(PanelCharges, 0, NumPanels*sizeof(cdouble));

  cdouble IW = II*Omega;

  /*--------------------------------------------------------------*/
  /*- contributions of interior edges  ---------------------------*/
  /*--------------------------------------------------------------*/
  int BFIndex, ns, ne; 
  RWGSurface *S;
  RWGEdge *E;
  int PIOffset;
  for(BFIndex=0, ns=0; ns<G->NumSurfaces; ns++)
   { 
     S=G->Surfaces[ns];
     PIOffset=G->PanelIndexOffset[ns];

     for(ne=0; ne<S->NumEdges; ne++, BFIndex++)
      { 
        E=S->Edges[ne];
        PanelCharges[PIOffset + E->iPPanel] += KN->GetEntry(BFIndex) * E->Length  / IW;
        PanelCharges[PIOffset + E->iMPanel] -= KN->GetEntry(BFIndex) * E->Length  / IW;
      };

   };

  /*--------------------------------------------------------------*/
  /*- contributions of driven ports ------------------------------*/
  /*--------------------------------------------------------------*/
  int nPort, nPanel;
  RWGPort *Port;
  cdouble PortCurrent;
  for(nPort=0; nPort<NumPorts; nPort++)
   { 
     PortCurrent=PortCurrents[nPort];
     if (PortCurrent==0.0) continue;
     Port=Ports[nPort];

     S=Ports[nPort]->PSurface;
     PIOffset=G->PanelIndexOffset[S->Index];
     for(nPanel=0; nPanel<Port->NumPEdges; nPanel++)
      PanelCharges[PIOffset + Port->PPanelIndices[nPanel]] 
       -= Port->PLengths[nPanel]*PortCurrent/(IW*Port->PPerimeter);

     S=Ports[nPort]->MSurface;
     PIOffset=G->PanelIndexOffset[S->Index];
     for(nPanel=0; nPanel<Port->NumMEdges; nPanel++)
      PanelCharges[PIOffset + Port->MPanelIndices[nPanel]]
       += Port->MLengths[nPanel]*PortCurrent/(IW*Port->MPerimeter);
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetPortVoltages(RWGGeometry *G, HVector *KN,
                     RWGPort **Ports, int NumPorts, cdouble *PortCurrents,
                     cdouble Omega, cdouble *PortVoltages)
{
  cdouble IK = II*Omega;

  memset(PortVoltages,0,NumPorts*sizeof(cdouble)); 

  /***************************************************************/
  /* PanelCharges is a statically-maintained vector that stores  */
  /* the total charge on each panel in the geometry.             */
  /***************************************************************/
  static int NumPanels=0;
  static cdouble *PanelCharges=0;

  if (NumPanels!=G->TotalPanels)
   { NumPanels=G->TotalPanels;
     if (PanelCharges==0) free(PanelCharges);
     PanelCharges=(cdouble *)mallocEC(NumPanels * sizeof(cdouble));
   };

  /***************************************************************/
  /* do a first pass through all interior edges and all driven   */
  /* ports to compute the total charge on each panel.            */
  /***************************************************************/
  ComputePanelCharges(G, KN, Ports, NumPorts, PortCurrents, Omega, PanelCharges);

  /***************************************************************/
  /* ok, now that we know the charge on all panels, go through   */
  /* and compute the contribution of each panel to the potential */
  /* difference between the positive and negative sides of each  */
  /* port.                                                       */
  /***************************************************************/
  bool UseAveragePotentials=true;
  if( UseAveragePotentials )
   {
     cdouble VP, VM;
     Log("   scalar potential contribution");
     RWGPort *Port;
     RWGSurface *DestSurface, *SourceSurface;
     RWGPanel *DestPanel, *SourcePanel;
     int iDestPanel, iSourcePanel, PIOffset;
     int iQDest;
     double *EV[2];  // edge vertices 
     double *PV[3];  // panel vertices
     for(int nPort=0; nPort<NumPorts; nPort++)
      { 
        Port=Ports[nPort];

        /*--------------------------------------------------------------*/
        /*- average the scalar potential over all positive port edges   */
        /*--------------------------------------------------------------*/
        VP=0.0;
        DestSurface=Port->PSurface;
        for(int npe=0; npe<Port->NumPEdges; npe++)
         { iDestPanel = Port->PPanelIndices[npe];
           DestPanel  = DestSurface->Panels[iDestPanel];
           iQDest     = Port->PPaneliQs[npe];
           EV[0]      = DestSurface->Vertices + 3*(DestPanel->VI[ (iQDest+1)%3 ]);
           EV[1]      = DestSurface->Vertices + 3*(DestPanel->VI[ (iQDest+2)%3 ]);
           for(int ns=0; ns<G->NumSurfaces; ns++)
            for(SourceSurface=G->Surfaces[ns], PIOffset=G->PanelIndexOffset[ns], 
                iSourcePanel=0; iSourcePanel<SourceSurface->NumPanels; iSourcePanel++)
             {
               SourcePanel = SourceSurface->Panels[iSourcePanel];
               PV[0] = SourceSurface->Vertices + 3*(SourcePanel->VI[0]);
               PV[1] = SourceSurface->Vertices + 3*(SourcePanel->VI[1]);
               PV[2] = SourceSurface->Vertices + 3*(SourcePanel->VI[2]);
               VP += 2.0*ZVAC*PanelCharges[ iSourcePanel + PIOffset ] 
                      * GetEdgePanelInteraction(PV, EV, Omega);
             };
         };
        VP /= (double)(Port->NumPEdges);

        /*--------------------------------------------------------------*/
        /*- average the scalar potential over all negative port edges   */
        /*--------------------------------------------------------------*/
        VM=0.0;
        DestSurface=Port->MSurface;
        for(int npe=0; npe<Port->NumMEdges; npe++)
         { iDestPanel = Port->MPanelIndices[npe];
           DestPanel  = DestSurface->Panels[iDestPanel];
           iQDest     = Port->MPaneliQs[npe];
           EV[0]      = DestSurface->Vertices + 3*(DestPanel->VI[ (iQDest+1)%3 ]);
           EV[1]      = DestSurface->Vertices + 3*(DestPanel->VI[ (iQDest+2)%3 ]);
           for(int ns=0; ns<G->NumSurfaces; ns++)
            for(SourceSurface=G->Surfaces[ns], PIOffset=G->PanelIndexOffset[ns], 
                iSourcePanel=0; iSourcePanel<SourceSurface->NumPanels; iSourcePanel++)
             {
               SourcePanel = SourceSurface->Panels[iSourcePanel];
               PV[0] = SourceSurface->Vertices + 3*(SourcePanel->VI[0]);
               PV[1] = SourceSurface->Vertices + 3*(SourcePanel->VI[1]);
               PV[2] = SourceSurface->Vertices + 3*(SourcePanel->VI[2]);
               VM += 2.0*ZVAC*PanelCharges[ iSourcePanel + PIOffset ]
                      * GetEdgePanelInteraction(PV, EV, Omega);
             };
         };
        VM /= (double)(Port->NumMEdges);

        PortVoltages[nPort] = VP - VM;

      }; // for(int nPort=0; nPort<NumPorts; nPort++)

     // for(int nPort=0; nPort<NumPorts; nPort++)
   }; //if( UseAveragePotentials==false )

  /***************************************************************/
  /* next get the contribution to the port voltage from the line */
  /* integral of iwA                                             */
  /***************************************************************/
  cdouble iwAI;
  Log("   vector potential contribution");
  for(int nPort=0; nPort<NumPorts; nPort++)
   { iwAI=iwAIntegral(G, KN, Ports, NumPorts, PortCurrents, IK, 
                      Ports[nPort]->PRefPoint, Ports[nPort]->MRefPoint,
                      abs(PortVoltages[nPort]));
     PortVoltages[nPort] += iwAI;
   }; 

}
