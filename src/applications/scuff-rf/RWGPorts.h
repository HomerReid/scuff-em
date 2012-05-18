/*
 * RWGPorts.h   -- header file for libRWG handling of ports in 
 *              -- RF/microwave structures 
 *                                          
 * homer reid  -- 3/2011 -- 9/2011
 */

#ifndef RWGPORTS_H
#define RWGPORTS_H

#include <stdio.h>
#include <math.h>
#include <stdarg.h>

#include "libscuff.h"
#include <libhmat.h>
#include <libMatProp.h>

#include <libhrutil.h>

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct RWGPort
 { 
   RWGObject *PObject;
   int NumPEdges; 
   int *PPanelIndices;
   int *PPaneliQs;
   double *PLengths;
   double PPerimeter;
   double PRefPoint[3]; 

   RWGObject *MObject;
   int NumMEdges;
   int *MPanelIndices;
   int *MPaneliQs;
   double *MLengths;
   double MPerimeter;
   double MRefPoint[3];

 } RWGPort;

/***************************************************************/
/***************************************************************/
/***************************************************************/
RWGPort **ParsePortFile(RWGGeometry *G, const char *PortFileName, int *NumPorts);

void AddPortContributionsToRHS(RWGGeometry *G,
                               RWGPort **Ports, int NumPorts, cdouble *PortCurrents,
                               cdouble Omega, HVector *KN);

void GetPortVoltages(RWGGeometry *G, HVector *KN,
                     RWGPort **Ports, int NumPorts, cdouble *PortCurrents,
                     cdouble Omega, cdouble *PortVoltages);

void PlotPorts(const char *GPFileName, RWGPort **Ports, int NumPorts);
void DrawGMSHCircle(const char *PPFile, const char *Name,
                    double *X0, double Theta, double Phi, double Radius);


void ZToS(HMatrix *ZMatrix, HMatrix *SMatrix, double ZCharacteristic);
void ZToS(HMatrix *ZMatrix, HMatrix *SMatrix);
void SToZ(HMatrix *SMatrix, HMatrix *ZMatrix, double ZCharacteristic);
void SToZ(HMatrix *SMatrix, HMatrix *ZMatrix);

#endif 
