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
   RWGSurface *PSurface;
   int NumPEdges; 
   int *PPanelIndices;
   int *PPaneliQs;
   double *PLengths;
   double PPerimeter;
   double PRefPoint[3]; 

   RWGSurface *MSurface;
   int NumMEdges;
   int *MPanelIndices;
   int *MPaneliQs;
   double *MLengths;
   double MPerimeter;
   double MRefPoint[3];

 } RWGPort;

/***************************************************************/
/* this is a somewhat hacky feature implemented to allow a     */
/* calling program to request that only certain contributions  */
/* to the port voltages be included.                           */
/***************************************************************/
void SetContribOnly(const char *ContribOnly);

/***************************************************************/
/***************************************************************/
/***************************************************************/
RWGPort **ParsePortFile(RWGGeometry *G, const char *PortFileName, int *NumPorts);

void AddPortContributionsToRHS(RWGGeometry *G,
                               RWGPort **Ports, int NumPorts, cdouble *PortCurrents,
                               cdouble Omega, HVector *KN);

void AddPortContributionsToPSD(RWGGeometry *G,
                               RWGPort **Ports, int NumPorts, cdouble *PortCurrents,
                               cdouble Omega, HMatrix *PSD);

void GetPortVoltages(RWGGeometry *G, HVector *KN,
                     RWGPort **Ports, int NumPorts, cdouble *PortCurrents,
                     cdouble Omega, cdouble *PortVoltages);

void PlotPortsinGNUPLOT(const char *GPFileName, RWGPort **Ports, int NumPorts);
void PlotPortsInGMSH(RWGPort **Ports, int NumPorts, const char *format, ...);
void DrawGMSHCircle(const char *PPFile, const char *Name,
                    double *X0, double Theta, double Phi, double Radius);


void ZToS(HMatrix *ZMatrix, HMatrix *SMatrix, double ZCharacteristic);
void ZToS(HMatrix *ZMatrix, HMatrix *SMatrix);
void SToZ(HMatrix *SMatrix, HMatrix *ZMatrix, double ZCharacteristic);
void SToZ(HMatrix *SMatrix, HMatrix *ZMatrix);

#endif 
