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
 * RFSolver.h   -- header file for libRFSolver library providing
 *              -- extensions to SCUFF-EM for modeling of RF systems
 *                                          
 * homer reid   -- 3/2011 -- 3/2018
 */

#ifndef RFSOLVER_H
#define RFSOLVER_H

#include <stdio.h>
#include <math.h>
#include <stdarg.h>

#include "libscuff.h"
#include <libhmat.h>
#include <libMatProp.h>

#include <libhrutil.h>

namespace scuff {

#define _PLUS  0
#define _MINUS 1
#define NUMPOLARITIES 2

//components of potential/field vector
#define _PF_AX    0
#define _PF_AY    1
#define _PF_AZ    2
#define _PF_PHI   3
#define _PF_DXPHI 4
#define _PF_DYPHI 5
#define _PF_DZPHI 6
#define NPFC      7

/***************************************************************/
/* Data structures used to describe ports on RF devices        */
/***************************************************************/
typedef struct RWGPortEdge
 { 
   int ns;       // RWGSurface on which port edge lies
   int ne;       // (encoded) index of exterior edge
   int nPort;    // RWGPort to which this edge belongs
   int Pol;      // Polarity: \pm 1 for positive/negative port edge

   RWGPortEdge(int _ns, int _ne, int _nPort, int _Pol): ns(_ns), ne(_ne), nPort(_nPort), Pol(_Pol) {}

 } RWGPortEdge;

typedef struct RWGPort
 { vector<RWGPortEdge *> PortEdges[NUMPOLARITIES];
   double RefPoint[NUMPOLARITIES][3];
   double Perimeter[NUMPOLARITIES];
 } RWGPort;

typedef struct RWGPortList
 { vector<RWGPort *>     Ports;
   vector<RWGPortEdge *> PortEdges;
   double                RMinMax[6]; // bounding box
 } RWGPortList;

/***************************************************************/
/* Routines for handling ports in SEI scattering problems      */
/***************************************************************/
RWGPortList *ParsePortFile(RWGGeometry *G, const char *PortFileName);

HMatrix *GetPortBFInteractionMatrix(RWGGeometry *G, RWGPortList *PortList,
                                    cdouble Omega, HMatrix *PBFIMatrix=0,
                                    HMatrix *PPIMatrix=0);

void GetPortContributionToRHS(RWGGeometry *G, RWGPortList *PortList,
                              cdouble *PortCurrents, cdouble Omega, HVector *KN);

void AddPortContributionsToPSD(RWGGeometry *G,
                               RWGPortList *PortList, cdouble *PortCurrents,
                               cdouble Omega, HMatrix *PSD);

void PlotPortsInGMSH(RWGGeometry *G, RWGPortList *PortList, const char *format, ...);
void DrawGMSHCircle(const char *PPFile, const char *Name,
                    double *X0, double Theta, double Phi, double Radius);

/***************************************************************/
/* Z-to-S parameter conversion *********************************/
/***************************************************************/
void ZToS(HMatrix *ZMatrix, HMatrix *SMatrix, double ZCharacteristic);
void ZToS(HMatrix *ZMatrix, HMatrix *SMatrix);
void SToZ(HMatrix *SMatrix, HMatrix *ZMatrix, double ZCharacteristic);
void SToZ(HMatrix *SMatrix, HMatrix *ZMatrix);

/***************************************************************/
/* Routines for handling MOI (metal on insulator, i.e. thin    */
/* metal traces on the surface of a dielectric substrate,      */
/* possibly with ground plane)                                 */
/***************************************************************/
int Get1BFMOIFields(RWGGeometry *G, LayeredSubstrate *Substrate,
                    int ns, int ne, cdouble Omega, double *XDest,
                    cdouble *PFVector,
                    int Order=-1, bool Subtract=true, 
                    cdouble *SingularTerms=0);

int GetMOIMatrixElement(RWGGeometry *G, LayeredSubstrate *Substrate,
                        int nsa, int nea, int nsb, int neb,
                        cdouble Omega, cdouble *ME,
                        int Order=-1, bool Subtract=true, cdouble *Terms=0);

void AssembleMOIMatrixBlock(RWGGeometry *G, LayeredSubstrate *Substrate,
                            int nsa, int nsb, cdouble Omega, HMatrix *M);

void AssembleMOIMatrix(RWGGeometry *G, LayeredSubstrate *Substrate,
                       cdouble Omega, HMatrix *M);

void GetMOIRPFMatrices(RWGGeometry *G, LayeredSubstrate *Substrate, RWGPortList *PortList,
                       cdouble Omega, HMatrix *XMatrix,
                       HMatrix **pBFRPFMatrix, HMatrix **pPortRPFMatrix=0);

HMatrix *GetMOIFields(RWGGeometry *G, LayeredSubstrate *S, cdouble Omega,
                      HMatrix *XMatrix, HVector *KN, RWGPortList *PortList,
                      cdouble *PortCurrents, HMatrix **PFContributions=0);

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ProcessEPFile(RWGGeometry *G, HVector *KN, cdouble Omega,
                   RWGPortList *PortList, cdouble *PortCurrents,
                   char *EPFile, char *FileBase);

void ProcessFVMesh(RWGGeometry *G, HVector *KN, cdouble Omega,
                   RWGPortList *PortList, cdouble *PortCurrents,
                   char *FVMesh, char *FileBase);

void ComputeSZParms(RWGGeometry *G, RWGPortList *PortList, cdouble Omega,
                    HMatrix *M, char *FileBase, bool SParms);

} // namespace scuff 
#endif // #ifndef RFSOLVER
