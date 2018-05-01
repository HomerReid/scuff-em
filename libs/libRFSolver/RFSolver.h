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
 * RFSolver.h   -- header file for libRFSolver, a library that extends
 *              -- the SCUFF-EM core library for modeling of RF and
 *              -- microwave devices
 *                                          
 * homer reid   -- 3/2011 -- 4/2018
 */

#ifndef RFSOLVER_H
#define RFSOLVER_H

#include <stdio.h>
#include <math.h>
#include <stdarg.h>

#include "libscuff.h"
#include <libhmat.h>
#include <libMatProp.h>
#include "EquivalentEdgePairs.h"

#include <libhrutil.h>

#include <map>

namespace scuff {

#define FREQ2OMEGA (2.0*M_PI/300.0)
#define OMEGA2FREQ (1/(FREQ2OMEGA))

#define _PLUS  0
#define _MINUS 1
#define NUMPOLARITIES 2

// components of potential/field vector
#define _PF_AX    0
#define _PF_AY    1
#define _PF_AZ    2
#define _PF_PHI   3
#define _PF_DXPHI 4
#define _PF_DYPHI 5
#define _PF_DZPHI 6
#define NPFC      7

//
#define CONTRIBUTION_BF   1    // basis functions
#define CONTRIBUTION_PORT 2    // ports
#define CONTRIBUTION_IWA  4    // iwA term in E-field
#define CONTRIBUTION_DPHI 8    // -grad(Phi) term in E-field
#define CONTRIBUTION_ALL  0xFF

/***************************************************************/
/* Data structures used to describe ports on RF devices        */
/***************************************************************/
typedef struct RWGPortEdge
 { 
   int ns;       // RWGSurface on which port edge lies
   int ne;       // (encoded) index of exterior edge
   int nPort;    // RWGPort to which this edge belongs
   int Pol;      // Polarity: \pm 1 if this edge belongs to the positive/negative terminal of port #nport
   double Sign;  // \pm 1 if the port current flows in same/opposite direction of usual RWG current

   RWGPortEdge(int _ns, int _ne, int _nPort, int _Pol, double _Sign):
    ns(_ns), ne(_ne), nPort(_nPort), Pol(_Pol), Sign(_Sign) {}

 } RWGPortEdge;

typedef vector<RWGPortEdge *> RWGPortEdgeList;

typedef struct RWGPort
 { RWGPortEdgeList PortEdges[NUMPOLARITIES];
   double Perimeter[NUMPOLARITIES];
 } RWGPort;

typedef struct RWGPortList
 { vector<RWGPort *> Ports;
   RWGPortEdgeList   PortEdges;
   double            RMinMax[6]; // bounding box
 } RWGPortList;

/***************************************************************/
/* an "RFSolver" is a SCUFF-EM geometry, plus a list of ports, */
/* plus various internally-cached data (operating frequency,   */
/***************************************************************/
class RFSolver
 {
public:
    ////////////////////////////////////////////////////
    // API methods
    ////////////////////////////////////////////////////

    // constructors and geometry-definition routines
    RFSolver(const char *scuffgeoFileName, const char *portFileName);
    RFSolver(const char *GDSIIFileName);
    void InitSolver();

    ~RFSolver();

    // functions designed to allow python users to define
    // layered substrates
    void SetSubstratePermittivity(cdouble Epsilon);
    void SetSubstrateThickness(double h);
    void AddGroundPlane(double zGP);
    void AddSubstrateLayer(double zInterface, cdouble Epsilon, cdouble Mu=1.0);
    void SetSubstrateFile(const char *SubstrateFile);
    void InitializeSubstrate();

    void PlotGeometry(const char *PPFormat, ...);
    void PlotGeometry();

    // system assembly
    void AssembleSystemMatrix(double Freq);
    void Solve(cdouble *PortCurrents);
    void Solve(cdouble PortCurrent, int WhichPort);

    // low-level post-processing
    HMatrix *GetFields(HMatrix *XMatrix, HMatrix *PFMatrix=0);
    void GetFields(double X[3], cdouble PF[NPFC]);
    // alternative implementation of GetFields that uses RPF ("reduced potential/field") matrices
    HMatrix *GetFieldsViaRPFMatrices(HMatrix *XMatrix);
    HMatrix *GetPanelSourceDensities(HMatrix *PSDMatrix=0);

    // high-level post-processing
    void ProcessEPFile(char *EPFile, char *OutFileName=0);
    HMatrix *ProcessFVMesh(char *FVMesh, char *FVMeshTransFile, char *OutFileBase=0);
    HMatrix *GetZMatrix(HMatrix *ZMatrix=0, HMatrix **pZTerms=0);

// private:

    ////////////////////////////////////////////////////
    // internal methods
    ////////////////////////////////////////////////////
    void EvalSourceDistribution(const double X[3], cdouble iwSigmaK[4]);
    void AddMinusIdVTermsToZMatrix(HMatrix *KMatrix, HMatrix *ZMatrix);
    void AssemblePortBFInteractionMatrix();
    void UpdateSystemMatrix();

    ////////////////////////////////////////////////////
    // internal data fields
    ////////////////////////////////////////////////////

    // info on the geometry
    RWGGeometry *G;
    RWGPortList *PortList;
    EquivalentEdgePairTable *EEPTable;
    int NumPorts;
    char *FileBase;

    // internally stored variables
    cdouble Omega;          // set by most recent call to AssembleSystemMatrix
    HMatrix *M;             // LU-factorized SIE matrix, set by AssembleSystemMatrix
    bool MClean;            // false if we need to recompute M at this frequency 
    HMatrix *PBFIMatrix;    // port <--> basis-function interaction matrix
    HMatrix *PPIMatrix;     // port <--> port interaction matrix
    bool PBFIClean;         // false if PBFIMatrix/PPIMatrix need to be recomputed at this frequency
    cdouble *PortCurrents;  // set by most recent call to SolveSystem()
    HVector *KN;            // set by most recent call to SolveSystem()

    // caching of system matrix blocks when geometrical transformations are present
    bool DisableSystemBlockCache;
    cdouble OmegaCache;
    HMatrix **TBlocks; 
    HMatrix **UBlocks;

    int RetainContributions; // used to retain/exclude specific contributions to quantities for debugging

    ////////////////////////////////////////////////////
    // stuff to facilitate python-driven sessions
    ////////////////////////////////////////////////////
    std::map<double, char *> SubstrateLayers;
    char *SubstrateFile;
    bool SubstrateInitialized;

 }; // class RFSolver;

/***************************************************************/
/***************************************************************/
/***************************************************************/
RWGPortList *ParsePortFile(RWGGeometry *G, const char *PortFileName);

/***************************************************************/
/* Routines for handling MOI (metal on insulator, i.e. thin    */
/* metal traces on the surface of a dielectric substrate,      */
/* possibly with ground plane)                                 */
/***************************************************************/
int Get1BFMOIFields(RWGGeometry *G, int ns, int ne,
                    cdouble Omega, double *XDest, cdouble *PFVector,
                    int Order=-1, bool Subtract=true, 
                    cdouble *SingularTerms=0);

int GetMOIMatrixElement(RWGGeometry *G, int nsa, int nea, int nsb, int neb,
                        cdouble Omega, cdouble *ME,
                        int Order=-1, bool Subtract=true, cdouble *Terms=0);

void AssembleMOIMatrixBlock(RWGGeometry *G, int nsa, int nsb,
                            cdouble Omega, HMatrix *Block, int OffsetA=0, int OffsetB=0,
                            EquivalentEdgePairTable *EEPTable=0);

void AssembleMOIMatrix(RWGGeometry *G, cdouble Omega, HMatrix *M, EquivalentEdgePairTable *EEPTable=0);

/***************************************************************/
/* utility routines for converting Z <--> S parameters *********/
/***************************************************************/
HMatrix *Z2S(HMatrix *Z, HMatrix *S=0, double ZCharacteristic=50.0);
HMatrix *S2Z(HMatrix *S, HMatrix *Z=0, double ZCharacteristic=50.0);

} // namespace scuff 
#endif // #ifndef RFSOLVER
