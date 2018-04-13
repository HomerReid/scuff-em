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
 * scuff-neq   -- a standalone code within the scuff-em suite
 *             -- for implementing the fluctuating-surface-current
 *             -- approach to nonequilibrium phenomena (more 
 *             -- specifically, for computing heat radiation, 
 *             -- heat transfer, and nonequilibrium casimir forces) 
 *
 * homer reid  -- 5/2012
 */
#ifndef SCUFFNEQ_H
#define SCUFFNEQ_H

#include <libhrutil.h>
#include <libhmat.h>
#include <libscuff.h>

using namespace scuff;

// these are 'quantity flags' and 'quantity indices' used to 
// differentiate the various quantities that may be computed
// (power flux and i-directed momentum flux for i=x,y,z)

#define QFLAG_PABS      1
#define QFLAG_PRAD      2
#define QFLAG_XFORCE    4
#define QFLAG_YFORCE    8
#define QFLAG_ZFORCE   16
#define QFLAG_XTORQUE  32
#define QFLAG_YTORQUE  64
#define QFLAG_ZTORQUE 128
#define QFLAG_ALL     255

#define QINDEX_PABS    0
#define QINDEX_PRAD    1
#define QINDEX_XFORCE  2
#define QINDEX_YFORCE  3
#define QINDEX_ZFORCE  4
#define QINDEX_XTORQUE 5
#define QINDEX_YTORQUE 6
#define QINDEX_ZTORQUE 7

#define NUMSRFLUX      12

// quadrature methods 
#define QMETHOD_ADAPTIVE 0
#define QMETHOD_CLIFF    1
#define QMETHOD_TRAPSIMP 2

#define MAXPFTMETHODS 6

/****************************************************************/
/* SNEQData ('scuff-neq data') is a structure that contains all */
/* information needed to run computations a given frequency.    */
/****************************************************************/
typedef struct SNEQData
 {
   /*--------------------------------------------------------------*/
   /*- information on the geometry and geometrical transforms     -*/
   /*--------------------------------------------------------------*/
   RWGGeometry *G;
   GTCList GTCs;
   int NumTransformations;
   char *WriteCache;

   /*--------------------------------------------------------------*/
   /*- information on the calculations requested by the user       */
   /*--------------------------------------------------------------*/
   HMatrix *SRXMatrix; // eval points for spatially-resolved fluxes
   HMatrix *SRFMatrix; // storage for spatially-resolved fluxes
   int NX;             // number of spatially-resolved points
   int NumSRQs;        // number of spatially-resolved quantities

   /*--------------------------------------------------------------*/
   /*- options for computing power, force, torque -----------------*/
   /*--------------------------------------------------------------*/
   PFTOptions PFTOpts;
   HMatrix *PFTMatrix;
   HMatrix *PFTByRegion, **RegionRegionPFT;
   int NumPFTMethods;
   int PFTMethods[MAXPFTMETHODS];
   char *SIFluxFileNames[MAXPFTMETHODS];

   /*--------------------------------------------------------------*/
   /* storage for the BEM matrix and its subblocks.                */
   // Note: Buffer[0..N] are pointers into an internally-allocated */
   // chunk of memory used as a workspace in the GetFlux() routine */
   /*--------------------------------------------------------------*/
   HMatrix *M;        // BEM matrix 
   HMatrix *DRMatrix; // dressed Rytov matrix for a given source body
   HMatrix **TInt;    // TInt[ns], TExt[ns] = interior and exterior
   HMatrix **TExt;    // contributions to BEM block for surface #ns
   HMatrix **U;       // U[nb] = // off-diagonal U-matrix block #nb 

   /*--------------------------------------------------------------*/
   /*- miscellaneous other options                                -*/
   /*--------------------------------------------------------------*/
   char *FileBase;
   HVector *DSIOmegaPoints;
   bool PlotFlux;       // generate flux plots
   bool OmitSelfTerms;
   int SourceOnly, DestOnly;

 } SNEQData;

/*--------------------------------------------------------------*/
/*- in CreateSNEQData.cc ---------------------------------------*/
/*--------------------------------------------------------------*/
SNEQData *CreateSNEQData(char *GeoFile, char *TransFile,
                         int *PFTMethods, int NumPFTMethods,
                         char *EPFile, char *pFileBase);

/*--------------------------------------------------------------*/
/*- in GetFlux.cc ----------------------------------------------*/
/*--------------------------------------------------------------*/
void WriteFlux(SNEQData *SNEQD, cdouble Omega, double *kBloch=0);

#endif
