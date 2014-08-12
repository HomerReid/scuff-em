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

#define QFLAG_POWER    1
#define QFLAG_XFORCE   2
#define QFLAG_YFORCE   4
#define QFLAG_ZFORCE   8
#define QFLAG_XTORQUE 16
#define QFLAG_YTORQUE 32
#define QFLAG_ZTORQUE 64

#define QINDEX_POWER   0
#define QINDEX_XFORCE  1
#define QINDEX_YFORCE  2
#define QINDEX_ZFORCE  3
#define QINDEX_XTORQUE 4
#define QINDEX_YTORQUE 5
#define QINDEX_ZTORQUE 6

#define MAXQUANTITIES 7

#define NUMSIPFT MAXQUANTITIES

// quadrature methods 
#define QMETHOD_ADAPTIVE 0
#define QMETHOD_CLIFF    1
#define QMETHOD_TRAPSIMP 2

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
   GTComplex **GTCList;
   int NumTransformations;
   char *WriteCache;

   /*--------------------------------------------------------------*/
   /*- information on the calculations requested by the user       */
   /*--------------------------------------------------------------*/
   int QuantityFlags;
   int NQ;
   int NSNQ;
   int NTNSNQ;
   int PlotFlux;
   bool SymGPower;

   /*--------------------------------------------------------------*/
   /* storage for the BEM matrix and its subblocks                 */
   /*--------------------------------------------------------------*/
   HMatrix *W;        // BEM matrix 
   HMatrix **T;       // T[ns] = T-matrix block for surface #ns
   HMatrix **TSelf;   //
   HMatrix **U;       // U[ns*NS + nsp] = // U-matrix block for surfaces #ns, #nsp

   // Buffer[0..N] are pointers into an internally-allocated
   // chunk of memory used as a workspace in the GetTrace() routine.
   void *Buffer[MAXQUANTITIES+1];

   /*--------------------------------------------------------------*/
   /* storage for sparse PFT matrices                              */
   /*--------------------------------------------------------------*/
   // NeedMatrix[nm] = 1 if we need overlap matrix of type #nm
   // SArray[ ns*8 + nm ] = overlap matrix of type #nm for
   // surface #ns (note 8 = SCUFF_NUM_OMATRICES)
   bool NeedMatrix[SCUFF_NUM_OMATRICES]; 
   SMatrix ***SArray;

   // radius and number of quadrature points per dimension
   // for surface-integral PFT
   double SIRadius;
   int SINumPoints;

   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   bool *OmegaConverged;

   /*--------------------------------------------------------------*/
   /*- miscellaneous other options                                -*/
   /*--------------------------------------------------------------*/
   double RelTol, AbsTol;  // integration tolerances
   char *FileBase;
   bool UseExistingData;

 } SNEQData;

/*--------------------------------------------------------------*/
/*- in CreateSNEQData.cc ---------------------------------------*/
/*--------------------------------------------------------------*/
SNEQData *CreateSNEQData(char *GeoFile, char *TransFile,
                         int QuantityFlags, char *FileBase);

/*--------------------------------------------------------------*/
/*- in GetFlux.cc ----------------------------------------------*/
/*--------------------------------------------------------------*/
int GetIndex(SNEQData *SNEQD, int nt, int nss, int nsd, int nq);
void GetFlux(SNEQData *SNEQD, cdouble Omega, double *kBloch, double *Flux);
void GetFlux(SNEQData *SNEQD, cdouble Omega, double *Flux);

/*--------------------------------------------------------------*/
/*- in Quadrature.cc          ----------------------------------*/
/*--------------------------------------------------------------*/
void WriteDataToOutputFile(SNEQData *SNEQD, double *I, double *E);

void GetOmegaIntegral_Adaptive(SNEQData *SNEQD,
                               double OmegaMin, double OmegaMax,
                               double *TObjects, double TEnvironment,
                               double *I, double *E);

void GetOmegaIntegral_TrapSimp(SNEQData *SNEQD,
                               double OmegaMin, double OmegaMax,
                               double *TObjects, double TEnvironment,
                               int Intervals, double *I, double *E);

void GetOmegaIntegral_Cliff(SNEQData *SNEQD,
                            double OmegaMin, double OmegaMax,
                            double *TObjects, double TEnvironment,
                            double *I, double *E);

/*--------------------------------------------------------------*/
/*- in SIPFT.cc ------------------------------------------------*/
/*--------------------------------------------------------------*/
void GetSIPFTMatrices(RWGGeometry *G, int WhichSurface,
                      RWGSurface *BS, double R, int NumPoints,
                      cdouble Omega, bool NeedMatrix[NUMSIPFT],
                      HMatrix *MSIPFT[NUMSIPFT]);

/*--------------------------------------------------------------*/
/*- in PlotFlux.cc ---------------------------------------------*/
/*--------------------------------------------------------------*/
void CreateFluxPlot(SNEQData *SNEQD, cdouble Omega, char *Tag);

#endif
