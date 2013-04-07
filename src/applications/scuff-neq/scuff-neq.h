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

/****************************************************************/
/* SNEQData ('scuff-neq data') is a structure that contains all */
/* information needed to run computations a given frequency.    */
/****************************************************************/
typedef struct SNEQData
 {
   RWGGeometry *G;
   char *WriteCache;

   int PlotFlux;

   GTComplex **GTCList;
   int NumTransformations;

   int QuantityFlags;
   int NQ;
   int NSNQ;
   int NTNSNQ;

   // HMatrix structures for the BEM matrix and its subblocks
   HMatrix *W;        // BEM matrix 
   HMatrix **T;       // T[ns] = T-matrix block for surface #ns
   HMatrix **TSelf;   //
   HMatrix **U;       // U[ns*NS + nsp] = // U-matrix block for surfaces #ns, #nsp
   void *Buffer[3];

   // the NMth slot in this array of flags is 1 iff we will need
   // to compute the NMth type of overlap matrix
   bool NeedMatrix[SCUFF_NUM_OMATRICES]; 

   // SMatrix structures for overlap matrices.
   // SArray[ ns*8 + nm ] = overlap matrix of type #nm for surface #ns 
   // (here 8 = SCUFF_NUM_OMATRICES)
   SMatrix ***SArray;
   
   bool UseSGJFormalism;

   char *FileBase;

 } SNEQData;

SNEQData *CreateSNEQData(char *GeoFile, char *TransFile, 
                         int WhichQuantities, int PlotFlux, 
                         char *FileBase, bool UseSGJFormalism);

int GetIndex(SNEQData *SNEQD, int nt, int nss, int nsd, int nq);
void GetFlux(SNEQData *SNEQD, cdouble Omega, double *Flux);

void EvaluateFrequencyIntegral(SNEQData *SNEQD,
                               double OmegaMin, double OmegaMax,
                               double *TObjects, double TEnvironment,
                               double AbsTol, double RelTol,
                               double *I, double *E);

void EvaluateFrequencyIntegral2(SNEQData *SNEQD,
                                double OmegaMin, double OmegaMax,
                                double *TObjects, double TEnvironment,
                                int Intervals, double *I, double *E);

void SpeedTest(char *Greeting);
#endif
