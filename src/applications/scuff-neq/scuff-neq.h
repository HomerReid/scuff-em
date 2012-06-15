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

#define QFLAG_POWER  1
#define QFLAG_XFORCE 2
#define QFLAG_YFORCE 4
#define QFLAG_ZFORCE 8

#define QINDEX_POWER  0
#define QINDEX_XFORCE 1
#define QINDEX_YFORCE 2
#define QINDEX_ZFORCE 3

#define MAXQUANTITIES 4

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
   int NONQ;
   int NTNONQ;

   // HMatrix structures for the BEM matrix and its subblocks
   HMatrix *W;        // BEM matrix 
   HMatrix **T;       // T[no] = T-matrix block for object #no
   HMatrix **U;       // U[no*NO + nop] = // U-matrix block for objects #no, #nop

   // SMatrix structures for overlap matrices
   // note: nq=0,1,2,3 for O^{PF}, O^{xMF}, O^{yMF}, O^{zMF}, 
   // (the matrices for power, x-force, y-force, and z-force.)
   SMatrix **OMatrices; // OMatrices[ no*4 + nq ] = overlap matrix for object no, quantity nq

   // frequency-resolved output files for each object 
   char **ByOmegaFileNames;

 } SNEQData;

SNEQData *CreateSNEQData(char *GeoFile, char *TransFile, int WhichQuantities, int PlotFlux);

void GetFrequencyIntegrand(SNEQData *SNEQD, cdouble Omega, double *FI);

void EvaluateFrequencyIntegral(SNEQData *SNEQD, 
                               double OmegaMin, double OmegaMax,
                               double *TObjects, double TEnvironment, 
                               double AbsTol, double RelTol,
                               double *I, double *E);

#endif
