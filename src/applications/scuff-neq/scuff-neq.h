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
   char *ByOmegaFile;
   char *WriteCache;

   int PlotFlux;

   GTComplex **GTCList;
   int NumTransformations;

   int QuantityFlags;
   int NQ;
   int NQNO;
   int NTNQNO;

   // HMatrix structures for the BEM matrix and its subblocks
   HMatrix *W;        // BEM matrix 
   HMatrix **T;       // T[no] = T-matrix block for object #no
   HMatrix **U;       // U[no*NO + nop] = // U-matrix block for objects #no, #nop

   // SMatrix structures for overlap matrices
   // note: nq=0,1,2,3 for O^{PF}, O^{xMF}, O^{yMF}, O^{zMF}, 
   // (the matrices for power, x-force, y-force, and z-force.)
   SMatrix **OMatrices // OMatrices[ no*4 + nq ] = nqth overlap matrix for object no

   // frequency-resolved output files for each object 
   char **ByOmegaFileNames;
   //  FILE **ByOmegaFiles;

 } SNEQData;

SNEQData *CreateSNEQData(char *GeoFile, char *TransFile, char *ByOmegaFile, 
                         int WhichQuantities, int PlotFlux);

void GetFrequencyIntegrand(SNEQData *SNEQD, cdouble Omega, double *FI);

#endif
