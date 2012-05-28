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

#define QUANTITY_POWER  1
#define QUANTITY_XFORCE 2
#define QUANTITY_YFORCE 4
#define QUANTITY_ZFORCE 8

/****************************************************************/
/* SNEQData ('scuff-neq data') is a structure that contains all */
/* information needed to run computations a given frequency.    */
/****************************************************************/
typedef struct SNEQData
 {
   RWGGeometry *G;
   char *ByOmegaFile;

   GTComplex **GTCList;
   int NumTransformations;

   int WhichQuantities;
   int NumQuantities;     
   int NTNQ;

   char *WriteCache;

   HMatrix *W;        // BEM matrix 
   HMatrix **T;       // T[no] = T-matrix block for object #no
   HMatrix **U;       // U[no*NO + nop] = // U-matrix block for objects #no, #nop

   SMatrix **OPF;     // OCross[no] = power-flux matrix for object #no
   SMatrix **OiMF;    // OCross[3*no + i] = i-momentum flux matrix for object #no

 } SNEQData;

SNEQData *CreateSNEQData(char *GeoFile, char *TransFile, char *ByOmegaFile, 
                         int WhichQuantities, int PlotFlux);

void GetFrequencyIntegrand(SNEQData *SNEQD, cdouble Omega, double *FI);

#endif
