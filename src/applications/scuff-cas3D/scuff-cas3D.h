/*
 * scuff-cas3D.h  -- header file for scuff-cas3D.cc
 *
 * homer reid  -- 2/2012
 */
#ifndef SCUFFCAS3D_H
#define SCUFFCAS3D_H

#include <libhrutil.h>
#include <libhmat.h>
#include <libscuff.h>

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
// SC3Data ('scuff-cas3D data') is a structure that contains all 
// information needed to compute the contribution of a single
// imaginary frequency to the Casimir quantities.
typedef struct SC3Data
 {
   RWGGeometry *G;
   char *ByXiFile;

   int N, N1;
   HMatrix **TBlocks, **UBlocks, **dUBlocks, *M, *dM;

   GTComplex **GTCList;
   int NumTransformations;

   char *WriteCache;
   int nThread;

 } SC3Data;

SC3Data *CreateSC3Data(char *GeoFile, char *TransFile, 
                       char *ByOmegaFile, int nThread);

void GetFrequencyIntegrand(SC3Data *SC3D, double Xi, double *FI);

#endif // #define SCUFFCAS3D_H
