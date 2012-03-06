/*
 * scuff-heat  -- a standalone code within the scuff-em suite
 *             -- for solving heat-transfer and radiation problems
 *
 * homer reid  -- 2/2012
 */
#ifndef SCUFFHEAT_H
#define SCUFFHEAT_H

#include <libhrutil.h>
#include <libhmat.h>
#include <libscuff.h>

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
// SHData ('scuff-heat data') is a structure that contains all 
// information needed to compute the heat radiation or transfer
// at a given frequency.
typedef struct SHData
 {
   RWGGeometry *G;
   char *ByOmegaFile;

   int N1, N2;
   HMatrix **TSelf, **TMedium, **UMedium;
   HMatrix *SymG1, *SymG2;
   HMatrix *W, *W21, *W21SymG1, *W21DSymG2;

   HVector *DV;
   int PlotFlux;

   GTComplex **GTCList;
   int NumTransformations;

   int nThread;

 } SHData;

SHData *CreateSHData(char *GeoFile, char *TransFile, int PlotFlux,
                     char *ByOmegaFile, int nThread);

void GetFrequencyIntegrand(SHData *SHD, cdouble Omega, double *FI);

#endif
