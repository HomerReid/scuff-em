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

typedef struct ScuffHeatData
 {
   RWGGeometry *G;
   char *ByOmegaFile;
   HMatrix *M0, *M1, *M2;
   int nThread;

 } ScuffHeatData;

void GetFrequencyIntegrand(ScuffHeatData *SHD, cdouble Omega, double *I);

#endif
