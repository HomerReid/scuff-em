/*
 * scuff-scatter.h -- a standalone code within the scuff-em suite
 *                 -- for solving scattering problems
 *
 * homer reid      -- 6/2011--2/2012
 */
#ifndef SCUFFSCATTER_H
#define SCUFFSCATTER_H

#include <libhrutil.h>
#include <libhmat.h>
#include "libIncField.h"
#include "libscuff.h"

using namespace scuff;


/***************************************************************/
/* data structure containing everything needed to execute a    */
/* scattering calculation                                      */
/***************************************************************/
typedef struct SSData
 {
   RWGGeometry *G;
   HMatrix *M;
   HVector *RHS, *KN;
   cdouble Omega;
   IncField *IF;
   double PowerRadius;
   int nThread;
 } SSData;
 

/***************************************************************/
/* these are the 'output modules' that compute and process the */
/* scattered fields in various ways.                           */
/***************************************************************/
void GetPower(SSData *SSD, char *PowerFile);
void GetMoments(SSData *SSD, char *MomentFile);
void ProcessEPFile(SSData *SSData, char *EPFileName);
void CreateFluxPlot(SSData *SSData, char *MeshFileName);

#endif
