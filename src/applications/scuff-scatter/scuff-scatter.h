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

#ifdef SCUFF
	#include "libscuff.h"
#else // RWG
	#include "libRWG.h"
#endif


/***************************************************************/
/* data structure containing everything needed to pass data    */
/* to the output modules                                       */
/***************************************************************/
typedef struct SSData
 {
   RWGGeometry *G;
   RWGGeometry *G0;
   HVector *KN;
   cdouble Omega;
   void *opIFD;
   int nThread;
 } SSData;
 

/***************************************************************/
/* these are the 'output modules' that compute and process the */
/* scattered fields in various ways.                           */
/***************************************************************/
void GetPower(SSData *SSD, double R, double *PScat, double *PTot);
void GetPower_BF(SSData *SSD, double R, double *PScat, double *PTot);
void ProcessEPFile(SSData *SSData, char *EPFileName, char *ObjectLabel);
void CreateFluxPlot(SSData *SSData, char *MeshFileName);

#endif
