/*
 * scuff-test-EEIs.h  -- header for scuff-test-EEI standalone test program
 *                       for computing edge-edge interactions 
 * 
 * homer reid         -- 11/2005 -- 11/2011
 */
#ifndef SCUFF_TEST_EEIS_H
#define SCUFF_TEST_EEIS_H

#include "libscuff.h" 
#include "libscuffInternals.h" 

using namespace scuff;

/******************************************************************/
/* routine for brute-force evaluation of panel-panel interactions */
/******************************************************************/
void GetPPIs_BruteForce(GetPPIArgStruct *Args, int PlotFits);

#endif
