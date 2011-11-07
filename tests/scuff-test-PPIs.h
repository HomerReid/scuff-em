/*
 * scuff-test-PPIs.h  -- header for scuff-test-PPI standalone test program
 *                       for computing panel-panel integrals
 * 
 * homer reid         -- 11/2005 -- 11/2011
 */
#ifndef SCUFF_TEST_PPIS_H
#define SCUFF_TEST_PPIS_H

#include "libscuff.h" 
#include "libscuffInternals.h" 

/***************************************************************/
/* routine for brute-force evaluation of panel-panel integrals */
/***************************************************************/
void GetPPIs_BruteForce(GetPPIArgStruct *Args);

#endif
