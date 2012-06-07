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
 * scuff-caspol.h    -- header file for scuff-caspol
 *
 * homer reid        -- 1/2012
 *
 */
#ifndef SCUFF_CASPOL_H
#define SCUFF_CASPOL_H

#include <stdio.h>
#include <math.h>
#include <stdarg.h>

#include <libhrutil.h>
#include <libhmat.h>
#include <libMDInterp.h>
#include <libscuff.h>

//using namespace scuff;

// default tolerances for numerical summation and integration
#define ABSTOL 1.0e-8
#define RELTOL 1.0e-3

// the minimum frequency at which we do calculations
#define XIMIN 1.0e-6

/***************************************************************/
/* PolModel is a simple class used to describe the frequency-  */
/* dependent polarizability of atoms and molecules. it exports */
/* just a single class method, GetPolarizability(), which takes*/
/* an imaginary frequency as an input and computes values of   */
/* the polarizability tensor at that frequency.                */
/***************************************************************/
class PolModel
 {
public:
   // constructor (the parameter is a filename)
   PolModel(const char *PolFile);

   // routine to compute the polarizability 
   void GetPolarizability(double Xi, double *Alpha);

   // implementation-dependent data
   Interp1D *PolInterp;

 };

using namespace scuff;

/***************************************************************/
/* SCPData ('scuff-caspol-data') is the basic structure passed */
/* around among the various routines; it contains everything   */
/* needed to compute the casimir-polder potential              */
/***************************************************************/
typedef struct SCPData
 {
   RWGGeometry *G;
   HMatrix *M;
   HVector *KN;

   PolModel *PM;

   HMatrix *EPList;

   int nThread;
   double AbsTol, RelTol;
   double XiMin;
   FILE *ByXiFile;

 } SCPData; 

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetCPIntegrand(SCPData *SCP, double Xi, double *U);
void EvaluateMatsubaraSum(SCPData *SCP, double Temperature, double *U);
void EvaluateFrequencyIntegral(SCPData *SCP, double *U);

#endif
