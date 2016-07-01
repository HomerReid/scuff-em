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
 * BZIntegration.h  -- header for BZIntegration.cc
 *
 * homer reid       -- 3/2015
 */
#ifndef BZINTEGRATION_H
#define BZINTEGRATION_H

#include "libhrutil.h"
#include "libhmat.h"
#include "libTriInt.h"

// values for BZIMethod
#define BZI_DEFAULT  0
#define BZI_CC       1
#define BZI_TC       2
#define BZI_RADIAL   3

// default parameters for adaptive integration
#define DEF_BZIMAXEVALS 1000    // max # brillouin-zone samples
#define DEF_BZIRELTOL   1.0e-2  // relative tolerance
#define DEF_BZIABSTOL   1.0e-10 // absolute tolerance

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef void (*BZIFunction)(void *UserData,
                            cdouble Omega, double *kBloch,
                            double *BZIntegrand);

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct GetBZIArgStruct
{
  // information on the Brillouin-zone integrand function
  BZIFunction BZIFunc;
  void *UserData;
  int FDim;            // number of doubles in the integrand vector
  int SymmetryFactor;  // either 1, 2, 4, or 8

  // information on the lattice geometry
  // RLBasis = "reciprocal lattice basis"
  //         = 3 x D matrix, columns = reciprocal lattice vectors
  // RLBasis[nc][nd] = ncth cartesian component of ndth RL vector
  HMatrix *RLBasis;
  double BZVolume;

  // information on the cubature scheme
  int BZIMethod;
  int Order;      // cubature order for fixed cubature schemes
  int MaxEvals;   // max # integrand samples for adaptive schemes
  double RelTol;  // relative tolerance for adaptive schemes
  double AbsTol;  // adaptive tolerance for adaptive schemes
  
  // fields used internally that may be ignored by the caller
  double kRho;
  cdouble Omega;
  int BufSize;
  double *DataBuffer[3]; // internally allocated

  // return values 
  int NumCalls;       // actual # integrand samples (return value)
  double *BZIError;   // error (for adaptive schemes)
  
} GetBZIArgStruct;

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetBZIntegral(GetBZIArgStruct *Args, cdouble Omega,
                   double *BZIntegral);

GetBZIArgStruct *InitBZIArgs(int argc=0, char **argv=0);
void UpdateBZIArgs(GetBZIArgStruct *BZIArgs, HMatrix *RLBasis,
                   double RLVolume);

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *GetRLBasis(HMatrix *LBasis,
                    double *pLVolume,
                    double *pRLVolume);

double SetRLBasis(HMatrix *LBasis,
                  HMatrix *RBasis,
                  double *pLVolume,
                  double *pRLVolume);

#endif //#ifndef BZINTEGRATION_H
