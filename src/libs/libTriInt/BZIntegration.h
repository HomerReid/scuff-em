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

#define BZI_DEFAULT  0
#define BZI_CC       1
#define BZI_TC       2
#define BZI_ADAPTIVE 3

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
  int FDim;         // number of doubles in the integrand vector
  bool BZSymmetric; // true if f(kx,ky) = f(ky,kx)

  // information on the cubature scheme
  int BZIMethod;
  int Order;      // cubature order for fixed cubature schemes
  int MaxEvals;   // max # integrand samples for adaptive schemes
  double RelTol;  // relative tolerance for adaptive schemes
  double AbsTol;  // adaptive tolerance for adaptive schemes
  bool Reduced;   // get full BZ integral from reduced BZ integral
                  //  times 2 (1D) or 4 (2D)

  // RLBasis = "reciprocal lattice basis"
  //         = 3 x D matrix, columns = reciprocal lattice vectors
  // RLBasis[nc][nd] = ncth cartesian component of ndth RL vector
  HMatrix *RLBasis;
  double BZVolume;
  
  // fields used internally that may be ignored by the caller
  int BZIErrorSize;
  cdouble Omega;

  // return values 
  int NumCalls;     // actual # integrand samples (return value)
  double *BZIError; // internally allocated; only non-null for adaptive schemes
  
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
