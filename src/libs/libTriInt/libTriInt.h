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
 * libTriInt.h  -- header for TriInt library providing routines for
 *              -- fixed-order and adaptive cubature over a triangle 
 *
 * homer reid   -- 12/2008
 */

#ifndef LIBTRIINT_H
#define LIBTRIINT_H 

#include <libhrutil.h>
#include <libhmat.h>
#include <libSGJC.h>

/***************************************************************/
/* Type definition for user's integrand routine.               */
/*                                                             */
/*  On input, X[0..2] are cartesian coordinates of evaluation  */
/*  point (guaranteed to be in the interior of the triangle    */
/*  specified by the user in the call to TriInt) and UserData  */
/*  is the void pointer supplied by the user in the call to    */
/*  TriInt.                                                    */
/*  On output, F[0...nFun-1] must be the values of the         */
/*  vector at the point X, where nFun is the parameter passed  */
/*  TriInt.                                                    */ 
/***************************************************************/
typedef void (*TriIntFun)(double *X, void *UserData, double *F);

/***************************************************************/
/* Prototypes for user-callable functions.                     */
/***************************************************************/
/* routines that are actually useful */

double *GetTCR(int Order, int *NumPts);


int TriIntFixed(TriIntFun F, int nFun, void *UserData,
                double *V1, double *V2, double *V3, 
                int Order, double *Result);

void TriIntEmbedded(TriIntFun F, int nFun, void *UserData,
                    double *V1, double *V2, double *V3,
                    double *I, double *Error);

/***************************************************************/
/***************************************************************/
/***************************************************************/
void *CreateDCUTRIWorkspace(int numfun, int maxpts);
int DCUTRI(void *opW, double **Vertices, TriIntFun func,
           void *parms, double epsabs, double epsrel,
           double *I, double *abserr);

/***************************************************************/
/* Return a pointer to an internal buffer containing the (1D)  */
/* Clenshaw-Curtis rule with numpts points. Only odd values of */
/* NumPts between 5 and 99 are supported. (Otherwise, the      */
/* return value is NULL.)                                      */
/***************************************************************/
double *GetCCRule(int NumPts);

/***************************************************************/
/***************************************************************/
/***************************************************************/
double *GetLebedevRule(int NumPts);

/***************************************************************/
/* embedded clenshaw-curtis cubature in 2D *********************/
/***************************************************************/
void ECC2D(int p, double xMin[2], double xMax[2],
           integrand Integrand, void *UserData, int nFun,
           bool xySymmetric,
           double *AllValues, double *OuterValues,
           double *Result, double *Error);

void ECC(int p, double xMin, double xMax,
         integrand Integrand, void *UserData, int nFun,
         double *AllValues, double *OuterValues,
         double *Result, double *Error);

/***************************************************************/
/***************************************************************/
/***************************************************************/
int CCCubature(int Order, unsigned fdim, integrand f, void *fdata,
	       unsigned dim, const double *xmin, const double *xmax, 
	       size_t maxEval, double reqAbsError, double reqRelError,
               error_norm norm, double *Integral, double *Error);

int RRCubature(int Order, int *Orders, 
               int FDim, integrand f, void *UserData,
	       int IDim, const double *Lower, const double *Upper,
               double *Integral, double *Error);

/***************************************************************/
/* CliffFunction is almost identical to integrand in libSGJC,  */
/* but it takes an additional argument consisting of a vector  */
/* of bools; if Skip[nf]==true then the evaluation of integrand*/
/* component nf is unnecessary and may be skipped. (Useful for */
/* integrating a vector-valued function in which the integrals */
/* of some components converge more rapidly than others.       */
/***************************************************************/
typedef int (*CliffFunction) (unsigned ndim, const double *x, void *params,
                              unsigned fdim, const bool *Skip, double *fval);

int IntegrateCliffFunction(CliffFunction fCliff, void *UserData, int nFun,
                           double xMin, double xMax, double xCliff,
                           double AbsTol, double RelTol, 
                           double *Integral, double *Error,
                           char *LogFileName=0);

/***************************************************************/
/* Routine for adaptive lattice summation to a desired error   */
/* tolerance, Sum = \sum_{L} F(L)                              */
/* where F is a user-supplied vector-valued summand            */
/* and L ranges over all points in a D-dimensional lattice.    */
/*                                                             */
/* note that SummandFunction should ACCUMULATE the             */
/* contribution of lattice point U to the sum, i.e. it should  */
/* implement something like                                    */
/*  Sum[ns] += f_{ns}[U]                                       */
/* and not                                                     */
/*  Sum[ns] = f_{ns}[U].                                       */
/***************************************************************/
typedef void (*SummandFunction)(double *L, void *UserData, double *Sum);

int GetLatticeSum(SummandFunction Summand, void *UserData, int nSum,
                  HMatrix *LBasis, double *Sum,
                  double AbsTol=0.0, double RelTol=1.0e-2, 
                  int MaxCells=1000);

#endif 
