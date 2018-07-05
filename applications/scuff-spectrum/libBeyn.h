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
 * libBeyn.h -- header file for libBeyn, a simple implementation of
 *           -- Beyn's algorithm for nonlinear eigenproblems
 *           -- 
 *           -- This is packaged together with SCUFF-EM, but
 *           -- it is really a standalone independent entity
 *           -- for general-purpose use in solving nonlinear
 *           -- eigenproblems.
 */


#ifndef LIBBEYN_H
#define LIBBEYN_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>

#include <libhmat.h>

/***************************************************************/
/* prototype for user-supplied function passed to BeynMethod.  */
/* The user's function should replace VHat with                */
/*  Inverse[ M(z) ] * VHat.                                    */
/***************************************************************/
typedef void (*BeynFunction)(cdouble z, void *UserData, HMatrix *VHat, HMatrix *MVHat);

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct BeynSolver
{
   int M;   // dimension of matrices
   int L;   // number of columns of VHat matrix

   HVector *Eigenvalues, *EVErrors, *Residuals;
   HMatrix *Eigenvectors;
   HMatrix *A0, *A1, *A0Coarse, *A1Coarse, *MInvVHat;
   HMatrix *VHat;
   HVector *Sigma;
   cdouble *Workspace;

 } BeynSolver;

// constructor, destructor
BeynSolver *CreateBeynSolver(int M, int L);
void DestroyBeynSolver(BeynSolver *Solver);

// reset the random matrix VHat used in the Beyn algorithm
// 
void ReRandomize(BeynSolver *Solver, unsigned int RandSeed=0);

// for both of the following routines,
// the return value is the number of eigenvalues found,
// and the eigenvalues and eigenvectors are stored in the
// Lambda and Eigenvectors fields of the BeynSolver structure

// Beyn method for circular contour of radius R,
// centered at z0, using N quadrature points
int BeynSolve(BeynSolver *Solver,
              BeynFunction UserFunction, void *UserData,
              cdouble z0, double R, int N=25);

// Beyn method for elliptical contour of horizontal, vertical
// radii Rx, Ry, centered at z0, using N quadrature points
int BeynSolve(BeynSolver *Solver,
              BeynFunction UserFunction, void *UserData,
              cdouble z0, double Rx, double Ry, int N=25);

#endif
