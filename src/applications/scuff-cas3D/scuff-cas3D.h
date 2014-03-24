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
 * scuff-cas3D.h  -- header file for 3D casimir module in scuff-EM suite
 *
 * homer reid     -- 2/2012
 */
#ifndef SCUFFCAS3D_H
#define SCUFFCAS3D_H

#include <libhrutil.h>
#include <libhmat.h>
#include <libscuff.h>

using namespace scuff;

/***************************************************************/
/* brillouin zone integration schemes **************************/
/***************************************************************/
#define BZIMETHOD_ADAPTIVE	0
#define BZIMETHOD_MP3		1
#define BZIMETHOD_MP6		2
#define BZIMETHOD_MP10		3
#define BZIMETHOD_MP15		4
#define BZIMETHOD_DEFAULT	BZIMETHOD_ADAPTIVE

#define QUANTITY_ENERGY  1
#define QUANTITY_XFORCE  2
#define QUANTITY_YFORCE  4
#define QUANTITY_ZFORCE  8
#define QUANTITY_TORQUE1 16
#define QUANTITY_TORQUE2 32
#define QUANTITY_TORQUE3 64

#define PREAMBLE_OUT   0
#define PREAMBLE_BYXI  1
#define PREAMBLE_BYXIK 2

/******************************************************************/
/* SC3Data ('scuff-cas3D data') is a structure that contains all  */
/* information needed to compute the contribution of a single     */
/* imaginary frequency and kBloch value to the Casimir quantities.*/
/******************************************************************/
typedef struct SC3Data
 {
   RWGGeometry *G;

   int N, N1;
   HMatrix **TBlocks, **UBlocks, **dUBlocks, *M, *dM;
   int *ipiv;
   HVector *MInfLUDiagonal;

   int WhichQuantities;
   int NumQuantities;

   int NumTorqueAxes;      // this number is in the range 0--3
   double TorqueAxes[9];   // [0,1,2] = x,y,z comps of 1st torque axis; [3,4,5] = 2nd axis, etc
   double GammaMatrix[27];

   double RLBasisVectors[2][2]; // reciprocal-lattice basis vectors

   GTComplex **GTCList;
   int NumTransformations;

   int NTNQ;

   int *Converged;

   double Xi;

   int BZIMethod;
   char *ByXiFileName, *ByXikBlochFileName;
   char *WriteCache;

   // these quantities are used to set limits for adaptive frequency integrations
   int MaxXiPoints, MaxkBlochPoints;
   double AbsTol, RelTol;

   // 20130427 alternative energy calculation
   bool NewEnergyMethod;
   HMatrix *MM1MInf;

   bool UseExistingData;

 } SC3Data;

SC3Data *CreateSC3Data(RWGGeometry *G, char *TransFile,
                       int WhichQuantities, int NumQuantities,
                       int NumTorqueAxes, double TorqueAxes[9],
                       bool NewEnergyMethod);

void WriteFilePreamble(FILE *f, SC3Data *SC3D, int PreambleType);

/***************************************************************/
/* The total Casimir energy (or force, or torque) is an        */
/* integral, over both the brillouin zone and the positive     */
/* imaginary frequency axis, of an integrand F(\xi, kBloch):   */
/*                                                             */
/*  \int_{imag freq} \int_{brillouin zone} F(\xi, kBloch)      */
/*                                                             */
/* The innermost integrand here, F(\xi, kBloch), is returned   */
/* by the routine GetCasimirIntegrand.                         */
/*                                                             */
/* The outer integrand, i.e. the quantity \int_{BZ} F(\xi, k), */
/* is returned by the routine GetXiIntegrand. In other words,  */
/* GetFrequencyIntegrand computes the integral over the        */
/* BZ of the quantity returned by GetCasimirIntegrand.         */
/*                                                             */
/* Note that for non-periodic geometries there is no Brillouin */
/* zone and the two routines return the same thing.            */
/*                                                             */
/* Note: 'EFT' stands for 'energy, force, and torque.'         */
/***************************************************************/
void GetCasimirIntegrand(SC3Data *SC3D, double Xi, double *kBloch, double *EFT);
void GetXiIntegrand(SC3Data *SC3D, double Xi, double *EFT);
void GetXiIntegral(SC3Data *SC3D, double *EFT, double *Error);
void GetXiIntegral2(SC3Data *SC3D, int NumIntervals, double *I, double *E);
void GetMatsubaraSum(SC3Data *SC3D, double Temperature, double *EFT, double *Error);

#endif // #define SCUFFCAS3D_H
