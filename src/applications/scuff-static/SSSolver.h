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
 * SSSolver.h   -- header file for 'scuff static solver' class
 *
 *
 * homer reid   -- 5/2013
 */

#ifndef SSGEOMETRY_H
#define SSGEOMETRY_H

#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <complex>
#include <cmath>

#include <libhrutil.h>
#include <libhmat.h>
#include <libscuff.h>

#include "StaticSubstrate.h"

#ifdef HAVE_CONFIG_H
	#include <config.h>
#endif

namespace scuff {

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- Class definition for SSSolver                              -*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
typedef void (StaticField)(double *x, void *UserData, double PhiE[4]);

enum SurfType     { PEC = 0, DIELECTRIC=1, LAMBDASURFACE=2 };
enum IntegralType { PHIINTEGRAL = 0, ENORMALINTEGRAL=1 };

/****************************************************************/
/****************************************************************/
/***************************************************************/
class SSSolver   
 { 
   /*--------------------------------------------------------------*/ 
   /*- public class methods ---------------------------------------*/ 
   /*--------------------------------------------------------------*/ 
  public:  

   /* constructor / destructor */
   SSSolver(const char *GeoFileName,
            const char *SubstrateFile=0,
            int pLogLevel = SCUFF_NOLOGGING);
   ~SSSolver();

   /* routines for allocating, and then filling in, the electrostatic BEM matrix */
   HMatrix *AllocateBEMMatrix();
   HMatrix *AssembleBEMMatrix(HMatrix *M = NULL);

   /* lower-level routine for assembling individual BEM matrix blocks */
   void AssembleBEMMatrixBlock(int nsa, int nsb,
                               HMatrix *M, int RowOffset=0, int ColOffset=0);

   /* routines for allocating, and then filling in, the RHS vector */
   HVector *AllocateRHSVector();
   HVector *AssembleRHSVector(double *Potentials, StaticField *SF, 
                              void *UserData, HVector *RHS = NULL);

   /* routine for calculating electric dipole moment */
   HMatrix *GetCartesianMoments(HVector *Sigma, HMatrix *Moments);
   HVector *GetSphericalMoments(HVector *Sigma, int lMax, HVector *Moments);
   HVector *GetSphericalMoments(HVector *Sigma, int WhichSurface, 
                                int lMax, HVector *Moments);

   /* compute fields */
   HMatrix *GetFields(StaticField *SF, void *UserData, HVector *Sigma, HMatrix *X, HMatrix *PhiE);

   /* visualization */
   void PlotChargeDensity(HVector *Sigma, const char *format, ...);

   /*--------------------------------------------------------------------*/ 
   /*- class methods intended for internal use only, i.e. which          */ 
   /*- would be private if we cared about the public/private distinction */
   /*--------------------------------------------------------------------*/ 

   double GetPPI(RWGSurface *Sa, int npa, RWGSurface *Sb, int npb, int WhichIntegral);
   void GetPhiE(int ns, int np, double *X, double PhiE[4]);

   /*--------------------------------------------------------------------*/ 
   /*- class data fields intended for internal use only, i.e. which -----*/ 
   /*- would be private if we cared about the public/private distinction-*/ 
   /*--------------------------------------------------------------------*/ 
   RWGGeometry *G;

   SubstrateData *Substrate;

   char *TransformLabel;
   char *FileBase;

   /*--------------------------------------------------------------*/
   /*- helper functions for contributions of layered dielectric    */
   /*- substrates                                                  */
   /*--------------------------------------------------------------*/
   void AddSubstrateContributionsToBEMMatrixBlock(int nsa, int nsb,
                                                  HMatrix *M,
                                                  int RowOffset,
                                                  int ColOffset);
   void AddSubstratePhiE(int ns, int np, double *X, double PhiE[4]);

 };

}

#endif // #ifndef SSGEOMETRY_H
