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
 * GBarAccelerator.h --
 *                   --
 * Homer Reid        -- 6/2014
 */

#ifndef GBARACCELERATOR_H
#define GBARACCELERATOR_H

#include "libscuff.h"
#include "libMDInterp.h"

namespace scuff {

/***************************************************************/
/* routine for computing the periodic green's function via     */
/* ewald summation                                             */
/***************************************************************/
void GBarVDEwald(double *R, cdouble k, double *kBloch,
                 double (*LBV)[3], int LDim,
                 double E, bool ExcludeInnerCells, cdouble *GBarVD);

/***************************************************************/
/* interpolation-based acceleration of periodic GF evaluation  */
/***************************************************************/
typedef struct GBarAccelerator
 {
   cdouble k;
   double *kBloch;
   bool ExcludeInnerCells;

   int LDim;
   double LBV[3][3];
   double RhoMin, RhoMax, LMax;

   bool ForceFullEwald;

   Interp2D *I2D;
   Interp3D *I3D;

 } GBarAccelerator;

/***************************************************************/
/***************************************************************/
/***************************************************************/
GBarAccelerator *CreateGBarAccelerator(HMatrix *LBasis,
                                       double RhoMin, double RhoMax,
                                       cdouble k, double *kBloch,
                                       double RelTol, bool ExcludeInnerCells);

void DestroyGBarAccelerator(GBarAccelerator *GBA);

cdouble GetGBar(double R[3], GBarAccelerator *GBA,
                cdouble *dGBar=0, cdouble *ddGBar=0,
                bool ForceFullEwald=false);

} // namespace scuff 

#endif // #ifndef GBARACCELERATOR_H
