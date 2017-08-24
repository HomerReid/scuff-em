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
 * libSubstrate.h -- header file for C++ library that implements
 *                -- implicit handling of layered material substrates
 *                -- (both electrostatic and full-wave cases)
 *
 * homer reid  -- 3/2017 - 9/2017
 */

#ifndef LIBSUBSTRATE_H
#define LIBSUBSTRATE_H

#include <stdio.h>
#include <math.h>
#include <stdarg.h>

#include <libhrutil.h>
#include <libhmat.h>
#include <libMatProp.h>
#include <libMDInterp.h>

/***************************************************************/
/* data structure for layered material substrate               */
/***************************************************************/
class LayeredSubstrate
 {
public:
   LayeredSubstrate(const char *SubstrateFile);
   ~LayeredSubstrate();

   // full-wave case: get the contribution of the substrate
   //                 to the 6x6 dyadic Green's function
   //                 giving the E,H fields at XD due to 
   //                 J, M currents at XS
   void GetDeltaG(cdouble Omega, double XD[3], double XS[3],
                  cdouble DeltaG[6][6]);

   // electrostatic case: get the contribution of the substrate
   //                     to the potential and electrostatic field
   //                     at XD due to a point source at XS 
   void GetDeltaPhiE(double XD[3], double XS[3],
                     double PhiE[4], double *pG0Correction=0);

   void InitStaticAccelerator1D(double RhoMin, double RhoMax, double z);

// private:

// internal ("private") class methods
   double GetStaticG0Correction(double z);
   double GetStaticG0Correction(double zD, double zS);
   void GetSigmaTwiddle(double zS, double q, double *SigmaTwiddle);
   void UpdateCachedEpsMu(cdouble Omega);
   void GetqIntegral(double RhoMag, double zD, double zS,
                     double qIntegral[3]);

// internal ("private") class data
   char *ErrMsg;
   int NumLayers;
   MatProp *MPMedium;   // material properties of ambient (uppermost) medium
   MatProp **MPLayer;   // MPLayer[n] = properties of layer #(n+1)
   cdouble EpsMedium;   // 
   cdouble  *EpsLayer;  // EpsLayer[n] = permittivity of layer #(n+1)
   cdouble  *MuLayer;   // MuLayer[n]  = permeability of layer #(n+1)
   cdouble OmegaCache;  // frequency at which EpsMedium, EpsLayer, MuLayer were cached
   double *zLayer;      // z[n] = z-coordinate of upper surface of layer #(n+1)
   double zGP;
   int qMaxEval;        // convergence parameters for q integration
   double qAbsTol;
   double qRelTol;
   int PPIOrder;
   int PhiEOrder;
   int WhichIntegral;

   Interp1D *I1D;
   double I1DRhoMin, I1DRhoMax, I1DZ;
 };

void AddPhiE0(double XDest[3], double xs, double ys, double zs, double Q, double PhiE[4]);

#endif // LIBSUBSTRATE_H
