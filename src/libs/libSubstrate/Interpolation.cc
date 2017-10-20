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
 * libSubstrate/Interpolation.cc -- accelerate calculations of
 *                               -- substrate Green's functions
 *                               -- using interpolation tables
 *
 * homer reid             -- 3/2017
 *
 */
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>

#include "libhrutil.h"
#include "libMDInterp.h"
#include "libSGJC.h"
#include "libSpherical.h"
#include "libSubstrate.h"

/***************************************************************/
/* entry point with the proper prototype for passage to the    */
/* Interp1D() initialization routine.                          */
/***************************************************************/
typedef struct fInterpData
 {
   LayeredSubstrate *Substrate;
   double zD, zS;
 } fInterpData;

void fInterp1DStatic(double Rho, void *UserData, double *fInterp)
{
  fInterpData *fID   = (fInterpData *)UserData;
  LayeredSubstrate *Substrate = fID->Substrate;
  double zD                   = fID->zD;
  double zS                   = fID->zS;

  double qI[3], dqI[3];
  if (Rho==0.0)
   { double DeltaRho=1.0e-5;
     Substrate->GetqIntegral(Rho,          zD, zS, qI);
     Substrate->GetqIntegral(Rho+DeltaRho, zD, zS, dqI);
     VecPlusEquals(dqI, -1.0, qI, 3);
     VecScale(dqI, 1.0/DeltaRho, 3);
   }
  else
   { double DeltaRho=1.0e-5 * fabs(Rho);
     Substrate->GetqIntegral(Rho+DeltaRho, zD, zS, dqI);
     Substrate->GetqIntegral(Rho-DeltaRho, zD, zS, qI);
     VecPlusEquals(dqI, -1.0, qI, 3);
     VecScale(dqI, 1.0/(2.0*DeltaRho), 3);
     Substrate->GetqIntegral(Rho,          zD, zS, qI);
   };
 
  fInterp[0] =  qI[0];
  fInterp[1] = dqI[0];
  fInterp[2] =  qI[1];
  fInterp[3] = dqI[1];
  fInterp[4] =  qI[2];
  fInterp[5] = dqI[2];

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void LayeredSubstrate::InitStaticAccelerator1D(double RhoMin,
                                               double RhoMax,
                                               double z)
{
  if (      I1D 
       &&  (I1DOmega   ==0.0     )
       &&  (I1DRhoMin  <= RhoMin )
       &&  (I1DRhoMax  >= RhoMax )
       &&  (I1DZ       == z      )
     )
   { Log(" reusing existing substrate interpolation table");
     return;
   };

  I1DOmega  = 0.0;
  I1DRhoMin = fmin(RhoMin, I1DRhoMin);
  I1DRhoMax = fmax(RhoMax, I1DRhoMax);
  
  if (I1D) delete I1D;

  Log(" (re)allocating substrate I1D(%g,%g,%g)...",I1DRhoMin,I1DRhoMax,z);

  int NGrid=100;
  char *s=getenv("SCUFF_SUBSTRATE_NGRID");
  if (s)
  sscanf(s,"%i",&NGrid);
   
  struct fInterpData MyfID, *fID=&MyfID;
  fID->Substrate = this;
  fID->zD        = fID->zS=z;
  I1D= new Interp1D(RhoMin, RhoMax, NGrid, 3, fInterp1DStatic,
                    (void *)fID);
  I1DZ      = z;

  Log(" ...done with substrate I1D");

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void LayeredSubstrate::InitAccelerator1D(cdouble Omega, double RhoMin, double RhoMax, double z, bool EEOnly)
{
  if (      I1D 
       &&  (I1DOmega   ==0.0     )
       &&  (I1DRhoMin  <= RhoMin )
       &&  (I1DRhoMax  >= RhoMax )
       &&  (I1DZ       == z      )
     )
   { Log(" reusing existing substrate interpolation table");
     return;
   };

  I1DOmega  = 0.0;
  I1DRhoMin = fmin(RhoMin, I1DRhoMin);
  I1DRhoMax = fmax(RhoMax, I1DRhoMax);
  
  if (I1D) delete I1D;

  Log(" (re)allocating substrate I1D(%g,%g,%g)...",I1DRhoMin,I1DRhoMax,z);

  int NGrid=100;
  char *s=getenv("SCUFF_SUBSTRATE_NGRID");
  if (s)
  sscanf(s,"%i",&NGrid);
   
  struct fInterpData MyfID, *fID=&MyfID;
  fID->Substrate = this;
  fID->zD        = fID->zS=z;
  I1D= new Interp1D(RhoMin, RhoMax, NGrid, 3, fInterp1DStatic,
                    (void *)fID);
  I1DZ      = z;

  Log(" ...done with substrate I1D");

}
