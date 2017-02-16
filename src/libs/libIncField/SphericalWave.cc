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
 * SphericalWave .cc   -- spherical wave implementation of IncField
 *
 * homer reid          -- 11/2009 -- 2/2012
 */

#include <string.h>
#include <math.h>

#include <libhrutil.h>
#include <libIncField.h>
#include <libSpherical.h>

#define II cdouble(0.0,1.0)

/**********************************************************************/
/* constructor and field-setting routines *****************************/
/**********************************************************************/
SphericalWave::SphericalWave(int NewL, int NewM, int NewType)
{
  L=NewL;
  M=NewM;
  Type=NewType;
}

void SphericalWave::SetL(int NewL) { L=NewL; }
void SphericalWave::SetM(int NewM) { M=NewM; }
void SphericalWave::SetType(int NewType) { Type=NewType; }

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
void SphericalWave::GetFields(const double X[3], cdouble EHC[6])
{

  cdouble K=sqrt(Eps*Mu) * Omega;
  cdouble Z=ZVAC*sqrt(Mu/Eps);
  
  // convert the evaluation point to spherical coordinates 
  double XX[3]; 
  memcpy(XX, X, 3*sizeof(double));
  double r, Theta, Phi;
  CoordinateC2S(XX, &r, &Theta, &Phi);

  // get the M and N vector spherical harmonics 
  cdouble MVec[3], NVec[3];
  GetMNlm(L, M, K, r, Theta, Phi, LS_REGULAR, MVec, NVec);
  
  // set the spherical components of E and H to the 
  // proper linear combinations of the M and N functions
  cdouble EHS[6]; // 'E,H spherical'
  if (Type==SW_MAGNETIC)
   { 
     EHS[0] = MVec[0];
     EHS[1] = MVec[1];
     EHS[2] = MVec[2];
     EHS[3] = -NVec[0] / Z;
     EHS[4] = -NVec[1] / Z;
     EHS[5] = -NVec[2] / Z;
   }
  else // Type==SW_ELECTRIC
   { 
     EHS[0] = NVec[0];
     EHS[1] = NVec[1];
     EHS[2] = NVec[2];
     EHS[3] = MVec[0] / Z;
     EHS[4] = MVec[1] / Z;
     EHS[5] = MVec[2] / Z;
   };

  // convert the spherical components of E and H to
  // cartesian components 
  VectorS2C(Theta, Phi, EHS+0, EHC+0);
  VectorS2C(Theta, Phi, EHS+3, EHC+3);

}
