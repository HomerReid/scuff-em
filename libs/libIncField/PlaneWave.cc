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
 * PlaneWave.cc   -- plane wave implementation of IncField
 *
 * homer reid     -- 11/2009 -- 2/2012
 */

#include <string.h>
#include <math.h>
#include <stdio.h>

#include "libIncField.h"
#define II cdouble(0.0,1.0)

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
PlaneWave::PlaneWave(const cdouble pE0[3], const double pnHat[3],
		     const char *Label)
{
  memcpy(E0, pE0, 3*sizeof(cdouble));
  memcpy(nHat, pnHat, 3*sizeof(double));
  SetRegionLabel(Label);
}

PlaneWave::PlaneWave()
{  E0[0]=1.0;    E0[1]=0.0;    E0[2]=0.0;
 nHat[0]=0.0;  nHat[1]=0.0;  nHat[2]=1.0;
  SetRegionLabel(0);
}

PlaneWave::PlaneWave(const PlaneWave &PW) : IncField(PW)
{ 

  // Copy PW attributes
  memcpy(E0, PW.E0, 3*sizeof(cdouble));
  memcpy(nHat, PW.nHat, 3*sizeof(double));
  //SetRegionLabel(PW.RegionLabel);
}

PlaneWave::~PlaneWave()
{ 
  // no malloc'ed data to free
}

void PlaneWave::SetE0(cdouble pE0[3])
 { memcpy(E0, pE0, 3*sizeof(cdouble)); }

void PlaneWave::SetnHat(double pnHat[3])
 { memcpy(nHat, pnHat, 3*sizeof(double)); }

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
void PlaneWave::GetFields(const double X[3], cdouble EH[6])
{

  cdouble K=sqrt(Eps*Mu) * Omega;
  cdouble Z=ZVAC*sqrt(Mu/Eps);
  cdouble ExpFac=exp(II*K*(nHat[0]*X[0] + nHat[1]*X[1] + nHat[2]*X[2]));

  EH[0] = E0[0] * ExpFac;
  EH[1] = E0[1] * ExpFac;
  EH[2] = E0[2] * ExpFac;

  /* H = (nHat \cross E) / Z */
  EH[3] = (nHat[1]*EH[2] - nHat[2]*EH[1]) / Z;
  EH[4] = (nHat[2]*EH[0] - nHat[0]*EH[2]) / Z ;
  EH[5] = (nHat[0]*EH[1] - nHat[1]*EH[0]) / Z;
  
} 

/***************************************************************/
/* overrides the default implementation of this method in      */
/* the base class                                              */
/***************************************************************/
void PlaneWave::GetFieldGradients(const double X[3], cdouble dEH[3][6])
{
  cdouble EH[6];
  GetFields(X, EH);
  cdouble K=sqrt(Eps*Mu) * Omega;
  for(int i=0; i<3; i++)
   for(int j=0; j<6; j++)
    dEH[i][j] = II*K*nHat[i]*EH[j];
}
