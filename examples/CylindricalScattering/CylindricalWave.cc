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
 * CylindricalWave .cc   -- spherical wave implementation of IncField
 *
 * homer reid          -- 11/2009 -- 2/2012
 */

#include <string.h>
#include <math.h>

#include <libhrutil.h>
#include <libIncField.h>
#include <libSpherical.h>
#include "CylindricalWave.h"

#define II cdouble(0.0,1.0)

/**********************************************************************/
/* RPZ = 'rho, phi, z' ************************************************/
/**********************************************************************/
void CoordinateCar2Cyl(double XYZ[3], double RPZ[3])
{
  RPZ[0] = sqrt(XYZ[0]*XYZ[0] + XYZ[1]*XYZ[1]);
  RPZ[1] = atan2(XYZ[1], XYZ[0]);
  RPZ[2] = XYZ[2];
}

void VectorCyl2Car(double Phi, cdouble VCyl[3], cdouble VCar[3])
{ 
  double CosPhi=cos(Phi), SinPhi=sin(Phi);

  VCar[0] = CosPhi*VCyl[0] - SinPhi*VCyl[1];
  VCar[1] = SinPhi*VCyl[0] + CosPhi*VCyl[1];
  VCar[2] = VCyl[2];
}

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
void GetCylMN(cdouble k0, int Nu, double kz, int WaveType,
              double RPZ[3], cdouble MVec[3], cdouble NVec[3])
{
  cdouble kRho2 = k0*k0 - kz*kz;
  cdouble kRho  = sqrt(kRho2);
  cdouble Zeta  = kRho*RPZ[0]; 

  char WhichFunction='J';
  switch(WaveType)
   { case CW_REGULAR:  WhichFunction='J'; break;
     case CW_OUTGOING: WhichFunction='O'; break;
     case CW_INCOMING: WhichFunction='T'; break;
     default: ErrExit("unknown wave type in GetCylMN");
   };
   
  cdouble RNu, NuRNuOverZeta, RNuPrime, RArray[3];
  double Workspace[12];
  if (Nu==0)
   { int NuMin=0;
     int NumNu=2;
     AmosBessel(WhichFunction, Zeta, NuMin, NumNu, false, RArray, Workspace);
     RNu           = RArray[0];
     RNuPrime      = -1.0*RArray[1];
     NuRNuOverZeta = 0.0;
   }
  else
   { int NuMin=Nu-1;
     int NumNu=3;
     AmosBessel(WhichFunction, Zeta, NuMin, NumNu, false, RArray, Workspace);
     RNu          = RArray[1];
     RNuPrime     = 0.5*(RArray[0] - RArray[2]);
     if (RPZ[0]==0.0)
      NuRNuOverZeta = (Nu==1 && WaveType==CW_REGULAR) ? 0.5 : 0.0;
     else
      NuRNuOverZeta = ((double)Nu)*RNu/Zeta;
   };

  cdouble ExpFac = exp( II*(Nu*RPZ[1] + kz*RPZ[2]) );

  MVec[0] = II*NuRNuOverZeta * ExpFac;
  MVec[1] = -RNuPrime * ExpFac;
  MVec[2] = 0.0;

  NVec[0] = -1.0*kz*RNuPrime*ExpFac/k0;
  NVec[1] = -II*kz*NuRNuOverZeta*ExpFac/k0;
  NVec[2] =  II*kRho*RNu*ExpFac/k0;

}

/**********************************************************************/
/* constructor and field-setting routines *****************************/
/**********************************************************************/
CylindricalWave::CylindricalWave(int NewNu, double Newkz, 
                                 int NewPolarization, int NewType)
{
  Nu=NewNu;
  kz=Newkz;
  Polarization=NewPolarization;
  WaveType=NewType;
}

void CylindricalWave::SetNu(int NewNu) 
 { Nu=NewNu; }
void CylindricalWave::SetKz(double Newkz) 
 { kz=Newkz; }
void CylindricalWave::SetWaveType(int Type) 
 { WaveType=Type; }
void CylindricalWave::SetPolarization(int NewPolarization) 
 { Polarization=NewPolarization; }

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
void CylindricalWave::GetFields(const double X[3], cdouble EH[6])
{
  cdouble k0=sqrt(Eps*Mu) * Omega;
  cdouble Z=ZVAC*sqrt(Mu/Eps);
  
  // convert the evaluation point to cylindrical coordinates
  double XYZ[3], RPZ[3];
  memcpy(XYZ, X, 3*sizeof(double));
  CoordinateCar2Cyl(XYZ, RPZ);

  // get the M and N vector cylindrical wave functions
  cdouble MVec[3], NVec[3];
  GetCylMN(k0, Nu, kz, WaveType, RPZ, MVec, NVec);

  // set the spherical components of E and H to the
  // proper linear combinations of the M and N functions
  cdouble EHCyl[6]; // 'E,H cylindrical'
  if (Polarization==CW_TE2Z)
   { 
     EHCyl[0] = MVec[0];
     EHCyl[1] = MVec[1];
     EHCyl[2] = MVec[2];
     EHCyl[3] = -NVec[0] / Z;
     EHCyl[4] = -NVec[1] / Z;
     EHCyl[5] = -NVec[2] / Z;
   }
  else // Type==CW_TM2Z
   { 
     EHCyl[0] = NVec[0];
     EHCyl[1] = NVec[1];
     EHCyl[2] = NVec[2];
     EHCyl[3] = MVec[0] / Z;
     EHCyl[4] = MVec[1] / Z;
     EHCyl[5] = MVec[2] / Z;
   };

  // convert the cylindrical components of E and H to
  // cartesian components
  cdouble EHCar[6];
  double Phi = RPZ[1];
  VectorCyl2Car(Phi, EHCyl+0, EHCar+0);
  VectorCyl2Car(Phi, EHCyl+3, EHCar+3);

  // rotate to a new coordinate system in which the 
  // cylinder axis is the X axis; we do this because  
  // the X axis is always the axis of periodicity for 
  // SCUFF-EM geometries of infinite extent in one
  // spatial direction 
  EH[0*3 + 0] =      EHCar[0*3 + 2];
  EH[0*3 + 1] =      EHCar[0*3 + 1];
  EH[0*3 + 2] = -1.0*EHCar[0*3 + 0];

  EH[1*3 + 0] =      EHCar[1*3 + 2];
  EH[1*3 + 1] =      EHCar[1*3 + 1];
  EH[1*3 + 2] = -1.0*EHCar[1*3 + 0];

}
