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
 * CylindricalWave.h -- definition for the CylindricalWave implementation
 *                 -- of the IncField class 
 *
 * Homer Reid      -- 7/2012
 */

#ifndef CYLINDRICALWAVE_H
#define CYLINDRICALWAVE_H 

#include <libhrutil.h>
#include <libIncField.h>

// 'tranverse-electric to z' or 'transverse-magnetic to z'
// TE2Z means the E-field has no z-component
// TM2Z means the H-field has no z-component
#define CW_TE2Z     0
#define CW_TM2Z     1

#define CW_REGULAR  0
#define CW_OUTGOING 1
#define CW_INCOMING 2

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
class CylindricalWave : public IncField
 { 
 public:
   int Nu;           // spherical wave indices
   double kz;        // z wavenumber
   int Polarization; // CW_TE2Z or CW_TM2Z
   int WaveType;     // CW_REGULAR, CW_OUTGOING, CW_INCOMING

   CylindricalWave(int Nu=0, double kz=0.0, int Polarization=CW_TE2Z,
                   int WaveType=CW_INCOMING);

   void SetNu(int Nu);
   void SetKz(double kz);
   void SetPolarization(int Polarization);
   void SetWaveType(int Type);

   void GetFields(const double X[3], cdouble EH[6]);

 };

#endif // #ifndef CYLINDRICALWAVE_H
