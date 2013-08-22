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
 * PolModel.cc   -- implementation of a simple class for describing
 *               -- the frequency-dependent polarizability of atoms
 *               -- and molecules
 *
 * homer reid    -- 1/2012
 *
 */

#include "scuff-caspol.h"
#include <libhmat.h>

/***************************************************************/
/* PolModel class constructor                                  */
/***************************************************************/
PolModel::PolModel(const char *Atom)
{
  ErrMsg=0;

  Name=strdup(Atom);

  NumPoints = NUMPOINTS_DPB;

  XiPoints = XiPoints_DPB;

  if ( !strcasecmp(Atom,"H")  || !strcasecmp(Atom,"Hydrogen") )
   PolPoints = PolPoints_DPB_Hydrogen;
  else if ( !strcasecmp(Atom,"Li") || !strcasecmp(Atom,"Lithium") )
   PolPoints = PolPoints_DPB_Lithium;
  else if ( !strcasecmp(Atom,"Na") || !strcasecmp(Atom,"Sodium") )
   PolPoints = PolPoints_DPB_Sodium;
  else if ( !strcasecmp(Atom,"K")  || !strcasecmp(Atom,"Potassium") )
   PolPoints = PolPoints_DPB_Potassium;
  else if ( !strcasecmp(Atom,"Rb") || !strcasecmp(Atom,"Rubidium") )
   PolPoints = PolPoints_DPB_Rubidium;
  else if ( !strcasecmp(Atom,"Cs") || !strcasecmp(Atom,"Cesium") )
   PolPoints = PolPoints_DPB_Cesium;
  else if ( !strcasecmp(Atom,"Fr") || !strcasecmp(Atom,"Francium") )
   PolPoints = PolPoints_DPB_Francium;
  else
   { ErrMsg=strdup("unknown atom type %s",Atom);
     return;
   };

   // initialize the interpolator
   PolInterp = new Interp1D(XiPoints, PolPoints, NumPoints, 1);

}

/***************************************************************/
/* routine to compute polarizability tensor at a given         */
/* frequency.                                                  */
/*                                                             */
/* on input, Xi is the imaginary frequency in units of         */
/* 3e14 rad/sec.                                               */
/*                                                             */
/* Alpha points to a preallocated 3x3 real-valued HMatrix.     */
/*                                                             */
/* On return, Alpha is filled in with the components of the    */
/* polarizability at this imaginary frequency.                 */
/***************************************************************/
void PolModel::GetPolarizability(double Xi, HMatrix *Alpha)
 {
   // the constant here convertes Xi from my units, 
   // in which '1' == 3e14 rad/sec, to atomic units, 
   // in which '1' == 2.598e+17 rad/sec
   double AlphaDiag;
   PolInterp->Evaluate(0.00115493*Xi , &AlphaDiag);

   Alpha->Zero();
   Alpha->SetEntry(0,0,AlphaDiag);
   Alpha->SetEntry(1,1,AlphaDiag);
   Alpha->SetEntry(2,2,AlphaDiag);

 }

/***************************************************************/
/* Note: DPB stands for "Derevianko, Porsev, Babb."            */
/***************************************************************/
#define NUMPOINTS_DPB 50;
double XiPoints_DPB[NUMPOINTS_DPB] =
{0.00178065, 0.00937463, 0.0230068, 0.042631, 0.0681817,
 0.0995817, 0.136748, 0.179602, 0.228071, 0.282107,
 0.341686, 0.406823, 0.477579, 0.554072, 0.636488,
 0.725091, 0.820235, 0.922382, 1.03212, 1.15017,
 1.27743, 1.41501, 1.56425, 1.72679, 1.90461,
 2.10016, 2.31644, 2.55713, 2.82683, 3.13128,
 3.47776, 3.87553, 4.3366, 4.87665, 5.51655,
 6.28448, 7.21927, 8.37558, 9.83228, 11.7066,
 14.179, 17.5384, 22.2715, 29.2508, 40.168,
 58.6667, 93.8285, 173.862, 426.684, 2246.37
};

double wPoints_DPB[NUMPOINTS_DPB] =
{0.00456886, 0.0106185, 0.0166378, 0.0225996, 0.0284889, 
 0.0342971, 0.0400227, 0.0456718, 0.0512592, 0.0568084, 
 0.0623519, 0.0679318, 0.0735998, 0.0794176, 0.085458, 
 0.0918056, 0.0985589, 0.105832, 0.113758, 0.122493, 
 0.132222, 0.143164, 0.155587, 0.169813, 0.18624, 
 0.205361, 0.227799, 0.254343, 0.286005, 0.324106, 
 0.370381, 0.427154, 0.497572, 0.585976, 0.698465, 
 0.843785, 1.03477, 1.29076, 1.64181, 2.13626, 
 2.85525, 3.94176, 5.66354, 8.56095, 13.8343, 
 24.5132, 49.7405, 125.732, 483.298, 5763.83
};

double PolPoints_DPB_Hydrogen[NUMPOINTS_DPB] = 
{
};
