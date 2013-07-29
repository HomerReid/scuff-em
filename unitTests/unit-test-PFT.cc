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
 * unit-test-PFT.cc -- SCUFF-EM unit tests for PFT calculations
 * 
 * homer reid       -- 11/2005 -- 10/2011
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>
#include "libscuff.h"
#include "libscuffInternals.h"

using namespace scuff;

#define II cdouble (0.0,1.0)

/***************************************************************/
/* tables of power and force data for spheres of radius R=1um  */
/* and various material composition.                           */
/*                                                             */
/* in each case, powers and forces are for the case of         */
/* illumination by a z-traveling plane wave with E-field       */
/* polarized in the +x direction.                              */
/*                                                             */
/* data table entries:                                         */
/*  (1) omega (angular frequency, units 3e14 rad/sec)          */
/*  (2) efficiency for power absorption                        */
/*  (3) efficiency for power scattering                        */
/*  (4) efficiency for z-force                                 */
/*                                                             */
/* where efficiencies are cross-sections divided by \pi R^2.   */
/*                                                             */
/* cross-sections for {absorbed, scattered} power are the      */
/* total power {absorbed by, scattered from} the object        */
/* divided by incident power flux. the cross-section for the   */
/* force is the z-force on the object divided by the incident  */
/* momentum flux.                                              */
/***************************************************************/
double GoldPFT[4]
 =

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{ 
  SetLogFileName("scuff-unit-tests.log");
  G->SetLogLevel(SCUFF_VERBOSELOGGING);
  Log("SCUFF-EM PFT unit tests running on %s",GetHostName());

  int FailedCases=0;

  /***************************************************************/ 
  /***************************************************************/ 
  /***************************************************************/ 
  RWGGeometry *G = new RWGGeometry("PECSphere.scuffgeo");
  for(Omega=0.001; Omega<=10.0; Omega*=sqrt(10.0)
  
  /***************************************************************/ 
  /***************************************************************/ 
  /***************************************************************/ 
  if (FailedCases>0)
   abort();

  printf("All tests successfully passed.\n");
  return 0;

}
