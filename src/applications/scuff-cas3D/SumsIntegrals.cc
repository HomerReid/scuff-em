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
 * SumsIntegrals.cc -- routines for evaluating matsubara sums,
 *                  -- imaginary-frequency integrals, and 
 *                  -- brillouin-zone integrations
 *
 * homer reid       -- 2/2012
 *
 */

#include "scuff-cas3D.h"
#include "libscuffInternals.h"

/***************************************************************/
/* Evaluate the contribution of a single Xi, kBloch point to   */
/* the Casimir quantities.                                     */
/***************************************************************/
void GetCasimirIntegrand(SC3Data *SC3D, double Xi, double *kBloch, double *EFT)
{
}

/***************************************************************/
/* Integrate over the positive imaginary frequency axis to get */
/* the total Casimir quantities at zero temperature.           */
/***************************************************************/
void GetXiIntegral(SC3Data *SC3D, double *EFT)
{
}

/***************************************************************/
/* Evaluate the Matsubara sum to get the total Casimir         */
/* quantities at temperature Temperature degrees Kelvin.       */
/***************************************************************/
void GetMatsubaraSum(SC3Data *SC3D, double Temperature, double *EFT)
{
}
