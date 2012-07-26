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
 * scuff-test-EEIs.h  -- header for scuff-test-EEI standalone test program
 *                       for computing edge-edge interactions 
 * 
 * homer reid         -- 11/2005 -- 11/2011
 */
#ifndef SCUFF_TEST_EEIS_H
#define SCUFF_TEST_EEIS_H

#include "libscuff.h" 
#include "libscuffInternals.h" 

using namespace scuff;

/******************************************************************/
/* routine for brute-force evaluation of panel-panel interactions */
/******************************************************************/
void GetPPIs_BruteForce(GetPPIArgStruct *Args, int PlotFits);

#endif
