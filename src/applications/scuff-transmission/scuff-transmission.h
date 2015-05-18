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
 * scuff-transmission.cc -- a command-line tool to compute transmission
 *                       -- through thin films and metamaterial arrays
 *
 * homer reid  -- 5/2012
 */
#ifndef SCUFFTRANSMISSION_H
#define SCUFFTRANSMISSION_H

#include <libhrutil.h>
#include <libhmat.h>
#include <libscuff.h>

using namespace scuff;

// incident-field polarizations 
#define POL_TE  0
#define POL_TM  1
#define NUMPOLS 2

// source/destination regions for plane waves 
#define REGION_UPPER 0
#define REGION_LOWER 1
#define NUMREGIONS   2

// in GetFlux.cc
void GetFlux(RWGGeometry *G, IncField *IF, HVector *KN,
             cdouble Omega, double *kBloch, int NQPoints,
             double ZAbove, double ZBelow, double *Flux);

// in GetAmplitudes.cc
void GetPlaneWaveAmplitudes(RWGGeometry *G, HVector *KN,
                            cdouble Omega, double *kBloch,
                            int WhichRegion, bool IsUpper,
                            cdouble TETM[2], bool WriteByKNFile=false);


#endif // SCUFFTRANSMISSION_H
