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
 * scuff-heat  -- a standalone code within the scuff-em suite
 *             -- for solving heat-transfer and radiation problems
 *
 * homer reid  -- 2/2012
 */
#ifndef SCUFFHEAT_H
#define SCUFFHEAT_H

#include <libhrutil.h>
#include <libhmat.h>
#include <libscuff.h>

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
// SHData ('scuff-heat data') is a structure that contains all 
// information needed to compute the heat radiation or transfer
// at a given frequency.
typedef struct SHData
 {
   RWGGeometry *G;
   char *ByOmegaFile;

   int N1, N2;
   HMatrix **TSelf, **TMedium, **UMedium;
   HMatrix *SymG1, *SymG2;
   HMatrix *W, *W21, *W21SymG1, *W21DSymG2;
   HMatrix *Scratch;

   HVector *DV;
   int PlotFlux;

   GTComplex **GTCList;
   int NumTransformations;

   char *WriteCache;
   int nThread;

 } SHData;

SHData *CreateSHData(char *GeoFile, char *TransFile, int PlotFlux,
                     char *ByOmegaFile, int nThread);

void GetFrequencyIntegrand(SHData *SHD, cdouble Omega, double *FI);

#endif
