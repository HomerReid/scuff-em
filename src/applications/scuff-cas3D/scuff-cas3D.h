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
 * scuff-cas3D.h  -- header file for scuff-cas3D.cc
 *
 * homer reid  -- 2/2012
 */
#ifndef SCUFFCAS3D_H
#define SCUFFCAS3D_H

#include <libhrutil.h>
#include <libhmat.h>
#include <libscuff.h>

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
// SC3Data ('scuff-cas3D data') is a structure that contains all 
// information needed to compute the contribution of a single
// imaginary frequency to the Casimir quantities.
typedef struct SC3Data
 {
   RWGGeometry *G;
   char *ByXiFile;

   int N, N1;
   HMatrix **TBlocks, **UBlocks, **dUBlocks, *M, *dM;

   GTComplex **GTCList;
   int NumTransformations;

   char *WriteCache;
   int nThread;

 } SC3Data;

SC3Data *CreateSC3Data(char *GeoFile, char *TransFile, 
                       char *ByOmegaFile, int nThread);

void GetFrequencyIntegrand(SC3Data *SC3D, double Xi, double *FI);

#endif // #define SCUFFCAS3D_H
