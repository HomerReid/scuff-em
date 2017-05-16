
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
 * scuff-static.h -- a standalone code within the scuff-em suite
 *                -- for solving electrostatics problems
 *
 * homer reid     -- 6/2011--2/2017
 */
#ifndef SCUFFSTATIC_H
#define SCUFFSTATIC_H

#include <libhrutil.h>
#include <libhmat.h>
#include "libIncField.h"
#include <libscuff.h>
#include <SSSolver.h>
#include <cmatheval.h>

using namespace scuff;

/***************************************************************/
/* routines in OutputModules.cc ********************************/
/***************************************************************/
void WritePolarizabilities(SSSolver *SSS, HMatrix *M,
                           HVector *Sigma, char *FileName);

void WriteCapacitanceMatrix(SSSolver *SSS, HMatrix *M,
                            HVector *Sigma, char *CapFile);

void WriteCMatrix(SSSolver *SSS, HMatrix *M,
                  HVector *Sigma, int lMax,
                  char *TextFileName, char *HDF5FileName);

void Solve(SSSolver *SSS, HMatrix *M, HVector *Sigma,
           char *PotFile, char *PhiExt, int ConstFieldDirection);

void WriteFields(SSSolver *SSS, HVector *Sigma,
                 char *PhiExt, int ConstFieldDirection,
                 char **EPFiles, int nEPFiles);

#endif
