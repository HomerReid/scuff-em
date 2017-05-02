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
 * StaticSubstrate.h -- header file for implicit handling of 
 *                      layered dielectric substrates in scuff-static
 *
 * homer reid        -- 3/2017
 */

#ifndef STATICSUBSTRATE_H
#define STATICSUBSTRATE_H

#include <stdio.h>
#include <math.h>
#include <stdarg.h>

#include <libhrutil.h>
#include <libhmat.h>
#include <libMatProp.h>

#define MAXLAYER 10

/***************************************************************/
/* data structure for layered dielectric substrate             */
/***************************************************************/
typedef struct SubstrateData
 { int NumLayers;
   MatProp *MPMedium;
   cdouble EpsMedium;
   MatProp **MPLayer;
   cdouble  *EpsLayer;
   double *zLayer;
   double zGP;
   int qMaxEval;
   double qAbsTol;
   double qRelTol;
   int PPIOrder;
   int PhiEOrder;
   //int PPIRelTol;
   //int PPIAbsTol;
   int WhichIntegral;
 } SubstrateData;

/***************************************************************/
/***************************************************************/
/***************************************************************/
SubstrateData *CreateSubstrateData(const char *FileName, char **pErrMsg);
void DestroySubstrateData(SubstrateData *SD);

void GetStaticSubstrateGFCorrection(SubstrateData *SD,
                                    double XD[3], double XS[3],
                                    double PhiE[4],
                                    bool RetainSameLayerContributions=true);

#endif // STATICSUBSTRATE_H
