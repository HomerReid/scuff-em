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
 *
 * libSGJC.h -- a couple of quick wrappers around the SGJ cubature 
 *           -- routines, provided for backward compatibility
 */

#include <stdlib.h>
#include "libSGJC.h"

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct NewIntegrandData 
 { oldintegrand f;
   void *fdata;
 } NewIntegrandData;

int MyNewIntegrand(unsigned ndim, const double *x, void *params,
                   unsigned fdim, double *fval)
{
  NewIntegrandData *NID=(NewIntegrandData *)params;
  NID->f(ndim, x, NID->fdata, fdim, fval);
  return 0;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int adapt_integrate(unsigned fdim, oldintegrand f, void *fdata,
		    unsigned dim, const double *xmin, const double *xmax, 
		    unsigned maxEval, double reqAbsError, double reqRelError, 
		    double *val, double *err)
{
  NewIntegrandData MyNID, *NID=&MyNID;

  MyNID.f=f;
  MyNID.fdata=fdata;

  return hcubature(fdim, MyNewIntegrand, (void *)NID, dim, xmin, xmax, 
                   maxEval, reqAbsError, reqRelError, ERROR_INDIVIDUAL, 
                   val, err);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int adapt_integrate_log(unsigned fdim, oldintegrand f, void *fdata, 
                        unsigned dim, const double *xmin, const double *xmax, 
                        unsigned maxEval, double reqAbsError, double reqRelError, 
                        double *val, double *err, const char *LogFileName, int LogInterval)
{
  (void )LogFileName; // unused
  (void )LogInterval; // unused

  return adapt_integrate(fdim, f, fdata, dim, xmin, xmax, 
                         maxEval, reqAbsError, reqRelError, val, err);
}
