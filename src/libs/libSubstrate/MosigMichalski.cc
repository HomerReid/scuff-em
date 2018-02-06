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
 * MosigMichalski.cc
 *
 * homer reid  -- 1/2018
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

#include <libhrutil.h>
#include <libhmat.h>
#include <libSGJC.h>
#include <vector>

#define MAXDIM 10

typedef std::vector<int> ivector;
typedef std::vector<double> dvector;
typedef std::vector<cdouble> zvector;

#define II cdouble(0.0,1.0)

/***************************************************************/
/* compute Omega_k, i.e. the remainder-ratio estimate used at  */
/*  step k of the Mosig-Michalski iteration.                   */
/* If OmegaVector is nonzero, we assume it has been prefilled  */
/* with remainder values, so just return OmegaVector[k].       */
/* Otherwise we use the 't', 'u', or 'v' estimators defined in */
/* the Michalski-Mosig paper, as selected by the tuv           */
/***************************************************************/
cdouble GetOmega(int k, cdouble *uValues, int Stride, int Offset,
                 double *xValues, char tuv, cdouble *OmegaValues)
{ 
  if (OmegaValues)
   return OmegaValues[k*Stride + Offset];

  cdouble uk = uValues[k*Stride + Offset], ukm1=uValues[(k-1)*Stride + Offset];
  double xk  = xValues[k], xkm1 = xValues[k-1];
  switch(tuv)
   { case 't': return uk/ukm1;
     case 'v': return ukm1*uk/(ukm1-uk);
     case 'u':
     default:  return xk*uk/(xkm1*ukm1);
   };
}

/***************************************************************/
/* Carry out the Kth-order Mosig-Michalski iteration to compute*/
/* S_0^(K) for a single series function, i.e. a single         */
/* integrand component.                                        */
/*                                                             */
/* uValues[ k*Stride + Offset ] = u[k]                         */
/* xValues[ k ]                 = x[k]                         */
/*                                                             */
/* uValues must be populated up to slot k=K                    */
/* xValues must be populated up to slot k=K+1                  */
/*                                                             */
/* If Workspace is non-null it must point to a buffer of       */
/* length at least K^2 (this could be reduced to N with more   */
/* complicated coding).                                        */
/***************************************************************/
cdouble MosigMichalski(int K, double *xValues,
                       cdouble *uValues, int Stride, int Offset,
                       double Mu,
                       cdouble *Workspace, cdouble *OmegaValues, char tuv)
{
  bool OwnsWorkspace = (Workspace==0);
  if (OwnsWorkspace)
   Workspace = new cdouble[K*K];

  HMatrix SMatrix(K,K,LHM_COMPLEX,Workspace);
  SMatrix.SetEntry(0,0,uValues[0*Stride + Offset]);

  // outer loop 
  for(int k=1; k<K; k++)
   {
     SMatrix.SetEntry(0, k, SMatrix.GetEntry(0,k-1) + uValues[k*Stride+Offset]);

     // inner loop to compute the $k$th counterdiagonal of the S matrix
     for(int m=1, n=k-1; m<=k; m++, n--)
      { double xp=xValues[n+1], x=xValues[n];
        cdouble Omega = GetOmega(k, uValues, Stride, Offset, xValues, tuv, OmegaValues);
        cdouble Eta = Omega/(1.0 + Mu*(m-1)*(xp-x)/x);
        cdouble Sp=SMatrix.GetEntry(m-1, n+1);
        cdouble Sm=SMatrix.GetEntry(m-1, n);
        SMatrix.SetEntry(m, n, (Sp - Eta*Sm)/(1.0-Eta));
      };

   };
  cdouble RetVal = SMatrix.GetEntry(K-1,0);
  if (OwnsWorkspace) delete Workspace;
  return RetVal;
}
