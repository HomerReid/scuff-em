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
 * GetDyadicGF.cc  -- compute the scattering parts of the electric and
 *                    magnetic dyadic Green's functions
 *
 * homer reid      -- 5/2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <ctype.h>

#include <libhrutil.h>

#include "libscuff.h"
#include "libscuffInternals.h"

#define II cdouble(0.0,1.0)

namespace scuff {

/***************************************************************/
/* This routine compute the scattering parts of the electric   */
/* and magnetic dyadic green's functions (DGFs) at a point X.  */
/*                                                             */
/* The calculation proceeds by (a) placing a point dipole      */
/* source at X, (b) solving a scattering problem with the      */
/* incident field emanating from this point source, and then   */
/* (c) evaluating the scattered fields back at X.              */
/*                                                             */
/* If the original point dipole source was an electric dipole  */
/* pointing in the i direction, then the scattered E-field at  */
/* X gives us the ith column of the electric DGF.              */
/*                                                             */
/* If the original point dipole source was a *magnetic* dipole */
/* pointing in the i direction, then the scattered H-field at  */
/* X gives us the ith column of the magnetic DGF.              */
/*                                                             */
/* Inputs:                                                     */
/*                                                             */
/*  XEval:  cartesian coordinates of evaluation point          */
/*  XSource: cartesian coordinates of source points            */
/*                                                             */
/*  Omega:  angular frequency                                  */
/*                                                             */
/*  kBloch: 1D or 2D Bloch vector; set to NULL for             */
/*          non-periodic geometries                            */
/*                                                             */
/*  M:     the LU-factorized BEM matrix. Note that the caller  */
/*         is responsible for calling AssembleBEMMatrix() and  */
/*         LUFactorize() before calling this routine.          */
/*                                                             */
/*  KN:    an HVector used internally within this routine that */
/*         that must have been obtained from a previous call   */
/*         to AllocateRHSVector().                             */
/*                                                             */
/* Outputs:                                                    */
/*                                                             */
/*  On return, GE[i][j] and GM[i][j] are respectively the      */
/*  i,j components of the electric and magnetic DGF tensors,   */
/*  with 'Scat' and 'Tot' indicating the scattering parts      */
/*  and the full DGFs.                                         */
/*                                                             */
/***************************************************************/
void RWGGeometry::GetDyadicGFs(double XEval[3], double XSource[3],
                               cdouble Omega, double *kBloch,
                               HMatrix *M, HVector *KN,
                               cdouble GEScat[3][3], cdouble GMScat[3][3],
                               cdouble GETot[3][3], cdouble GMTot[3][3])
{
  if (M==0 || M->NR != TotalBFs || M->NC!=M->NR )
   ErrExit("%s:%i: invalid M matrix passed to GetDyadicGFs()",__FILE__,__LINE__);

  if (KN==0 || KN->N != TotalBFs )
   ErrExit("%s:%i: invalid K vector passed to GetDyadicGFs()",__FILE__,__LINE__);

  if ( (LDim>0 && kBloch==0) || (LDim==0 && kBloch!=0) )
   ErrExit("%s:%i: incorrect kBloch specification",__FILE__,__LINE__);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble P[3]={1.0, 0.0, 0.0};
  PointSource PS(XSource, P);

  cdouble EH[6];
  cdouble EpsRel, MuRel;
  int nr=GetRegionIndex(XSource);
  RegionMPs[nr]->GetEpsMu(Omega, &EpsRel, &MuRel);
  //cdouble IKZ = II*Omega*Mu*ZVAC;
  cdouble k2 = EpsRel*MuRel*Omega*Omega;
  cdouble Z2 = ZVAC*ZVAC*MuRel/EpsRel;

  double r=VecDistance(XSource, XEval);

  for(int Mu=0; Mu<3; Mu++)
   { 
     // set point source to point in the ith direction
     memset(P, 0, 3*sizeof(cdouble));
     P[Mu]=1.0;
     PS.SetP(P);

     // solve the scattering problem for an electric point source 
     PS.SetType(LIF_ELECTRIC_DIPOLE);
     AssembleRHSVector(Omega, kBloch, &PS, KN);
     M->LUSolve(KN);
     GetFields(0, KN, Omega, kBloch, XEval, EH);
     GEScat[0][Mu]=EH[0] / k2;
     GEScat[1][Mu]=EH[1] / k2;
     GEScat[2][Mu]=EH[2] / k2;
     
     // add the incident fields unless XSource=XEval
     if (r>0.0)
      { PS.GetFields(XEval, EH);
        GETot[0][Mu] = GEScat[0][Mu] + EH[0]/k2;
        GETot[1][Mu] = GEScat[1][Mu] + EH[1]/k2;
        GETot[2][Mu] = GEScat[2][Mu] + EH[2]/k2;
      };

     // solve the scattering problem for a magnetic point source 
     PS.SetType(LIF_MAGNETIC_DIPOLE);
     AssembleRHSVector(Omega, kBloch, &PS, KN);
     M->LUSolve(KN);
     GetFields(0, KN, Omega, kBloch, XEval, EH);
     GMScat[0][Mu]=EH[3] * Z2/k2;
     GMScat[1][Mu]=EH[4] * Z2/k2;
     GMScat[2][Mu]=EH[5] * Z2/k2;

     // add the incident fields unless XSource=XEval
     if (r>0.0)
      { PS.GetFields(XEval, EH);
        GMTot[0][Mu] = GMScat[0][Mu] + EH[3]*Z2/k2;
        GMTot[1][Mu] = GMScat[1][Mu] + EH[4]*Z2/k2;
        GMTot[2][Mu] = GMScat[2][Mu] + EH[5]*Z2/k2;
      };
   };
}

/***************************************************************/
/* get just the scattering part of the ``diagonal'' element of */
/* the DGFs (in which the source and destination points        */
/* coincide), working at a fixed frequency and Bloch vector.   */
/* Useful for LDOS computations and Casimir-Polder forces.     */
/***************************************************************/
void RWGGeometry::GetDyadicGFs(double X[3], cdouble Omega, double *kBloch,
                               HMatrix *M, HVector *KN,
                               cdouble GEScat[3][3], cdouble GMScat[3][3])
{
  cdouble Dummy[3][3];
  GetDyadicGFs(X, X, Omega, kBloch, M, KN, GEScat, GMScat, Dummy, Dummy);
}

/***************************************************************/
/* get the scattering part of the diagonal DGF elements,       */
/* working at a specific frequency but not a specific kBloch.  */
/* For non-periodic geometries this is just a single call to   */
/* GetDyadicGFs() above, but for periodic geometries it        */
/* involves an integral over the Brillouin zone.               */
/***************************************************************/
#if 0
void RWGGeometry::GetDyadicGFs(double X[3], cdouble Omega,
                               HMatrix *M, HVector *KN,
                               cdouble GEScat[3][3], cdouble GMScat[3][3])
{
  if (LDim==0) // non-periodic case
   { GetDyadicGFs(X, Omega, kBloch, M, KN, GEScat, GMScat);
     return;
   };
}
#endif

} // namespace scuff
