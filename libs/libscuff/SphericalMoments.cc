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
 * SphericalMoments.cc -- libscuff code for computing spherical
 *                     -- multipole moments
 * 
 * homer reid          -- 11/2005 -- 12/2011
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>
#include <libSpherical.h>
#include <libSGJC.h>

#include "libscuff.h"

namespace scuff {

/***************************************************************/
/* integrand routine passed to adaptive integrator *************/
/***************************************************************/
typedef struct GEMNMIData
 {
   double *V0, A[3], B[3], *Q;
   double *R0;
   cdouble K;
   int lMax;
 } GEMNMIData;

void GetEdgeMNMomentsIntegrand(unsigned ndim, const double *x, void *params,
			       unsigned fdim, double *fval)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  GEMNMIData *D=(GEMNMIData *)params;
  double *V0 = D->V0;
  double *A  = D->A;
  double *B  = D->B;
  double *Q  = D->Q;
  double *R0 = D->R0;
  cdouble K  = D->K;
  int lMax   = D->lMax;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double u=x[0];
  double v=u*x[1];
  double R[3], F[3], RmR0[3];
  int Mu;
  for(Mu=0; Mu<3; Mu++)
   { R[Mu] = V0[Mu] + u*A[Mu] + v*B[Mu]; 
     F[Mu] = R[Mu] - Q[Mu];
     RmR0[Mu] = R[Mu] - R0[Mu];
   };

  double r, Theta, Phi;
  double FS[3]; 
  CoordinateC2S(RmR0,&r,&Theta,&Phi);
  VectorC2S(Theta, Phi, F, FS);
  
  /* the jacobian of the coordinate transformation is u, which */
  /* i put into FS since it is a factor in all integrand components */
  FS[0]*=u;
  FS[1]*=u;
  FS[2]*=u;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int Alpha, NAlpha=(lMax+1)*(lMax+1);
  cdouble M[3*NAlpha], N[3*NAlpha];
  GetMNlmArray(lMax, K, r, Theta, Phi, LS_INTERIOR, M, N);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble *zf=(cdouble *)fval;
  for(Alpha=0; Alpha<NAlpha; Alpha++)
   { zf[2*Alpha+0] = M[3*Alpha+0]*FS[0] + M[3*Alpha+1]*FS[1] + M[3*Alpha+2]*FS[2];
     zf[2*Alpha+1] = N[3*Alpha+0]*FS[0] + N[3*Alpha+1]*FS[1] + N[3*Alpha+2]*FS[2];
   };
  
}


/***************************************************************/
/* Compute the spherical multipole moments of the electric     */
/* current distribution described by a single RWG basis        */
/* function.                                                   */
/*                                                             */
/* R0 is the origin about which moments are computed. R0 may   */
/* be set to NULL, in which case the origin is taken to be the */
/* centroid of the edge (defined as the midpoint of the common */
/* edge shared by the panel pair).                             */
/*                                                             */
/* AM and AN must point to buffers with enough room to store   */
/* (lMax+1)*(lMax+1) cdoubles.                                 */
/***************************************************************/
void GetEdgeMNMoments(RWGObject *O, int ne, cdouble K, int lMax, 
                      double *R0, cdouble *AM, cdouble *AN)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  RWGEdge *E = O->Edges[ne];
  double *QP = O->Vertices + 3*E->iQP;
  double *V1 = O->Vertices + 3*E->iV1;
  double *V2 = O->Vertices + 3*E->iV2;
  double *QM = E->iQM==-1 ? 0 : O->Vertices + 3*E->iQM;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  GEMNMIData MyGEMNMIData, *D=&MyGEMNMIData;
  D->K    = K;
  D->lMax = lMax;
  D->R0   = R0 ? R0 : E->Centroid;

  /***************************************************************/ 
  /* (lMax+1)*(lMax+1) multipoles, times 2 (M,N multipoles),     */
  /* times 2 (cdouble-valued integrand vector)                   */
  /***************************************************************/
  int fdim=2*2*(lMax+1)*(lMax+1);
  new double FP[fdim]; 
  new double FM[fdim]; 
  new double Error[fdim];

  double Lower[2]={0.0, 0.0};
  double Upper[2]={1.0, 1.0};
   
  /***************************************************************/
  /* contribution of positive panel ******************************/
  /***************************************************************/
  D->V0 = QP;
  VecSub(V1, QP, D->A);
  VecSub(V2, V1, D->B);
  D->Q  = QP;
  
  adapt_integrate(fdim, GetEdgeMNMomentsIntegrand, (void *)D, 2,
		  Lower, Upper, 0, ABSTOL, RELTOL, FP, Error);
   

  /***************************************************************/
  /* contribution of negative panel if present *******************/
  /***************************************************************/
  if (QM)
   { D->V0 = QM;
     VecSub(V1, QM, D->A);
     VecSub(V2, V1, D->B);
     D->Q  = QM;
     adapt_integrate(fdim, GetEdgeMNMomentsIntegrand, (void *)D, 2,
                     Lower, Upper, 0, ABSTOL, RELTOL, FM, Error);
   }
  else
   memset(FM, fdim, 0*sizeof(double)); 

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double PreFac=E->Length;
  int Alpha, NAlpha=(lMax+1)*(lMax+1);
  for (Alpha=0; Alpha<NAlpha; Alpha++)
   { AM[Alpha] = PreFac*(FP[2*Alpha+0] - FM[2*Alpha+0]);
     AN[Alpha] = PreFac*(FP[2*Alpha+1] - FM[2*Alpha+1]);
   };

  delete[] FP;
  delete[] FM;
  delete[] Error;

}

} // namespace scuff
