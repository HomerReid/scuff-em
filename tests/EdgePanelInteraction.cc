/* Copyright (C) 2005-2011 M. T. Homer Reid *
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
 * EdgePanelInteractions.cc -- libscuff routines for evaluating the interactions
 *                          -- between a single RWG panel and the line charge
 *                          -- associated with a half-RWG basis function     
 * 
 * homer reid               -- 6/2012
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "libhrutil.h"
#include "libTriInt.h"

#include "libscuff.h"
#include "libscuffInternals.h"

namespace scuff {

#define II cdouble(0.0,1.0)

/*--------------------------------------------------------------*/
/*- PART 1: Routine to compute edge-panel interaction using     */
/*-         fixed-order numerical cubature.                     */
/*-                                                             */
/*- V[0][0..2] = cartesian coordinates of panel vertex 1        */
/*- V[1][0..2] = cartesian coordinates of panel vertex 2        */
/*- V[2][0..2] = cartesian coordinates of panel vertex 3        */
/*- L1[0..2]   = cartesian coordinates of edge vertex 1         */
/*- L2[0..2]   = cartesian coordinates of edge vertex 2         */
/*--------------------------------------------------------------*/
cdouble GetEPI_Cubature(double **V, double **L, cdouble K)
                      
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double *V0, A[3], B[3], *V0P, DL[3];
  V0=V[0];
  VecSub(V[1], V[0], A);
  VecSub(V[2], V[0], B);
  V0P=L[0];
  VecSub(L[1],L[0],DL);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double *TCR;  // triangle cubature rule 
  int NumTCPts;
  if (HighOrder)
   TCR=GetTCR(20, &NumTCPts);
  else
   TCR=GetTCR(4, &NumTCPts);

  double *QR;   // quadrature rule 
  int NumQPts;

  /***************************************************************/
  /* outer loop (cubature over triangle) *************************/
  /***************************************************************/
  int np, ncp, npp, ncpp;
  int Mu, Nu;
  double u, v, w, up, wp;
  double X[3], XP[3], R[3];
  double r;
  cdouble Phi;
  cdouble GOuter, GInner;
  GOuter=0.0;
  for(np=ncp=0; np<NumTCPts; np++) 
   { 
     u=TCR[ncp++]; v=TCR[ncp++]; w=TCR[ncp++];

     /***************************************************************/
     /* set X and F=X-Q *********************************************/
     /***************************************************************/
     for(Mu=0; Mu<3; Mu++)
      X[Mu] = V0[Mu] + u*A[Mu] + v*B[Mu];

     /***************************************************************/
     /* inner loop (quadrature over edge)                           */
     /***************************************************************/
     GInner=0.0;
     for(npp=ncpp=0; npp<NumQPts; npp++)
      { 
        up=QR[ncpp++]; wp=QR[ncpp++];

        /***************************************************************/ 
        /* set XP and FP=XP-QP *****************************************/
        /***************************************************************/
        for(Mu=0; Mu<3; Mu++)
         { XP[Mu] = V0P[Mu] + up*DL[Mu];
           R[Mu] = X[Mu] - XP[Mu];
         };
      
        /***************************************************************/
        /* inner integrand  ********************************************/
        /***************************************************************/
        r=VecNorm(R);
        Phi = exp(II*k*r)/(4.0*M_PI*r);
        GInner += wp*Phi;

      }; /* for(npp=ncpp=0; npp<NumPts; npp++) */

     /*--------------------------------------------------------------*/
     /*- accumulate contributions to outer integral                  */
     /*--------------------------------------------------------------*/
     GOuter+=w*HInner[0];

   }; // for(np=ncp=0; np<nPts; np++) 

  return GOuter;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct EPITDData
{ 
  double *A, *B, *L;
  cdouble K;

} EPITDData;

void EPITDIntegrand(unsigned ndim, const double *xy, void *params,
		    unsigned fdim, double *fval);
{
  (void) ndim;
  (void) fdim;

  EPITDData *EPITDD = (EPITDData *)params;

  double *A = EPITDD->A;
  double *B = EPITDD->B;
  double *L = EPITDD->L;
  cdouble IK = II*EPITDD->K;

  double x=xy[0];
  double y=xy[1];

  int p;
  double RR1[3], RR2[3], RR3[3];
  for(p=0; p<3; p++)
   { RR1[p] = x*A[p] + x*y*B[p]       + L[p];
     RR2[p] =   A[p] + (1.0-x)*y*B[p] + x*L[p];
     RR3[p] =   A[p] + (1.0-x*y)*B[p] + x*L[p];
   };
  
  double R1 = sqrt( RR1[0]*RR1[0] + RR1[1]*RR1[1] + RR1[2]*RR1[2] );
  double R2 = sqrt( RR2[0]*RR2[0] + RR2[1]*RR2[1] + RR2[2]*RR2[2] );
  double R3 = sqrt( RR3[0]*RR3[0] + RR3[1]*RR3[1] + RR3[2]*RR3[2] );
  
  cdouble *zfval=(cdouble *)fval;

  zfval[0] =        x*exp(IK*R1)*ExpRel(2,-IK*R1) / (8.0*M_PI*R1)
             +(1.0-x)*exp(IK*R2)*ExpRel(2,-IK*R2) / (8.0*M_PI*R2)
                   +x*exp(IK*R3)*ExpRel(2,-IK*R3) / (8.0*M_PI*R3);
  
}


/***************************************************************/
/***************************************************************/
/***************************************************************/
cdouble GetEPI_TaylorDuffy(double *V0, double *A, double *B, double *DL, cdouble K)
{
  
  EPITDData;
 
  double Lower[2] = {0.0, 0.0 };
  double Upper[2] = {1.0, 1.0 };
  adapt_integrate(2, EPITDIntegrand, (void *)Data, 2, Lower, Upper,
		  MAXEVALS, ABSTOL, RELTOL, (double *)&G, Error);

  return G;

}

/***************************************************************/
/* PV[0][0..2] = cartesian coordinates of panel vertex 1       */
/* PV[1][0..2] = cartesian coordinates of panel vertex 2       */
/* PV[2][0..2] = cartesian coordinates of panel vertex 3       */
/*                                                             */
/* EVL[0][0..2] = cartesian coordinates of edge vertex 1       */
/* EV[1][0..2] = cartesian coordinates of edge vertex 2        */
/*                                                             */
/* Note:  What this routine computes is simply the unadorned   */
/*        integral                                             */
/*                                                             */
/*        \int_{panel} dx \int_{edge} dy G(x-y)                */
/*                                                             */
/*       where G(r)=e^{ikr}/(4*pi*r).                          */
/*                                                             */
/*       This means that the calling routine must multiply the */
/*       return value by a prefactor of                        */
/*                                                             */
/*         -l_\alpha l_\beta/K^2                               */
/*                                                             */
/*       to get the value of the <f|G|h> inner product.        */
/***************************************************************/
cdouble GetEdgePanelInteraction(double **PV, double **EV, cdouble K)
{ 
  
  int EIEV[2]; // 'edge index of equal vertex' 
  int PIEV[2]; // 'panel index of equal vertex'
  int i, j, ncv=0; 
  for (i=0; i<2; i++)
   for (j=0; j<3; j++)
    if ( VecEqualFloat(L[i], V[j]) )
     { if (ncv==2) ErrExit("%s:%i: internal error \n",__FILE__,__LINE__);
       EIEV[ncv]=i;
       PIEV[ncv]=j;
       ncv++;
     };

  if ( ncv == 0 )
   return GetEPI_Cubature(PV, EV, K);

  double A[3], B[3], DL[3];

  int i=PIEV[0], ip1=(i+1)%3, ip2=(i+2)%3;
  int j=EIEV[0], jp1=(j+1)%2;
  VecSub(PV[ip1], PV[i],   A);
  VecSub(PV[ip2], PV[ip1], B);
  VecSub(EV[jp1], EV[j],   DL);
  
  return GetEPI_TaylorDuffy(PV[0], A, B, DL, K);
  
}
