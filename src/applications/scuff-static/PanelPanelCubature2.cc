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
 * PanelCubature2.cc -- libscuff routines for computing numerical 
 *                   -- cubatures over panels, pairs of panels,
 *                   -- basis functions, or pairs of basis functions
 * 
 * homer reid        -- 11/2005 -- 5/2013
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>
#include <libTriInt.h>

#include "libscuff.h"
#include "libscuffInternals.h"
#include "libSGJC.h"

using namespace scuff;

typedef void (*PC2Function)(double X[3], void *UserData, double Weight, double *Integral);

typedef void (*PPC2Function)(double XA[3], double XB[3], void *UserData,
                             double Weight, double *Integral);

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetPanelCubature2(RWGSurface *S, int np, PC2Function Integrand,
                       void *UserData, int IDim, int Order,
                       double *Integral)
{
  int NumPts;
  double *TCR=GetTCR(Order,&NumPts);
  if (TCR==0) 
   ErrExit("invalid cubature order in GetPanelPanelCubature2");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGPanel *P  = S->Panels[np];
  double *Q    = S->Vertices + 3*(P->VI[0]);
  double *V1   = S->Vertices + 3*(P->VI[1]);
  double *V2   = S->Vertices + 3*(P->VI[2]);

  double L1[3], L2[3];
  VecSub(V1, Q, L1);
  VecSub(V2, Q, L2);
  double Area = P->Area;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  memset(Integral, 0, IDim*sizeof(double));
  for(int n=0, ncp=0; n<NumPts; n++)
   { 
     double u=TCR[ncp++], v=TCR[ncp++], w=TCR[ncp++], X[3];
     for(int Mu=0; Mu<3; Mu++)
      X[Mu] = Q[Mu] + u*L1[Mu] + v*L2[Mu];

     Integrand(X,UserData,2.0*Area*w,Integral);
   };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetPanelPanelCubature2(RWGSurface *SA, int npA,
                            RWGSurface *SB, int npB,
                            PPC2Function Integrand,
                            void *UserData, int IDim,
                            int Order, double *Integral)
{
  int NumPts;
  double *TCR=GetTCR(Order,&NumPts);
  if (TCR==0) 
   ErrExit("invalid cubature order in GetPanelPanelCubature2");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGPanel *PA  = SA->Panels[npA];
  double *QA    = SA->Vertices + 3*(PA->VI[0]);
  double *V1A   = SA->Vertices + 3*(PA->VI[1]);
  double *V2A   = SA->Vertices + 3*(PA->VI[2]);

  double L1A[3], L2A[3];
  VecSub(V1A, QA, L1A);
  VecSub(V2A, QA, L2A);
  double AreaA = PA->Area;

  RWGPanel *PB  = SB->Panels[npB];
  double *QB    = SB->Vertices + 3*(PB->VI[0]);
  double *V1B   = SB->Vertices + 3*(PB->VI[1]);
  double *V2B   = SB->Vertices + 3*(PB->VI[2]);

  double L1B[3], L2B[3];
  VecSub(V1B, QB, L1B);
  VecSub(V2B, QB, L2B);
  double AreaB = PB->Area;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  memset(Integral, 0, IDim*sizeof(double));
  for(int nA=0, ncpA=0; nA<NumPts; nA++)
   { 
     double uA=TCR[ncpA++], vA=TCR[ncpA++], wA=TCR[ncpA++], XA[3];
     for(int Mu=0; Mu<3; Mu++)
      XA[Mu] = QA[Mu] + uA*L1A[Mu] + vA*L2A[Mu];

     for(int nB=0, ncpB=0; nB<NumPts; nB++)
      { 
        double uB=TCR[ncpB++], vB=TCR[ncpB++], wB=TCR[ncpB++], XB[3];
        for(int Mu=0; Mu<3; Mu++)
         XB[Mu] = QB[Mu] + uB*L1B[Mu] + vB*L2B[Mu];

        Integrand(XA,XB,UserData,4.0*AreaA*AreaB*wA*wB,Integral);
      };
   };
}
