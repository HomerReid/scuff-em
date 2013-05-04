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
 * PanelCubature .cc -- libscuff routines for computing numerical 
 *                   -- cubatures over panels or over pairs of panels
 * 
 * homer reid        -- 11/2005 -- 4/2013
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
#include "PanelCubature.h"

namespace scuff {

/***************************************************************/
/* data structure passed to GPCIntegrand ***********************/
/***************************************************************/
typedef struct GPCIData
 { 
   PCFunction Integrand;
   void *UserData;

   double *V0, A[3], B[3];

   double Area;
   double *nHat;

   HVector *KN;
   bool IsPEC;
   int Offset;
   int EdgeIndices[3];
   double *Q[3];
   double RWGPreFac[3];
 
 } GPCIData;

/***************************************************************/
/* integrand passed to adaptive cubature routine for cubature  */
/* over a single panel                                         */
/***************************************************************/
void GPCIntegrand(unsigned ndim, const double *uv, void *params,
                  unsigned fdim, double *fval)
{
   /*--------------------------------------------------------------*/
   /*- unpack fields from GPCIData --------------------------------*/
   /*--------------------------------------------------------------*/
   GPCIData *Data = (GPCIData *)params;
   double *V0     = Data->V0;
   double *A      = Data->A;
   double *B      = Data->B;
   double Area    = Data->Area;

   /*--------------------------------------------------------------*/
   /*- get the evaluation point from the standard-triangle coords  */
   /*--------------------------------------------------------------*/
   double u = uv[0];
   double v = u*uv[1];
   double Jacobian= 2.0 * Area * u;

   double X[3];
   X[0] = V0[0] + u*A[0] + v*B[0];
   X[1] = V0[1] + u*A[1] + v*B[1];
   X[2] = V0[2] + u*A[2] + v*B[2];

   /*--------------------------------------------------------------*/
   /*- get the surface currents at the evaluation point           -*/
   /*--------------------------------------------------------------*/
   cdouble K[3], N[3];
   K[0]=K[1]=K[2]=N[0]=N[1]=N[2]=0.0;
   HVector *KN = Data->KN;
   if (KN)
    { int IsPEC  = Data->IsPEC;
      int Offset = Data->Offset;
      for(int nce=0; nce<3; nce++)
       {  
         int ne           = Data->EdgeIndices[nce];
         if (ne<0) continue; // this happens if the edge is an exterior edge

         double *Q        = Data->Q[nce];
         double RWGPreFac = Data->RWGPreFac[nce];

         cdouble KAlpha = (IsPEC ? KN->GetEntry(Offset + ne) : KN->GetEntry(Offset + 2*ne) );
         cdouble NAlpha = (IsPEC ? 0.0 : -1.0*ZVAC*KN->GetEntry(Offset + 2*ne+1) );

         K[0] += KAlpha * RWGPreFac * (X[0]-Q[0]);
         K[1] += KAlpha * RWGPreFac * (X[1]-Q[1]);
         K[2] += KAlpha * RWGPreFac * (X[2]-Q[2]);
         N[0] += NAlpha * RWGPreFac * (X[0]-Q[0]);
         N[1] += NAlpha * RWGPreFac * (X[1]-Q[1]);
         N[2] += NAlpha * RWGPreFac * (X[2]-Q[2]);
       };
    };

   /*--------------------------------------------------------------*/
   /*- fill in the PCData structure passed to the user's function, */
   /*- then call user's function                                   */
   /*--------------------------------------------------------------*/
   PCData MyPCData, *PCD=&MyPCData;
   PCD->nHat = Data->nHat; 
   PCD->K    = K;
   PCD->N    = N;

   Data->Integrand(X, PCD, Data->UserData, fval);

   /*--------------------------------------------------------------*/
   /*- put in Jacobian factors ------------------------------------*/
   /*--------------------------------------------------------------*/
   for(int n=0; n<ndim; n++)
    fval[n]*=Jacobian;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetPanelCubature(RWGGeometry *G, int ns, int np,
                      PCFunction Integrand, void *UserData, int IDim,
                      int MaxEvals, double RelTol, double AbsTol,
                      cdouble Omega, HVector *KN,
                      double *Result)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGSurface *S = G->Surfaces[ns];
  RWGPanel *P   = S->Panels[np];
  double *V0    = S->Vertices + 3*(P->VI[0]);
  double *V1    = S->Vertices + 3*(P->VI[1]);
  double *V2    = S->Vertices + 3*(P->VI[2]);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  struct GPCIData MyGPCIData, *Data=&MyGPCIData;
  Data->V0 = V0;
  VecSub(V1, V0, Data->A);
  VecSub(V2, V1, Data->B);

  Data->nHat      = P->ZHat;
  Data->Area      = P->Area;

  Data->Integrand = Integrand;
  Data->UserData  = UserData;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Data->KN = KN;
  if (KN)
   { 
     Data->IsPEC  = S->IsPEC;
     Data->Offset = G->BFIndexOffset[ns];
     for(int nce=0; nce<3; nce++)
      { 
        int ne = P->EI[nce];
        Data->EdgeIndices[nce] = ne;
        if (ne<0) continue; // this happens if the edge is an exterior edge

        Data->Q[nce] = S->Vertices + 3*(P->VI[nce]);

        RWGEdge *E = S->Edges[ne];
        double Sign = ( np == E->iMPanel ? -1.0 : 1.0 );
        Data->RWGPreFac[nce] = Sign * E->Length / (2.0*P->Area);
      };
   };
 
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double *Error  = new double[IDim];
  double Lower[2]={0.0, 0.0};
  double Upper[2]={1.0, 1.0};
  adapt_integrate(IDim, GPCIntegrand, (void *)Data,
		  2, Lower, Upper, MaxEvals, AbsTol, RelTol, Result, Error);

  delete[] Error;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetPanelPanelCubature(RWGGeometry *G, int ns1, int np1, int ns2, int np2,
                           PCFunction *Integrand, void *UserData, int IDim,
                           int Order, double RelTol, double AbsTol,
                           cdouble Omega, HVector *KN,
                           double *Result)
{
}

} // namespace scuff
