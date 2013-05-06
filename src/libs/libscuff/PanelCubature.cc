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

   cdouble Omega;

   int NumContributingEdges;
   cdouble KAlpha[3], NAlpha[3];
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
   GPCIData *Data           = (GPCIData *)params;
   double *V0               = Data->V0;
   double *A                = Data->A;
   double *B                = Data->B;
   double Area              = Data->Area;
   int NumContributingEdges = Data->NumContributingEdges;

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
   cdouble K[3], DivK, N[3], DivN;
   K[0]=K[1]=K[2]=DivK=N[0]=N[1]=N[2]=DivN=0.0;
   for(int nce=0; nce<NumContributingEdges; nce++)
    { 
      cdouble KAlpha   = Data->KAlpha[nce];
      cdouble NAlpha   = Data->NAlpha[nce];
      double *Q        = Data->Q[nce];
      double RWGPreFac = Data->RWGPreFac[nce];

      K[0] += KAlpha * RWGPreFac * (X[0]-Q[0]);
      K[1] += KAlpha * RWGPreFac * (X[1]-Q[1]);
      K[2] += KAlpha * RWGPreFac * (X[2]-Q[2]);
      DivK += 2.0*KAlpha * RWGPreFac;
      N[0] += NAlpha * RWGPreFac * (X[0]-Q[0]);
      N[1] += NAlpha * RWGPreFac * (X[1]-Q[1]);
      N[2] += NAlpha * RWGPreFac * (X[2]-Q[2]);
      DivN += 2.0*NAlpha * RWGPreFac;
    };

   /*--------------------------------------------------------------*/
   /*- fill in the PCData structure passed to the user's function, */
   /*- then call user's function                                   */
   /*--------------------------------------------------------------*/
   PCData MyPCData, *PCD=&MyPCData;
   PCD->nHat  = Data->nHat; 
   PCD->Omega = Data->Omega;
   PCD->K     = K;
   PCD->DivK  = DivK;
   PCD->N     = N;
   PCD->DivN  = DivN;

   Data->Integrand(X, PCD, Data->UserData, fval);

   /*--------------------------------------------------------------*/
   /*- put in Jacobian factors ------------------------------------*/
   /*--------------------------------------------------------------*/
   for(int n=0; n<fdim; n++)
    fval[n]*=Jacobian;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetPanelCubature(RWGGeometry *G, int ns, int np,
                      PCFunction Integrand, void *UserData, int IDim,
                      int MaxEvals, double RelTol, double AbsTol,
                      cdouble Omega, HVector *KN, int iQ, double BFSign,
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

  Data->Omega     = Omega;
  Data->nHat      = P->ZHat;
  Data->Area      = P->Area;

  Data->Integrand = Integrand;
  Data->UserData  = UserData;

  /*--------------------------------------------------------------*/
  /*- set up the calculation that the GPCIntegrand routine will  -*/
  /*- perform to compute the surface currents at each integration-*/
  /*- point.                                                     -*/
  /*- If the user gave us a non-NULL KN vector, then the surface -*/
  /*- current is computed by summing the contributions of up to  -*/
  /*- 3 RWG basis functions. Otherwise, if KN==NULL but the user -*/
  /*- specified a reasonable value for iQ, then the surface      -*/
  /*- current is taken to be that of the RWG basis function      -*/
  /*- assocated to edge #iQ, populated with unit strength, and   -*/
  /*- signed with sign BFSign.                                   -*/
  /*--------------------------------------------------------------*/
  bool IsPEC = S->IsPEC;
  int Offset = G->BFIndexOffset[ns];
  double Sign;
  Data->NumContributingEdges = 0;
  if (!KN && (0<=iQ && iQ<=2) ) // single unit-strength basis function contributes
   { Sign=BFSign;
     Data->KAlpha[0]    = 1.0; 
     Data->NAlpha[0]    = 0.0;
     Data->Q[0]         = S->Vertices + 3*(P->VI[iQ]);
     Data->RWGPreFac[0] = Sign*S->Edges[P->EI[iQ]]->Length / (2.0*P->Area);
     Data->NumContributingEdges++;
   }
  else // up to three basis functions contribute 
   { 
     for(int nce=0; nce<3; nce++)
      { 
        int ne = P->EI[nce];
        if (ne<0 || ne>=S->NumEdges)
         continue; // this can happen for exterior edges

        RWGEdge *E = S->Edges[ne];
        Sign = ( np == E->iMPanel ? -1.0 : 1.0 );

        int dnce = Data->NumContributingEdges;
        Data->KAlpha[dnce]    = IsPEC ? KN->GetEntry( Offset + ne ) : KN->GetEntry(Offset+2*ne);
        Data->NAlpha[dnce]    = IsPEC ? 0.0 : -ZVAC*KN->GetEntry(Offset+2*ne+1);
        Data->Q[dnce]         = S->Vertices + 3*(P->VI[nce]);
        Data->RWGPreFac[dnce] = Sign * E->Length / (2.0*P->Area);
        Data->NumContributingEdges++;
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
void GetBFCubature(RWGGeometry *G, int ns, int ne,
                   PCFunction Integrand, void *UserData, int IDim,
                   int MaxEvals, double RelTol, double AbsTol,
                   cdouble Omega, HVector *KN,
                   double *Result)
{
  double *ResultP = new double[IDim];
  double *ResultM = new double[IDim];

  RWGEdge *E = G->Surfaces[ns]->Edges[ne];

  GetPanelCubature(G, ns, E->iPPanel, Integrand, UserData, 
                   IDim, MaxEvals, RelTol, AbsTol, Omega, KN, 
                   E->PIndex, +1.0, ResultP);

  GetPanelCubature(G, ns, E->iMPanel, Integrand, UserData, 
                   IDim, MaxEvals, RelTol, AbsTol, Omega, KN,
                   E->MIndex, +1.0, ResultM);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  for(int n=0; n<IDim; n++)
   Result[n] = ResultP[n] + ResultM[n];

  delete[] ResultP;
  delete[] ResultM;
  
}

/***************************************************************/
/* data structure passed to GPPCIntegrand **********************/
/***************************************************************/
typedef struct GPPCIData
{
   PPCFunction Integrand;
   void *UserData;

   double *V01, A1[3], B1[3];
   double V02[3], A2[3], B2[3];
   double Area1, Area2;
   double *nHat1, *nHat2;

   cdouble Omega;
   HVector *KN;
   bool IsPEC1, IsPEC2;
   int Offset1, Offset2;
   int EdgeIndices1[3], EdgeIndices2[3];
   double *Q1[3], *Q2[3];
   double RWGPreFac1[3], RWGPreFac2[3];

} GPPCIData;

/***************************************************************/
/* integrand passed to adaptive cubature routine for cubature  */
/* over a pair of panels                                       */
/***************************************************************/
void GPPCIntegrand(unsigned ndim, const double *uv, void *params,
                   unsigned fdim, double *fval)
{
   /*--------------------------------------------------------------*/
   /*- unpack fields from GPPCIData -------------------------------*/
   /*--------------------------------------------------------------*/
   GPPCIData *Data = (GPPCIData *)params;

   double *V01    = Data->V01;
   double *A1     = Data->A1;
   double *B1     = Data->B1;
   double Area1   = Data->Area1;

   double *V02    = Data->V02;
   double *A2     = Data->A2;
   double *B2     = Data->B2;
   double Area2   = Data->Area2;

   /*--------------------------------------------------------------*/
   /*- get the evaluation points from the standard-triangle coords */
   /*--------------------------------------------------------------*/
   double u1 = uv[0];
   double v1 = u1*uv[1];
   double u2 = uv[2];
   double v2 = u2*uv[3];
   double Jacobian= 4.0 * Area1 * Area2 * u1 * u2;

   double X1[3];
   X1[0] = V01[0] + u1*A1[0] + v1*B1[0];
   X1[1] = V01[1] + u1*A1[1] + v1*B1[1];
   X1[2] = V01[2] + u1*A1[2] + v1*B1[2];

   double X2[3];
   X2[0] = V02[0] + u2*A2[0] + v2*B2[0];
   X2[1] = V02[1] + u2*A2[1] + v2*B2[1];
   X2[2] = V02[2] + u2*A2[2] + v2*B2[2];

   /*--------------------------------------------------------------*/
   /*- get the surface currents at the evaluation point           -*/
   /*--------------------------------------------------------------*/
   cdouble K1[3], N1[3], DivK1, DivN1, K2[3], N2[3], DivK2, DivN2;
   K1[0]=K1[1]=K1[2]=DivK1=N1[0]=N1[1]=N1[2]=DivN1=0.0;
   K2[0]=K2[1]=K2[2]=DivK2=N2[0]=N2[1]=N2[2]=DivN2=0.0;
   HVector *KN = Data->KN;
   if (KN)
    { 
      int IsPEC1  = Data->IsPEC1;
      int IsPEC2  = Data->IsPEC2;
      int Offset1 = Data->Offset1;
      int Offset2 = Data->Offset2;
      for(int nce=0; nce<3; nce++)
       {  
         double *Q;
         double RWGPreFac;
         cdouble KAlpha;
         cdouble NAlpha;
         int ne;
  
         ne = Data->EdgeIndices1[nce];
         if (ne>=0)
          { Q         = Data->Q1[nce];
            RWGPreFac = Data->RWGPreFac1[nce];
            KAlpha = IsPEC1 ? KN->GetEntry(Offset1 +   ne)
                            : KN->GetEntry(Offset1 + 2*ne);
            NAlpha = IsPEC1 ? 0.0
                            : -1.0*ZVAC*KN->GetEntry(Offset1 + 2*ne+1);

            K1[0] += KAlpha * RWGPreFac * (X1[0]-Q[0]);
            K1[1] += KAlpha * RWGPreFac * (X1[1]-Q[1]);
            K1[2] += KAlpha * RWGPreFac * (X1[2]-Q[2]);
            DivK1 += 2.0 * KAlpha * RWGPreFac;
            N1[0] += NAlpha * RWGPreFac * (X1[0]-Q[0]);
            N1[1] += NAlpha * RWGPreFac * (X1[1]-Q[1]);
            N1[2] += NAlpha * RWGPreFac * (X1[2]-Q[2]);
            DivN1 += 2.0 * NAlpha * RWGPreFac;
          };
  
         ne = Data->EdgeIndices2[nce];
         if (ne>=0)
          { Q         = Data->Q2[nce];
            RWGPreFac = Data->RWGPreFac2[nce];
            KAlpha = IsPEC2 ? KN->GetEntry(Offset2 +   ne)
                            : KN->GetEntry(Offset2 + 2*ne);
            NAlpha = IsPEC2 ? 0.0 
                            : -1.0*ZVAC*KN->GetEntry(Offset2 + 2*ne+1);

            K2[0] += KAlpha * RWGPreFac * (X2[0]-Q[0]);
            K2[1] += KAlpha * RWGPreFac * (X2[1]-Q[1]);
            K2[2] += KAlpha * RWGPreFac * (X2[2]-Q[2]);
            DivK2 += 2.0 * KAlpha * RWGPreFac;
            N2[0] += NAlpha * RWGPreFac * (X2[0]-Q[0]);
            N2[1] += NAlpha * RWGPreFac * (X2[1]-Q[1]);
            N2[2] += NAlpha * RWGPreFac * (X2[2]-Q[2]);
            DivN2 += 2.0 * NAlpha * RWGPreFac;
         };
       };  // for(int nce=0; nce<3; nce++)
    }; // if (KN)

   /*--------------------------------------------------------------*/
   /*- fill in the PCData structure passed to the user's function, */
   /*- then call user's function                                   */
   /*--------------------------------------------------------------*/
   PPCData MyPPCData, *PPCD=&MyPPCData;
   PPCD->Omega = Data->Omega;
   PPCD->nHat1 = Data->nHat1; 
   PPCD->nHat2 = Data->nHat2;
   PPCD->K1    = K1;
   PPCD->N1    = N1;
   PPCD->DivK1 = DivK1;
   PPCD->DivN1 = DivN1;
   PPCD->K2    = K2;
   PPCD->N2    = N2;
   PPCD->DivK2 = DivK2;
   PPCD->DivN2 = DivN2;

   Data->Integrand(X1, X2, PPCD, Data->UserData, fval);

   /*--------------------------------------------------------------*/
   /*- put in Jacobian factors ------------------------------------*/
   /*--------------------------------------------------------------*/
   for(int n=0; n<fdim; n++)
    fval[n]*=Jacobian;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetPanelPanelCubature(RWGGeometry *G, int ns1, int np1, int ns2, int np2,
                           double *Displacement,
                           PPCFunction Integrand, void *UserData, int IDim,
                           int MaxEvals, double RelTol, double AbsTol,
                           cdouble Omega, HVector *KN,
                           double *Result)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGSurface *S1 = G->Surfaces[ns1];
  RWGPanel *P1   = S1->Panels[np1];
  double *V01    = S1->Vertices + 3*(P1->VI[0]);
  double *V11    = S1->Vertices + 3*(P1->VI[1]);
  double *V21    = S1->Vertices + 3*(P1->VI[2]);

  RWGSurface *S2 = G->Surfaces[ns2];
  RWGPanel *P2   = S2->Panels[np2];
  double *V02    = S2->Vertices + 3*(P2->VI[0]);
  double *V12    = S2->Vertices + 3*(P2->VI[1]);
  double *V22    = S2->Vertices + 3*(P2->VI[2]);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  struct GPPCIData MyGPPCIData, *Data=&MyGPPCIData;
  Data->V01 = V01;
  VecSub(V11, V01, Data->A1);
  VecSub(V21, V11, Data->B1);

  Data->Omega     = Omega;

  Data->nHat1     = P1->ZHat;
  Data->Area1     = P1->Area;

  if (Displacement)
   { Data->V02[0] = V02[0] + Displacement[0];
     Data->V02[1] = V02[1] + Displacement[1];
     Data->V02[2] = V02[2] + Displacement[2];
   }
  else 
   { Data->V02[0] = V02[0];
     Data->V02[1] = V02[1];
     Data->V02[2] = V02[2];
   };
  VecSub(V12, V02, Data->A2);
  VecSub(V22, V12, Data->B2);

  Data->nHat2     = P2->ZHat;
  Data->Area2     = P2->Area;

  Data->Integrand = Integrand;
  Data->UserData  = UserData;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Data->KN = KN;
  if (KN)
   { 
     Data->IsPEC1  = S1->IsPEC;
     Data->Offset1 = G->BFIndexOffset[ns1];

     Data->IsPEC2  = S2->IsPEC;
     Data->Offset2 = G->BFIndexOffset[ns2];

     for(int nce=0; nce<3; nce++)
      { 
        RWGEdge *E;
        double Sign;
        int ne;
    
        ne = P1->EI[nce];
        Data->EdgeIndices1[nce] = ne;
        if (ne>=0) 
         { 
           Data->Q1[nce] = S1->Vertices + 3*(P1->VI[nce]);
           E = S1->Edges[ne];
           Sign = ( np1 == E->iMPanel ? -1.0 : 1.0 );
           Data->RWGPreFac1[nce] = Sign * E->Length / (2.0*P1->Area);
         };

        ne = P2->EI[nce];
        Data->EdgeIndices2[nce] = ne;
        if (ne>=0) 
         { 
           Data->Q2[nce] = S2->Vertices + 3*(P2->VI[nce]);
           E = S2->Edges[ne];
           Sign = ( np2 == E->iMPanel ? -1.0 : 1.0 );
           Data->RWGPreFac2[nce] = Sign * E->Length / (2.0*P2->Area);
         };

      };
   };
 
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double *Error  = new double[IDim];
  double Lower[4]={0.0, 0.0, 0.0, 0.0};
  double Upper[4]={1.0, 1.0, 1.0, 1.0};
  adapt_integrate(IDim, GPPCIntegrand, (void *)Data,
		  4, Lower, Upper, MaxEvals, AbsTol, RelTol, Result, Error);

  delete[] Error;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetBFBFCubature(RWGGeometry *G, int ns1, int ne1, int ns2, int ne2,
                     double *Displacement, 
                     PPCFunction Integrand, void *UserData, int IDim,
                     int MaxEvals, double RelTol, double AbsTol,
                     cdouble Omega, HVector *KN,
                     double *Result)
{
  double *ResultPP = new double[IDim];
  double *ResultPM = new double[IDim];
  double *ResultMP = new double[IDim];
  double *ResultMM = new double[IDim];

  RWGEdge *E1 = G->Surfaces[ns1]->Edges[ne1];
  RWGEdge *E2 = G->Surfaces[ns2]->Edges[ne2];

  int np1, np2;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  np1 = E1->iPPanel;
  np2 = E2->iPPanel;
  GetPanelPanelCubature(G, ns1, np1, ns2, np2, Displacement,
                        Integrand, UserData, 
                        IDim, MaxEvals, RelTol, AbsTol, Omega, KN, 
                        ResultPP);

  np1 = E1->iPPanel;
  np2 = E2->iMPanel;
  GetPanelPanelCubature(G, ns1, np1, ns2, np2, Displacement,
                        Integrand, UserData, 
                        IDim, MaxEvals, RelTol, AbsTol, Omega, KN, 
                        ResultPM);

  np1 = E1->iMPanel;
  np2 = E2->iPPanel;
  GetPanelPanelCubature(G, ns1, np1, ns2, np2, Displacement,
                        Integrand, UserData, 
                        IDim, MaxEvals, RelTol, AbsTol, Omega, KN, 
                        ResultMP);
  np1 = E1->iMPanel;
  np2 = E2->iMPanel;
  GetPanelPanelCubature(G, ns1, np1, ns2, np2, Displacement,
                        Integrand, UserData, 
                        IDim, MaxEvals, RelTol, AbsTol, Omega, KN, 
                        ResultMM);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  for(int n=0; n<IDim; n++)
   Result[n] = ResultPP[n] + ResultPM[n] + ResultMP[n] + ResultMM[n];

  delete[] ResultPP;
  delete[] ResultPM;
  delete[] ResultMP;
  delete[] ResultMM;
  
}

} // namespace scuff
