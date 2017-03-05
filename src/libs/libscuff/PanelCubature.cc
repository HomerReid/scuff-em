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

   bool UseSquareMapping;

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
void GPCIntegrand(unsigned ndim, const double *XiEta, void *params,
                  unsigned fdim, double *fval)
{
  (void) ndim; // unused 

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
   /*- get the evaluation points from the standard-triangle coords */
   /*--------------------------------------------------------------*/
   double Xi, Eta, Jacobian;
   if (Data->UseSquareMapping)
    { Xi  = XiEta[0];
      Eta = Xi*XiEta[1];
      Jacobian= 2.0 * Area * Xi;
    }
   else
    { Xi  = XiEta[0];
      Eta = XiEta[1];
      Jacobian = 2.0 * Area;
    };

   double X[3];
   X[0] = V0[0] + Xi*A[0] + Eta*B[0];
   X[1] = V0[1] + Xi*A[1] + Eta*B[1];
   X[2] = V0[2] + Xi*A[2] + Eta*B[2];

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
   for(unsigned n=0; n<fdim; n++)
    fval[n]*=Jacobian;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetPanelCubature(RWGGeometry *G, int ns, int np,
                      PCFunction Integrand, void *UserData, int IDim,
                      int MaxEvals, double RelTol, double AbsTol,
                      cdouble Omega, HVector *KN, 
                      int iQ, double BFSign,
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
  /* evaluate the two-dimensional integral by fixed-order         */
  /* cubature or by adaptive cubature                             */
  /*--------------------------------------------------------------*/
  if (MaxEvals==6 || MaxEvals==21 || MaxEvals==78)
   { Data->UseSquareMapping=false;
     int NumPts, Order = (MaxEvals==6) ? 4 : (MaxEvals==21) ? 9 : 20;
     double uv[2], *TCR=GetTCR(Order,&NumPts);
     memset(Result, 0, IDim*sizeof(double));
     double *DeltaResult=new double[IDim];
     for(int ncp=0; ncp<NumPts; ncp++)
      { double u1 = TCR[3*ncp+0];
        double v1 = TCR[3*ncp+1];
        double w  = TCR[3*ncp+2];
        uv[0] = u1+v1; uv[1] = v1;
        GPCIntegrand(2, uv, (void *)Data, IDim, DeltaResult);
        for(int n=0; n<IDim; n++)
         Result[n] += w*DeltaResult[n];
      };
     delete[] DeltaResult;
   }
  else
   { 
     Data->UseSquareMapping=true;
     double *Error  = new double[IDim];
     double Lower[2]={0.0, 0.0};
     double Upper[2]={1.0, 1.0};
     adapt_integrate(IDim, GPCIntegrand, (void *)Data, 2, 
                     Lower, Upper, MaxEvals, AbsTol, RelTol, 
                     Result, Error);
     delete[] Error;
   };
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

  RWGEdge *E = G->Surfaces[ns]->Edges[ne];

  GetPanelCubature(G, ns, E->iPPanel, Integrand, UserData, 
                   IDim, MaxEvals, RelTol, AbsTol, Omega, KN, 
                   E->PIndex, +1.0, Result);

  if (E->iMPanel!=-1)
   { double *ResultM = new double[IDim];
     GetPanelCubature(G, ns, E->iMPanel, Integrand, UserData, 
                      IDim, MaxEvals, RelTol, AbsTol, Omega, KN,
                      E->MIndex, -1.0, ResultM);
     VecPlusEquals(Result, 1.0, ResultM, IDim);
     delete[] ResultM;
   };

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

   bool UseSquareMapping;

   cdouble Omega;

   int NumContributingEdges1;
   cdouble KAlpha1[3], NAlpha1[3];
   double *Q1[3];
   double RWGPreFac1[3];

   int NumContributingEdges2;
   cdouble KAlpha2[3], NAlpha2[3];
   double *Q2[3];
   double RWGPreFac2[3];

} GPPCIData;

/***************************************************************/
/* integrand passed to adaptive cubature routine for cubature  */
/* over a pair of panels                                       */
/***************************************************************/
void GPPCIntegrand(unsigned ndim, const double *XiEta, void *params,
                   unsigned fdim, double *fval)
{
  (void)ndim;

   /*--------------------------------------------------------------*/
   /*- unpack fields from GPPCIData -------------------------------*/
   /*--------------------------------------------------------------*/
   GPPCIData *Data = (GPPCIData *)params;

   double *V01    = Data->V01;
   double *A1     = Data->A1;
   double *B1     = Data->B1;
   double Area1   = Data->Area1;
   int NumContributingEdges1 = Data->NumContributingEdges1;

   double *V02    = Data->V02;
   double *A2     = Data->A2;
   double *B2     = Data->B2;
   double Area2   = Data->Area2;
   int NumContributingEdges2 = Data->NumContributingEdges2;

   /*--------------------------------------------------------------*/
   /*- get the evaluation points from the standard-triangle coords */
   /*--------------------------------------------------------------*/
   double Xi1, Eta1, Xi2, Eta2, Jacobian;
   if (Data->UseSquareMapping)
    { Xi1  = XiEta[0];
      Eta1 = Xi1*XiEta[1];
      Xi2  = XiEta[2];
      Eta2 = Xi2*XiEta[3];
      Jacobian= 4.0 * Area1 * Area2 * Xi1 * Xi2;
    }
   else
    { Xi1  = XiEta[0];
      Eta1 = XiEta[1];
      Xi2  = XiEta[2];
      Eta2 = XiEta[3];
      Jacobian = 4.0 * Area1 * Area2;
    };

   double X1[3];
   X1[0] = V01[0] + Xi1*A1[0] + Eta1*B1[0];
   X1[1] = V01[1] + Xi1*A1[1] + Eta1*B1[1];
   X1[2] = V01[2] + Xi1*A1[2] + Eta1*B1[2];

   double X2[3];
   X2[0] = V02[0] + Xi2*A2[0] + Eta2*B2[0];
   X2[1] = V02[1] + Xi2*A2[1] + Eta2*B2[1];
   X2[2] = V02[2] + Xi2*A2[2] + Eta2*B2[2];

   /*--------------------------------------------------------------*/
   /*- get the surface currents at the evaluation point           -*/
   /*--------------------------------------------------------------*/
   cdouble K1[3], N1[3], DivK1, DivN1, K2[3], N2[3], DivK2, DivN2;
   K1[0]=K1[1]=K1[2]=DivK1=N1[0]=N1[1]=N1[2]=DivN1=0.0;
   K2[0]=K2[1]=K2[2]=DivK2=N2[0]=N2[1]=N2[2]=DivN2=0.0;
   for(int nce=0; nce<NumContributingEdges1; nce++)
    { 
      cdouble KAlpha   = Data->KAlpha1[nce];
      cdouble NAlpha   = Data->NAlpha1[nce];
      double *Q        = Data->Q1[nce];
      double RWGPreFac = Data->RWGPreFac1[nce];

      K1[0] += KAlpha * RWGPreFac * (X1[0]-Q[0]);
      K1[1] += KAlpha * RWGPreFac * (X1[1]-Q[1]);
      K1[2] += KAlpha * RWGPreFac * (X1[2]-Q[2]);
      DivK1 += 2.0*KAlpha * RWGPreFac;
      N1[0] += NAlpha * RWGPreFac * (X1[0]-Q[0]);
      N1[1] += NAlpha * RWGPreFac * (X1[1]-Q[1]);
      N1[2] += NAlpha * RWGPreFac * (X1[2]-Q[2]);
      DivN1 += 2.0*NAlpha * RWGPreFac;
    };
   for(int nce=0; nce<NumContributingEdges2; nce++)
    { 
      cdouble KAlpha   = Data->KAlpha2[nce];
      cdouble NAlpha   = Data->NAlpha2[nce];
      double *Q        = Data->Q2[nce];
      double RWGPreFac = Data->RWGPreFac2[nce];

      K2[0] += KAlpha * RWGPreFac * (X2[0]-Q[0]);
      K2[1] += KAlpha * RWGPreFac * (X2[1]-Q[1]);
      K2[2] += KAlpha * RWGPreFac * (X2[2]-Q[2]);
      DivK2 += 2.0*KAlpha * RWGPreFac;
      N2[0] += NAlpha * RWGPreFac * (X2[0]-Q[0]);
      N2[1] += NAlpha * RWGPreFac * (X2[1]-Q[1]);
      N2[2] += NAlpha * RWGPreFac * (X2[2]-Q[2]);
      DivN2 += 2.0*NAlpha * RWGPreFac;
    };

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
   for(unsigned n=0; n<fdim; n++)
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
                           int iQ1, double BFSign1, int iQ2, double BFSign2,
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
  /*- set up the calculation that the GPPCIntegrand routine will -*/
  /*- perform to compute the surface currents at each integration-*/
  /*- point in the first panel.                                  -*/
  /*--------------------------------------------------------------*/
  bool IsPEC = S1->IsPEC;
  int Offset = G->BFIndexOffset[ns1];
  double Sign;
  Data->NumContributingEdges1 = 0;
  if (!KN && (0<=iQ1 && iQ1<=2) ) // single unit-strength basis function contributes
   { 
     Sign=BFSign1;
     Data->NumContributingEdges1 = 1;
     Data->KAlpha1[0]            = 1.0; 
     Data->NAlpha1[0]            = 0.0;
     Data->Q1[0]                 = S1->Vertices + 3*(P1->VI[iQ1]);
     Data->RWGPreFac1[0]         = Sign*S1->Edges[P1->EI[iQ1]]->Length / (2.0*P1->Area);
   }
  else // up to three basis functions contribute 
   { 
     for(int nce=0; nce<3; nce++)
      { 
        int ne = P1->EI[nce];
        if (ne<0 || ne>=S1->NumEdges)
         continue; // this can happen for exterior edges

        RWGEdge *E = S1->Edges[ne];
        Sign = ( np1 == E->iMPanel ? -1.0 : 1.0 );

        int dnce = Data->NumContributingEdges1;
        Data->KAlpha1[dnce]    = IsPEC ? KN->GetEntry( Offset + ne ) : KN->GetEntry(Offset+2*ne);
        Data->NAlpha1[dnce]    = IsPEC ? 0.0 : -ZVAC*KN->GetEntry(Offset+2*ne+1);
        Data->Q1[dnce]         = S1->Vertices + 3*(P1->VI[nce]);
        Data->RWGPreFac1[dnce] = Sign * E->Length / (2.0*P1->Area);
        Data->NumContributingEdges1++;
      };
   };

  /*--------------------------------------------------------------*/
  /*- set up the calculation that the GPPCIntegrand routine will -*/
  /*- perform to compute the surface currents at each integration-*/
  /*- point in the second panel.                                 -*/
  /*--------------------------------------------------------------*/
  IsPEC = S2->IsPEC;
  Offset = G->BFIndexOffset[ns2];
  Data->NumContributingEdges2 = 0;
  if (!KN && (0<=iQ2 && iQ2<=2) ) // single unit-strength basis function contributes
   {
     Sign=BFSign2;
     Data->NumContributingEdges2 = 1;
     Data->KAlpha2[0]            = 1.0; 
     Data->NAlpha2[0]            = 0.0;
     Data->Q2[0]                 = S2->Vertices + 3*(P2->VI[iQ2]);
     Data->RWGPreFac2[0]         = Sign*S2->Edges[P2->EI[iQ2]]->Length / (2.0*P2->Area);
   }
  else // up to three basis functions contribute 
   { 
     for(int nce=0; nce<3; nce++)
      { 
        int ne = P2->EI[nce];
        if (ne<0 || ne>=S2->NumEdges)
         continue; // this can happen for exterior edges

        RWGEdge *E = S2->Edges[ne];
        Sign = ( np2 == E->iMPanel ? -1.0 : 1.0 );

        int dnce = Data->NumContributingEdges2;
        Data->KAlpha2[dnce]    = IsPEC ? KN->GetEntry( Offset + ne ) : KN->GetEntry(Offset+2*ne);
        Data->NAlpha2[dnce]    = IsPEC ? 0.0 : -ZVAC*KN->GetEntry(Offset+2*ne+1);
        Data->Q2[dnce]         = S2->Vertices + 3*(P2->VI[nce]);
        Data->RWGPreFac2[dnce] = Sign * E->Length / (2.0*P2->Area);
        Data->NumContributingEdges2++;
      };
   };
 
  /*--------------------------------------------------------------*/
  /* evaluate the four-dimensional integral by fixed-order        */
  /* cubature or by adaptive cubature                             */
  /*--------------------------------------------------------------*/
  if (MaxEvals==36 || MaxEvals==441)
   { Data->UseSquareMapping=false;
     int NumPts, Order = (MaxEvals==36) ? 4 : 9;
     double uv[4], *TCR=GetTCR(Order,&NumPts);
     memset(Result, 0, IDim*sizeof(double));
     double *DeltaResult=new double[IDim];
     for(int np=0; np<NumPts; np++)
      for(int npp=0; npp<NumPts; npp++)
       { double u1=TCR[3*np+0];  double v1=TCR[3*np+1];  double w=TCR[3*np+2];
         double u2=TCR[3*npp+0]; double v2=TCR[3*npp+1]; double wp=TCR[3*npp+2];
         uv[0] = u1+v1; uv[1] = v1; uv[2] = u2+v2; uv[3] = v2;
         GPPCIntegrand(4, uv, (void *)Data, IDim, DeltaResult);
         for(int n=0; n<IDim; n++)
          Result[n] += w*wp*DeltaResult[n];
       };
     delete[] DeltaResult;
   }
  else
   { Data->UseSquareMapping=true;
     double *Error  = new double[IDim];
     double Lower[4]={0.0, 0.0, 0.0, 0.0};
     double Upper[4]={1.0, 1.0, 1.0, 1.0};
     adapt_integrate(IDim, GPPCIntegrand, (void *)Data,
                     4, Lower, Upper, MaxEvals, AbsTol, RelTol, Result, Error);
     delete[] Error;
   };
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
                        E1->PIndex, +1.0, E2->PIndex, +1.0, ResultPP);

  np1 = E1->iPPanel;
  np2 = E2->iMPanel;
  GetPanelPanelCubature(G, ns1, np1, ns2, np2, Displacement,
                        Integrand, UserData, 
                        IDim, MaxEvals, RelTol, AbsTol, Omega, KN, 
                        E1->PIndex, +1.0, E2->MIndex, -1.0, ResultPM);

  np1 = E1->iMPanel;
  np2 = E2->iPPanel;
  GetPanelPanelCubature(G, ns1, np1, ns2, np2, Displacement,
                        Integrand, UserData, 
                        IDim, MaxEvals, RelTol, AbsTol, Omega, KN, 
                        E1->MIndex, -1.0, E2->PIndex, +1.0, ResultMP);

  np1 = E1->iMPanel;
  np2 = E2->iMPanel;
  GetPanelPanelCubature(G, ns1, np1, ns2, np2, Displacement,
                        Integrand, UserData, 
                        IDim, MaxEvals, RelTol, AbsTol, Omega, KN, 
                        E1->MIndex, -1.0, E2->MIndex, -1.0, ResultMM);
                        
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

/***************************************************************/
/* 20151118 streamlined implementation of GetBFCubature and    */
/* GetBFBFCubature that are 2x and 4x faster respectively.     */
/*                                                             */
/* In these versions, the user's integrand function must       */
/* *accumulate* contributions to Integral with weight Weight.  */
/***************************************************************/
void GetBFCubature2(RWGGeometry *G, int ns, int ne,
                    PCFunction2 Integrand, void *UserData, int IDim,
                    int Order, double *Integral, int PanelOnly)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGSurface *S = G->Surfaces[ns];
  RWGEdge *E    = S->Edges[ne];
  double Length = E->Length;
  double *QP    = S->Vertices + 3*(E->iQP);
  double *V1    = S->Vertices + 3*(E->iV1);
  double *V2    = S->Vertices + 3*(E->iV2);
  double *QM    = (E->iQM==-1) ? 0 : S->Vertices + 3*(E->iQM);

  double A[2][3], B[2][3], Area[2], *Q[2];
  VecSub(V1, QP, A[0]);
  VecSub(V2, QP, B[0]);
  Q[0]=QP;
  Area[0] = S->Panels[E->iPPanel]->Area;
  if (QM)
   { VecSub(V1, QM, A[1]);
     VecSub(V2, QM, B[1]);
     Q[1]=QM;
     Area[1] = S->Panels[E->iMPanel]->Area;
   };

  int NumPM=QM ? 2 : 1;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int NumPts;
  double *TCR=GetTCR(Order,&NumPts);
  if (TCR==0) ErrExit("invalid cubature order in GetBFCubature2");
  memset(Integral, 0, IDim*sizeof(double));
  for(int np=0, ncp=0; np<NumPts; np++)
   { 
     double u=TCR[ncp++];
     double v=TCR[ncp++];
     double w=TCR[ncp++];

     for(int PM=0; PM<NumPM; PM++)
      { 
        if (PanelOnly!=-1 && PM!=PanelOnly) 
         continue;

        double b[3], X[3], Sign=(PM==0) ? 1.0 : -1.0;
        double PreFac=Sign*Length/(2.0*Area[PM]);
        for(int Mu=0; Mu<3; Mu++)
         { b[Mu] = u*A[PM][Mu] + v*B[PM][Mu];
           X[Mu] = b[Mu]       + Q[PM][Mu];
           b[Mu] *= PreFac;
         };

        Integrand(X,b,2.0*PreFac,UserData,2.0*Area[PM]*w,Integral);
      };
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetBFBFCubature2(RWGSurface *S, int ne,
                      RWGSurface *SP, int neP,
                      PPCFunction2 Integrand,
                      void *UserData, int IDim,
                      int Order, double *Integral,
                      int PanelOnlyA, int PanelOnlyB)
{
  int NumPts;
  double *TCR=GetTCR(Order,&NumPts);
  if (TCR==0) ErrExit("invalid cubature order in GetBFCubature2");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGEdge *E    = S->Edges[ne];
  double Length = E->Length;
  double *QP    = S->Vertices + 3*(E->iQP);
  double *V1    = S->Vertices + 3*(E->iV1);
  double *V2    = S->Vertices + 3*(E->iV2);
  double *QM    = (E->iQM==-1) ? 0 : S->Vertices + 3*(E->iQM);

  double A[2][3], B[2][3], Area[2], *Q[2];
  VecSub(V1, QP, A[0]);
  VecSub(V2, QP, B[0]);
  Q[0]=QP;
  Area[0] = S->Panels[E->iPPanel]->Area;
  if (QM)
   { VecSub(V1, QM, A[1]);
     VecSub(V2, QM, B[1]);
     Q[1]=QM;
     Area[1] = S->Panels[E->iMPanel]->Area;
   };

  int NumPM=QM ? 2 : 1;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGEdge *EP    = SP->Edges[neP];
  double LengthP = EP->Length;
  double *QPP    = SP->Vertices + 3*(EP->iQP);
  double *V1P    = SP->Vertices + 3*(EP->iV1);
  double *V2P    = SP->Vertices + 3*(EP->iV2);
  double *QMP    = (EP->iQM==-1) ? 0 : SP->Vertices + 3*(EP->iQM);

  double AP[2][3], BP[2][3], AreaP[2], *QPArray[2];
  VecSub(V1P, QPP, AP[0]);
  VecSub(V2P, QPP, BP[0]);
  QPArray[0]=QPP;
  AreaP[0] = SP->Panels[EP->iPPanel]->Area;
  if (QMP)
   { VecSub(V1P, QMP, AP[1]);
     VecSub(V2P, QMP, BP[1]);
     QPArray[1]=QMP;
     AreaP[1] = SP->Panels[EP->iMPanel]->Area;
   };

  int NumPMP=QMP ? 2 : 1;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  memset(Integral, 0, IDim*sizeof(double));
  for(int np=0, ncp=0; np<NumPts; np++)
   { 
     double u=TCR[ncp++];
     double v=TCR[ncp++];
     double w=TCR[ncp++];

     for(int PM=0; PM<NumPM; PM++)
      { 
        if (PanelOnlyA!=-1 && PanelOnlyA!=PM) continue;

        double b[3], X[3], Sign=(PM==0) ? 1.0 : -1.0;
        double PreFac=Sign*Length/(2.0*Area[PM]);
        for(int Mu=0; Mu<3; Mu++)
         { b[Mu] = u*A[PM][Mu] + v*B[PM][Mu];
           X[Mu] = b[Mu]       + Q[PM][Mu];
           b[Mu] *= PreFac;
         };

        for(int npP=0, ncpP=0; npP<NumPts; npP++)
         { 
           double uP=TCR[ncpP++];
           double vP=TCR[ncpP++];
           double wP=TCR[ncpP++];

           for(int PMP=0; PMP<NumPMP; PMP++)
            { 
              if (PanelOnlyB!=-1 && PanelOnlyB!=PMP) continue;
              double bP[3], XP[3], SignP=(PMP==0) ? 1.0 : -1.0;
              double PreFacP=SignP*LengthP/(2.0*AreaP[PMP]);
              for(int Mu=0; Mu<3; Mu++)
               { bP[Mu] = uP*AP[PMP][Mu] + vP*BP[PMP][Mu];
                 XP[Mu] = bP[Mu]         + QPArray[PMP][Mu];
                 bP[Mu] *= PreFacP;
               };
  

              Integrand(X,b,2.0*PreFac,XP,bP,2.0*PreFacP,UserData,
                        4.0*Area[PM]*AreaP[PMP]*w*wP,Integral);
            };
         };        
      };
   };
}

void GetBFBFCubature2(RWGGeometry *G,
                      int ns, int ne, int nsP, int neP,
                      PPCFunction2 Integrand,
                      void *UserData, int IDim,
                      int Order, double *Integral,
                      int PanelOnlyA, int PanelOnlyB)
{
  GetBFBFCubature2(G->Surfaces[ns], ne, G->Surfaces[nsP], neP,
                   Integrand, UserData, IDim, Order, Integral,
                   PanelOnlyA, PanelOnlyB);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct BFBFCTDData
{
  PPCFunction2 *UserFunction;
  void *UserData; 
  int NTA, NTB, ncv;
  double LA, AreaA[2], L1A[2][3], L2A[3], *QA[2];
  double LB, AreaB[2], L1B[2][3], L2B[3], *QB[2];
 
} BFBFCTDData;

/***************************************************************/
/***************************************************************/
/***************************************************************/
void wyToXiEta(int ncv, int WhichRegion, const double wyyy[4], double XiEtaJ[5])
{
  double w=wyyy[0], y1=wyyy[1], y2=wyyy[2], y3=wyyy[3];
  double Xi1=0.0, Xi2=0.0, Eta1=0.0, Eta2=0.0, J=0.0;
  if (ncv>=2)
   { 
     double u1=0.0, u2=0.0, Xi1Min=0.0, Xi1Max=0.0;
     switch(WhichRegion)
      { case 0: u1  = -w*y1;              u2  = -w*y1*y2;
                Xi2 = w*(1.0-y1+y1*y2);
                Xi1Min=Xi2+u2-u1; Xi1Max=1.0;
                break;

        case 1: u1  = w*y1;               u2  = w*y1*y2;
                Xi2 = w*(1.0-y1);
                Xi1Min=Xi2; Xi1Max=1.0-u1;
                break;

        case 2: u1  = -w*y1*y2;           u2  = w*y1*(1.0-y2);
                Xi2 = w*(1.0-y1);
                Xi1Min=Xi2+u2-u1; Xi1Max=1.0;
                break;

        case 3: u1  = w*y1*y2;            u2  = -w*y1*(1.0-y2);
                Xi2 = w*(1.0-y1*y2);
                Xi1Min=Xi2; Xi1Max=1.0-u1;
                break;

        case 4: u1  = -w*y1*y2;           u2  = -w*y1;
                Xi2 = w;
                Xi1Min=Xi2; Xi1Max=1.0;
                break;

        case 5: u1  = w*y1*y2;            u2  = w*y1;
                Xi2 = w*(1.0-y1);
                Xi1Min=Xi2+u2-u1; Xi1Max=1.0-u1;
                break;
      };
     Xi1=Xi1Min + (Xi1Max-Xi1Min)*y3;
     J=w*w*y1*(Xi1Max-Xi1Min);
     Eta1 = u1 + Xi1;
     Eta2 = u2 + Xi2;
   }
  else // (ncv==1 || ncv==0)
   { 
     switch(WhichRegion)
      { case 0: Xi1=w;    Xi2=w*y1;   Eta1=w*y2;   Eta2=w*y2*y3;
                break;
        case 1: Eta1=w;  Eta2=w*y1;    Xi1=w*y2;    Xi2=w*y2*y3;
                break;
      };
     J=w*w*w*y2;
   };

  XiEtaJ[0] = Xi1;
  XiEtaJ[1] = Xi2;
  XiEtaJ[2] = Eta1;
  XiEtaJ[3] = Eta2;
  XiEtaJ[4] = J;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int BFBFCTDIntegrand(unsigned ndim, const double *wyyy, void *pBFBFCTDData,
                  unsigned fdim, double *fval)
{ 
  (void )ndim; // unused 

  memset(fval, 0, fdim*sizeof(double));

  BFBFCTDData *Data             = (BFBFCTDData *)pBFBFCTDData;
  void *UserData             = Data->UserData;
  PPCFunction2 *UserFunction = Data->UserFunction;
  double LA                  = Data->LA;
  double *AreaA              = Data->AreaA;
  int NTA                    = Data->NTA;
  double LB                  = Data->LB;
  double *AreaB              = Data->AreaB;
  int NTB                    = Data->NTB;
  int ncv                    = Data->ncv;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NumRegions = (ncv>=2) ? 6 : 2;
  for(int nr=0; nr<NumRegions; nr++)
   { 
     double XiEtaJ[5];
     wyToXiEta(ncv, nr, wyyy, XiEtaJ);
     double Xi1A = XiEtaJ[0];
     double Xi2A = XiEtaJ[1];
     double Xi1B = XiEtaJ[2];
     double Xi2B = XiEtaJ[3];
     double J    = XiEtaJ[4];
     if (J==0.0) continue;

     for(int ntA=0; ntA<NTA; ntA++)
      for(int ntB=0; ntB<NTB; ntB++)
       { 
         double DivbA = (ntA==0 ? 1.0 : -1.0 ) * LA/AreaA[ntA];
         double DivbB = (ntB==0 ? 1.0 : -1.0 ) * LB/AreaB[ntB];

         double bA[3], xA[3], bB[3], xB[3];

         VecLinComb(Xi1A, Data->L1A[ntA], Xi2A, Data->L2A, bA);
         VecLinComb(1.0, bA, 1.0, Data->QA[ntA], xA);
         VecScale(bA, DivbA / 2.0);

         VecLinComb(Xi1B, Data->L1B[ntB], Xi2B, Data->L2B, bB);
         VecLinComb(1.0, bB, 1.0, Data->QB[ntB], xB);
         VecScale(bB, DivbB / 2.0);

         double Weight = 4.0 * AreaA[ntA] * AreaB[ntB] * J;
         UserFunction(xA, bA, DivbA, xB, bB, DivbB, UserData, Weight, fval);
       };
   };

  return 0;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetBFBFCubatureTD(RWGSurface *SA, int neA,
                       RWGSurface *SB, int neB,
                       PPCFunction2 Integrand, void *UserData,
                       int IDim, double *Integral, double *Error,
                       int MaxEvals, double AbsTol, double RelTol)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGEdge    *EA  = SA->Edges[neA];
  double     *QPA = SA->Vertices + 3*(EA->iQP);
  double     *V1A = SA->Vertices + 3*(EA->iV1);
  double     *V2A = SA->Vertices + 3*(EA->iV2);
  double     *QMA = (EA->iQM==-1) ? 0 : SA->Vertices + 3*(EA->iQM);

  RWGEdge    *EB=SB->Edges[neB];
  double     *QPB = SB->Vertices + 3*(EB->iQP);
  double     *V1B = SB->Vertices + 3*(EB->iV1);
  double     *V2B = SB->Vertices + 3*(EB->iV2);
  double     *QMB = (EB->iQM==-1) ? 0 : SB->Vertices + 3*(EB->iQM);
  
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  BFBFCTDData MyData, *Data = &MyData;
  Data->UserFunction     = Integrand;
  Data->UserData         = UserData;

  Data->LA       = EA->Length;
  Data->NTA      = QMA ? 2 : 1;
  Data->AreaA[0] = SA->Panels[EA->iPPanel]->Area;
  Data->QA[0]    = QPA;
  VecSub(V1A, QPA, Data->L1A[0]);
  VecSub(V2A, V1A, Data->L2A);
  if (QMA)
   { Data->QA[1]    = QMA;
     Data->AreaA[1] = SA->Panels[EA->iMPanel]->Area;
     VecSub(V1A, QMA, Data->L1A[1]);
   };

  Data->LB       = EB->Length;
  Data->NTB      = QMB ? 2 : 1;
  Data->AreaB[0] = SB->Panels[EB->iPPanel]->Area;
  Data->QB[0]    = QPB;
  VecSub(V1B, QPB, Data->L1B[0]);
  VecSub(V2B, V1B, Data->L2B);
  if (QMB)
   { Data->QB[1]    = QMB;
     Data->AreaB[1] = SB->Panels[EB->iMPanel]->Area;
     VecSub(V1B, QMB, Data->L1B[1]);
   };

  Data->ncv = AssessBFPair(SA, neA, SB, neB);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double Lower[4]={0.0, 0.0, 0.0, 0.0};
  double Upper[4]={1.0, 1.0, 1.0, 1.0};
  hcubature(IDim, BFBFCTDIntegrand, (void *)Data, 4, Lower, Upper,
            MaxEvals, AbsTol, RelTol, ERROR_INDIVIDUAL,
            Integral, Error);
}

void GetBFBFCubatureTD(RWGGeometry *G,
                       int nsA, int neA, int nsB, int neB,
                       PPCFunction2 Integrand, void *UserData,
                       int IDim, double *Integral, double *Error,
                       int MaxEvals, double AbsTol, double RelTol)
{
  GetBFBFCubatureTD(G->Surfaces[nsA], neA, G->Surfaces[nsB], neB,  
                    Integrand, UserData, IDim, Integral, Error,
                    MaxEvals, AbsTol, RelTol);
}

} // namespace scuff
