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
#include "TaylorDuffy.h"

namespace scuff {

#define II cdouble(0.0,1.0)

/**********************************************************************/
/*a 'panel-panel cubature' function accepts input parameters giving  */
/**********************************************************************/
typedef void PCFunction(double *x, double *b, double Divb,
                        void *UserData, double *Result);

typedef void PPCFunction(double *x, double *b, double Divb,
                         double *xp, double *bp, double Divbp,
		   	 void *UserData, double *Result);

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
typedef struct PCIData 
 {
   double *V0, double A[3], double B[3], double *Q;
   double PreFac;
   PCFunction Integrand;
   void *UserData;

 } PCIData;

void PCIntegrand(unsigned ndim, const double *x, void *params,
                 unsigned fdim, double *fval)
{ 
  (void *)ndim; // unused

  /*--------------------------------------------------------------*/
  /*- unpack fields from input data structure --------------------*/
  /*--------------------------------------------------------------*/
  PCIData *Data         = (PCIData *)params;
  double *V0            = Data->V0;
  double *A             = Data->A;
  double *B             = Data->B;
  double *Q             = Data->Q;
  double *PreFac        = Data->PreFac;
  PCFunction *Integrand = Data->Integrand;
  void *UserData        = Data->UserData;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double u=x[0];
  double v=u*x[1];
  double Jacobian=u;

  double x[3], double f[3];

  x[0] = V0[0] + u*A[0] + v*B[0];
  x[1] = V0[1] + u*A[1] + v*B[1];
  x[2] = V0[2] + u*A[2] + v*B[2];

  f[0] = PreFac*(x[0] - Q[0]);
  f[1] = PreFac*(x[1] - Q[1]);
  f[2] = PreFac*(x[2] - Q[2]);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Integrand(x, f, UserData, fVal);

  for(int nf=0; nf<fdim; nf++)
   fVal[nf] *= Jacobian;

}

/**********************************************************************/
/* evaluate a numerical cubature over a single panel.                 */
/**********************************************************************/
void GetPanelCubature(RWGSurface *S, int np, int iQ,
                      PCFunction *Integrand, void *UserData, int IDim,
                      int Order, double RelTol, double AbsTol, 
                      double *Result);
{
  /*--------------------------------------------------------------*/
  /*- unpack panel data ------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGPanel *P   = S->Panels[np];
  double *V0    = S->Vertices + 3*(P->VI[ (iQ+0) % 3 ]);
  double *V1    = S->Vertices + 3*(P->VI[ (iQ+1) % 3 ]);
  double *V2    = S->Vertices + 3*(P->VI[ (iQ+2) % 3 ]);
  double Length = VecDistance(V2,V1);
  double Area   = P->Area;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (Order==0)
   { PCIData MyPCIData, *Data=&MyPCIData;

     Data->V0        = S->Vertices + 3*P->VI[0];
     Data->A         = S->Vertices + 3*P->VI[0];
     Data->B         = S->Vertices + 3*P->VI[0];
     Data->Q         = S->Vertices + 3*P->VI[0];
     Data->RWGPreFac = Length / (2.0*Area);
     Data->Integrand = Integrand;
     Data->UserData  = UserData;

     double Lower[2]={0.0,0.0};
     double Lower[2]={1.0,1.0};
     adapt_integrate

   }

  VecSub(Va[1], Va[0], A);
  VecSub(Va[2], Va[0], B);
  Q=Qa;

  double *V0P, AP[3], BP[3], *QP;
  V0P=Vb[0];
  VecSub(Vb[1], Vb[0], AP);
  VecSub(Vb[2], Vb[0], BP);
  QP=Qb;

  /***************************************************************/
  /***************************************************************/
  double *TCR;
  int NumPts;
  switch(Order)
   { 
     case 4:   TCR=GetTCR(4, &NumPts);
     case 20:  TCR=GetTCR(20, &NumPts);
     default:  TCR=0;
   };

  /***************************************************************/
  /* outer loop **************************************************/
  /***************************************************************/
  int np, ncp, npp, ncpp;
  int Mu;
  double u, v, w, up, vp, wp;
  double X[3], F[3], XP[3], FP[3], R[3];
  cdouble k = Args->k;
  memset(H,0,2*sizeof(cdouble));
  if (GradH) memset(GradH,0,6*sizeof(cdouble));
  if (dHdT) memset(dHdT,0,2*NumTorqueAxes*sizeof(cdouble));
  for(np=ncp=0; np<NumPts; np++) 
   { 
     u=TCR[ncp++]; v=TCR[ncp++]; w=TCR[ncp++];

     /***************************************************************/
     /* set X and F=X-Q *********************************************/
     /***************************************************************/
     for(Mu=0; Mu<3; Mu++)
      { X[Mu] = V0[Mu] + u*A[Mu] + v*B[Mu];
        F[Mu] = X[Mu] - Q[Mu];
      };

     /***************************************************************/
     /* inner loop to calculate value of inner integrand ************/
     /***************************************************************/
     memset(HInner,0,2*sizeof(cdouble));
     if (GradH) memset(GradHInner,0,6*sizeof(cdouble));
     if (dHdT) memset(dHdTInner,0,2*NumTorqueAxes*sizeof(cdouble));
     for(npp=ncpp=0; npp<NumPts; npp++)
      { 
        up=TCR[ncpp++]; vp=TCR[ncpp++]; wp=TCR[ncpp++];

        /***************************************************************/ 
        /* set XP and FP=XP-QP *****************************************/
        /***************************************************************/
        for(Mu=0; Mu<3; Mu++)
         { XP[Mu] = V0P[Mu] + up*AP[Mu] + vp*BP[Mu];
           FP[Mu] = XP[Mu] - QP[Mu];
           R[Mu] = X[Mu] - XP[Mu];
         };
      
        if ( Args->GInterp )
         AssembleInnerPPIIntegrand_Interp(wp, k, R, F, FP, 
                                          Args->GInterp, 
                                          NumTorqueAxes, GammaMatrix, 
                                          HInner, GradHInner, dHdTInner);
        else
         AssembleInnerPPIIntegrand_NoInterp(wp, k, R, X, F, FP, 
                                            DeSingularize,
                                            NumTorqueAxes, GammaMatrix, 
                                            HInner, GradHInner, dHdTInner);

      }; /* for(npp=ncpp=0; npp<NumPts; npp++) */

     /*--------------------------------------------------------------*/
     /*- accumulate contributions to outer integral                  */
     /*--------------------------------------------------------------*/
     H[0]+=w*HInner[0];
     H[1]+=w*HInner[1];
     if (GradH)
      for(Mu=0; Mu<6; Mu++)
       GradH[Mu]+=w*GradHInner[Mu];
     if (dHdT)
      for(Mu=0; Mu<2*NumTorqueAxes; Mu++)
       dHdT[Mu]+=w*dHdTInner[Mu];

   }; // for(np=ncp=0; np<nPts; np++) 

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  memcpy(Args->H, H, 2*sizeof(cdouble));
  if (GradH) memcpy(Args->GradH, GradH, 6*sizeof(cdouble));
  if (dHdT) memcpy(Args->dHdT, dHdT, 2*NumTorqueAxes*sizeof(cdouble));

}

} // namespace scuff
