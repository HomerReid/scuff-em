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
 * GetPanelPotentials.cc
 *
 * homer reid  -- 9/2011
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <libhrutil.h>
#include <libhmat.h>
#include <libscuff.h>
#include <libSGJC.h>

#define ABSTOL 1.0e-8
#define RELTOL 1.0e-4

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
int WhichCase;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

#define II cdouble(0.0,1.0)

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct PAPVData
 { 
   double A[3];   // V2 - V1 
   double B[3];   // V3 - V2 
   cdouble IK;
 } PAPVData;

/***************************************************************/
/* integrand components:                                       */
/*                                                             */
/* 0: \int_0^1 du exp(X*Alpha) / X                             */
/* 1: \int_0^1 du u*exp(X*Alpha) / X                           */
/* 2: alpha*\int_0^1 du u*exp(X*Alpha) / X                     */
/*                                                             */
/* where X = ik|A + Alpha*B|                                   */
/*                                                             */
/***************************************************************/
static void PAPVIntegrand(unsigned ndim, const double *x, void *params, unsigned fdim, double *fval)
{ 
  (void) ndim;
  (void) fdim;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  PAPVData *D=(PAPVData *)params;
  double *A=D->A;
  double *B=D->B;
  cdouble IK=D->IK;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double Alpha, ApBAlpha[3], X;
  cdouble IKX, ExpFac;

  Alpha=x[0];
  ApBAlpha[0] = A[0] + Alpha*B[0];
  ApBAlpha[1] = A[1] + Alpha*B[1];
  ApBAlpha[2] = A[2] + Alpha*B[2];
  X = VecNorm(ApBAlpha);
  IKX = IK*X;
  ExpFac=exp(IKX);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble *zf=(cdouble *)fval;
  zf[0] = (ExpFac - 1.0) / (IKX*X);
  zf[1] = (1.0 + (IKX-1.0)*ExpFac) / (IKX*IKX*X);
  zf[2] = Alpha*zf[1];
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
double GetAngle(double *X, double *V1, double *V2)
{  
  double A[3], B[3];
 
  VecSub(V1, X, A);
  VecSub(V2, X, B);

  return acos( VecDot(A, B) / (VecNorm(A)*VecNorm(B)) );
}


/***************************************************************/
/* compute the contributions to the scalar & vector potentials */
/* from a single panel at an observation point lying at a      */
/* panel vertex.                                               */
/*                                                             */
/* inputs:                                                     */
/*  V1=X is the observation point (a panel vertex)             */
/*  V2 and V3 are the other panel vertices                     */
/*  Q is the RWG current source vertex                         */
/*  L is the panel edge length                                 */
/*  IK = i*sqrt{Eps*Mu}*Omega  (for real frequencies)          */
/*        -sqrt{Eps*Mu}*Xi     (for imag frequencies)          */
/*                                                             */
/* outputs:                                                    */
/*  PhiA[0]    = scalar potential                              */
/*  PhiA[1..3] = cartesian components of vector potential * iw */
/***************************************************************/
void GetPotentialsAtPanelVertex(double *V1, double *V2, double *V3,
                                double *Q, double PreFac, cdouble IK, 
                                cdouble *PhiA)
{
  /***************************************************************/
  /* compute the three integrals over Alpha **********************/
  /***************************************************************/
  PAPVData MyPAPVData, *D=&(MyPAPVData);
  VecSub(V2, V1, D->A);
  VecSub(V3, V2, D->B);
  D->IK=IK;

  double Lower=0.0;
  double Upper=1.0;
  cdouble Integral[3]; 
  double Error[6];

  adapt_integrate(6, PAPVIntegrand, (void *)D, 1, &Lower, &Upper,
		  0, ABSTOL, RELTOL, (double *)Integral, Error);
  
  /***************************************************************/
  /* multiply the integrals by 2*A (A=panel area), the jacobian  */
  /* of the transformation to (u,v) variables.                   */
  /* ordinarily we skip this step because the 2*A factor is      */
  /* cancelled by the corresponding factor in the denominator    */
  /* of the RWG basis function prefactor; however, here we can't */
  /* do this because the triangle over which we are integrating  */
  /* maybe a subtriangle of the original triangle in which case  */
  /* it does not have the same area.                             */
  /***************************************************************/
  double A[3], B[3], Z[3], TwoA;
  VecSub(V2, V1, A);
  VecSub(V3, V1, B);
  VecCross(A, B, Z);
  TwoA=VecNorm(Z);
  Integral[0] *= TwoA;
  Integral[1] *= TwoA;
  Integral[2] *= TwoA;

  /***************************************************************/
  /* stamp the Alpha integrals into their appropriate places in  */ 
  /* final expressions for Phi and A.                            */ 
  /***************************************************************/
  PhiA[0] = 2.0*ZVAC * PreFac * Integral[0] / (4.0*M_PI*IK);

  int Mu;
  for(Mu=0; Mu<3; Mu++)
   PhiA[Mu+1] = IK*ZVAC*PreFac * 
                    ( (V1[Mu] - Q[Mu] ) * Integral[0] 
                     +(V2[Mu] - V1[Mu]) * Integral[1] 
                     +(V3[Mu] - V2[Mu]) * Integral[2]
                    ) / (4.0*M_PI);
  
}  

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct GPPData 
 {
   double *V1, A[3], B[3], *Q, *X;
   cdouble IK;
 
 } GPPData;

static void GPPIntegrand(unsigned ndim, const double *x, void *params, 
                         unsigned fdim, double *fval)
{
  (void )ndim; 
  (void )fdim;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  GPPData *GPPD=(GPPData *)params;
  double *V1 = GPPD->V1;
  double *A  = GPPD->A;
  double *B  = GPPD->B;
  double *Q  = GPPD->Q;
  double *X  = GPPD->X;
  cdouble IK = GPPD->IK;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double u=x[0], v=u*x[1];
  double XP[3], XPmQ[3];
  double r;
  int Mu;
  cdouble ExpFac;

  for(r=0.0, Mu=0; Mu<3; Mu++)
   { XP[Mu]    = V1[Mu] + u*A[Mu]  + v*B[Mu];
     XPmQ[Mu]  = XP[Mu] - Q[Mu];
     r        += (X[Mu]-XP[Mu])*(X[Mu]-XP[Mu]);
   };
  r=sqrt(r);
  ExpFac=exp( IK * r ) / (4.0*M_PI*r);
  
  cdouble *zfval=(cdouble *)fval;
  zfval[0] = u*ExpFac;
  zfval[1] = u*XPmQ[0] * ExpFac;
  zfval[2] = u*XPmQ[1] * ExpFac;
  zfval[3] = u*XPmQ[2] * ExpFac;
 
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetPanelPotentials(RWGObject *O, int np, int iQ, cdouble IK,
                        double *X, cdouble *PhiA)
{
  RWGPanel *P=O->Panels[np];
  cdouble PhiA1[4], PhiA2[4], PhiA3[4];
  int nv, i, Mu;
  double Length, PreFac;
  double Angle[3], AngleSum;
  
  double *V[3], *Q;

  V[0] = O->Vertices + 3*P->VI[0];
  V[1] = O->Vertices + 3*P->VI[1];
  V[2] = O->Vertices + 3*P->VI[2];
     Q = O->Vertices + 3*P->VI[iQ];

  Length=VecDistance( O->Vertices + 3*P->VI[ (iQ+1)%3 ], 
                      O->Vertices + 3*P->VI[ (iQ+2)%3 ]);

  PreFac = Length / (2.0*P->Area);

  /***************************************************************/
  /* first check whether or not the observation point lies at    */
  /* a panel vertex                                              */
  /***************************************************************/
  for(nv=0; nv<3; nv++)
   if ( VecDistance(X, V[nv]) < 1.0e-6*P->Radius )
    { GetPotentialsAtPanelVertex( V[nv], V[ (nv+1)%3 ], V[(nv+2)%3 ], Q,
                                  PreFac, IK, PhiA);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
WhichCase=1;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
      return; 
    };

  /***************************************************************/
  /* next check if the point lies on a panel edge ****************/
  /***************************************************************/
  for(AngleSum=0.0, i=0; i<3; i++)
   { 
     Angle[i]=GetAngle(X, V[i], V[(i+1)%3]);

     if ( fabs(Angle[i]-M_PI) < 1.0e-6 )
      { GetPotentialsAtPanelVertex( X, V[i],       V[(i+2)%3], Q, PreFac, IK, PhiA1);
        GetPotentialsAtPanelVertex( X, V[(i+1)%3], V[(i+2)%3], Q, PreFac, IK, PhiA2);
        for(Mu=0; Mu<4; Mu++) 
         PhiA[Mu]=PhiA1[Mu]+PhiA2[Mu];
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
WhichCase=2;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
        return;
      };

     AngleSum+=Angle[i];
   };

  /***************************************************************/
  /* next check if the point lies inside the panel ***************/
  /***************************************************************/
  if ( fabs(AngleSum-2.0*M_PI) < 1.0e-6 )
   { GetPotentialsAtPanelVertex(X, V[0], V[1], Q, PreFac, IK, PhiA1);
     GetPotentialsAtPanelVertex(X, V[0], V[2], Q, PreFac, IK, PhiA2);
     GetPotentialsAtPanelVertex(X, V[1], V[2], Q, PreFac, IK, PhiA3);
     for(Mu=0; Mu<4; Mu++) 
      PhiA[Mu]=PhiA1[Mu]+PhiA2[Mu]+PhiA3[Mu];
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
WhichCase=3;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
     return;
   };

  /***************************************************************/
  /* and if none of the above then the point lies outside the    */
  /* panel so ordinary cubature is fine.                         */
  /***************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
WhichCase=0;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  struct GPPData MyGPPData, *GPPD=&MyGPPData;
  GPPD->V1=V[0];
  VecSub(V[1], V[0], GPPD->A);
  VecSub(V[2], V[1], GPPD->B);
  GPPD->Q=Q;
  GPPD->X=X;
  GPPD->IK=IK;
  double Lower[2]={0.0, 0.0};
  double Upper[2]={1.0, 1.0};
  double Error[8];
  adapt_integrate(8, GPPIntegrand, (void *)GPPD, 2, Lower, Upper,
		  0, ABSTOL, RELTOL, (double *)PhiA, Error);

  PhiA[0] *= 2.0 * ZVAC * Length / (IK);
  PhiA[1] *= IK * ZVAC * Length;
  PhiA[2] *= IK * ZVAC * Length;
  PhiA[3] *= IK * ZVAC * Length;
  
}

/***************************************************************/
/* this is an alternate interface with a couple of differences */
/* from the above:                                             */
/*  (1) it only returns the scalar potential                   */
/*  (2) the value returned is normalized differently from the  */
/*      above: here the quantity returned is the scalar        */
/*      potential due to a unit CHARGE on the panel (as        */
/*      compared to the above, which is the potential due to   */
/*      the charge density corresponding to a unit surface-    */
/*      current density injected across one of the panel edges).*/
/*                                                             */
/***************************************************************/
cdouble GetPanelPotential(RWGObject *O, int np, cdouble IK, double *X)
{

  /* do the above computation assuming the edge in question */
  /* is the edge opposite panel vertex 0                    */
  cdouble PhiA[4];
  GetPanelPotentials(O, np, 0, IK, X, PhiA);

  /* get the length of the panel edge opposite vertex 0 */
  RWGPanel *P=O->Panels[np];
  double Length=VecDistance( O->Vertices + 3*P->VI[1], 
                             O->Vertices + 3*P->VI[2]
                           );

  /* return the appropriately renormalized scalar potential */
  return PhiA[0] * IK / Length;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
static void GPP2Integrand(unsigned ndim, const double *x, void *params, 
                          unsigned fdim, double *fval)
{
  (void )ndim;
  (void )fdim;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  GPPData *GPPD=(GPPData *)params;
  double *V1 = GPPD->V1;
  double *A  = GPPD->A;
  double *B  = GPPD->B;
  double *Q  = GPPD->Q;
  double *X  = GPPD->X;
  cdouble IK = GPPD->IK;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double u=x[0], v=u*x[1];
  double XP[3], XmXP[3], XPmQ[3];
  double r;
  int Mu;
  cdouble Phi, Psi;

  for(r=0.0, Mu=0; Mu<3; Mu++)
   { XP[Mu]    = V1[Mu] + u*A[Mu]  + v*B[Mu];
     XPmQ[Mu]  = XP[Mu] - Q[Mu];
     XmXP[Mu]  = X[Mu] - XP[Mu];
     r        += XmXP[Mu]*XmXP[Mu];
   };
  r=sqrt(r);
  Phi=exp( IK * r ) / (4.0*M_PI*r);
  Psi= (IK - 1.0/r) * Phi / r;
  
  cdouble *zfval=(cdouble *)fval;

  /***************************************************************/
  /* components of A (vector potential) (actually A/\mu_0)       */
  /***************************************************************/
  zfval[0] = u * XPmQ[0] * Phi;
  zfval[1] = u * XPmQ[1] * Phi;
  zfval[2] = u * XPmQ[2] * Phi;

  /***************************************************************/
  /* components of curl A / \mu_0 ********************************/
  /***************************************************************/
  double XxF[3];
  VecCross(XmXP, XPmQ, XxF);
  zfval[3] = u * XxF[0] * Psi;
  zfval[4] = u * XxF[1] * Psi;
  zfval[5] = u * XxF[2] * Psi;

  /***************************************************************/
  /* components of grad phi / (IZoK) *****************************/
  /* (i.e. multiply the returned value by IZoK to get grad phi)  */
  /***************************************************************/
  zfval[6] = -u * 2.0 * XmXP[0] * Psi;
  zfval[7] = -u * 2.0 * XmXP[1] * Psi;
  zfval[8] = -u * 2.0 * XmXP[2] * Psi;
 
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetPanelPotentials2(RWGObject *O, int np, int iQ, 
                         cdouble Omega, double *X, 
                         cdouble *A, cdouble *CurlA, cdouble *GradPhi)
{
  RWGPanel *P=O->Panels[np];
  double Length;
  double *V[3], *Q;

  V[0] = O->Vertices + 3*P->VI[0];
  V[1] = O->Vertices + 3*P->VI[1];
  V[2] = O->Vertices + 3*P->VI[2];
     Q = O->Vertices + 3*P->VI[iQ];

  Length=VecDistance( O->Vertices + 3*P->VI[ (iQ+1)%3 ], 
                      O->Vertices + 3*P->VI[ (iQ+2)%3 ]);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  struct GPPData MyGPPData, *GPPD=&MyGPPData;
  GPPD->V1=V[0];
  VecSub(V[1], V[0], GPPD->A);
  VecSub(V[2], V[1], GPPD->B);
  GPPD->Q=Q;
  GPPD->X=X;
  GPPD->IK=II*Omega;

  double Lower[2]={0.0, 0.0};
  double Upper[2]={1.0, 1.0};
  cdouble Potentials[9], Error[9];
  adapt_integrate(18, GPP2Integrand, (void *)GPPD, 2, Lower, Upper,
		   0, ABSTOL, RELTOL, (double *)Potentials, (double *)Error);

  A[0] = Length * Potentials[0];
  A[1] = Length * Potentials[1];
  A[2] = Length * Potentials[2];

  CurlA[0] = Length * Potentials[3];
  CurlA[1] = Length * Potentials[4];
  CurlA[2] = Length * Potentials[5];

  GradPhi[0] = Length * Potentials[6];
  GradPhi[1] = Length * Potentials[7];
  GradPhi[2] = Length * Potentials[8];
  
}
