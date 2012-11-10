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
 * TaylorDuffy.cc: an implementation of a generalized Taylor-Duffy method
 *                 for evaluating singular panel-panel integrals 
 *
 * homer reid      5/2009 -- 2/2012
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <complex>

#include <libhrutil.h>
#include <libSGJC.h>

#include "libscuff.h"
#include "libscuffInternals.h"
#include "TaylorDuffy.h"

#ifdef USE_PTHREAD
#  include <pthread.h>
#endif

namespace scuff {

#define DEFABSTOL 0.0
//#define DEFRELTOL 1.0e-6
#define DEFRELTOL 1.0e-10
#define MAXFEVALS 10000
#define INTERVALS (MAXFEVALS/15)

#define II cdouble(0,1)

// defined at end for warning-disabling purposes
static double X_CT(TMWorkspace *TMW, int i, double x);
static double X_CE(TMWorkspace *TMW, int i, double x1, double x2);
static double X_CV(TMWorkspace *TMW, int i, double x1, double x2, double x3);

/***************************************************************/
/* integrands for x, x1x2, x1x2x3 integrals                    */
/***************************************************************/
static void xIntegrand(unsigned ndim, const double *x, void *parms,
                       unsigned nfun, double *f)
{ 
  (void) ndim; (void) nfun; // unused
  int i, Alpha, AlphaMin, AlphaMax;
  double X; 
  cdouble SiAlpha[7][5], Sum, *zf=(cdouble *)f;

  TMWorkspace *TMW=(TMWorkspace *)parms;
  TMW->nCalls++;

  TMW->SiAlphaFunc(x, TMW, TM_COMMONTRIANGLE, &AlphaMin, &AlphaMax, SiAlpha);

  for(Sum=0.0, i=1; i<=3; i++)
   { 
     X=X_CT(TMW, i, x[0]);
     for(Alpha=AlphaMin; Alpha<=AlphaMax; Alpha++)
      Sum+=TMW->InFunc(Alpha+1,TMW->GParam,X) * SiAlpha[i][Alpha];
   };
  
  *zf = Sum;
}

static void x1x2Integrand(unsigned ndim, const double *x, void *parms, unsigned nfun, double *f)
{ 
  (void) ndim; (void) nfun; // unused
  int i, Alpha, AlphaMin, AlphaMax;
  double X;
  cdouble SiAlpha[7][5], Sum, *zf=(cdouble *)f;

  TMWorkspace *TMW=(TMWorkspace *)parms;
  TMW->nCalls++;

  TMW->SiAlphaFunc(x, TMW, TM_COMMONEDGE, &AlphaMin, &AlphaMax, SiAlpha);
  for(Sum=0.0, i=1; i<=6; i++)
   { 
     X=X_CE(TMW, i, x[0], x[1]);

     for(Alpha=AlphaMin; Alpha<=AlphaMax; Alpha++)
      Sum+=TMW->InFunc(Alpha+2,TMW->GParam,X) * SiAlpha[i][Alpha];
   };

  zf[0]=x[0]*Sum;
}

static void x1x2x3Integrand(unsigned ndim, const double *x, void *parms, unsigned nfun, double *f)
{ 
  (void) ndim; (void) nfun; // unused
  int i, Alpha, AlphaMin, AlphaMax;
  double X; 
  cdouble SiAlpha[7][5], Sum, *zf=(cdouble *)f;

  TMWorkspace *TMW=(TMWorkspace *)parms;
  TMW->nCalls++;
  
  TMW->SiAlphaFunc(x, TMW, TM_COMMONVERTEX, &AlphaMin, &AlphaMax, SiAlpha);

  for(Sum=0.0, i=1; i<=2; i++)
   { 
     X=X_CV(TMW, i, x[0], x[1], x[2]);
     for(Alpha=AlphaMin; Alpha<=AlphaMax; Alpha++)
      Sum+=TMW->InFunc(Alpha+3,TMW->GParam,X) * SiAlpha[i][Alpha];
   };

  zf[0]=x[1]*Sum;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
cdouble TaylorDuffy(TaylorDuffyArgStruct *Args)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int WhichCase   = Args->WhichCase;
  int WhichG      = Args->WhichG;
  int WhichH      = Args->WhichH;
  cdouble GParam  = Args->GParam;
  double *V1      = Args->V1;
  double *V2      = Args->V2;
  double *V3      = Args->V3;
  double *V2P     = Args->V2P;
  double *V3P     = Args->V3P;
  double *Q       = Args->Q;
  double *QP      = Args->QP;
  double AbsTol   = Args->AbsTol;
  double RelTol   = Args->RelTol;
   
  TMWorkspace MyTMW, *TMW=&MyTMW;

  /***************************************************************/
  /* 1. choose In function based on which kernel was requested   */
  /***************************************************************/
  switch(WhichG)
   { case TM_RP:             
       TMW->InFunc=In_RP;
       break;
     case TM_EIKR_OVER_R: 
       TMW->InFunc=In_EIKROverR;
       break;
     case TM_GRADEIKR_OVER_R:
       TMW->InFunc=In_GradEIKROverR;
       break;
     default:
       ErrExit("%s:%i: unknown G function %i",__FILE__,__LINE__,WhichG);
   };
  TMW->GParam=GParam;

  /***************************************************************/
  /* 2. choose SiAlpha function based on which h-function was    */
  /*    requested                                                */
  /***************************************************************/
  switch(WhichH)
   { case TM_ONE:
       TMW->SiAlphaFunc=SiAlpha_One;
       break;
     case TM_DOT:
       TMW->SiAlphaFunc=SiAlpha_Dot;
       break;
     case TM_DOTPLUS:
       // note: DOTPLUS is equivalent to 'Dot' + F * 'One' 
       //  where 'F' = -4.0/(GParam*GParam);
       TMW->SiAlphaFunc=SiAlpha_DotPlus;
       break;
     case TM_CROSS:
       TMW->SiAlphaFunc=SiAlpha_Cross;
       break;
     case TM_ENORMAL:
       TMW->SiAlphaFunc=SiAlpha_ENormal;
       break;
   };

  /***************************************************************/
  /* 3. fill in geometric parameters.                            */
  /***************************************************************/
   /* the manual claims that certain parameters are not referenced    */
   /* in certain cases, which means the users might pass NULL for     */
   /* those parameters, which could cause core dumps unless we do the */
   /* following                                                       */
  if (WhichCase==TM_COMMONTRIANGLE)
   { V2P=V2; V3P=V3; } 
  else if (WhichCase==TM_COMMONEDGE)
   V2P=V2;

  double A[3], AP[3], B[3], BP[3], L[3], D[3], DP[3], QmQP[3], QxQP[3], TV[3];

  VecSub(V2,V1,A);
  VecSub(V3,V2,B);
  VecSub(V2P,V1,AP);
  VecSub(V3P,V2P,BP);
  VecSub(BP,B,L);

  TMW->A2    = VecDot(A,A);
  TMW->B2    = VecDot(B,B);
  TMW->AP2   = VecDot(AP,AP);
  TMW->BP2   = VecDot(BP,BP);
  TMW->L2    = VecDot(L,L);
  TMW->AdB   = VecDot(A,B);
  TMW->AdAP  = VecDot(A,AP);
  TMW->AdBP  = VecDot(A,BP);
  TMW->AdL   = VecDot(A,L);
  TMW->BdAP  = VecDot(B,AP);
  TMW->BdBP  = VecDot(B,BP);
  TMW->APdBP = VecDot(AP,BP);
  TMW->BPdL  = VecDot(BP,L);

  if (WhichH==TM_DOT || WhichH==TM_DOTPLUS )
   {
     VecSub(V1,Q,D);
     VecSub(V1,QP,DP);
     TMW->AdD   = VecDot(A,D);
     TMW->AdDP  = VecDot(A,DP);
     TMW->BdDP  = VecDot(B,DP);
     TMW->APdD  = VecDot(AP,D);
     TMW->BPdD  = VecDot(BP,D);
     TMW->DdDP  = VecDot(D,DP);
   }
  else if (WhichH==TM_CROSS)
   { 
     VecSub(Q,QP,QmQP);
     VecCross(Q,QP,QxQP);

     TMW->AdQxQP     = VecDot(A,  QxQP );
     TMW->APdQxQP    = VecDot(AP, QxQP );
     TMW->BdQxQP     = VecDot(B,  QxQP );
     TMW->BPdQxQP    = VecDot(BP, QxQP );
     TMW->LdQxQP     = VecDot(L,  QxQP );
     TMW->V1xAdQmQP  = VecDot( VecCross(V1,A,TV),  QmQP );
     TMW->V1xAPdQmQP = VecDot( VecCross(V1,AP,TV), QmQP );
     TMW->V1xBdQmQP  = VecDot( VecCross(V1,B,TV),  QmQP );
     TMW->V1xBPdQmQP = VecDot( VecCross(V1,BP,TV), QmQP );
     TMW->AxAPdQmQP  = VecDot( VecCross(A,AP,TV),  QmQP );
     TMW->AxBdQmQP   = VecDot( VecCross(A,B,TV),   QmQP );
     TMW->AxBPdQmQP  = VecDot( VecCross(A,BP,TV),  QmQP );
     TMW->BxAPdQmQP  = VecDot( VecCross(B,AP,TV),  QmQP );
     TMW->BxBPdQmQP  = VecDot( VecCross(B,BP,TV),  QmQP );

   }
  else if (WhichH==TM_ENORMAL)
   {
     /* note that the ENormal choice of 'h' function depends   */
     /* on a vector 'ZHat', which we pass via the 'Q' argument */
     /* which is otherwise not used for this h function        */
     TMW->APdZHat = VecDot(AP, Q);
     TMW->BPdZHat = VecDot(BP, Q);
   };

  /***************************************************************/
  /* 4. evaluate the integral over x, x1x2, or x1x2x3            */
  /***************************************************************/
  static double Lower[3]={0.0, 0.0, 0.0};
  static double Upper[3]={1.0, 1.0, 1.0};
  cdouble Result, Error;

  TMW->nCalls=0;

  switch(WhichCase)
   { 
     case TM_COMMONTRIANGLE:
       adapt_integrate(2, xIntegrand, (void *)TMW, 1, Lower, Upper,
                       0, AbsTol, RelTol, (double *)&Result, (double *)&Error);

       break;

     case TM_COMMONEDGE:
       adapt_integrate(2, x1x2Integrand, (void *)TMW, 2, Lower, Upper,
                       0, AbsTol, RelTol, (double *)&Result, (double *)&Error);
       break;

     case TM_COMMONVERTEX:

       adapt_integrate(2, x1x2x3Integrand, (void *)TMW, 3, Lower, Upper,
                       0, AbsTol, RelTol, (double *)&Result, (double *)&Error);

/*
       adapt_integrate_log(2, x1x2x3Integrand, (void *)TMW, 3, Lower, Upper,
                           0, AbsTol, RelTol, 
                           (double *)&Result, (double *)&Error,
                           "/tmp/SGJC.log",15);
*/

       break;
   };

  return Result/(4.0*M_PI);

}

/***************************************************************/
/* quick factorial function needed below                       */
/***************************************************************/
static inline double factorial(int n) 
{ 
  static double FactTable[7]={1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 620.0}; 

  if (n<=6)
   return FactTable[n];
  else
   return ((double)n)*factorial(n-1);
}

/***************************************************************/
/* this is what the GSL calls the 'relative exponential'       */
/* function. it is equal to exp(Z), minus the first n terms in */
/* the taylor expansion for exp(Z), divided by the (n+1)th     */
/* term in the taylor expansion for exp(Z). (this latter       */
/* normalization ensures that ExpRel(n,0) = 1 for all n.)      */
/***************************************************************/
#define EXPRELTOL  1.0e-8
#define EXPRELTOL2 EXPRELTOL*EXPRELTOL
cdouble ExpRel(int n, cdouble Z)
{
  int m;

  /*--------------------------------------------------------------*/
  /*- purely real case                                           -*/
  /*--------------------------------------------------------------*/
  if ( imag(Z)==0.0 )
   { 
     double Term, Sum;
     double X=real(Z);

     //////////////////////////////////////////////////
     // small-Z expansion
     //////////////////////////////////////////////////
     if ( fabs(X) < 0.1 )
      { Sum=1.0;
        for(Term=1.0, m=1; m<100; m++)
         { Term*=X/((double)(m+n));
           Sum+=Term;
           if ( fabs(Term) < EXPRELTOL*fabs(Sum) )
            break;
         };
        return Sum;
      }
     else
      { Sum=exp(X);
        for(Term=1.0, m=0; m<n; m++)
         { 
          Sum-=Term;
           Term*=X/((double)(m+1));
         };
        return Sum / Term;
      };
   }
  else
   { 
     cdouble Term, Sum;
     /*--------------------------------------------------------------*/
     /*- small-Z expansion                                          -*/
     /*--------------------------------------------------------------*/
     if ( abs(Z) < 0.1 )
      { for(Sum=Term=1.0, m=1; m<100; m++)
         { Term*=Z/((double)(m+n));
           Sum+=Term;
           if ( norm(Term) < EXPRELTOL2*norm(Sum) )
            break;
         };
        return Sum;
      }
     else
      {
        Sum=exp(Z);
        for(Term=1.0, m=0; m<n; m++)
         { 
           Sum-=Term;
           Term*=Z/((double)(m+1));
         };
        return Sum / Term;
      };
   };

} 

/***************************************************************/
/* the remainder of this file contains implementations of the  */
/* In and Si functions for the different types of kernels and  */
/* the different types of h-functions.                         */
/*                                                             */
/* note In(n,X) == \int_0^1 w^n g(wX) dw                       */ 
/***************************************************************/

// In function for g(r) = r^P
cdouble In_RP(int n, cdouble GParam, double X)
{ double RetVal; 
  double P = real(GParam);
  RetVal = pow(X, P) / ( ((double) n) + P + 1.0);
  return RetVal;
}

// In function for g(r) = exp(-K*r)/r
cdouble In_EMKROverR(int n, double K, double R)
{
  double KR=K*R;

  if ( abs(KR) < 0.1 )
   return exp(-KR)*ExpRel(n,KR) / ((double)(n)*R);
  else 
   { cdouble Sum, Term;
     int m;
     for(Sum=0.0, Term=1.0, m=0; m<n; m++)
      { 
        Sum+=Term;
        Term *= KR/((double)(m+1));
      };
     return (1.0 - exp(-K*R)*Sum) / (n*R*Term);
   };
}

// In function for g(r) = exp(I*K*r)/r
cdouble In_EIKROverR(int n, cdouble K, double R)
{ 
  cdouble KR=K*R;

  if ( (imag(KR)) > 40.0  )
   return factorial(n-1) / (R*pow(-II*K*R,n));
  else if ( real(K)==0.0 )
   return In_EMKROverR(n, imag(K), R);
  else if ( abs(KR) < 0.1 )
   return exp(II*K*R)*ExpRel(n,-II*K*R) / ((double)(n)*R);
  else 
   { cdouble Sum, Term, miKR=-II*KR;
     int m;
     for(Sum=0.0, Term=1.0, m=0; m<n; m++)
      { 
        Sum+=Term;
        Term *= miKR/((double)(m+1));
      };
     return (1.0 - exp(II*K*R)*Sum) / (n*R*Term);
   };
}

// In function for g(r) = (-KR-1) * exp(-K*r) / r^3
cdouble In_GradEMKROverR(int n, double K, double R)
{
  double KR=K*R;
  double PreFac=((double)(n-1))/((double)(n-2));

  if ( abs(KR) < 0.1 )
   return exp(-KR)*( 1.0 - PreFac*ExpRel(n-2,KR)) / (R*R*R); 
  else 
   { cdouble Sum, Term;
     int m;
     for(Sum=0.0, Term=1.0, m=0; m<n-2; m++)
      { 
        Sum+=Term;
        Term *= KR/((double)(m+1));
      };
     Sum+=Term/PreFac;
     return PreFac*(exp(-KR)*Sum - 1.0) / (Term*R*R*R);
   };
}

// In function for g(r) = (ikr-1.0)*exp(ikr)/r^3
cdouble In_GradEIKROverR(int n, cdouble K, double R)
{ 

  double PreFac=((double)(n-1))/((double)(n-2));
  cdouble KR=K*R;

  if ( (imag(KR)) > 40.0  )
   return PreFac*factorial(n-2) / (R*R*R*pow(-II*KR,n-2));
  else if ( real(K)==0.0 )
   return In_GradEMKROverR(n, imag(K), R);
  else if ( abs(KR) < 0.1 )
   return exp(II*KR)*( 1.0 - PreFac*ExpRel(n-2,-II*KR)) / (R*R*R); 
  else 
   { cdouble Sum, Term, miKR=-II*KR;
     int m;
     for(Sum=0.0, Term=1.0, m=0; m<n-2; m++)
      { 
        Sum+=Term;
        Term *= miKR/((double)(m+1));
      };
     Sum+=Term/PreFac;
     return PreFac*(exp(II*KR)*Sum - 1.0) / (Term*R*R*R);
   };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SiAlpha_One(const double *xVec, TMWorkspace *TMW, int WhichCase,
                 int *AlphaMin, int *AlphaMax, cdouble S[7][5])
{ 
  (void) xVec; (void) TMW; // unused
  switch(WhichCase)
   { 
     case TM_COMMONTRIANGLE: 
       *AlphaMin=0;
       *AlphaMax=2; 
       S[1][0]=S[2][0]=S[3][0]=1.0; 
       S[1][1]=S[2][1]=S[3][1]=-2.0;
       S[1][2]=S[2][2]=S[3][2]=1.0;
       return; 

     case TM_COMMONEDGE: 
       *AlphaMin=0;
       *AlphaMax=1; 
       S[1][0]=S[2][0]=S[3][0]=S[4][0]=S[5][0]=S[6][0]=1.0; 
       S[1][1]=S[2][1]=S[3][1]=S[4][1]=S[5][1]=S[6][1]=-1.0; 
       return;

     case TM_COMMONVERTEX: 
       *AlphaMin=0;
       *AlphaMax=0; 
       S[1][0]=S[2][0]=1.0;
       return;
   };

} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SiAlpha_Dot(const double *xVec, TMWorkspace *TMW, int WhichCase,
                 int *AlphaMin, int *AlphaMax, cdouble S[7][5])
{
  double x, x1, x1_2, x2, x2_2, x3;
  double A2, B2, AdAP, AdB, AdBP, AdD, BdD, AdDP, APdD, BdAP, BdBP, BdDP, BPdD, DdDP;

  switch(WhichCase)
   { 
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     case TM_COMMONTRIANGLE: 
       *AlphaMin=0;
       *AlphaMax=4; 
       x=xVec[0];

       A2=TMW->A2;
       B2=TMW->B2;
       AdB=TMW->AdB;
       AdD=TMW->APdD;
       AdDP=TMW->AdDP;
       BdD=TMW->BPdD;
       BdDP=TMW->BdDP;
       DdDP=TMW->DdDP;

       S[1][0] = (B2+3*A2+6*DdDP+2*BdDP+2*BdD+4*AdDP+4*AdD+3*AdB)/6.0;
       S[1][1] = ((2*B2+3*BdDP+3*BdD+4*AdB)*x-4*B2-8*A2-12*DdDP-6*BdDP-6*BdD-9*AdDP
                                              -9*AdD-10*AdB)/6.0;
       S[1][2] = (-B2-BdDP-BdD-2*AdB)*x+B2+A2+DdDP+BdDP+BdD+AdDP+AdD+2*AdB;
       S[1][3] = ((6*B2+3*BdDP+3*BdD+12*AdB)*x-4*B2-2*BdDP-2*BdD-AdDP-AdD-6*AdB)/6.0;
       S[1][4] = -((2*B2+4*AdB)*x-B2+A2-AdB)/6.0;
       
       S[2][0] = (B2+3*A2+6*DdDP+2*BdDP+2*BdD+4*AdDP+4*AdD+3*AdB)/6.0;
       S[2][1] = -((2*B2+4*A2+3*BdDP+3*BdD+3*AdDP+3*AdD+6*AdB)*x
                   +2*B2+4*A2+12*DdDP+3*BdDP+3*BdD+6*AdDP+6*AdD+4*AdB)/6.0;
       S[2][2] = (B2+A2+BdDP+BdD+AdDP+AdD+2*AdB)*x+DdDP;
       S[2][3] = -((6*B2+3*BdDP+3*BdD+3*AdDP+3*AdD+6*AdB)*x 
                   -2*B2-BdDP-BdD-2*AdDP-2*AdD)/6.0;
       S[2][4] = ((2*B2-2*A2)*x-B2+A2+AdB)/6;
       
       S[3][0] = (B2+3*A2+6*DdDP+2*BdDP+2*BdD+4*AdDP+4*AdD+3*AdB)/6;
       S[3][1] = -(  (4*A2+3*AdDP+3*AdD+2*AdB)*x
                   + 2*B2+4*A2+12*DdDP+3*BdDP+3*BdD+6*AdDP+6*AdD+4*AdB)/6.0;
       S[3][2] = (A2+AdDP+AdD)*x+DdDP;
       S[3][3] = -((3*AdDP+3*AdD-6*AdB)*x-2*B2-BdDP-BdD-2*AdDP-2*AdD)/6.0;
       S[3][4] = -((2*A2+4*AdB)*x+B2-A2-AdB)/6.0;

       return;

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     case TM_COMMONEDGE: 
       *AlphaMin=0;
       *AlphaMax=3; 

       x1=xVec[0];
       x1_2=x1*x1;
       x2=xVec[1];
       x2_2=x2*x2;

       A2=TMW->A2;
       AdB=TMW->AdB;
       AdBP=TMW->AdBP;
       AdD=TMW->AdD;
       AdDP=TMW->AdDP;
       BdBP=TMW->BdBP;
       BdDP=TMW->BdDP;
       BPdD=TMW->BPdD;
       DdDP=TMW->DdDP;

       S[1][0] = (2*A2+6*DdDP+3*AdDP+3*AdD)/6;
       S[1][1] = ((2*BdDP+AdB)*x1*x2+(-A2-2*BdDP-2*BPdD-2*AdD-AdBP-AdB)*x1
                   -2*DdDP+2*BdDP+2*BPdD+AdBP+AdB)/2;
       S[1][2] = -(((2*BdBP+2*AdB)*x1_2+(2*BdDP-2*BdBP)*x1)*x2
                 +(-2*BdBP-2*AdB)*x1_2+(-2*BdDP+4*BdBP-2*BPdD-2*AdD+2*AdB)*x1
                   +2*BdDP -2*BdBP+2*BPdD+AdDP+AdD)/2;
       S[1][3] = (((6*BdBP+6*AdB)*x1_2+(-6*BdBP-3*AdB)*x1)*x2
                 +(-6*BdBP-6*AdB)*x1_2+(3*A2+12*BdBP+3*AdBP+9*AdB)*x1
                  -2*A2-6*BdBP -3*AdBP-3*AdB)/6;
        
       S[2][0] = (2*A2+6*DdDP+3*AdDP+3*AdD)/6;
       S[2][1] = ((2*BPdD+AdBP)*x1*x2+(-A2-2*BdDP-2*BPdD-2*AdDP-AdBP-AdB)*x1
                  -2*DdDP+2*BdDP+2*BPdD+AdBP+AdB)/2;
       S[2][2] = -(((2*BdBP+2*AdBP)*x1_2+(-2*BdBP+2*BPdD)*x1)*x2
                   +(-2*BdBP-2*AdBP)*x1_2+(-2*BdDP+4*BdBP-2*BPdD-2*AdDP
                     +2*AdBP)*x1+2*BdDP-2*BdBP+2*BPdD+AdDP+AdD)/2;
       S[2][3] = (((6*BdBP+6*AdBP)*x1_2+(-6*BdBP-3*AdBP)*x1)*x2
                  +(-6*BdBP-6*AdBP)*x1_2+(3*A2+12*BdBP+9*AdBP+3*AdB)*x1
                    -2*A2-6*BdBP-3*AdBP-3*AdB)/6;
        
        
       S[3][0] = (2*A2+6*DdDP+3*AdDP+3*AdD)/6;
       S[3][1] = -((A2+2*BPdD+2*AdD+AdBP)*x1*x2+(2*BdDP+AdB)*x1+2*DdDP
                    -2*BdDP-2*BPdD -AdBP-AdB)/2;
       S[3][2] = (((2*BdBP+2*AdB)*x1_2+(-2*BdBP+2*BPdD+2*AdD-2*AdB)*x1)*x2
                   +(2*BdDP-2*BdBP)*x1-2*BdDP+2*BdBP-2*BPdD-AdDP-AdD)/2;
       S[3][3] = -(((6*BdBP+6*AdB)*x1_2+(-3*A2-6*BdBP-3*AdBP-6*AdB)*x1)*x2
                  +(-6*BdBP-3*AdB)*x1+2*A2+6*BdBP+3*AdBP+3*AdB)/6;
        
       S[4][0] = (2*A2+6*DdDP+3*AdDP+3*AdD)/6;
       S[4][1] = -((A2+2*BdDP+2*AdDP+AdB)*x1*x2+(2*BPdD+AdBP)*x1+2*DdDP
                   -2*BdDP-2*BPdD-AdBP-AdB)/2;
       S[4][2] = (((2*BdBP+2*AdBP)*x1_2+(2*BdDP-2*BdBP+2*AdDP-2*AdBP)*x1)*x2
               +(-2*BdBP+2*BPdD)*x1-2*BdDP+2*BdBP-2*BPdD-AdDP-AdD)/2;
       S[4][3] = -(((6*BdBP+6*AdBP)*x1_2+(-3*A2-6*BdBP-6*AdBP-3*AdB)*x1)*x2
               +(-6*BdBP-3*AdBP)*x1+2*A2+6*BdBP+3*AdBP+3*AdB)/6;
        
        
       S[5][0] = (2*A2+6*DdDP+3*AdDP+3*AdD)/6;
       S[5][1] = -((A2+2*AdD)*x1*x2+(2*BPdD+AdBP)*x1+2*DdDP-2*BdDP-2*BPdD
                   -AdBP-AdB)/2;
       S[5][2] = ((2*AdD-2*AdB)*x1*x2+(-2*BdBP+2*BPdD)*x1-2*BdDP+2*BdBP
                  -2*BPdD-AdDP-AdD)/2; 
       S[5][3] = ((3*A2+6*AdB)*x1*x2+(6*BdBP+3*AdBP)*x1-2*A2-6*BdBP-3*AdBP
                  -3*AdB)/6;
       
       
       S[6][0] = (2*A2+6*DdDP+3*AdDP+3*AdD)/6;
       S[6][1] = -((A2+2*AdDP)*x1*x2+(2*BdDP+AdB)*x1+2*DdDP-2*BdDP-2*BPdD
                   -AdBP-AdB)/2;
       S[6][2] = ((2*AdDP-2*AdBP)*x1*x2+(2*BdDP-2*BdBP)*x1-2*BdDP+2*BdBP
                   -2*BPdD -AdDP-AdD)/2;
       S[6][3] = ((3*A2+6*AdBP)*x1*x2+(6*BdBP+3*AdB)*x1-2*A2-6*BdBP-3*AdBP
                   -3*AdB)/6;
       
       return;
        
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     case TM_COMMONVERTEX: 
       *AlphaMin=0;
       *AlphaMax=2;

       x1=xVec[0];
       x2=xVec[1];
       x3=xVec[2];

       A2=TMW->A2;
       AdAP=TMW->AdAP;
       AdB=TMW->AdB;
       AdBP=TMW->AdBP;
       AdD=TMW->AdD;
       AdDP=TMW->AdDP;
       BdAP=TMW->BdAP;
       BdBP=TMW->BdBP;
       BdDP=TMW->BdDP;
       BPdD=TMW->BPdD;
       APdD=TMW->APdD;
       BPdD=TMW->BPdD;
       DdDP=TMW->DdDP;

       S[1][0] = DdDP;
       S[1][1] = BPdD*x2*x3+APdD*x2+BdDP*x1+AdDP;
       S[1][2] = (BdBP*x1+AdBP)*x2*x3+(BdAP*x1+AdAP)*x2;

       S[2][0] = DdDP;
       S[2][1] = BdDP*x2*x3+AdDP*x2+BPdD*x1+APdD;
       S[2][2] = (BdBP*x1+BdAP)*x2*x3+(AdBP*x1+AdAP)*x2;
              
       return;

   }; // switch(WhichCase)
}

/***************************************************************/
/* this is equivalent to 'Dot' - 4/(K*K) * 'One'               */
/* where K is the parameter that enters into the kernela       */
/* e^{i*K*r}/(4*pi*r)                                          */
/***************************************************************/
void SiAlpha_DotPlus(const double *xVec, TMWorkspace *TMW, int WhichCase,
                     int *AlphaMin, int *AlphaMax, cdouble S[7][5])
{

   /***************************************************************/
   /* get the 'Dot' contributions *********************************/
   /***************************************************************/
   SiAlpha_Dot(xVec, TMW, WhichCase, AlphaMin, AlphaMax, S);

   /***************************************************************/
   /* add in the 'One' contributions ******************************/
   /***************************************************************/
   cdouble Factor = -4.0/(TMW->GParam*TMW->GParam);
   switch(WhichCase)
     { 
       case TM_COMMONTRIANGLE: 
         S[1][0] += Factor;
         S[2][0] += Factor;
         S[3][0] += Factor;

         S[1][1] += -2.0*Factor;
         S[2][1] += -2.0*Factor;
         S[3][1] += -2.0*Factor;

         S[1][2] += Factor;
         S[2][2] += Factor;
         S[3][2] += Factor;
         break;

       case TM_COMMONEDGE: 
         S[1][0] += Factor;
         S[2][0] += Factor;
         S[3][0] += Factor;
         S[4][0] += Factor;
         S[5][0] += Factor;
         S[6][0] += Factor;

         S[1][1] -= Factor;
         S[2][1] -= Factor;
         S[3][1] -= Factor;
         S[4][1] -= Factor;
         S[5][1] -= Factor;
         S[6][1] -= Factor;
         break;

       case TM_COMMONVERTEX: 
         S[1][0] += Factor;
         S[2][0] += Factor;
         break;

     }; //switch(WhichCase)

}
               
/***************************************************************/
/***************************************************************/
/***************************************************************/
void SiAlpha_Cross(const double *xVec, TMWorkspace *TMW, int WhichCase,
                   int *AlphaMin, int *AlphaMax, cdouble S[7][5])
{
  double x1, x1_2, x2, x3;

  double AdQxQP, APdQxQP, BdQxQP, BPdQxQP, LdQxQP; 
  double V1xAdQmQP, V1xBdQmQP, V1xAPdQmQP, V1xBPdQmQP; 
  double AxAPdQmQP, BxAPdQmQP, AxBdQmQP, AxBPdQmQP, BxBPdQmQP;

  switch(WhichCase)
   { 
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     case TM_COMMONTRIANGLE: 
       *AlphaMin=0; 
       *AlphaMax=-1;
       return;

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     case TM_COMMONEDGE: 
       *AlphaMin=1; 
       *AlphaMax=3;

       x1=xVec[0];
       x1_2=x1*x1;
       x2=xVec[1];

       AdQxQP     = TMW->AdQxQP;
       BPdQxQP    = TMW->BPdQxQP;
       LdQxQP     = TMW->LdQxQP;
       V1xAdQmQP  = TMW->V1xAdQmQP;
       V1xBdQmQP  = TMW->V1xBdQmQP;
       V1xBPdQmQP = TMW->V1xBPdQmQP;
       AxBdQmQP   = TMW->AxBdQmQP;
       AxBPdQmQP  = TMW->AxBPdQmQP;
       BxBPdQmQP  = TMW->BxBPdQmQP;

       S[1][1] = -((2*V1xBdQmQP+2*LdQxQP-2*BPdQxQP+AxBdQmQP)*x1*x2
                  +(-2*V1xBdQmQP+2*V1xBPdQmQP+2*V1xAdQmQP-2*LdQxQP-AxBdQmQP
                    +AxBPdQmQP-2*AdQxQP)*x1+2*V1xBdQmQP-2*V1xBPdQmQP
                    +2*LdQxQP+AxBdQmQP-AxBPdQmQP)/2;
       S[1][2] = ((-BxBPdQmQP+AxBdQmQP)*x1_2+(V1xBdQmQP+LdQxQP+BxBPdQmQP-BPdQxQP)*x1)
                  *x2+(BxBPdQmQP-AxBdQmQP)*x1_2
                 +(-V1xBdQmQP+V1xBPdQmQP+V1xAdQmQP-LdQxQP-2*BxBPdQmQP+AxBdQmQP-AdQxQP)
                  *x1+V1xBdQmQP-V1xBPdQmQP+LdQxQP+BxBPdQmQP;
       S[1][3] = (((2*BxBPdQmQP-2*AxBdQmQP)*x1_2+(-2*BxBPdQmQP+AxBdQmQP)*x1)*x2
                 +(-2*BxBPdQmQP+2*AxBdQmQP)*x1_2+
                  (4*BxBPdQmQP-3*AxBdQmQP+AxBPdQmQP)*x1
                  -2*BxBPdQmQP+AxBdQmQP-AxBPdQmQP)/2;
       
       S[2][1] = ((2*V1xBPdQmQP-2*BPdQxQP+AxBPdQmQP)*x1*x2
               +(2*V1xBdQmQP-2*V1xBPdQmQP+2*V1xAdQmQP+2*LdQxQP+AxBdQmQP-AxBPdQmQP
                            -2*AdQxQP)
                *x1-2*V1xBdQmQP+2*V1xBPdQmQP-2*LdQxQP-AxBdQmQP+AxBPdQmQP)/2;
       S[2][2] = ((-BxBPdQmQP-AxBPdQmQP)*x1_2+(-V1xBPdQmQP+BxBPdQmQP+BPdQxQP)*x1)*x2
               +(BxBPdQmQP+AxBPdQmQP)*x1_2
               +(-V1xBdQmQP+V1xBPdQmQP-V1xAdQmQP-LdQxQP-2*BxBPdQmQP-AxBPdQmQP+AdQxQP)
                *x1+V1xBdQmQP-V1xBPdQmQP+LdQxQP+BxBPdQmQP;
       S[2][3] = (((2*BxBPdQmQP+2*AxBPdQmQP)*x1_2+(-2*BxBPdQmQP-AxBPdQmQP)*x1)*x2
               +(-2*BxBPdQmQP-2*AxBPdQmQP)*x1_2+(4*BxBPdQmQP-AxBdQmQP+3*AxBPdQmQP)*x1
               -2*BxBPdQmQP+AxBdQmQP-AxBPdQmQP)/2;
       
       S[3][1] = -((2*V1xBPdQmQP+2*V1xAdQmQP-2*BPdQxQP+AxBPdQmQP-2*AdQxQP)*x1*x2
               +(-2*V1xBdQmQP-2*LdQxQP+2*BPdQxQP-AxBdQmQP)*x1+2*V1xBdQmQP
               -2*V1xBPdQmQP+2*LdQxQP+AxBdQmQP-AxBPdQmQP)/2;
       S[3][2] = ((BxBPdQmQP-AxBdQmQP)*x1_2+(V1xBPdQmQP+V1xAdQmQP-BxBPdQmQP-BPdQxQP
                                                       +AxBdQmQP-AdQxQP)*x1)*x2
                 +(-V1xBdQmQP-LdQxQP-BxBPdQmQP+BPdQxQP)*x1+V1xBdQmQP-V1xBPdQmQP+LdQxQP
                 +BxBPdQmQP;
       S[3][3] = -(((2*BxBPdQmQP-2*AxBdQmQP)*x1_2
               +(-2*BxBPdQmQP+2*AxBdQmQP-AxBPdQmQP)*x1)*x2
               +(-2*BxBPdQmQP+AxBdQmQP)*x1+2*BxBPdQmQP-AxBdQmQP+AxBPdQmQP)/2;
       
       S[4][1] = ((2*V1xBdQmQP+2*V1xAdQmQP+2*LdQxQP-2*BPdQxQP+AxBdQmQP-2*AdQxQP)*x1*x2
               +(-2*V1xBPdQmQP+2*BPdQxQP-AxBPdQmQP)*x1-2*V1xBdQmQP+2*V1xBPdQmQP
               -2*LdQxQP-AxBdQmQP+AxBPdQmQP)/2;
       S[4][2] = ((BxBPdQmQP+AxBPdQmQP)*x1_2+(-V1xBdQmQP-V1xAdQmQP-LdQxQP-BxBPdQmQP
                                                        +BPdQxQP-AxBPdQmQP+AdQxQP)*x1)*x2
               +(V1xBPdQmQP-BxBPdQmQP-BPdQxQP)*x1+V1xBdQmQP-V1xBPdQmQP+LdQxQP+BxBPdQmQP;
       S[4][3] = -(  ((2*BxBPdQmQP+2*AxBPdQmQP)*x1_2 
                      +(-2*BxBPdQmQP+AxBdQmQP-2*AxBPdQmQP)*x1)*x2
                    +(-2*BxBPdQmQP-AxBPdQmQP)*x1+2*BxBPdQmQP-AxBdQmQP+AxBPdQmQP)/2;
       
       S[5][1] = -((2*V1xAdQmQP-2*AdQxQP)*x1*x2+(2*V1xBPdQmQP-2*BPdQxQP+AxBPdQmQP)*x1
                                               +2*V1xBdQmQP-2*V1xBPdQmQP+2*LdQxQP
                                               +AxBdQmQP-AxBPdQmQP)/2;
       S[5][2] = (V1xAdQmQP+AxBdQmQP-AdQxQP)*x1*x2
               +(V1xBPdQmQP-BxBPdQmQP-BPdQxQP)*x1+V1xBdQmQP-V1xBPdQmQP+LdQxQP
               +BxBPdQmQP;
       S[5][3] = -(2*AxBdQmQP*x1*x2+(-2*BxBPdQmQP-AxBPdQmQP)*x1+2*BxBPdQmQP-AxBdQmQP
                                   +AxBPdQmQP)/2;
       
       S[6][1] = ( (2*V1xAdQmQP-2*AdQxQP)*x1*x2
                  +(2*V1xBdQmQP+2*LdQxQP-2*BPdQxQP +AxBdQmQP)*x1
                  - 2*V1xBdQmQP+2*V1xBPdQmQP-2*LdQxQP-AxBdQmQP+AxBPdQmQP)/2;

       S[6][2] = (-V1xAdQmQP-AxBPdQmQP+AdQxQP)*x1*x2
                 +(-V1xBdQmQP-LdQxQP-BxBPdQmQP+BPdQxQP)*x1+V1xBdQmQP-V1xBPdQmQP+LdQxQP
                 +BxBPdQmQP;
       S[6][3] = (2*AxBPdQmQP*x1*x2+(2*BxBPdQmQP-AxBdQmQP)*x1-2*BxBPdQmQP+AxBdQmQP
                                   -AxBPdQmQP)/2;

       return;

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     case TM_COMMONVERTEX: 

       *AlphaMin=1; 
       *AlphaMax=2;

       x1=xVec[0];
       x2=xVec[1];
       x3=xVec[2];

       AdQxQP     = TMW->AdQxQP;
       APdQxQP    = TMW->APdQxQP;
       BdQxQP     = TMW->BdQxQP;
       BPdQxQP    = TMW->BPdQxQP;
       V1xAdQmQP  = TMW->V1xAdQmQP;
       V1xAPdQmQP = TMW->V1xAPdQmQP;
       V1xBdQmQP  = TMW->V1xBdQmQP;
       V1xBPdQmQP = TMW->V1xBPdQmQP;
       AxAPdQmQP  = TMW->AxAPdQmQP;
       AxBPdQmQP  = TMW->AxBPdQmQP;
       BxAPdQmQP  = TMW->BxAPdQmQP;
       BxBPdQmQP  = TMW->BxBPdQmQP;

       S[1][1] = (V1xBPdQmQP-BPdQxQP)*x2*x3+(V1xAPdQmQP-APdQxQP)*x2
                  +(-V1xBdQmQP+BdQxQP)*x1-V1xAdQmQP+AdQxQP;
       S[1][2] = (BxBPdQmQP*x1+AxBPdQmQP)*x2*x3+(BxAPdQmQP*x1+AxAPdQmQP)*x2;

       S[2][1] = (-V1xBdQmQP+BdQxQP)*x2*x3+(-V1xAdQmQP+AdQxQP)*x2
                  +(V1xBPdQmQP-BPdQxQP)*x1+V1xAPdQmQP-APdQxQP;
       S[2][2] = (BxBPdQmQP*x1+BxAPdQmQP)*x2*x3+(AxBPdQmQP*x1+AxAPdQmQP)*x2;

       return;

   };
}
               
/***************************************************************/
/***************************************************************/
/***************************************************************/
void SiAlpha_ENormal(const double *xVec, TMWorkspace *TMW, int WhichCase,
                     int *AlphaMin, int *AlphaMax, cdouble S[7][5])
{
  double x1, x2, x3;
  double APdZHat, BPdZHat;

  switch(WhichCase)
   { 
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     case TM_COMMONTRIANGLE: 
       *AlphaMin=0; 
       *AlphaMax=-1;
       return;

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     case TM_COMMONEDGE: 
       *AlphaMin=1; 
       *AlphaMax=2;

       x1=xVec[0];
       x2=xVec[1];

       BPdZHat  = TMW->BPdZHat;  

       S[1][1] = -(BPdZHat*(1.0 - x1));
       S[1][2] = -(BPdZHat*(-1.0 + x1));
       S[2][1] = -(BPdZHat*(1.0 - x1 + x1*x2));
       S[2][2] = -(BPdZHat*(-1.0 + x1 - x1*x2));
       S[3][1] = -(BPdZHat*(1.0 - x1*x2));
       S[3][2] = -(BPdZHat*(-1.0 + x1*x2));
       S[4][1] = -(BPdZHat*(1.0 - x1));
       S[4][2] = -(BPdZHat*(-1.0 + x1));
       S[5][1] = -(BPdZHat*(1.0 - x1));
       S[5][2] = -(BPdZHat*(-1.0 + x1));
       S[6][1] = -BPdZHat;
       S[6][2] = BPdZHat;
       return;
       
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     case TM_COMMONVERTEX: 

       *AlphaMin=1; 
       *AlphaMax=1;

       x1=xVec[0];
       x2=xVec[1];
       x3=xVec[2];

       APdZHat = TMW->APdZHat;
       BPdZHat = TMW->BPdZHat;

       S[1][1] = -APdZHat*x2 - BPdZHat*x2*x3;
       S[2][1] = -APdZHat - BPdZHat*x1;

       return;

   };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void InitTaylorDuffyArgs(TaylorDuffyArgStruct *Args)
{ 
  Args->Q=0;
  Args->QP=0;
  Args->AbsTol=DEFABSTOL;
  Args->RelTol=DEFRELTOL;
}

/***************************************************************/
/***************************************************************/
/* X functions *************************************************/
/***************************************************************/
/***************************************************************/

// Disable spurious gcc warnings for variables initialized in case statements.
// (These functions are at the end of the file so that we only disable the 
//  warnings for these three functions.)
#pragma GCC diagnostic ignored "-Wuninitialized"

// common triangle
static double X_CT(TMWorkspace *TMW, int i, double x)
{
  double u1, u2;

  switch(i)
   { case 1: u1=1.0;  u2=x;     break;
     case 2: u1=x;    u2=(x-1); break;
     case 3: u1=x;    u2=1.0;   break;
   };
  return sqrt( TMW->A2*u1*u1 + 2.0*TMW->AdB*u1*u2 + TMW->B2*u2*u2);
}

// common edge
static double X_CE(TMWorkspace *TMW, int i, double x1, double x2)
{ 
  double u1, u2, xi2;

  switch(i)
   { case 1: u1=-x1;    u2=-x1*x2;       xi2=(1.0-x1+x1*x2);   break;
     case 2: u1=x1;     u2=x1*x2;        xi2=(1.0-x1);         break;
     case 3: u1=-x1*x2; u2=x1*(1.0-x2);  xi2=(1.0-x1);         break;
     case 4: u1=x1*x2;  u2=-x1*(1.0-x2); xi2=(1.0-x1*x2);      break;
     case 5: u1=-x1*x2; u2=-x1;          xi2=1.0;              break;
     case 6: u1=x1*x2;  u2=x1;           xi2=(1.0-x1);         break;
   };
  
  return sqrt(  u1*u1*TMW->A2 + u2*u2*TMW->BP2 + xi2*xi2*TMW->L2 
              + 2.0*(u1*(u2*TMW->AdBP + xi2*TMW->AdL) + u2*xi2*TMW->BPdL) );
}

// common vertex 
static double X_CV(TMWorkspace *TMW, int i, double x1, double x2, double x3)
{ 
  double xi1, xi2, eta1, eta2;

  switch(i)
   { case 1: xi1=1.0; xi2=x1; eta1=x2; eta2=x2*x3; break;
     case 2: xi1=x2; xi2=x2*x3; eta1=1.0; eta2=x1; break;
   };
    
  return sqrt(     xi1*xi1*TMW->A2    +   xi2*xi2*TMW->B2 
               + eta1*eta1*TMW->AP2   + eta2*eta2*TMW->BP2
               + 2.0*xi1*(xi2*TMW->AdB - eta1*TMW->AdAP - eta2*TMW->AdBP)
               - 2.0*xi2*(eta1*TMW->BdAP + eta2*TMW->BdBP) 
               + 2.0*eta1*eta2*TMW->APdBP
             );
}

} // namespace scuff
