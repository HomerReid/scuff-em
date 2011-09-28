/* 
 * TaylorMaster.cc: implementation of master formulae for Taylor's methods
 *                  for evaluating panel-panel integrals 
 *
 * homer reid       5/18/2009
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <gsl/gsl_sf.h>

#include <complex>

#include <libhrutil.h>
#include <libSGJC.h>

#include "libscuff.h"
#include "TaylorMaster.h"

#define ABSTOL 0.0
#define RELTOL 1.0e-8
#define MAXFEVALS 10000
#define INTERVALS (MAXFEVALS/15)

#define II cdouble(0,1)

/***************************************************************/
/***************************************************************/
/* X functions *************************************************/
/***************************************************************/
/***************************************************************/
// common triangle
double X_CT(TMWorkspace *TMW, int i, double x)
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
double X_CE(TMWorkspace *TMW, int i, double x1, double x2)
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
double X_CV(TMWorkspace *TMW, int i, double x1, double x2, double x3)
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

/***************************************************************/
/* integrands for x, x1x2, x1x2x3 integrals                    */
/***************************************************************/
static void xIntegrand(unsigned ndim, const double *x, void *parms,
                       unsigned nfun, double *f)
{ 
  int i, Alpha, AlphaMin, AlphaMax;
  double X, SiAlpha[7][5];
  cdouble Sum, *zf=(cdouble *)f;

  TMWorkspace *TMW=(TMWorkspace *)parms;

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
  int i, Alpha, AlphaMin, AlphaMax;
  double X, SiAlpha[7][5];
  cdouble Sum, *zf=(cdouble *)f;

  TMWorkspace *TMW=(TMWorkspace *)parms;

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
  int i, Alpha, AlphaMin, AlphaMax;
  double X, SiAlpha[7][5];
  cdouble Sum, *zf=(cdouble *)f;

  TMWorkspace *TMW=(TMWorkspace *)parms;
  
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
cdouble TaylorMaster(int WhichCase, int WhichG, int WhichH,
                     cdouble GParam, double HParam,
                     double *V1, double *V2, double *V3,
                     double *V2P, double *V3P, 
                     double *Q, double *QP, double RefVal)
{
  TMWorkspace MyTMW, *TMW=&MyTMW;

  /***************************************************************/
  /* 1. choose In function based on which kernel was requested   */
  /***************************************************************/
  switch(WhichG)
   { case TM_RP:             
       TMW->InFunc=In_RP;
       break;
     case TM_EIKR_OVER_R: 
       if ( real(GParam) == 0.0 )
        TMW->InFunc=In_EMKROverR;   
       else
        TMW->InFunc=In_EIKROverR;
       break;
     case TM_GRADEIKR_OVER_R:
       if ( real(GParam) == 0. )
        TMW->InFunc=In_GradEMKROverR;   
       else
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
       TMW->SiAlphaFunc=SiAlpha_DotPlus;
       TMW->HParam=HParam;
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
  switch(WhichCase)
   { 
     case TM_COMMONTRIANGLE:
       adapt_integrate(2, xIntegrand, (void *)TMW, 1, Lower, Upper,
                       0, RELTOL*fabs(RefVal), RELTOL,
                       (double *)&Result, (double *)&Error);

       break;

     case TM_COMMONEDGE:
       adapt_integrate(2, x1x2Integrand, (void *)TMW, 2, Lower, Upper,
                       0, RELTOL*fabs(RefVal), RELTOL,
                       (double *)&Result, (double *)&Error);
       break;

     case TM_COMMONVERTEX:
       adapt_integrate(2, x1x2x3Integrand, (void *)TMW, 3, Lower, Upper,
                       0, RELTOL*fabs(RefVal), RELTOL, 
                       (double *)&Result, (double *)&Error);

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
/* complex version of gsl's 'relative exponential' function.   */
/* this is similar to but different from the cExpRel() function*/
/* defined in LFunctions.cc.                                   */
/***************************************************************/
#define EXPRELTOL  1.0e-8
#define EXPRELTOL2 EXPRELTOL*EXPRELTOL
cdouble cExpRel2(int n, cdouble x)
{
  int m;
  cdouble Term, Sum;
  double mag2Term, mag2Sum;

  Sum=1.0;
  for(Term=1.0, m=1; m<100; m++)
   { Term*=x/((double)(m+n));
     Sum+=Term;
     mag2Term=norm(Term);
     mag2Sum=norm(Sum);
     if ( mag2Term < EXPRELTOL2*mag2Sum )
      break;
   };
  return Sum;
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

// In function for g(r) = e^{-Kappa*r}/r 
cdouble In_EMKROverR(int n, cdouble GParam, double X)
{ 
  double RetVal, Kappa = imag(GParam);
  double KX=Kappa*X;

  /* 20100610 note: if KX>40.0, then 
      exp(-KX)*gsl_sf_exprel_n(n,KX) 
     agrees with
      n! / (KX)^n
     to 12 digits (tested this for n up to 5).
     this modification essentially allows us to 
     go up to arbitrarily large values of kappa.
  */
  if ( KX > 40.0  )
   RetVal = factorial(n-1) / (X * pow(KX,n));
  else
   RetVal = exp(-KX)*gsl_sf_exprel_n(n,KX) / (X*(double)(n));
  
  return RetVal;
}

// g(r) = -(kr+1.0)*exp(-kr)/r^3
cdouble In_GradEMKROverR(int n, cdouble GParam, double R)
{ 
  double RetVal, Kappa = imag(GParam);
  double KR=Kappa*R;
  RetVal=-Kappa*exp(-KR)/(R*R)
         *(  gsl_sf_exprel_n(n-1,KR)/(    (double)(n-1) )
           + gsl_sf_exprel_n(n-2,KR)/( KR*(double)(n-2) )
          );
  return RetVal;
}

// g(r) = exp(I*K*r)/r
cdouble In_EIKROverR(int n, cdouble K, double R)
{ 
  cdouble KR=K*R;

  /* 20100610 note: if imag(KR) > 40.0, then
      \int_0^1 dw w^n e^{w*I*K*R} / (w*R)
     agrees with
      (n-1)! / ( R* (-I*K*R)^(n) )
     to high accuracy.
  */

  if ( (imag(KR)) > 40.0  )
   return factorial(n-1) / (R*pow(-II*K*R,n));
  else
   return exp(II*K*R)*cExpRel2(n,-II*K*R) / ((double)(n)*R);
}

// g(r) = (ikr-1.0)*exp(ikr)/r^3
cdouble In_GradEIKROverR(int n, cdouble K, double R)
{ 
  cdouble KR=K*R;

  return II*K*exp(II*KR)/(R*R)
          *(      cExpRel2(n-1,-II*KR)/(       (double)(n-1) )
             -1.0*cExpRel2(n-2,-II*KR)/( II*KR*(double)(n-2) )
           );
}
/*
  return (1.0fi*Kappa*cexp(I*K*R)/(R*R))
        *( cExpRel2(n-1,-1.0*1.0fi*K*R)/((double)n-1)
          -cExpRel2(n-2,-1.0*1.0fi*K*R)/((double)n-2)
                     ) / (R*R);
}
  RetVal=-Kappa*exp(-KX)/(X*X)
         *(  gsl_sf_exprel_n(n-1,KX)/(    (double)(n-1) )
           + gsl_sf_exprel_n(n-2,KX)/( KX*(double)(n-2) )
          );
*/


/***************************************************************/
/***************************************************************/
/***************************************************************/
void SiAlpha_One(const double *xVec, TMWorkspace *TMW, int WhichCase,
                 int *AlphaMin, int *AlphaMax, double S[7][5])
{ 
  double x, x1, x2, x3;

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
                 int *AlphaMin, int *AlphaMax, double S[7][5])
{
  double x, x1, x1_2, x2, x2_2, x3;
  double A2, B2, AdAP, AdB, AdBP, AdD, AdDP, APdD, BdAP, BdBP, BdDP, BPdD, DdDP;

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
       AdDP=TMW->AdDP;
       BdDP=TMW->BdDP;

       S[1][0] = (B2+3*A2+2*BdDP+4*AdDP+3*AdB)/6;
       S[1][1] = ((2*B2+3*BdDP+4*AdB)*x-4*B2-8*A2-6*BdDP-9*AdDP-10*AdB)/6;
       S[1][2] = (-B2-BdDP-2*AdB)*x+B2+A2+BdDP+AdDP+2*AdB;
       S[1][3] = ((6*B2+3*BdDP+12*AdB)*x-4*B2-2*BdDP-AdDP-6*AdB)/6;
       S[1][4] = -((2*B2+4*AdB)*x-B2+A2-AdB)/6;
       
       S[2][0] = (B2+3*A2+2*BdDP+4*AdDP+3*AdB)/6;
       S[2][1] = -((2*B2+4*A2+3*BdDP+3*AdDP+6*AdB)*x
                    +2*B2+4*A2+3*BdDP+6*AdDP+4*AdB)/6;
       S[2][2] = (B2+A2+BdDP+AdDP+2*AdB)*x;
       S[2][3] = -((6*B2+3*BdDP+3*AdDP+6*AdB)*x-2*B2-BdDP-2*AdDP)/6;
       S[2][4] = ((2*B2-2*A2)*x-B2+A2+AdB)/6;
       
       S[3][0] = (B2+3*A2+2*BdDP+4*AdDP+3*AdB)/6;
       S[3][1] = -((4*A2+3*AdDP+2*AdB)*x+2*B2+4*A2+3*BdDP+6*AdDP+4*AdB)/6;
       S[3][2] = (A2+AdDP)*x;
       S[3][3] = -((3*AdDP-6*AdB)*x-2*B2-BdDP-2*AdDP)/6;
       S[3][4] = -((2*A2+4*AdB)*x+B2-A2-AdB)/6;
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
/* this is equivalent to 'dot' + Factor * 'One'                */
/***************************************************************/
void SiAlpha_DotPlus(const double *xVec, TMWorkspace *TMW, int WhichCase,
                     int *AlphaMin, int *AlphaMax, double S[7][5])
{
   double SOne[7][5];
   int AlphaMinSOne, AlphaMaxSOne;
   double Factor = TMW->HParam;

   /***************************************************************/
   /* get the 'Dot' contributions *********************************/
   /***************************************************************/
   SiAlpha_Dot(xVec, TMW, WhichCase, AlphaMin, AlphaMax, S);

   /***************************************************************/
   /* add in the 'One' contributions ******************************/
   /***************************************************************/
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

};
               
/***************************************************************/
/***************************************************************/
/***************************************************************/
void SiAlpha_Cross(const double *xVec, TMWorkspace *TMW, int WhichCase,
                   int *AlphaMin, int *AlphaMax, double S[7][5])
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
       S[4][3] = -(((2*BxBPdQmQP+2*AxBPdQmQP)*x1_2+(-2*BxBPdQmQP+AxBdQmQP-2*AxBPdQmQP)*x1)
               *x2+(-2*BxBPdQmQP-AxBPdQmQP)*x1+2*BxBPdQmQP-AxBdQmQP+AxBPdQmQP)/2;
       
       S[5][1] = -((2*V1xAdQmQP-2*AdQxQP)*x1*x2+(2*V1xBPdQmQP-2*BPdQxQP+AxBPdQmQP)*x1
                                               +2*V1xBdQmQP-2*V1xBPdQmQP+2*LdQxQP
                                               +AxBdQmQP-AxBPdQmQP)/2;
       S[5][2] = (V1xAdQmQP+AxBdQmQP-AdQxQP)*x1*x2
               +(V1xBPdQmQP-BxBPdQmQP-BPdQxQP)*x1+V1xBdQmQP-V1xBPdQmQP+LdQxQP
               +BxBPdQmQP;
       S[5][3] = -(2*AxBdQmQP*x1*x2+(-2*BxBPdQmQP-AxBPdQmQP)*x1+2*BxBPdQmQP-AxBdQmQP
                                   +AxBPdQmQP)/2;
       
       S[6][1] = ((2*V1xAdQmQP-2*AdQxQP)*x1*x2+(2*V1xBdQmQP+2*LdQxQP-2*BPdQxQP
                                                           +AxBdQmQP)
                                               *x1-2*V1xBdQmQP+2*V1xBPdQmQP-2*LdQxQP
                                              -AxBdQmQP+AxBPdQmQP)/2;
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
                     int *AlphaMin, int *AlphaMax, double S[7][5])
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
