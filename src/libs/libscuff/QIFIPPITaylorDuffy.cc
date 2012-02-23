/* 
 * FIPPITaylorDuffy.cc: a specialized version of my master taylor-duffy 
 *                      routine that handles the special case of computing
 *                      the frequency-independent panel-panel integrals
 *
 * homer reid          11/2009 - 1/2012
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>
#include <libSGJC.h>

#include "libscuff.h"
#include "libscuffInternals.h"
#include "TaylorDuffy.h"

#define NUMGS 5
#define NUMHS 9
#define NUMFUNCS 42

#define ABSTOL 0.0
#define RELTOL 1.0e-6
#define MAXFEVALS 10000

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct FIPPITDWorkspace
 { 
   /* geometric data on triangles */
   double A[3], B[3], AP[3], BP[3]; 
   double V0xAP[3], V0xBP[3], AxV0[3], BxV0[3];
   double AxAP[3], AxBP[3], BxAP[3], BxBP[3];

   double A2, B2, AP2, BP2, L2;
   double AdB, AdAP, AdBP, AdL;
   double BdAP, BdBP, APdBP, BPdL;

   int Calls;

 } FIPPITDWorkspace;

/***************************************************************/
/***************************************************************/
/***************************************************************/
static void GetSCE(double x1, double x2, double SCE[NUMHS][7][4]);
static void GetSCE_RM3(FIPPITDWorkspace *W, 
                       double x1, double x2, double SCE[6][7][4]);

static void GetSCV(double x1, double x2, double x3, double SCV[NUMHS][3][3]);
static void GetSCV_RM3(FIPPITDWorkspace *W, 
                       double x1, double x2, double x3, double SCV[6][3][3]);

static void GetIn(double X, int AlphaMin, int AlphaMax, double I[NUMGS][6]);

/***************************************************************/
/***************************************************************/
/* X functions *************************************************/
/***************************************************************/
/***************************************************************/

// common edge
static double X_CE(FIPPITDWorkspace *W, int i, double x1, double x2)
{ 
  double u1, u2, xi2, r;

  switch(i)
   { case 1: u1=-x1;    u2=-x1*x2;       xi2=(1.0-x1+x1*x2);   break;
     case 2: u1=x1;     u2=x1*x2;        xi2=(1.0-x1);         break;
     case 3: u1=-x1*x2; u2=x1*(1.0-x2);  xi2=(1.0-x1);         break;
     case 4: u1=x1*x2;  u2=-x1*(1.0-x2); xi2=(1.0-x1*x2);      break;
     case 5: u1=-x1*x2; u2=-x1;          xi2=1.0;              break;
     case 6: u1=x1*x2;  u2=x1;           xi2=(1.0-x1);         break;
   };
  
  r=sqrt(   u1*u1*W->A2 + u2*u2*W->BP2 + xi2*xi2*W->L2 
          + 2.0*( u1*(u2*W->AdBP + xi2*W->AdL)+ u2*xi2*W->BPdL ) 
        );

  if (isnan(r)) 
   return 0.0;

  return r;

}

// common vertex 
static double X_CV(FIPPITDWorkspace *W, int i, double x1, double x2, double x3)
{ 
  double xi1, xi2, eta1, eta2, r;

  switch(i)
   { case 1: xi1=1.0; xi2=x1; eta1=x2; eta2=x2*x3; break;
     case 2: xi1=x2; xi2=x2*x3; eta1=1.0; eta2=x1; break;
   };
    
  r=sqrt(   xi1*xi1*W->A2 + xi2*xi2*W->B2
          + eta1*eta1*W->AP2 + eta2*eta2*W->BP2
          + 2.0*xi1*(xi2*W->AdB - eta1*W->AdAP - eta2*W->AdBP)
          - 2.0*xi2*(eta1*W->BdAP + eta2*W->BdBP) 
          + 2.0*eta1*eta2*W->APdBP
        );

  if (isnan(r)) 
   return 0.0;

  return r;

}

/***************************************************************/
/* integrands for x1x2, x1x2x3 integrals                       */
/* (there is no x integral in this case because we do not      */
/*  compute static panel-panel integrals for the common-       */
/*  triangle case)                                             */
/***************************************************************/
static void x1x2Integrand(unsigned ndim, const double *x, void *parms,
                          unsigned nfun, double *f)
{ 
  int i, Alpha, ri, ng, nh;
  double X, In[NUMGS][6], SCE[9][7][4], SCE_RM3[6][7][4];

  FIPPITDWorkspace *W=(FIPPITDWorkspace *)parms;
W->Calls++;
  
  GetSCE(x[0],x[1],SCE);
  GetSCE_RM3(W, x[0], x[1], SCE_RM3);

  memset(f,0,NUMFUNCS*sizeof(double));
  for(i=1; i<=6; i++)
   { 
     X=X_CE(W, i, x[0], x[1]);
     GetIn(X, 2, 5, In);

     ri=0;
     ng=0;
     for(nh=0; nh<6; nh++,ri++)
      for(Alpha=1; Alpha<=3; Alpha++)
       f[ri] += x[0] * In[ng][Alpha+2] * SCE_RM3[nh][i][Alpha];

     for(ng=1; ng<NUMGS; ng++)
      for(nh=0; nh<NUMHS; nh++,ri++)
       for(Alpha=0; Alpha<=3; Alpha++)
        f[ri] += x[0] * In[ng][Alpha+2] * SCE[nh][i][Alpha];
   };

}

static void x1x2x3Integrand(unsigned ndim, const double *x, void *parms,
                            unsigned nfun, double *f)
{ 
  int i, Alpha, ri, ng, nh;
  double X, In[NUMGS][6], SCV[9][3][3], SCV_RM3[6][3][3];

  FIPPITDWorkspace *W=(FIPPITDWorkspace *)parms;
W->Calls++;
  
  GetSCV(x[0],x[1],x[2],SCV);
  GetSCV_RM3(W, x[0], x[1], x[2], SCV_RM3);

  memset(f,0,NUMFUNCS*sizeof(double));
  for(i=1; i<=2; i++)
   { 
     X=X_CV(W, i, x[0], x[1], x[2]);
     GetIn(X, 3, 5, In);

     ri=0;
     ng=0;
     for(nh=0; nh<6; nh++,ri++)
      for(Alpha=1; Alpha<=2; Alpha++)
       f[ri] += x[1] * In[ng][Alpha+3] * SCV_RM3[nh][i][Alpha]; 

     for(ng=1; ng<NUMGS; ng++)
      for(nh=0; nh<NUMHS; nh++,ri++)
       for(Alpha=0; Alpha<=2; Alpha++)
        f[ri] += x[1] * In[ng][Alpha+3] * SCV[nh][i][Alpha]; 
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ComputeQIFIPPIData_TaylorDuffy(double *V1, double *V2, double *V3, 
                                    double *V2P, double *V3P, 
                                    QIFIPPIData *QIFD)
{
  FIPPITDWorkspace MyFIPPITDWorkspace, *W=&MyFIPPITDWorkspace;
  static double Lower[3]={0.0, 0.0, 0.0};
  static double Upper[3]={1.0, 1.0, 1.0};
  double Result[NUMFUNCS], Error[NUMFUNCS];
  double A[3], AP[3], B[3], BP[3], L[3];
  int rp;

  /***************************************************************/
  /* 1. compute geometric parameters.                            */
  /***************************************************************/
  VecSub(V2,V1,A);
  VecSub(V3,V2,B);
  VecSub(V2P,V1,AP);
  VecSub(V3P,V2P,BP);
  VecSub(BP,B,L);

  W->A2    = VecDot(A,A);
  W->B2    = VecDot(B,B);
  W->AP2   = VecDot(AP,AP);
  W->BP2   = VecDot(BP,BP);
  W->L2    = VecDot(L,L);
  W->AdB   = VecDot(A,B);
  W->AdAP  = VecDot(A,AP);
  W->AdBP  = VecDot(A,BP);
  W->AdL   = VecDot(A,L);
  W->BdAP  = VecDot(B,AP);
  W->BdBP  = VecDot(B,BP);
  W->APdBP = VecDot(AP,BP);
  W->BPdL  = VecDot(BP,L);

  memcpy(W->A,  A, 3*sizeof(double));
  memcpy(W->B,  B, 3*sizeof(double));
  memcpy(W->AP, AP, 3*sizeof(double));
  memcpy(W->BP, BP, 3*sizeof(double));
  VecCross(V1,AP,W->V0xAP);
  VecCross(V1,BP,W->V0xBP);
  VecCross(A,V1,W->AxV0);
  VecCross(B,V1,W->BxV0);
  VecCross(A,AP,W->AxAP);
  VecCross(B,AP,W->BxAP);
  VecCross(A,BP,W->AxBP);
  VecCross(B,BP,W->BxBP);

  /***************************************************************/
  /* 2. evaluate the integral over x1x2 or x1x2x3                */
  /***************************************************************/
  if( V2P==V2 ) // common-edge case 
   adapt_integrate(NUMFUNCS, x1x2Integrand, (void *)W, 2, Lower, Upper,
                   MAXFEVALS, ABSTOL, RELTOL, Result, Error);
  else // common-vertex case 
   adapt_integrate(NUMFUNCS, x1x2x3Integrand, (void *)W, 3, Lower, Upper,
                   MAXFEVALS, ABSTOL, RELTOL, Result, Error);

  /***************************************************************/
  /* 3. pack the integrals into the appropriate slots in the     */
  /*    QIFIPPIData structure.                                   */
  /***************************************************************/
  memcpy(QIFD->xMxpRM3,   Result+0, 3*sizeof(double));
  memcpy(QIFD->xXxpRM3,   Result+3, 3*sizeof(double));
  memcpy(QIFD->uvupvpRM1, Result+6, 9*sizeof(double));
  // skip the next 9 entries in the I vector ... 
  memcpy(QIFD->uvupvpR1,  Result+24, 9*sizeof(double));
  memcpy(QIFD->uvupvpR2,  Result+33, 9*sizeof(double));

}

/***************************************************************/
/* In routine **************************************************/
/***************************************************************/
static void GetIn(double X, int AlphaMin, int AlphaMax, double In[NUMGS][6])
{ 
  int A;
  double dA;

  for( A=AlphaMin; A<=AlphaMax; A++ )
   { dA=(double)A;
     In[0][A]= 1.0   / ( (dA-2.0) * X*X*X); 
     In[1][A]= 1.0   / ( dA * X );
     In[2][A]= 1.0   / ( dA + 1.0 );
     In[3][A]= X     / ( dA + 2.0 );
     In[4][A]= X*X   / ( dA + 3.0 );
   };
}

/***************************************************************/
/* content of this routine computer-generated by maxima        */
/***************************************************************/
void GetSCE(double x1, double x2, double SCE[9][7][4])
{
  SCE[0][1][0] = 1.0;
  SCE[0][1][1] = -1.0;
  SCE[0][1][2] = 0.0;
  SCE[0][1][3] = 0.0;
  SCE[0][2][0] = 1.0;
  SCE[0][2][1] = -1.0;
  SCE[0][2][2] = 0.0;
  SCE[0][2][3] = 0.0;
  SCE[0][3][0] = 1.0;
  SCE[0][3][1] = -1.0;
  SCE[0][3][2] = 0.0;
  SCE[0][3][3] = 0.0;
  SCE[0][4][0] = 1.0;
  SCE[0][4][1] = -1.0;
  SCE[0][4][2] = 0.0;
  SCE[0][4][3] = 0.0;
  SCE[0][5][0] = 1.0;
  SCE[0][5][1] = -1.0;
  SCE[0][5][2] = 0.0;
  SCE[0][5][3] = 0.0;
  SCE[0][6][0] = 1.0;
  SCE[0][6][1] = -1.0;
  SCE[0][6][2] = 0.0;
  SCE[0][6][3] = 0.0;
  SCE[1][1][0] = 1.0/2.0;
  SCE[1][1][1] = -x1;
  SCE[1][1][2] = (2.0*x1-1.0)/2.0;
  SCE[1][1][3] = 0.0;
  SCE[1][2][0] = 1.0/2.0;
  SCE[1][2][1] = 0.0;
  SCE[1][2][2] = -1.0/2.0;
  SCE[1][2][3] = 0.0;
  SCE[1][3][0] = 1.0/2.0;
  SCE[1][3][1] = -x1*x2;
  SCE[1][3][2] = (2.0*x1*x2-1.0)/2.0;
  SCE[1][3][3] = 0.0;
  SCE[1][4][0] = 1.0/2.0;
  SCE[1][4][1] = 0.0;
  SCE[1][4][2] = -1.0/2.0;
  SCE[1][4][3] = 0.0;
  SCE[1][5][0] = 1.0/2.0;
  SCE[1][5][1] = -x1*x2;
  SCE[1][5][2] = (2.0*x1*x2-1.0)/2.0;
  SCE[1][5][3] = 0.0;
  SCE[1][6][0] = 1.0/2.0;
  SCE[1][6][1] = 0.0;
  SCE[1][6][2] = -1.0/2.0;
  SCE[1][6][3] = 0.0;
  SCE[2][1][0] = 0.0;
  SCE[2][1][1] = -x1+1.0;
  SCE[2][1][2] = x1-1.0;
  SCE[2][1][3] = 0.0;
  SCE[2][2][0] = 0.0;
  SCE[2][2][1] = x1*x2-x1+1.0;
  SCE[2][2][2] = -x1*x2+x1-1.0;
  SCE[2][2][3] = 0.0;
  SCE[2][3][0] = 0.0;
  SCE[2][3][1] = -x1*x2+1.0;
  SCE[2][3][2] = x1*x2-1.0;
  SCE[2][3][3] = 0.0;
  SCE[2][4][0] = 0.0;
  SCE[2][4][1] = -x1+1.0;
  SCE[2][4][2] = x1-1.0;
  SCE[2][4][3] = 0.0;
  SCE[2][5][0] = 0.0;
  SCE[2][5][1] = -x1+1.0;
  SCE[2][5][2] = x1-1.0;
  SCE[2][5][3] = 0.0;
  SCE[2][6][0] = 0.0;
  SCE[2][6][1] = 1.0;
  SCE[2][6][2] = -1.0;
  SCE[2][6][3] = 0.0;
  SCE[3][1][0] = 1.0/2.0;
  SCE[3][1][1] = 0.0;
  SCE[3][1][2] = -1.0/2.0;
  SCE[3][1][3] = 0.0;
  SCE[3][2][0] = 1.0/2.0;
  SCE[3][2][1] = -x1;
  SCE[3][2][2] = (2.0*x1-1.0)/2.0;
  SCE[3][2][3] = 0.0;
  SCE[3][3][0] = 1.0/2.0;
  SCE[3][3][1] = 0.0;
  SCE[3][3][2] = -1.0/2.0;
  SCE[3][3][3] = 0.0;
  SCE[3][4][0] = 1.0/2.0;
  SCE[3][4][1] = -x1*x2;
  SCE[3][4][2] = (2.0*x1*x2-1.0)/2.0;
  SCE[3][4][3] = 0.0;
  SCE[3][5][0] = 1.0/2.0;
  SCE[3][5][1] = 0.0;
  SCE[3][5][2] = -1.0/2.0;
  SCE[3][5][3] = 0.0;
  SCE[3][6][0] = 1.0/2.0;
  SCE[3][6][1] = -x1*x2;
  SCE[3][6][2] = (2.0*x1*x2-1.0)/2.0;
  SCE[3][6][3] = 0.0;
  SCE[4][1][0] = 1.0/3.0;
  SCE[4][1][1] = -x1/2.0;
  SCE[4][1][2] = 0.0;
  SCE[4][1][3] = (3.0*x1-2.0)/6.0;
  SCE[4][2][0] = 1.0/3.0;
  SCE[4][2][1] = -x1/2.0;
  SCE[4][2][2] = 0.0;
  SCE[4][2][3] = (3.0*x1-2.0)/6.0;
  SCE[4][3][0] = 1.0/3.0;
  SCE[4][3][1] = -x1*x2/2.0;
  SCE[4][3][2] = 0.0;
  SCE[4][3][3] = (3.0*x1*x2-2.0)/6.0;
  SCE[4][4][0] = 1.0/3.0;
  SCE[4][4][1] = -x1*x2/2.0;
  SCE[4][4][2] = 0.0;
  SCE[4][4][3] = (3.0*x1*x2-2.0)/6.0;
  SCE[4][5][0] = 1.0/3.0;
  SCE[4][5][1] = -x1*x2/2.0;
  SCE[4][5][2] = 0.0;
  SCE[4][5][3] = (3.0*x1*x2-2.0)/6.0;
  SCE[4][6][0] = 1.0/3.0;
  SCE[4][6][1] = -x1*x2/2.0;
  SCE[4][6][2] = 0.0;
  SCE[4][6][3] = (3.0*x1*x2-2.0)/6.0;
  SCE[5][1][0] = 0.0;
  SCE[5][1][1] = -(x1-1.0)/2.0;
  SCE[5][1][2] = 0.0;
  SCE[5][1][3] = (x1-1.0)/2.0;
  SCE[5][2][0] = 0.0;
  SCE[5][2][1] = (x1*x2-x1+1.0)/2.0;
  SCE[5][2][2] = -x1*x1*x2+x1*x1-x1;
  SCE[5][2][3] = ((2.0*x1*x1-x1)*x2-2.0*x1*x1+3.0*x1-1.0)/2.0;
  SCE[5][3][0] = 0.0;
  SCE[5][3][1] = -(x1*x2-1.0)/2.0;
  SCE[5][3][2] = 0.0;
  SCE[5][3][3] = (x1*x2-1.0)/2.0;
  SCE[5][4][0] = 0.0;
  SCE[5][4][1] = -(x1-1.0)/2.0;
  SCE[5][4][2] = (x1*x1-x1)*x2;
  SCE[5][4][3] = -((2.0*x1*x1-2.0*x1)*x2-x1+1.0)/2.0;
  SCE[5][5][0] = 0.0;
  SCE[5][5][1] = -(x1-1.0)/2.0;
  SCE[5][5][2] = 0.0;
  SCE[5][5][3] = (x1-1.0)/2.0;
  SCE[5][6][0] = 0.0;
  SCE[5][6][1] = 1.0/2.0;
  SCE[5][6][2] = -x1*x2;
  SCE[5][6][3] = (2.0*x1*x2-1.0)/2.0;
  SCE[6][1][0] = 0.0;
  SCE[6][1][1] = x1*x2-x1+1.0;
  SCE[6][1][2] = -x1*x2+x1-1.0;
  SCE[6][1][3] = 0.0;
  SCE[6][2][0] = 0.0;
  SCE[6][2][1] = -x1+1.0;
  SCE[6][2][2] = x1-1.0;
  SCE[6][2][3] = 0.0;
  SCE[6][3][0] = 0.0;
  SCE[6][3][1] = -x1+1.0;
  SCE[6][3][2] = x1-1.0;
  SCE[6][3][3] = 0.0;
  SCE[6][4][0] = 0.0;
  SCE[6][4][1] = -x1*x2+1.0;
  SCE[6][4][2] = x1*x2-1.0;
  SCE[6][4][3] = 0.0;
  SCE[6][5][0] = 0.0;
  SCE[6][5][1] = 1.0;
  SCE[6][5][2] = -1.0;
  SCE[6][5][3] = 0.0;
  SCE[6][6][0] = 0.0;
  SCE[6][6][1] = -x1+1.0;
  SCE[6][6][2] = x1-1.0;
  SCE[6][6][3] = 0.0;
  SCE[7][1][0] = 0.0;
  SCE[7][1][1] = (x1*x2-x1+1.0)/2.0;
  SCE[7][1][2] = -x1*x1*x2+x1*x1-x1;
  SCE[7][1][3] = ((2.0*x1*x1-x1)*x2-2.0*x1*x1+3.0*x1-1.0)/2.0;
  SCE[7][2][0] = 0.0;
  SCE[7][2][1] = -(x1-1.0)/2.0;
  SCE[7][2][2] = 0.0;
  SCE[7][2][3] = (x1-1.0)/2.0;
  SCE[7][3][0] = 0.0;
  SCE[7][3][1] = -(x1-1.0)/2.0;
  SCE[7][3][2] = (x1*x1-x1)*x2;
  SCE[7][3][3] = -((2.0*x1*x1-2.0*x1)*x2-x1+1.0)/2.0;
  SCE[7][4][0] = 0.0;
  SCE[7][4][1] = -(x1*x2-1.0)/2.0;
  SCE[7][4][2] = 0.0;
  SCE[7][4][3] = (x1*x2-1.0)/2.0;
  SCE[7][5][0] = 0.0;
  SCE[7][5][1] = 1.0/2.0;
  SCE[7][5][2] = -x1*x2;
  SCE[7][5][3] = (2.0*x1*x2-1.0)/2.0;
  SCE[7][6][0] = 0.0;
  SCE[7][6][1] = -(x1-1.0)/2.0;
  SCE[7][6][2] = 0.0;
  SCE[7][6][3] = (x1-1.0)/2.0;
  SCE[8][1][0] = 0.0;
  SCE[8][1][1] = 0.0;
  SCE[8][1][2] = (-x1*x1+x1)*x2+x1*x1-2.0*x1+1.0;
  SCE[8][1][3] = (x1*x1-x1)*x2-x1*x1+2.0*x1-1.0;
  SCE[8][2][0] = 0.0;
  SCE[8][2][1] = 0.0;
  SCE[8][2][2] = (-x1*x1+x1)*x2+x1*x1-2.0*x1+1.0;
  SCE[8][2][3] = (x1*x1-x1)*x2-x1*x1+2.0*x1-1.0;
  SCE[8][3][0] = 0.0;
  SCE[8][3][1] = 0.0;
  SCE[8][3][2] = (x1*x1-x1)*x2-x1+1.0;
  SCE[8][3][3] = (-x1*x1+x1)*x2+x1-1.0;
  SCE[8][4][0] = 0.0;
  SCE[8][4][1] = 0.0;
  SCE[8][4][2] = (x1*x1-x1)*x2-x1+1.0;
  SCE[8][4][3] = (-x1*x1+x1)*x2+x1-1.0;
  SCE[8][5][0] = 0.0;
  SCE[8][5][1] = 0.0;
  SCE[8][5][2] = -x1+1.0;
  SCE[8][5][3] = x1-1.0;
  SCE[8][6][0] = 0.0;
  SCE[8][6][1] = 0.0;
  SCE[8][6][2] = -x1+1.0;
  SCE[8][6][3] = x1-1.0;
}

/***************************************************************/
/* content of this routine computer-generated by maxima        */
/***************************************************************/
void GetSCE_RM3(FIPPITDWorkspace *W, 
                double x1, double x2, double SCE[6][7][4])
{ 
  double *A       = W->A; 
  double *B       = W->B; 
  double *AP      = W->AP; 
  double *BP      = W->BP; 
  double *V0xAP   = W->V0xAP;
  double *V0xBP   = W->V0xBP;
  double *AxV0    = W->AxV0;
  double *BxV0    = W->BxV0;
  double *AxAP    = W->AxAP;
  double *AxBP    = W->AxBP;
  double *BxAP    = W->BxAP;
  double *BxBP    = W->BxBP;
  
  SCE[0][1][1] = B[0]*x1*x2+(BP[0]-B[0]+AP[0])*x1-BP[0]+B[0];
  SCE[0][1][2] = -(2.0*B[0]*x1*x2+(2.0*BP[0]-2.0*B[0]+2.0*AP[0])*x1-2.0*BP[0]+2.0*B[0]-AP[0] +A[0])/2.0;
  SCE[0][1][3] = 0.0;
  SCE[0][2][1] = -BP[0]*x1*x2+(BP[0]-B[0]-A[0])*x1-BP[0]+B[0];
  SCE[0][2][2] = (2.0*BP[0]*x1*x2+(-2.0*BP[0]+2.0*B[0]+2.0*A[0])*x1+2.0*BP[0]-2.0*B[0]+AP[0] -A[0])/2.0;
  SCE[0][2][3] = 0.0;
  SCE[0][3][1] = (BP[0]+AP[0])*x1*x2-B[0]*x1-BP[0]+B[0];
  SCE[0][3][2] = -((2.0*BP[0]+2.0*AP[0])*x1*x2-2.0*B[0]*x1-2.0*BP[0]+2.0*B[0]-AP[0]+A[0])/2.0;
  SCE[0][3][3] = 0.0;
  SCE[0][4][1] = (-B[0]-A[0])*x1*x2+BP[0]*x1-BP[0]+B[0];
  SCE[0][4][2] = ((2.0*B[0]+2.0*A[0])*x1*x2-2.0*BP[0]*x1+2.0*BP[0]-2.0*B[0]+AP[0]-A[0])/2.0;
  SCE[0][4][3] = 0.0;
  SCE[0][5][1] = AP[0]*x1*x2+BP[0]*x1-BP[0]+B[0];
  SCE[0][5][2] = -(2.0*AP[0]*x1*x2+2.0*BP[0]*x1-2.0*BP[0]+2.0*B[0]-AP[0]+A[0])/2.0;
  SCE[0][5][3] = 0.0;
  SCE[0][6][1] = -A[0]*x1*x2-B[0]*x1-BP[0]+B[0];
  SCE[0][6][2] = (2.0*A[0]*x1*x2+2.0*B[0]*x1+2.0*BP[0]-2.0*B[0]+AP[0]-A[0])/2.0;
  SCE[0][6][3] = 0.0;
  SCE[1][1][1] = B[1]*x1*x2+(BP[1]-B[1]+AP[1])*x1-BP[1]+B[1];
  SCE[1][1][2] = -(2.0*B[1]*x1*x2+(2.0*BP[1]-2.0*B[1]+2.0*AP[1])*x1-2.0*BP[1]+2.0*B[1]-AP[1] +A[1]) /2.0;
  SCE[1][1][3] = 0.0;
  SCE[1][2][1] = -BP[1]*x1*x2+(BP[1]-B[1]-A[1])*x1-BP[1]+B[1];
  SCE[1][2][2] = (2.0*BP[1]*x1*x2+(-2.0*BP[1]+2.0*B[1]+2.0*A[1])*x1+2.0*BP[1]-2.0*B[1]+AP[1] -A[1]) /2.0;
  SCE[1][2][3] = 0.0;
  SCE[1][3][1] = (BP[1]+AP[1])*x1*x2-B[1]*x1-BP[1]+B[1];
  SCE[1][3][2] = -((2.0*BP[1]+2.0*AP[1])*x1*x2-2.0*B[1]*x1-2.0*BP[1]+2.0*B[1]-AP[1]+A[1]) /2.0;
  SCE[1][3][3] = 0.0;
  SCE[1][4][1] = (-B[1]-A[1])*x1*x2+BP[1]*x1-BP[1]+B[1];
  SCE[1][4][2] = ((2.0*B[1]+2.0*A[1])*x1*x2-2.0*BP[1]*x1+2.0*BP[1]-2.0*B[1]+AP[1]-A[1])/2.0;
  SCE[1][4][3] = 0.0;
  SCE[1][5][1] = AP[1]*x1*x2+BP[1]*x1-BP[1]+B[1];
  SCE[1][5][2] = -(2.0*AP[1]*x1*x2+2.0*BP[1]*x1-2.0*BP[1]+2.0*B[1]-AP[1]+A[1])/2.0;
  SCE[1][5][3] = 0.0;
  SCE[1][6][1] = -A[1]*x1*x2-B[1]*x1-BP[1]+B[1];
  SCE[1][6][2] = (2.0*A[1]*x1*x2+2.0*B[1]*x1+2.0*BP[1]-2.0*B[1]+AP[1]-A[1])/2.0;
  SCE[1][6][3] = 0.0;
  SCE[2][1][1] = B[2]*x1*x2+(BP[2]-B[2]+AP[2])*x1-BP[2]+B[2];
  SCE[2][1][2] = -(2.0*B[2]*x1*x2+(2.0*BP[2]-2.0*B[2]+2.0*AP[2])*x1-2.0*BP[2]+2.0*B[2]-AP[2] +A[2]) /2.0;
  SCE[2][1][3] = 0.0;
  SCE[2][2][1] = -BP[2]*x1*x2+(BP[2]-B[2]-A[2])*x1-BP[2]+B[2];
  SCE[2][2][2] = (2.0*BP[2]*x1*x2+(-2.0*BP[2]+2.0*B[2]+2.0*A[2])*x1+2.0*BP[2]-2.0*B[2]+AP[2] -A[2]) /2.0;
  SCE[2][2][3] = 0.0;
  SCE[2][3][1] = (BP[2]+AP[2])*x1*x2-B[2]*x1-BP[2]+B[2];
  SCE[2][3][2] = -((2.0*BP[2]+2.0*AP[2])*x1*x2-2.0*B[2]*x1-2.0*BP[2]+2.0*B[2]-AP[2]+A[2]) /2.0;
  SCE[2][3][3] = 0.0;
  SCE[2][4][1] = (-B[2]-A[2])*x1*x2+BP[2]*x1-BP[2]+B[2];
  SCE[2][4][2] = ((2.0*B[2]+2.0*A[2])*x1*x2-2.0*BP[2]*x1+2.0*BP[2]-2.0*B[2]+AP[2]-A[2])/2.0;
  SCE[2][4][3] = 0.0;
  SCE[2][5][1] = AP[2]*x1*x2+BP[2]*x1-BP[2]+B[2];
  SCE[2][5][2] = -(2.0*AP[2]*x1*x2+2.0*BP[2]*x1-2.0*BP[2]+2.0*B[2]-AP[2]+A[2])/2.0;
  SCE[2][5][3] = 0.0;
  SCE[2][6][1] = -A[2]*x1*x2-B[2]*x1-BP[2]+B[2];
  SCE[2][6][2] = (2.0*A[2]*x1*x2+2.0*B[2]*x1+2.0*BP[2]-2.0*B[2]+AP[2]-A[2])/2.0;
  SCE[2][6][3] = 0.0;
  SCE[3][1][1] = ((2.0*BxV0[0]+BxAP[0])*x1*x2 +(-2.0*V0xBP[0]-2.0*V0xAP[0]-2.0*BxV0[0]-BxAP[0]-AxBP[0]-AxAP[0])*x1 +2.0*V0xBP[0]+2.0*BxV0[0]+BxAP[0]+AxBP[0]) /2.0;
  SCE[3][1][2] = -(((2.0*BxBP[0]+2.0*BxAP[0])*x1*x1+(2.0*BxV0[0]-2.0*BxBP[0])*x1)*x2 +(-2.0*BxBP[0]-2.0*BxAP[0])*x1*x1 +(-2.0*V0xBP[0]-2.0*V0xAP[0]-2.0*BxV0[0]+4*BxBP[0]+2.0*BxAP[0])*x1 +2.0*V0xBP[0]+V0xAP[0]+2.0*BxV0[0]-2.0*BxBP[0]+AxV0[0]) /2.0;
  SCE[3][1][3] = (((6.0*BxBP[0]+6.0*BxAP[0])*x1*x1+(-6.0*BxBP[0]-3.0*BxAP[0])*x1)*x2 +(-6.0*BxBP[0]-6.0*BxAP[0])*x1*x1 +(12.0*BxBP[0]+9*BxAP[0]+3.0*AxBP[0]+3.0*AxAP[0])*x1-6.0*BxBP[0] -3.0*BxAP[0]-3.0*AxBP[0]-2.0*AxAP[0]) /6.0;
  SCE[3][2][1] = ((2.0*V0xBP[0]+AxBP[0])*x1*x2 +(-2.0*V0xBP[0]-2.0*BxV0[0]-BxAP[0]-2.0*AxV0[0]-AxBP[0]-AxAP[0])*x1 +2.0*V0xBP[0]+2.0*BxV0[0]+BxAP[0]+AxBP[0]) /2.0;
  SCE[3][2][2] = -(((2.0*BxBP[0]+2.0*AxBP[0])*x1*x1+(2.0*V0xBP[0]-2.0*BxBP[0])*x1)*x2 +(-2.0*BxBP[0]-2.0*AxBP[0])*x1*x1 +(-2.0*V0xBP[0]-2.0*BxV0[0]+4*BxBP[0]-2.0*AxV0[0]+2.0*AxBP[0])*x1 +2.0*V0xBP[0]+V0xAP[0]+2.0*BxV0[0]-2.0*BxBP[0]+AxV0[0]) /2.0;
  SCE[3][2][3] = (((6.0*BxBP[0]+6.0*AxBP[0])*x1*x1+(-6.0*BxBP[0]-3.0*AxBP[0])*x1)*x2 +(-6.0*BxBP[0]-6.0*AxBP[0])*x1*x1 +(12.0*BxBP[0]+3.0*BxAP[0]+9*AxBP[0]+3.0*AxAP[0])*x1-6.0*BxBP[0] -3.0*BxAP[0]-3.0*AxBP[0]-2.0*AxAP[0]) /6.0;
  SCE[3][3][1] = -((2.0*V0xBP[0]+2.0*V0xAP[0]+AxBP[0]+AxAP[0])*x1*x2 +(2.0*BxV0[0]+BxAP[0])*x1-2.0*V0xBP[0]-2.0*BxV0[0]-BxAP[0]-AxBP[0]) /2.0;
  SCE[3][3][2] = (((2.0*BxBP[0]+2.0*BxAP[0])*x1*x1 +(2.0*V0xBP[0]+2.0*V0xAP[0]-2.0*BxBP[0]-2.0*BxAP[0])*x1) *x2 +(2.0*BxV0[0]-2.0*BxBP[0])*x1-2.0*V0xBP[0]-V0xAP[0]-2.0*BxV0[0]+2.0*BxBP[0] -AxV0[0]) /2.0;
  SCE[3][3][3] = -(((6.0*BxBP[0]+6.0*BxAP[0])*x1*x1 +(-6.0*BxBP[0]-6.0*BxAP[0]-3.0*AxBP[0]-3.0*AxAP[0])*x1) *x2 +(-6.0*BxBP[0]-3.0*BxAP[0])*x1+6.0*BxBP[0]+3.0*BxAP[0]+3.0*AxBP[0] +2.0*AxAP[0]) /6.0;
  SCE[3][4][1] = -((2.0*BxV0[0]+BxAP[0]+2.0*AxV0[0]+AxAP[0])*x1*x2 +(2.0*V0xBP[0]+AxBP[0])*x1-2.0*V0xBP[0]-2.0*BxV0[0]-BxAP[0]-AxBP[0]) /2.0;
  SCE[3][4][2] = (((2.0*BxBP[0]+2.0*AxBP[0])*x1*x1 +(2.0*BxV0[0]-2.0*BxBP[0]+2.0*AxV0[0]-2.0*AxBP[0])*x1) *x2 +(2.0*V0xBP[0]-2.0*BxBP[0])*x1-2.0*V0xBP[0]-V0xAP[0]-2.0*BxV0[0] +2.0*BxBP[0]-AxV0[0]) /2.0;
  SCE[3][4][3] = -(((6.0*BxBP[0]+6.0*AxBP[0])*x1*x1 +(-6.0*BxBP[0]-3.0*BxAP[0]-6.0*AxBP[0]-3.0*AxAP[0])*x1) *x2 +(-6.0*BxBP[0]-3.0*AxBP[0])*x1+6.0*BxBP[0]+3.0*BxAP[0]+3.0*AxBP[0] +2.0*AxAP[0]) /6.0;
  SCE[3][5][1] = -((2.0*V0xAP[0]+AxAP[0])*x1*x2+(2.0*V0xBP[0]+AxBP[0])*x1-2.0*V0xBP[0] -2.0*BxV0[0]-BxAP[0]-AxBP[0]) /2.0;
  SCE[3][5][2] = ((2.0*V0xAP[0]-2.0*BxAP[0])*x1*x2+(2.0*V0xBP[0]-2.0*BxBP[0])*x1 -2.0*V0xBP[0]-V0xAP[0]-2.0*BxV0[0] +2.0*BxBP[0]-AxV0[0]) /2.0;
  SCE[3][5][3] = ((6.0*BxAP[0]+3.0*AxAP[0])*x1*x2+(6.0*BxBP[0]+3.0*AxBP[0])*x1-6.0*BxBP[0] -3.0*BxAP[0]-3.0*AxBP[0]-2.0*AxAP[0]) /6.0;
  SCE[3][6][1] = -((2.0*AxV0[0]+AxAP[0])*x1*x2+(2.0*BxV0[0]+BxAP[0])*x1-2.0*V0xBP[0] -2.0*BxV0[0]-BxAP[0]-AxBP[0]) /2.0;
  SCE[3][6][2] = ((2.0*AxV0[0]-2.0*AxBP[0])*x1*x2+(2.0*BxV0[0]-2.0*BxBP[0])*x1 -2.0*V0xBP[0]-V0xAP[0]-2.0*BxV0[0] +2.0*BxBP[0]-AxV0[0]) /2.0;
  SCE[3][6][3] = ((6.0*AxBP[0]+3.0*AxAP[0])*x1*x2+(6.0*BxBP[0]+3.0*BxAP[0])*x1-6.0*BxBP[0] -3.0*BxAP[0]-3.0*AxBP[0]-2.0*AxAP[0]) /6.0;
  SCE[4][1][1] = ((2.0*BxV0[1]+BxAP[1])*x1*x2 +(-2.0*V0xBP[1]-2.0*V0xAP[1]-2.0*BxV0[1]-BxAP[1]-AxBP[1]-AxAP[1])*x1 +2.0*V0xBP[1]+2.0*BxV0[1]+BxAP[1]+AxBP[1]) /2.0;
  SCE[4][1][2] = -(((2.0*BxBP[1]+2.0*BxAP[1])*x1*x1+(2.0*BxV0[1]-2.0*BxBP[1])*x1)*x2 +(-2.0*BxBP[1]-2.0*BxAP[1])*x1*x1 +(-2.0*V0xBP[1]-2.0*V0xAP[1]-2.0*BxV0[1]+4*BxBP[1]+2.0*BxAP[1])*x1 +2.0*V0xBP[1]+V0xAP[1]+2.0*BxV0[1]-2.0*BxBP[1]+AxV0[1]) /2.0;
  SCE[4][1][3] = (((6.0*BxBP[1]+6.0*BxAP[1])*x1*x1+(-6.0*BxBP[1]-3.0*BxAP[1])*x1)*x2 +(-6.0*BxBP[1]-6.0*BxAP[1])*x1*x1 +(12.0*BxBP[1]+9*BxAP[1]+3.0*AxBP[1]+3.0*AxAP[1])*x1-6.0*BxBP[1] -3.0*BxAP[1]-3.0*AxBP[1]-2.0*AxAP[1]) /6.0;
  SCE[4][2][1] = ((2.0*V0xBP[1]+AxBP[1])*x1*x2 +(-2.0*V0xBP[1]-2.0*BxV0[1]-BxAP[1]-2.0*AxV0[1]-AxBP[1]-AxAP[1])*x1 +2.0*V0xBP[1]+2.0*BxV0[1]+BxAP[1]+AxBP[1]) /2.0;
  SCE[4][2][2] = -(((2.0*BxBP[1]+2.0*AxBP[1])*x1*x1+(2.0*V0xBP[1]-2.0*BxBP[1])*x1)*x2 +(-2.0*BxBP[1]-2.0*AxBP[1])*x1*x1 +(-2.0*V0xBP[1]-2.0*BxV0[1]+4*BxBP[1]-2.0*AxV0[1]+2.0*AxBP[1])*x1 +2.0*V0xBP[1]+V0xAP[1]+2.0*BxV0[1]-2.0*BxBP[1]+AxV0[1])/2.0;
  SCE[4][2][3] = (((6.0*BxBP[1]+6.0*AxBP[1])*x1*x1+(-6.0*BxBP[1]-3.0*AxBP[1])*x1)*x2 +(-6.0*BxBP[1]-6.0*AxBP[1])*x1*x1 +(12.0*BxBP[1]+3.0*BxAP[1]+9*AxBP[1]+3.0*AxAP[1])*x1-6.0*BxBP[1] -3.0*BxAP[1]-3.0*AxBP[1]-2.0*AxAP[1])/6.0;
  SCE[4][3][1] = -((2.0*V0xBP[1]+2.0*V0xAP[1]+AxBP[1]+AxAP[1])*x1*x2 +(2.0*BxV0[1]+BxAP[1])*x1-2.0*V0xBP[1]-2.0*BxV0[1]-BxAP[1]-AxBP[1]) /2.0;
  SCE[4][3][2] = (((2.0*BxBP[1]+2.0*BxAP[1])*x1*x1 +(2.0*V0xBP[1]+2.0*V0xAP[1]-2.0*BxBP[1]-2.0*BxAP[1])*x1) *x2 +(2.0*BxV0[1]-2.0*BxBP[1])*x1-2.0*V0xBP[1]-V0xAP[1]-2.0*BxV0[1]+2.0*BxBP[1] -AxV0[1]) /2.0;
  SCE[4][3][3] = -(((6.0*BxBP[1]+6.0*BxAP[1])*x1*x1 +(-6.0*BxBP[1]-6.0*BxAP[1]-3.0*AxBP[1]-3.0*AxAP[1])*x1) *x2 +(-6.0*BxBP[1]-3.0*BxAP[1])*x1+6.0*BxBP[1]+3.0*BxAP[1]+3.0*AxBP[1] +2.0*AxAP[1]) /6.0;
  SCE[4][4][1] = -((2.0*BxV0[1]+BxAP[1]+2.0*AxV0[1]+AxAP[1])*x1*x2 +(2.0*V0xBP[1]+AxBP[1])*x1-2.0*V0xBP[1]-2.0*BxV0[1]-BxAP[1]-AxBP[1]) /2.0;
  SCE[4][4][2] = (((2.0*BxBP[1]+2.0*AxBP[1])*x1*x1 +(2.0*BxV0[1]-2.0*BxBP[1]+2.0*AxV0[1]-2.0*AxBP[1])*x1) *x2 +(2.0*V0xBP[1]-2.0*BxBP[1])*x1-2.0*V0xBP[1]-V0xAP[1]-2.0*BxV0[1] +2.0*BxBP[1]-AxV0[1]) /2.0;
  SCE[4][4][3] = -(((6.0*BxBP[1]+6.0*AxBP[1])*x1*x1 +(-6.0*BxBP[1]-3.0*BxAP[1]-6.0*AxBP[1]-3.0*AxAP[1])*x1) *x2 +(-6.0*BxBP[1]-3.0*AxBP[1])*x1+6.0*BxBP[1]+3.0*BxAP[1]+3.0*AxBP[1] +2.0*AxAP[1]) /6.0;
  SCE[4][5][1] = -((2.0*V0xAP[1]+AxAP[1])*x1*x2+(2.0*V0xBP[1]+AxBP[1])*x1-2.0*V0xBP[1] -2.0*BxV0[1]-BxAP[1]-AxBP[1]) /2.0;
  SCE[4][5][2] = ((2.0*V0xAP[1]-2.0*BxAP[1])*x1*x2+(2.0*V0xBP[1]-2.0*BxBP[1])*x1 -2.0*V0xBP[1]-V0xAP[1]-2.0*BxV0[1] +2.0*BxBP[1]-AxV0[1]) /2.0;
  SCE[4][5][3] = ((6.0*BxAP[1]+3.0*AxAP[1])*x1*x2+(6.0*BxBP[1]+3.0*AxBP[1])*x1-6.0*BxBP[1] -3.0*BxAP[1]-3.0*AxBP[1]-2.0*AxAP[1]) /6.0;
  SCE[4][6][1] = -((2.0*AxV0[1]+AxAP[1])*x1*x2+(2.0*BxV0[1]+BxAP[1])*x1-2.0*V0xBP[1] -2.0*BxV0[1]-BxAP[1]-AxBP[1]) /2.0;
  SCE[4][6][2] = ((2.0*AxV0[1]-2.0*AxBP[1])*x1*x2+(2.0*BxV0[1]-2.0*BxBP[1])*x1 -2.0*V0xBP[1]-V0xAP[1]-2.0*BxV0[1] +2.0*BxBP[1]-AxV0[1]) /2.0;
  SCE[4][6][3] = ((6.0*AxBP[1]+3.0*AxAP[1])*x1*x2+(6.0*BxBP[1]+3.0*BxAP[1])*x1-6.0*BxBP[1] -3.0*BxAP[1]-3.0*AxBP[1]-2.0*AxAP[1]) /6.0;
  SCE[5][1][1] = ((2.0*BxV0[2]+BxAP[2])*x1*x2 +(-2.0*V0xBP[2]-2.0*V0xAP[2]-2.0*BxV0[2]-BxAP[2]-AxBP[2]-AxAP[2])*x1 +2.0*V0xBP[2]+2.0*BxV0[2]+BxAP[2]+AxBP[2]) /2.0;
  SCE[5][1][2] = -(((2.0*BxBP[2]+2.0*BxAP[2])*x1*x1+(2.0*BxV0[2]-2.0*BxBP[2])*x1)*x2 +(-2.0*BxBP[2]-2.0*BxAP[2])*x1*x1 +(-2.0*V0xBP[2]-2.0*V0xAP[2]-2.0*BxV0[2]+4*BxBP[2]+2.0*BxAP[2])*x1 +2.0*V0xBP[2]+V0xAP[2]+2.0*BxV0[2]-2.0*BxBP[2]+AxV0[2]) /2.0;
  SCE[5][1][3] = (((6.0*BxBP[2]+6.0*BxAP[2])*x1*x1+(-6.0*BxBP[2]-3.0*BxAP[2])*x1)*x2 +(-6.0*BxBP[2]-6.0*BxAP[2])*x1*x1 +(12.0*BxBP[2]+9*BxAP[2]+3.0*AxBP[2]+3.0*AxAP[2])*x1-6.0*BxBP[2] -3.0*BxAP[2]-3.0*AxBP[2]-2.0*AxAP[2]) /6.0;
  SCE[5][2][1] = ((2.0*V0xBP[2]+AxBP[2])*x1*x2 +(-2.0*V0xBP[2]-2.0*BxV0[2]-BxAP[2]-2.0*AxV0[2]-AxBP[2]-AxAP[2])*x1 +2.0*V0xBP[2]+2.0*BxV0[2]+BxAP[2]+AxBP[2]) /2.0;
  SCE[5][2][2] = -(((2.0*BxBP[2]+2.0*AxBP[2])*x1*x1+(2.0*V0xBP[2]-2.0*BxBP[2])*x1)*x2 +(-2.0*BxBP[2]-2.0*AxBP[2])*x1*x1 +(-2.0*V0xBP[2]-2.0*BxV0[2]+4*BxBP[2]-2.0*AxV0[2]+2.0*AxBP[2])*x1 +2.0*V0xBP[2]+V0xAP[2]+2.0*BxV0[2]-2.0*BxBP[2]+AxV0[2])/2.0;
  SCE[5][2][3] = (((6.0*BxBP[2]+6.0*AxBP[2])*x1*x1+(-6.0*BxBP[2]-3.0*AxBP[2])*x1)*x2 +(-6.0*BxBP[2]-6.0*AxBP[2])*x1*x1 +(12.0*BxBP[2]+3.0*BxAP[2]+9*AxBP[2]+3.0*AxAP[2])*x1-6.0*BxBP[2] -3.0*BxAP[2]-3.0*AxBP[2]-2.0*AxAP[2])/6.0;
  SCE[5][3][1] = -((2.0*V0xBP[2]+2.0*V0xAP[2]+AxBP[2]+AxAP[2])*x1*x2 +(2.0*BxV0[2]+BxAP[2])*x1-2.0*V0xBP[2]-2.0*BxV0[2]-BxAP[2]-AxBP[2]) /2.0;
  SCE[5][3][2] = (((2.0*BxBP[2]+2.0*BxAP[2])*x1*x1 +(2.0*V0xBP[2]+2.0*V0xAP[2]-2.0*BxBP[2]-2.0*BxAP[2])*x1) *x2 +(2.0*BxV0[2]-2.0*BxBP[2])*x1-2.0*V0xBP[2]-V0xAP[2]-2.0*BxV0[2]+2.0*BxBP[2] -AxV0[2]) /2.0;
  SCE[5][3][3] = -(((6.0*BxBP[2]+6.0*BxAP[2])*x1*x1 +(-6.0*BxBP[2]-6.0*BxAP[2]-3.0*AxBP[2]-3.0*AxAP[2])*x1) *x2 +(-6.0*BxBP[2]-3.0*BxAP[2])*x1+6.0*BxBP[2]+3.0*BxAP[2]+3.0*AxBP[2] +2.0*AxAP[2]) /6.0;
  SCE[5][4][1] = -((2.0*BxV0[2]+BxAP[2]+2.0*AxV0[2]+AxAP[2])*x1*x2 +(2.0*V0xBP[2]+AxBP[2])*x1-2.0*V0xBP[2]-2.0*BxV0[2]-BxAP[2]-AxBP[2]) /2.0;
  SCE[5][4][2] = (((2.0*BxBP[2]+2.0*AxBP[2])*x1*x1 +(2.0*BxV0[2]-2.0*BxBP[2]+2.0*AxV0[2]-2.0*AxBP[2])*x1) *x2 +(2.0*V0xBP[2]-2.0*BxBP[2])*x1-2.0*V0xBP[2]-V0xAP[2]-2.0*BxV0[2] +2.0*BxBP[2]-AxV0[2]) /2.0;
  SCE[5][4][3] = -(((6.0*BxBP[2]+6.0*AxBP[2])*x1*x1 +(-6.0*BxBP[2]-3.0*BxAP[2]-6.0*AxBP[2]-3.0*AxAP[2])*x1) *x2 +(-6.0*BxBP[2]-3.0*AxBP[2])*x1+6.0*BxBP[2]+3.0*BxAP[2]+3.0*AxBP[2] +2.0*AxAP[2]) /6.0;
  SCE[5][5][1] = -((2.0*V0xAP[2]+AxAP[2])*x1*x2+(2.0*V0xBP[2]+AxBP[2])*x1-2.0*V0xBP[2] -2.0*BxV0[2]-BxAP[2]-AxBP[2]) /2.0;
  SCE[5][5][2] = ((2.0*V0xAP[2]-2.0*BxAP[2])*x1*x2+(2.0*V0xBP[2]-2.0*BxBP[2])*x1 -2.0*V0xBP[2]-V0xAP[2]-2.0*BxV0[2] +2.0*BxBP[2]-AxV0[2]) /2.0;
  SCE[5][5][3] = ((6.0*BxAP[2]+3.0*AxAP[2])*x1*x2+(6.0*BxBP[2]+3.0*AxBP[2])*x1-6.0*BxBP[2] -3.0*BxAP[2]-3.0*AxBP[2]-2.0*AxAP[2]) /6.0;
  SCE[5][6][1] = -((2.0*AxV0[2]+AxAP[2])*x1*x2+(2.0*BxV0[2]+BxAP[2])*x1-2.0*V0xBP[2] -2.0*BxV0[2]-BxAP[2]-AxBP[2]) /2.0;
  SCE[5][6][2] = ((2.0*AxV0[2]-2.0*AxBP[2])*x1*x2+(2.0*BxV0[2]-2.0*BxBP[2])*x1 -2.0*V0xBP[2]-V0xAP[2]-2.0*BxV0[2] +2.0*BxBP[2]-AxV0[2]) /2.0;
  SCE[5][6][3] = ((6.0*AxBP[2]+3.0*AxAP[2])*x1*x2+(6.0*BxBP[2]+3.0*BxAP[2])*x1-6.0*BxBP[2] -3.0*BxAP[2]-3.0*AxBP[2]-2.0*AxAP[2]) /6.0;

  
}

/***************************************************************/
/* content of this routine computer-generated by maxima        */
/***************************************************************/
void GetSCV(double x1, double x2, double x3, double SCV[9][3][3])
{ 
  SCV[0][1][0] = 1.0;
  SCV[0][1][1] = 0.0;
  SCV[0][1][2] = 0.0;
  SCV[0][2][0] = 1.0;
  SCV[0][2][1] = 0.0;
  SCV[0][2][2] = 0.0;
  SCV[1][1][0] = 0.0;
  SCV[1][1][1] = x2;
  SCV[1][1][2] = 0.0;
  SCV[1][2][0] = 0.0;
  SCV[1][2][1] = 1.0;
  SCV[1][2][2] = 0.0;
  SCV[2][1][0] = 0.0;
  SCV[2][1][1] = x2*x3;
  SCV[2][1][2] = 0.0;
  SCV[2][2][0] = 0.0;
  SCV[2][2][1] = x1;
  SCV[2][2][2] = 0.0;
  SCV[3][1][0] = 0.0;
  SCV[3][1][1] = 1.0;
  SCV[3][1][2] = 0.0;
  SCV[3][2][0] = 0.0;
  SCV[3][2][1] = x2;
  SCV[3][2][2] = 0.0;
  SCV[4][1][0] = 0.0;
  SCV[4][1][1] = 0.0;
  SCV[4][1][2] = x2;
  SCV[4][2][0] = 0.0;
  SCV[4][2][1] = 0.0;
  SCV[4][2][2] = x2;
  SCV[5][1][0] = 0.0;
  SCV[5][1][1] = 0.0;
  SCV[5][1][2] = x2*x3;
  SCV[5][2][0] = 0.0;
  SCV[5][2][1] = 0.0;
  SCV[5][2][2] = x1*x2;
  SCV[6][1][0] = 0.0;
  SCV[6][1][1] = x1;
  SCV[6][1][2] = 0.0;
  SCV[6][2][0] = 0.0;
  SCV[6][2][1] = x2*x3;
  SCV[6][2][2] = 0.0;
  SCV[7][1][0] = 0.0;
  SCV[7][1][1] = 0.0;
  SCV[7][1][2] = x1*x2;
  SCV[7][2][0] = 0.0;
  SCV[7][2][1] = 0.0;
  SCV[7][2][2] = x2*x3;
  SCV[8][1][0] = 0.0;
  SCV[8][1][1] = 0.0;
  SCV[8][1][2] = x1*x2*x3;
  SCV[8][2][0] = 0.0;
  SCV[8][2][1] = 0.0;
  SCV[8][2][2] = x1*x2*x3;
}

/***************************************************************/
/* content of this routine computer-generated by maxima        */
/***************************************************************/
void GetSCV_RM3(FIPPITDWorkspace *W, 
                double x1, double x2, double x3, double SCV[6][3][3])
{ 
  double *A       = W->A; 
  double *B       = W->B; 
  double *AP      = W->AP; 
  double *BP      = W->BP; 
  double *V0xAP   = W->V0xAP;
  double *V0xBP   = W->V0xBP;
  double *AxV0    = W->AxV0;
  double *BxV0    = W->BxV0;
  double *AxAP    = W->AxAP;
  double *AxBP    = W->AxBP;
  double *BxAP    = W->BxAP;
  double *BxBP    = W->BxBP;

  SCV[0][1][1] = -BP[0]*x2*x3-AP[0]*x2+B[0]*x1+A[0];
  SCV[0][1][2] = 0.0;
  SCV[0][2][1] = B[0]*x2*x3+A[0]*x2-BP[0]*x1-AP[0];
  SCV[0][2][2] = 0.0;
  SCV[1][1][1] = -BP[1]*x2*x3-AP[1]*x2+B[1]*x1+A[1];
  SCV[1][1][2] = 0.0;
  SCV[1][2][1] = B[1]*x2*x3+A[1]*x2-BP[1]*x1-AP[1];
  SCV[1][2][2] = 0.0;
  SCV[2][1][1] = -BP[2]*x2*x3-AP[2]*x2+B[2]*x1+A[2];
  SCV[2][1][2] = 0.0;
  SCV[2][2][1] = B[2]*x2*x3+A[2]*x2-BP[2]*x1-AP[2];
  SCV[2][2][2] = 0.0;
  SCV[3][1][1] = V0xBP[0]*x2*x3+V0xAP[0]*x2+BxV0[0]*x1+AxV0[0];
  SCV[3][1][2] = (BxBP[0]*x1+AxBP[0])*x2*x3+(BxAP[0]*x1+AxAP[0])*x2;
  SCV[3][2][1] = BxV0[0]*x2*x3+AxV0[0]*x2+V0xBP[0]*x1+V0xAP[0];
  SCV[3][2][2] = (BxBP[0]*x1+BxAP[0])*x2*x3+(AxBP[0]*x1+AxAP[0])*x2;
  SCV[4][1][1] = V0xBP[1]*x2*x3+V0xAP[1]*x2+BxV0[1]*x1+AxV0[1];
  SCV[4][1][2] = (BxBP[1]*x1+AxBP[1])*x2*x3+(BxAP[1]*x1+AxAP[1])*x2;
  SCV[4][2][1] = BxV0[1]*x2*x3+AxV0[1]*x2+V0xBP[1]*x1+V0xAP[1];
  SCV[4][2][2] = (BxBP[1]*x1+BxAP[1])*x2*x3+(AxBP[1]*x1+AxAP[1])*x2;
  SCV[5][1][1] = V0xBP[2]*x2*x3+V0xAP[2]*x2+BxV0[2]*x1+AxV0[2];
  SCV[5][1][2] = (BxBP[2]*x1+AxBP[2])*x2*x3+(BxAP[2]*x1+AxAP[2])*x2;
  SCV[5][2][1] = BxV0[2]*x2*x3+AxV0[2]*x2+V0xBP[2]*x1+V0xAP[2];
  SCV[5][2][2] = (BxBP[2]*x1+BxAP[2])*x2*x3+(AxBP[2]*x1+AxAP[2])*x2;
}
