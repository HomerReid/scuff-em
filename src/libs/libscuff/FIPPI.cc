/*
 * FIPPI.cc     -- libscuff routines for working with frequency-independent
 *                 panel-panel integrals (FIPPIs)
 * 
 * homer reid   -- 11/2005 -- 11/2011
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>
#include <libSGJC.h>

#include "libscuff.h"
#include "libscuffInternals.h"

#define ABSTOL 1.0e-10
#define RELTOL 1.0e-6

#define _ONE 0
#define _UP  1
#define _VP  2
#define _U   3
#define _UUP 4
#define _UVP 5
#define _V   6
#define _VUP 7
#define _VVP 8

#define uvupvpR0_ONE (1.0/4.0)
#define uvupvpR0_UP  (1.0/6.0)
#define uvupvpR0_VP  (1.0/12.0)
#define uvupvpR0_U   (1.0/6.0)
#define uvupvpR0_UUP (1.0/9.0)
#define uvupvpR0_UVP (1.0/18.0)
#define uvupvpR0_V   (1.0/12.0)
#define uvupvpR0_VUP (1.0/18.0)
#define uvupvpR0_VVP (1.0/36.0)

void ComputeQIFIPPIData_TaylorDuffy(double *V1, double *V2, double *V3, 
                                    double *V2P, double *V3P, 
                                    QIFIPPIData *QIFD);

/***************************************************************/
/* 'vertex less than.' returns 1 if V1<V2, 0 otherwise.        */
/* vertices are sorted using a fairly obvious sorting scheme.  */
/***************************************************************/
static int VLT(double *V1, double *V2)
{
  double DV;

  DV=V1[0]-V2[0];
  if ( fabs(DV) > 1.0e-6*fabs(V1[0]) )
   return DV<0.0 ? 1 : 0;

  DV=V1[1]-V2[1];
  if ( fabs(DV) > 1.0e-6*fabs(V1[1]) )
   return DV<0.0 ? 1 : 0;
 
  DV=V1[2]-V2[2];
  return DV<0.0 ? 1 : 0;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int CanonicallyOrderVertices(double **Va, double *Qa,
                             double **Vb, double *Qb,
                             double **OVa, double **OQa,
                             double **OVb, double **OQb)
{ 
  /***************************************************************/
  /* if there are any common vertices, then we define the        */
  /* canonical ordering to be simply that computed by            */
  /* AssessPanelPair, albeit possibly with the panels swapped.   */
  /***************************************************************/
  OVa[0]=Va[0]; OVa[1]=Va[1]; OVa[2]=Va[2];
  OVb[0]=Vb[0]; OVb[1]=Vb[1]; OVb[2]=Vb[2];
  int ncv=AssessPanelPair(OVa, OVb);
  if ( ncv > 0 )
   { if ( VLT(OVa[0], OVb[0]) )
      { *OQa=Qa;
        *OQb=Qb;
        return 0;
      }
     else
      { double *TV;
        TV=OVa[0]; OVa[0]=OVb[0]; OVb[0]=TV;
        TV=OVa[1]; OVa[1]=OVb[1]; OVb[1]=TV;
        TV=OVa[2]; OVa[2]=OVb[2]; OVb[2]=TV;
        *OQa=Qb;
        *OQb=Qa;
        return 1;
      };
   };

  /***************************************************************/
  /* otherwise, sort the vertices in ascending order. ************/
  /***************************************************************/
  double *VMin, *VMed, *VMax;
  if ( VLT(Va[0], Va[1]) )
   { VMin=Va[0]; VMax=Va[1]; }
  else
   { VMin=Va[1]; VMax=Va[0]; }

  if ( VLT(Va[2], VMin) )
   { VMed=VMin; VMin=Va[2]; }
  else if ( VLT(Va[2], VMax) )
   { VMed=Va[2]; }
  else 
   { VMed=VMax; VMax=Va[2]; }

  double *VMinP, *VMedP, *VMaxP;
  if ( VLT(Vb[0], Vb[1]) )
   { VMinP=Vb[0]; VMaxP=Vb[1]; }
  else
   { VMinP=Vb[1]; VMaxP=Vb[0]; }

  if ( VLT(Vb[2], VMinP) )
   { VMedP=VMinP; VMinP=Vb[2]; }
  else if ( VLT(Vb[2], VMaxP) )
   { VMedP=Vb[2]; }
  else 
   { VMedP=VMaxP; VMaxP=Vb[2]; }

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if ( VLT(Va[0], Vb[0] ) )
   { 
     OVa[0]=VMin; 
     OVa[1]=VMed; 
     OVa[2]=VMax; 
     *OQa=Qa;

     OVb[0]=VMinP; 
     OVb[1]=VMedP; 
     OVb[2]=VMaxP; 
     *OQb=Qb;

     return 0;
   }
  else
   {
     OVa[0]=VMinP; 
     OVa[1]=VMedP; 
     OVa[2]=VMaxP; 
     *OQa=Qb;

     OVb[0]=VMin; 
     OVb[1]=VMed; 
     OVb[2]=VMax; 
     *OQb=Qa;

     return 1;
   };
  
}

/***************************************************************/
/* implementation of the FIPPIDT class *************************/
/***************************************************************/

/*--------------------------------------------------------------*/
/*- class constructor ------------------------------------------*/
/*--------------------------------------------------------------*/
FIPPIDataTable::FIPPIDataTable()
{
}

/*--------------------------------------------------------------*/
/*- class destructor  ------------------------------------------*/
/*--------------------------------------------------------------*/
FIPPIDataTable::~FIPPIDataTable()
{
} 

/*--------------------------------------------------------------*/
/*- routine for fetching a FIPPI data record from a FIPPIDT:    */
/*- we look through our table to see if we have a record for    */
/*- this panel pair, we return it if we do, and otherwise we    */
/*- compute a new FIPPI data record for this panel pair and     */
/*- add it to the table                                         */
/*- important note: the vertices are assumed to be canonically  */
/*- ordered on entry.                                           */
/*--------------------------------------------------------------*/
QIFIPPIData *FIPPIDataTable::GetQIFIPPIData(double **OVa, double **OVb)
{

  /***************************************************************/
  /* the search key is kinda stupid, just a string of 15 doubles */
  /* as follows:                                                 */
  /* 0--2    VMed  - VMin [0..2]                                 */
  /* 3--5    VMax  - VMin [0..2]                                 */
  /* 6--8    VMinP - VMin [0..2]                                 */
  /* 9--11   VMedP - VMin [0..2]                                 */
  /* 12--14  VMaxP - VMin [0..2]                                 */
  /***************************************************************/
  double Key[15];
  VecSub(OVa[1], OVa[0], Key+0 );
  VecSub(OVa[2], OVa[0], Key+3 );
  VecSub(OVb[0], OVa[0], Key+6 );
  VecSub(OVb[1], OVa[0], Key+9 );
  VecSub(OVb[2], OVa[0], Key+12);
  
}

/*--------------------------------------------------------------*/
/*- integrand routine used for evaluating FIPPIs by cubature.  -*/
/*- note: CFD = 'compute FIPPI data.'                          -*/
/*--------------------------------------------------------------*/
typedef struct CFDData
 {
   double *V0, A[3], B[3];
   double *V0P, AP[3], BP[3];
   int nCalls;
 } CFDData;

void CFDIntegrand4D(unsigned ndim, const double *x, void *params,
                    unsigned fdim, double *fval)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  CFDData *CFDD=(CFDData *)params;
  CFDD->nCalls++;

  double *V0=CFDD->V0;
  double *A=CFDD->A;
  double *B=CFDD->B;

  double *V0P=CFDD->V0P;
  double *AP=CFDD->AP;
  double *BP=CFDD->BP;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double u, v, up, vp, uup;
  u=x[0];
  v=u*x[1];
  up=x[2];
  vp=up*x[3];
  uup=u*up; // jacobian

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double X[3], XP[3], R[3], XxXP[3];
  double r, r2=0.0;
  int Mu;
  for (Mu=0; Mu<3; Mu++)
   { X[Mu]  = V0[Mu]  + u*A[Mu]   + v*B[Mu];
     XP[Mu] = V0P[Mu] + up*AP[Mu] + vp*BP[Mu];
     R[Mu]  = X[Mu] - XP[Mu];
     r2    += R[Mu]*R[Mu];
   };
  r=sqrt(r2);
  VecCross(X,XP,XxXP);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int nf=0;
  double oor, oor3;

  // we put the jacobian factors into these quantities for convenience
  oor=u*up/r;
  oor3=oor/r2;
  r*=u*up;
  r2*=u*up;

  fval[nf++] = R[0] * oor3;
  fval[nf++] = R[1] * oor3;
  fval[nf++] = R[2] * oor3;
  fval[nf++] = XxXP[0] * oor3;
  fval[nf++] = XxXP[1] * oor3;
  fval[nf++] = XxXP[2] * oor3;

  fval[nf++] = oor;
  fval[nf++] = up*oor;
  fval[nf++] = vp*oor;
  fval[nf++] = u*oor;
  fval[nf++] = u*up*oor;
  fval[nf++] = u*vp*oor;
  fval[nf++] = v*oor;
  fval[nf++] = v*up*oor;
  fval[nf++] = v*vp*oor;

  fval[nf++] = r;
  fval[nf++] = up*r;
  fval[nf++] = vp*r;
  fval[nf++] = u*r;
  fval[nf++] = u*up*r;
  fval[nf++] = u*vp*r;
  fval[nf++] = v*r;
  fval[nf++] = v*up*r;
  fval[nf++] = v*vp*r;

  fval[nf++] = r2;
  fval[nf++] = up*r2;
  fval[nf++] = vp*r2;
  fval[nf++] = u*r2;
  fval[nf++] = u*up*r2;
  fval[nf++] = u*vp*r2;
  fval[nf++] = v*r2;
  fval[nf++] = v*up*r2;
  fval[nf++] = v*vp*r2;

} 

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
void CFDIntegrand3D(unsigned ndim, const double *x, void *params,
                    unsigned fdim, double *fval)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  CFDData *CFDD=(CFDData *)params;
  CFDD->nCalls++;

  double *V0=CFDD->V0;
  double *A=CFDD->A;
  double *B=CFDD->B;

  double *V0P=CFDD->V0P;
  double *AP=CFDD->AP;
  double *BP=CFDD->BP;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double u, v, up;
  u=x[0];
  v=u*x[1];
  up=x[2];
  
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int Mu;
  double X[3], Y[3], Z[3], BPdY, Y2, a2, a, a3, vp0, vp02, vp03, vp04, b2; 
  a2=BPdY=Y2=0.0;
  for (Mu=0; Mu<3; Mu++)
   { 
     //Y[Mu]  = V0[Mu] - V0P[Mu] + u*A[Mu] + v*B[Mu] - up*AP[Mu];
     X[Mu]  = V0[Mu] + u*A[Mu] + v*B[Mu];
     Z[Mu]  = V0P[Mu] + up*AP[Mu];
     Y[Mu]  = X[Mu] - Z[Mu];
     a2    += BP[Mu]*BP[Mu];
     BPdY  += BP[Mu]*Y[Mu];
     Y2    += Y[Mu]*Y[Mu];
   };
  a=sqrt(a2);
  a3=a*a2;
  vp0=-BPdY/a2;
  vp02=vp0*vp0;
  vp03=vp02*vp0;
  vp04=vp03*vp0;
  b2=Y2/a2 - vp02;

  double S1, S2, LogFac, Sum, Sum3, Sum4;
  S1=sqrt( b2 + vp02 );
  S2=sqrt( b2 + (vp0+up)*(vp0+up) );
  LogFac=log( (S2 + (up+vp0)) / (S1 + vp0) );
  Sum=vp0+up;
  Sum3=Sum*Sum*Sum;
  Sum4=Sum3*Sum;

  // the RHSs here are the integrals evaluated in Section 10 of the memo.
  // note we put the jacobian factor (u) into these quantities for convenience.
  double OneRM3Int = u*( (up+vp0)/S2 - vp0/S1 ) / (a3*b2);
  double  vpRM3Int = -vp0*OneRM3Int + u*( 1.0/S1 - 1.0/S2 ) / a3;
  double OneRM1Int = u*LogFac/a;
  double  vpRM1Int = -vp0*OneRM1Int + u*(S2-S1)/a;
  double OneR1Int  = u*a*0.5*(b2*LogFac + (up+vp0)*S2 - vp0*S1);
  double  vpR1Int  = -vp0*OneR1Int + u*a*(S2*S2*S2 - S1*S1*S1) / 3.0;
  double OneR2Int  = u*a2*( (Sum3-vp03)/3.0 + up*b2 );
  double  vpR2Int  = -vp0*OneR2Int + u*a2*( (Sum4-vp04)/4.0 + up*(vp0 + 0.5*up)*b2 );

  if (    !isfinite(OneRM3Int) 
       || !isfinite(OneRM1Int)
       || !isfinite(OneR1Int)
       || !isfinite(OneR2Int) )
   { printf("Bawonkatage! b2=%e, vp0=%e, up=%e, S1=%e, S2=%e\n",b2,vp0,up,S1,S2);
   };
  
  double CPCT[3], CPLT[3]; // 'cross product constant term' and 'cross product linear term'
  VecCross(X,Z,CPCT);
  VecCross(X,BP,CPLT);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int nf=0;

  fval[nf++] = Y[0]*OneRM3Int - BP[0]*vpRM3Int;
  fval[nf++] = Y[1]*OneRM3Int - BP[1]*vpRM3Int;
  fval[nf++] = Y[2]*OneRM3Int - BP[2]*vpRM3Int;
  fval[nf++] = CPCT[0]*OneRM3Int + CPLT[0]*vpRM3Int;
  fval[nf++] = CPCT[1]*OneRM3Int + CPLT[1]*vpRM3Int;
  fval[nf++] = CPCT[2]*OneRM3Int + CPLT[2]*vpRM3Int;

  fval[nf++] = OneRM1Int;
  fval[nf++] = up*OneRM1Int;
  fval[nf++] = vpRM1Int;
  fval[nf++] = u*OneRM1Int;
  fval[nf++] = u*up*OneRM1Int;
  fval[nf++] = u*vpRM1Int;
  fval[nf++] = v*OneRM1Int;
  fval[nf++] = v*up*OneRM1Int;
  fval[nf++] = v*vpRM1Int;

  fval[nf++] = OneR1Int;
  fval[nf++] = up*OneR1Int;
  fval[nf++] = vpR1Int;
  fval[nf++] = u*OneR1Int;
  fval[nf++] = u*up*OneR1Int;
  fval[nf++] = u*vpR1Int;
  fval[nf++] = v*OneR1Int;
  fval[nf++] = v*up*OneR1Int;
  fval[nf++] = v*vpR1Int;

  fval[nf++] = OneR2Int;
  fval[nf++] = up*OneR2Int;
  fval[nf++] = vpR2Int;
  fval[nf++] = u*OneR2Int;
  fval[nf++] = u*up*OneR2Int;
  fval[nf++] = u*vpR2Int;
  fval[nf++] = v*OneR2Int;
  fval[nf++] = v*up*OneR2Int;
  fval[nf++] = v*vpR2Int;

} 

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
void ComputeQIFIPPIData_Cubature(double **Va, double **Vb, QIFIPPIData *QIFD)
{
  /*--------------------------------------------------------------*/
  /* fill in the data structure used to pass parameters to the    */
  /* integrand routine                                            */
  /*--------------------------------------------------------------*/
  CFDData MyCFDData, *CFDD=&MyCFDData;

  CFDD->V0 = Va[0];
  VecSub(Va[1], Va[0], CFDD->A);
  VecSub(Va[2], Va[1], CFDD->B);

  CFDD->V0P = Vb[0];
  VecSub(Vb[1], Vb[0], CFDD->AP);
  VecSub(Vb[2], Vb[1], CFDD->BP);

  /*--------------------------------------------------------------*/
  /*- evaluate the adaptive cubature over the pair of triangles  -*/
  /*--------------------------------------------------------------*/
  int fdim = 33;
  double F[fdim], E[fdim];

  CFDD->nCalls=0;
#ifdef CUBATURE4D
  double Lower[4]={0.0, 0.0, 0.0, 0.0};
  double Upper[4]={1.0, 1.0, 1.0, 1.0};
  adapt_integrate(fdim, CFDIntegrand4D, CFDD, 4, Lower, Upper,
                  0, ABSTOL, RELTOL, F, E);
#else
  double Lower[3]={0.0, 0.0, 0.0};
  double Upper[3]={1.0, 1.0, 1.0};
  CFDD->nCalls=0;
  adapt_integrate(fdim, CFDIntegrand3D, CFDD, 3, Lower, Upper,
                  0, ABSTOL, RELTOL, F, E);

  printf("FIPPI cubature: %i calls\n",CFDD->nCalls);
#endif
  printf("FIPPI cubature: %i calls\n",CFDD->nCalls);


  /*--------------------------------------------------------------*/
  /*- unpack the results into the output data record -------------*/
  /*--------------------------------------------------------------*/
  int nf=0;

  QIFD->xMxpRM3[0]   = F[nf++];
  QIFD->xMxpRM3[1]   = F[nf++];
  QIFD->xMxpRM3[2]   = F[nf++];
  QIFD->xXxpRM3[0]   = F[nf++];
  QIFD->xXxpRM3[1]   = F[nf++];
  QIFD->xXxpRM3[2]   = F[nf++];

  QIFD->uvupvpRM1[0] = F[nf++];
  QIFD->uvupvpRM1[1] = F[nf++];
  QIFD->uvupvpRM1[2] = F[nf++];
  QIFD->uvupvpRM1[3] = F[nf++];
  QIFD->uvupvpRM1[4] = F[nf++];
  QIFD->uvupvpRM1[5] = F[nf++];
  QIFD->uvupvpRM1[6] = F[nf++];
  QIFD->uvupvpRM1[7] = F[nf++];
  QIFD->uvupvpRM1[8] = F[nf++];

  QIFD->uvupvpR1[0] = F[nf++];
  QIFD->uvupvpR1[1] = F[nf++];
  QIFD->uvupvpR1[2] = F[nf++];
  QIFD->uvupvpR1[3] = F[nf++];
  QIFD->uvupvpR1[4] = F[nf++];
  QIFD->uvupvpR1[5] = F[nf++];
  QIFD->uvupvpR1[6] = F[nf++];
  QIFD->uvupvpR1[7] = F[nf++];
  QIFD->uvupvpR1[8] = F[nf++];

  QIFD->uvupvpR2[0] = F[nf++];
  QIFD->uvupvpR2[1] = F[nf++];
  QIFD->uvupvpR2[2] = F[nf++];
  QIFD->uvupvpR2[3] = F[nf++];
  QIFD->uvupvpR2[4] = F[nf++];
  QIFD->uvupvpR2[5] = F[nf++];
  QIFD->uvupvpR2[6] = F[nf++];
  QIFD->uvupvpR2[7] = F[nf++];
  QIFD->uvupvpR2[8] = F[nf++];

}

/*--------------------------------------------------------------*/
/*- routine for computing Q-independent FIPPIs.                 */
/*-                                                             */
/*- inputs:                                                     */
/*-                                                             */
/*-  Va[i][j] = jth cartesian coord of ith vertex of panel A    */
/*-  Vb       = similarly for panel B                           */
/*-  FD       = must point to a preallocated QIFIPPIData        */
/*--------------------------------------------------------------*/
void ComputeQIFIPPIData(double **Va, double **Vb, QIFIPPIData *QIFD)
{ 
  
  int ncv=AssessPanelPair(Va, Vb);

  /*--------------------------------------------------------------*/
  /*- if there are no common vertices, then use 4-dimensional    -*/
  /*- adaptive cubature over both triangles to compute the FIPPIs-*/
  /*--------------------------------------------------------------*/
  if (ncv==0)
   ComputeQIFIPPIData_Cubature(Va, Vb, QIFD);
  else
   ComputeQIFIPPIData_TaylorDuffy(Va[0], Va[1], Va[2], Vb[1], Vb[2], QIFD);
}

/*--------------------------------------------------------------*/
/*- routine for computing Q-dependent FIPPIs.                   */
/*--------------------------------------------------------------*/
void GetQDFIPPIData(double **Va, double *Qa, double **Vb, double *Qb, 
                    void *opFDT, QDFIPPIData *QDFD)
{
  int Flipped;
  double *OVa[3], *OQa, *OVb[3], *OQb;  // 'ordered vertices'
  QIFIPPIData MyQIFD, *QIFD;

  Flipped=CanonicallyOrderVertices(Va, Qa, Vb, Qb, OVa, &OQa, OVb, &OQb);

  /*--------------------------------------------------------------*/
  /*- get the Q-independent FIPPIs by looking them up in a table  */
  /*- if we have one or by computing them if we don't             */
  /*--------------------------------------------------------------*/
  if (opFDT)
   { 
     QIFD=((FIPPIDataTable *)opFDT)->GetQIFIPPIData(OVa, OVb);
   }
  else
   { 
     ComputeQIFIPPIData(OVa, OVb, &MyQIFD);
     QIFD=&MyQIFD;
   };
  
  /*--------------------------------------------------------------*/
  /*- now assemble the Q-dependent FIPPIs.                        */
  /*--------------------------------------------------------------*/
  double A[3],  B[3],  V0mQ[3];
  double AP[3], BP[3], V0PmQP[3]; 
  double Delta;

  VecSub(OVa[1], OVa[0], A);
  VecSub(OVa[2], OVa[1], B);
  VecSub(OVa[0], OQa, V0mQ);

  VecSub(OVb[1], OVb[0], AP);
  VecSub(OVb[2], OVb[1], BP);
  VecSub(OVb[0], OQb, V0PmQP);

  double VdVP = VecDot(V0mQ, V0PmQP);
  double VdAP = VecDot(V0mQ, AP);
  double VdBP = VecDot(V0mQ, BP);
  double AdVP = VecDot(A, V0PmQP);
  double AdAP = VecDot(A, AP);
  double AdBP = VecDot(A, BP);
  double BdVP = VecDot(B, V0PmQP);
  double BdAP = VecDot(B, AP);
  double BdBP = VecDot(B, BP);

  double V0QaxQb[3], QamQb[3], QaxQb[3], V0mV0P[3], Scratch[3];
  double *V0=OVa[0], *V0P=OVb[0];
  VecSub(V0, V0P, V0mV0P);
  VecSub(Qa, Qb, QamQb);
  VecCross(Qa, Qb, QaxQb);
  double hTimes0000 = VecDot(QaxQb, V0mV0P)  + VecDot(QamQb,  VecCross(V0,V0P, Scratch));
  double hTimes1000 = VecDot(QaxQb, A     )  + VecDot(QamQb,  VecCross(A,V0P,Scratch));
  double hTimes0100 = VecDot(QaxQb, B     )  + VecDot(QamQb,  VecCross(B,V0P,Scratch));
  double hTimes0010 = -VecDot(QaxQb, AP    ) + VecDot(QamQb,  VecCross(V0,AP,Scratch));
  double hTimes0001 = -VecDot(QaxQb, BP    ) + VecDot(QamQb,  VecCross(V0,BP,Scratch));
  double hTimes1010 = VecDot(QamQb,  VecCross(A,AP,Scratch));
  double hTimes1001 = VecDot(QamQb,  VecCross(A,BP,Scratch));
  double hTimes0110 = VecDot(QamQb,  VecCross(B,AP,Scratch));
  double hTimes0101 = VecDot(QamQb,  VecCross(B,BP,Scratch));

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  QDFD->hTimesRM3 = VecDot(QaxQb, QIFD->xMxpRM3) + VecDot(QamQb, QIFD->xXxpRM3);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  QDFD->hDotRM1 =   QIFD->uvupvpRM1[_ONE] * VdVP
                  + QIFD->uvupvpRM1[_UP ] * VdAP
                  + QIFD->uvupvpRM1[_VP ] * VdBP
                  + QIFD->uvupvpRM1[_U  ] * AdVP
                  + QIFD->uvupvpRM1[_UUP] * AdAP
                  + QIFD->uvupvpRM1[_UVP] * AdBP
                  + QIFD->uvupvpRM1[_V  ] * BdVP
                  + QIFD->uvupvpRM1[_VUP] * BdAP
                  + QIFD->uvupvpRM1[_VVP] * BdBP;

  QDFD->hNablaRM1 =  4.0*QIFD->uvupvpRM1[_ONE];

  QDFD->hTimesRM1 =   QIFD->uvupvpRM1[_ONE] * hTimes0000
                    + QIFD->uvupvpRM1[_U  ] * hTimes1000
                    + QIFD->uvupvpRM1[_V  ] * hTimes0100
                    + QIFD->uvupvpRM1[_UP ] * hTimes0010
                    + QIFD->uvupvpRM1[_VP ] * hTimes0001
                    + QIFD->uvupvpRM1[_UUP] * hTimes1010
                    + QIFD->uvupvpRM1[_UVP] * hTimes1001
                    + QIFD->uvupvpRM1[_VUP] * hTimes0110
                    + QIFD->uvupvpRM1[_VVP] * hTimes0101;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  QDFD->hDotR0  =   uvupvpR0_ONE * VdVP
                  + uvupvpR0_UP  * VdAP
                  + uvupvpR0_VP  * VdBP
                  + uvupvpR0_U   * AdVP
                  + uvupvpR0_UUP * AdAP
                  + uvupvpR0_UVP * AdBP
                  + uvupvpR0_V   * BdVP
                  + uvupvpR0_VUP * BdAP
                  + uvupvpR0_VVP * BdBP;

  QDFD->hNablaR0 = 1.0;

  QDFD->hTimesR0  =   uvupvpR0_ONE * hTimes0000
                    + uvupvpR0_U   * hTimes1000
                    + uvupvpR0_V   * hTimes0100
                    + uvupvpR0_UP  * hTimes0010
                    + uvupvpR0_VP  * hTimes0001
                    + uvupvpR0_UUP * hTimes1010
                    + uvupvpR0_UVP * hTimes1001
                    + uvupvpR0_VUP * hTimes0110
                    + uvupvpR0_VVP * hTimes0101;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  QDFD->hDotR1 =   QIFD->uvupvpR1[_ONE] * VdVP
                 + QIFD->uvupvpR1[_UP ] * VdAP
                 + QIFD->uvupvpR1[_VP ] * VdBP
                 + QIFD->uvupvpR1[_U  ] * AdVP
                 + QIFD->uvupvpR1[_UUP] * AdAP
                 + QIFD->uvupvpR1[_UVP] * AdBP
                 + QIFD->uvupvpR1[_V  ] * BdVP
                 + QIFD->uvupvpR1[_VUP] * BdAP
                 + QIFD->uvupvpR1[_VVP] * BdBP;

  QDFD->hNablaR1 =  4.0*QIFD->uvupvpR1[_ONE];

  QDFD->hTimesR1 =   QIFD->uvupvpR1[_ONE] * hTimes0000
                   + QIFD->uvupvpR1[_U  ] * hTimes1000
                   + QIFD->uvupvpR1[_V  ] * hTimes0100
                   + QIFD->uvupvpR1[_UP ] * hTimes0010
                   + QIFD->uvupvpR1[_VP ] * hTimes0001
                   + QIFD->uvupvpR1[_UUP] * hTimes1010
                   + QIFD->uvupvpR1[_UVP] * hTimes1001
                   + QIFD->uvupvpR1[_VUP] * hTimes0110
                   + QIFD->uvupvpR1[_VVP] * hTimes0101;

  /*--------------------------------------------------------------*/
  /*-- NOTE: the R2 integrals may actually be done analytically,  */
  /*--       but the calculation is so cumbersome that i do them  */
  /*--       numerically for now. in the future, maybe explore    */
  /*--       the costs and benefits of doing them analytically.   */
  /*--------------------------------------------------------------*/
  QDFD->hDotR2 =   QIFD->uvupvpR2[_ONE] * VdVP
                 + QIFD->uvupvpR2[_UP ] * VdAP
                 + QIFD->uvupvpR2[_VP ] * VdBP
                 + QIFD->uvupvpR2[_U  ] * AdVP
                 + QIFD->uvupvpR2[_UUP] * AdAP
                 + QIFD->uvupvpR2[_UVP] * AdBP
                 + QIFD->uvupvpR2[_V  ] * BdVP
                 + QIFD->uvupvpR2[_VUP] * BdAP
                 + QIFD->uvupvpR2[_VVP] * BdBP;

  QDFD->hNablaR2 =  4.0*QIFD->uvupvpR2[_ONE];

  if (Flipped)
   { QDFD->hTimesRM3 *= -1.0; 
     QDFD->hTimesRM1 *= -1.0; 
     QDFD->hTimesR1  *= -1.0; 
   };

}
