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
#include "TaylorMaster.h"

#define ABSTOL 1.0e-12
#define RELTOL 1.0e-8

#define _ONE 0
#define _U   1
#define _V   2
#define _UP  3
#define _UUP 4
#define _VUP 5
#define _VP  6
#define _UVP 7
#define _VVP 8


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
/*- 'vertex less than.' returns 1 if V1<V2, 0 otherwise.        */
/*- vertices are sorted using a fairly obvious sorting scheme.  */
/*--------------------------------------------------------------*/
int FIPPIDataTable::VLT(double *V1, double *V2)
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
/* create a unique search key for the given panel pair.        */
/* (note: search keys are unique up to rigid translation of    */
/*  both panels. if Pa,Pb and PaP, PbP are panel pairs that    */
/*  are equivalent under a rigid translation then we want them */
/*  to map to the same key since the FIPPIs are equal.         */
/***************************************************************/
void FIPPIDataTable::ComputeSearchKey(double **Va, double **Vb, double *Key)
{ 
  /***************************************************************/
  /* sort the vertices in ascending order ************************/
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
  /* the search key is kinda stupid, just a string of 15 doubles */
  /* as follows:                                                 */
  /* 0--2    VMed  - VMin [0..2]                                 */
  /* 3--5    VMax  - VMin [0..2]                                 */
  /* 6--8    VMinP - VMin [0..2]                                 */
  /* 9--11   VMedP - VMin [0..2]                                 */
  /* 12--14  VMaxP - VMin [0..2]                                 */
  /***************************************************************/
  VecSub(VMed,  VMin, Key+0 );
  VecSub(VMax,  VMin, Key+3 );
  VecSub(VMinP, VMin, Key+6 );
  VecSub(VMedP, VMin, Key+9 );
  VecSub(VMaxP, VMin, Key+12);
  
}

/*--------------------------------------------------------------*/
/*- routine for fetching a FIPPI data record from a FIPPIDT:    */
/*- we look through our table to see if we have a record for    */
/*- this panel pair, we return it if we do, and otherwise we    */
/*- compute a new FIPPI data record for this panel pair and     */
/*- add it to the table                                         */
/*--------------------------------------------------------------*/
QIFIPPIData *FIPPIDataTable::GetQIFIPPIData(double **Va, double **Vb)
{
  double Key[15];

  ComputeSearchKey(Va, Vb, Key);
}

void GetQIFIPPIData(double **Va, double **Vb)

/*--------------------------------------------------------------*/
/*- integrand routine used for evaluating FIPPIs by cubature.  -*/
/*- note: CFD = 'compute FIPPI data.'                          -*/
/*--------------------------------------------------------------*/
typedef struct CFDData
 {
   double *V0[3], A[3], B[3];
   double *V0P[3], AP[3], BP[3];
   int nCalls;
 } CFDData;

void CFDIntegrand(unsigned ndim, const double *x, void *params,
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
  fval[nf++] = u*oor;
  fval[nf++] = v*oor;
  fval[nf++] = up*oor;
  fval[nf++] = u*up*oor;
  fval[nf++] = v*up*oor;
  fval[nf++] = vp*oor;
  fval[nf++] = u*vp*oor;
  fval[nf++] = v*vp*oor;

  fval[nf++] = r;
  fval[nf++] = u*r;
  fval[nf++] = v*r;
  fval[nf++] = up*r;
  fval[nf++] = u*up*r;
  fval[nf++] = v*up*r;
  fval[nf++] = vp*r;
  fval[nf++] = u*vp*r;
  fval[nf++] = v*vp*r;

  fval[nf++] = r2;
  fval[nf++] = u*r2;
  fval[nf++] = v*r2;
  fval[nf++] = up*r2;
  fval[nf++] = u*up*r2;
  fval[nf++] = v*up*r2;
  fval[nf++] = vp*r2;
  fval[nf++] = u*vp*r2;
  fval[nf++] = v*vp*r2;

} 

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
void ComputeQIFIPPIData_Cubature(double **Va, double **Vb, QIFIPPIData *FD)
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
  double Lower[4]={0.0, 0.0, 0.0, 0.0};
  double Upper[4]={1.0, 1.0, 1.0, 1.0};
  int fdim = 33;
  double F[fdim], E[fdim];
CFDD->nCalls=0;
printf("%i \n",fdim);
  adapt_integrate(fdim, CFDIntegrand, CFDD, 4, Lower, Upper,
                  0, ABSTOL, RELTOL, F, E);
printf("FIPPI cubature: %i calls\n",CFDD->nCalls);

  /*--------------------------------------------------------------*/
  /*- unpack the results into the output data record -------------*/
  /*--------------------------------------------------------------*/
  int nf=0;

  FDR->xMxpRM3[0]   = F[nf++];
  FDR->xMxpRM3[1]   = F[nf++];
  FDR->xMxpRM3[2]   = F[nf++];
  FDR->xXxpRM3[3]   = F[nf++];
  FDR->xXxpRM3[4]   = F[nf++];
  FDR->xXxpRM3[5]   = F[nf++];

  FDR->uvupvpRM1[0] = F[nf++];
  FDR->uvupvpRM1[1] = F[nf++];
  FDR->uvupvpRM1[2] = F[nf++];
  FDR->uvupvpRM1[3] = F[nf++];
  FDR->uvupvpRM1[4] = F[nf++];
  FDR->uvupvpRM1[5] = F[nf++];
  FDR->uvupvpRM1[6] = F[nf++];
  FDR->uvupvpRM1[7] = F[nf++];
  FDR->uvupvpRM1[8] = F[nf++];

  FDR->uvupvpR1[0] = F[nf++];
  FDR->uvupvpR1[1] = F[nf++];
  FDR->uvupvpR1[2] = F[nf++];
  FDR->uvupvpR1[3] = F[nf++];
  FDR->uvupvpR1[4] = F[nf++];
  FDR->uvupvpR1[5] = F[nf++];
  FDR->uvupvpR1[6] = F[nf++];
  FDR->uvupvpR1[7] = F[nf++];
  FDR->uvupvpR1[8] = F[nf++];

  FDR->uvupvpR2[0] = F[nf++];
  FDR->uvupvpR2[1] = F[nf++];
  FDR->uvupvpR2[2] = F[nf++];
  FDR->uvupvpR2[3] = F[nf++];
  FDR->uvupvpR2[4] = F[nf++];
  FDR->uvupvpR2[5] = F[nf++];
  FDR->uvupvpR2[6] = F[nf++];
  FDR->uvupvpR2[7] = F[nf++];
  FDR->uvupvpR2[8] = F[nf++];

}

/*--------------------------------------------------------------*/
/*- note: on entry it is assumed that a call to AssessPanelPair */
/*- has already been made and that Va, Vb are in the order      */
/*- they were put in by that routine                            */
/*--------------------------------------------------------------*/
void ComputeQIFIPPIData_TaylorDuffy(int ncv, double **Va, double **Vb, 
                                    QIFIPPIData *FD)
{ 

  /*--------------------------------------------------------------*/
  /*- otherwise (there are common vertices) compute the FIPPIs   -*/
  /*- using the taylor-duffy method                              -*/
  /*--------------------------------------------------------------*/
  int WhichCase;
  switch (ncv)
   { case 1: WhichCase=TM_COMMONVERTEX;   break;
     case 2: WhichCase=TM_COMMONEDGE;     break;
     case 3: WhichCase=TM_COMMONTRIANGLE; break; // should never happen
   };

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
void ComputeQIFIPPIData(double **Va, double **Vb, FIPPIData *FD)
{ 
  int ncv=AssessPanelPair(Va, Vb);

  /*--------------------------------------------------------------*/
  /*- if there are no common vertices, then use 4-dimensional    -*/
  /*- adaptive cubature over both triangles to compute the FIPPIs-*/
  /*--------------------------------------------------------------*/
  if (ncv==0)
   ComputeQIFIPPIData_Cubature(Va, Vb, FD);
  else
   ComputeQIFIPPIData_TaylorDuffy(ncv, Va, Vb, FD);

}

/*--------------------------------------------------------------*/
/*- routine for computing Q-dependent FIPPIs.                   */
/*--------------------------------------------------------------*/
void GetQDFIPPIData(double **Va, double *Qa, double **Vb, double *Qb, 
                    void *opFDT, QDFIPPIData *QDFD)
{
  double *OVa[3], *OQa, *OVb[3], *OQb;  // 'ordered vertices'
  QIFIPPIData MyQIFD, *QIFD;

  /*--------------------------------------------------------------*/
  /*- get the Q-independent FIPPIs by looking them up in a table  */
  /*- if we have one or by computing them if we don't             */
  /*--------------------------------------------------------------*/
  if (opFDT)
   { 
     Flipped=CanonicallyOrderVertices(Va, Qa, Vb, Qb, OVa, &Qa, OVb, &Qb);
     QIFD=((FIPPIDataTable *)opFDT)->GetQIFIPPIData(OVa, OVb);
   }
  else
   { 
     ComputeQIFIPPIData(Va, Vb, &MyQIFD);
     QIFD=&MyQIFD;

     OVa[0]=Va[0];
     OVa[1]=Va[1];
     OVa[2]=Va[2];
     OQa=Qa;

     OVb[0]=Vb[0];
     OVb[1]=Vb[1];
     OVb[2]=Vb[2];
     OQb=Qb;

     Flipped=0;
   };
  
  double Sign = ( Flipped ? -1.0 : 1.0 );

  /*--------------------------------------------------------------*/
  /*- now assemble the Q-dependent FIPPIs.                        */
  /*--------------------------------------------------------------*/
  double A[3],  B[3],  V0mQ[3];
  double AP[3], BP[3], V0PmQP[3], 
  double Delta;

  VecSub(OVa[1], OVa[0], A);
  VecSub(OVa[2], OVa[1], B);
  VecSub(OVa[0], OQa, V0mQ[3]);

  VecSub(OVb[1], OVb[0], AP);
  VecSub(OVb[2], OVb[1], BP);
  VecSub(OVb[0], OQb, V0PmQP[3]);

  double VdVP, AdVP, BdVP, VdAP, AdAP, BdAP, VdBP, AdBP, BdBP;

  VecDot(V0mQ, V0PmQP, VdVP);
  VecDot(A,    V0PmQP, AdVP);
  VecDot(B,    V0PmQP, BdVP);
  VecDot(V0mQ, AP,     VdAP);
  VecDot(A,    AP,     AdAP);
  VecDot(B,    AP,     BdAP);
  VecDot(V0mQ, BP,     VdBP);
  VecDot(A,    BP,     AdBP);
  VecDot(B,    BP,     BdBP);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  QDFD->hDotRM1 =   QIFD->uvupvpRM1[_ONE] * VdVP
                  + QIFD->uvupvpRM1[_U  ] * AdVP
                  + QIFD->uvupvpRM1[_V  ] * BdVP
                  + QIFD->uvupvpRM1[_UP ] * VdAP
                  + QIFD->uvupvpRM1[_UUP] * AdAP
                  + QIFD->uvupvpRM1[_VUP] * BdAP
                  + QIFD->uvupvpRM1[_VP ] * VdBP
                  + QIFD->uvupvpRM1[_UVP] * AdBP
                  + QIFD->uvupvpRM1[_VVP] * BdBP;

  QDFD->hNablaRM1 =  4.0*QIFD->uvupvpRM1[_ONE];

  // QDFD->hTimesRM1 =; 

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  QDFD->hDotR0  =   VdVP / 4.0
                  + AdVP / 6.0
                  + BdVP / 12.0
                  + VdAP / 6.0
                  + AdAP / 9.0
                  + BdAP / 18.0
                  + VdBP / 12.0
                  + AdBP / 18.0
                  + BdBP / 36.0;

  QDFD->hNablaR0 = 1.0;

  // QDFD->hTimesR0=; 

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  QDFD->hDotR1 =   QIFD->uvupvpR1[_ONE] * VdVP
                 + QIFD->uvupvpR1[_U  ] * AdVP
                 + QIFD->uvupvpR1[_V  ] * BdVP
                 + QIFD->uvupvpR1[_UP ] * VdAP
                 + QIFD->uvupvpR1[_UUP] * AdAP
                 + QIFD->uvupvpR1[_VUP] * BdAP
                 + QIFD->uvupvpR1[_VP ] * VdBP
                 + QIFD->uvupvpR1[_UVP] * AdBP
                 + QIFD->uvupvpR1[_VVP] * BdBP;

  QDFD->hNablaRM1 =  4.0*QIFD->uvupvpR1[_ONE];

  // QDFD->hTimesRM1 =; 

  /*--------------------------------------------------------------*/
  /*-- NOTE: the R2 integrals may actually be done analytically,  */
  /*--       but the calculation is so cumbersome that i do them  */
  /*--       numerically for now. in the future, maybe explore    */
  /*--       the costs and benefits of doing them analytically.   */
  /*--------------------------------------------------------------*/
  QDFD->hDotR2 =   QIFD->uvupvpR2[_ONE] * VdVP
                 + QIFD->uvupvpR2[_U  ] * AdVP
                 + QIFD->uvupvpR2[_V  ] * BdVP
                 + QIFD->uvupvpR2[_UP ] * VdAP
                 + QIFD->uvupvpR2[_UUP] * AdAP
                 + QIFD->uvupvpR2[_VUP] * BdAP
                 + QIFD->uvupvpR2[_VP ] * VdBP
                 + QIFD->uvupvpR2[_UVP] * AdBP
                 + QIFD->uvupvpR2[_VVP] * BdBP;

  QDFD->hNablaRM2 =  4.0*QIFD->uvupvpR2[_ONE];

  // QDFD->hTimesRM1 =; 

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/

}
