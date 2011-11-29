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
int FIPPIDataRecord::VLT(double *V1, double *V2)
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
int FIPPIDataRecord::ComputeSearchKey(double **Va, double **Vb, double *Key)
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
FIPPIDataRecord *FIPPIDataTable::GetFIPPIDataRecord(double **Va, double **Vb,
                                                    int NeedDerivatives)
{
  double Key[15];

  ComputeSearchKey(Va, Vb, Key);
}

/*--------------------------------------------------------------*/
/*- integrand routine used for evaluating FIPPIs by cubature   -*/
/*- note: CFDR = 'compute FIPPI data record'                   -*/
/*--------------------------------------------------------------*/
typedef struct CFDRData
 {
   double *V0, A[3], B[3], *Q;
   double *V0P, AP[3], BP[3], *QP;
   int NeedDerivatives;
   int nCalls;
 } CFDRData;

void CFDRIntegrand(unsigned ndim, const double *x, void *params,
                   unsigned fdim, double *fval)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  CFDRData *CFDRD=(CFDRData *)params;
  CFDRD->nCalls++;

  double *V0=CFDRD->V0;
  double *A=CFDRD->A;
  double *B=CFDRD->B;
  double *Q=CFDRD->Q;

  double *V0P=CFDRD->V0P;
  double *AP=CFDRD->AP;
  double *BP=CFDRD->BP;
  double *QP=CFDRD->QP;

  int NeedDerivatives=CFDRD->NeedDerivatives;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double u=x[0];
  double v=u*x[1];
  double up=x[2];
  double vp=up*x[3];

  double F[3], X[3], FP[3], XP[3], FxFP[3], R[3];
  double r, r2=0.0, r3;
  int Mu;
  for (Mu=0; Mu<3; Mu++)
   { X[Mu]  = V0[Mu] + u*A[Mu] + v*B[Mu];
     F[Mu]  = X[Mu] - Q[Mu];
     XP[Mu] = V0P[Mu] + up*AP[Mu] + vp*BP[Mu];
     FP[Mu] = XP[Mu] - QP[Mu];
     R[Mu]  = X[Mu]-XP[Mu];
     r2    += R[Mu]*R[Mu];
   };
  r=sqrt(r2);
  r3=r*r2;

  /*--------------------------------------------------------------*/
  /*- assemble output vector -------------------------------------*/
  /*--------------------------------------------------------------*/
  double hDot=u*up*VecDot(F, FP);
  fval[ 0] = hDot / r;
  fval[ 1] = hDot;
  fval[ 2] = hDot*r;
  fval[ 3] = hDot*r2;

  double hNabla=u*up*4.0;
  fval[ 4] = hNabla / r;
  fval[ 5] = hNabla;
  fval[ 6] = hNabla*r;
  fval[ 7] = hNabla*r2;

  double hTimes=u*up*VecDot( VecCross(F, FP, FxFP), R );
  fval[ 8] = hTimes / r3;
  fval[ 9] = hTimes / r;
  fval[10] = hTimes;  
  fval[11] = hTimes * r;

  /*--------------------------------------------------------------*/
  /*- additional entries in the output vector that are only needed*/
  /*- for derivative calculations                                 */
  /*--------------------------------------------------------------*/
  if (NeedDerivatives)
   { 
     fval[12] = hDot   / r3;
     fval[13] = hNabla / r3;
     fval[14] = hTimes / (r2*r3);

     FxFP[0]*=u*up;
     FxFP[1]*=u*up;
     FxFP[2]*=u*up;

     fval[15] = FxFP[0] / r3;
     fval[16] = FxFP[1] / r3;
     fval[17] = FxFP[2] / r3;

     fval[18] = FxFP[0] / r;
     fval[19] = FxFP[1] / r;
     fval[20] = FxFP[2] / r;

     fval[21] = FxFP[0];
     fval[22] = FxFP[1];
     fval[23] = FxFP[2];

     fval[24] = FxFP[0] * r;
     fval[25] = FxFP[1] * r;
     fval[26] = FxFP[2] * r;
   };

} 

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
void ComputeFIPPIDataRecord_Cubature(double **Va, double *Qa,
                                     double **Vb, double *Qb,
                                     int NeedDerivatives,
                                     FIPPIDataRecord *FDR)
{
  FDR->HaveDerivatives=NeedDerivatives;
  
  /*--------------------------------------------------------------*/
  /* fill in the data structure used to pass parameters to the    */
  /* integrand routine                                            */
  /*--------------------------------------------------------------*/
  CFDRData MyCFDRData, *CFDRD=&MyCFDRData;

  CFDRD->V0 = Va[0];
  VecSub(Va[1], Va[0], CFDRD->A);
  VecSub(Va[2], Va[1], CFDRD->B);
  CFDRD->Q  = Qa;

  CFDRD->V0P = Vb[0];
  VecSub(Vb[1], Vb[0], CFDRD->AP);
  VecSub(Vb[2], Vb[1], CFDRD->BP);
  CFDRD->QP = Qb;

  CFDRD->NeedDerivatives = NeedDerivatives;

  /*--------------------------------------------------------------*/
  /*- evaluate the adaptive cubature                             -*/
  /*--------------------------------------------------------------*/
  double Lower[4]={0.0, 0.0, 0.0, 0.0};
  double Upper[4]={1.0, 1.0, 1.0, 1.0};
  int fdim = NeedDerivatives ? 27 : 12;
  double F[fdim], E[fdim];
CFDRD->nCalls=0;
  adapt_integrate(fdim, CFDRIntegrand, CFDRD, 4, Lower, Upper,
                  0, ABSTOL, RELTOL, F, E);
printf("FIPPI cubature: %i calls\n",CFDRD->nCalls);

  /*--------------------------------------------------------------*/
  /*- unpack the results into the output data record -------------*/
  /*--------------------------------------------------------------*/
  FDR->hDotRM1   = F[ 0];
  FDR->hDotR0    = F[ 1];
  FDR->hDotR1    = F[ 2];
  FDR->hDotR2    = F[ 3];

  FDR->hNablaRM1 = F[ 4];
  FDR->hNablaR0  = F[ 5];
  FDR->hNablaR1  = F[ 6];
  FDR->hNablaR2  = F[ 7];

  FDR->hTimesRM3 = F[ 8];
  FDR->hTimesRM1 = F[ 9];
  FDR->hTimesR0  = F[10];
  FDR->hTimesR1  = F[11];

  if (NeedDerivatives)
   { 
     FDR->hDotRM3    = F[12];
     FDR->hNablaRM3  = F[13];
     FDR->hTimesRM5  = F[14];
     
     FDR->dhTimesdRMuRM3[0] = F[15];
     FDR->dhTimesdRMuRM3[1] = F[16];
     FDR->dhTimesdRMuRM3[2] = F[17];

     FDR->dhTimesdRMuRM1[0] = F[18];
     FDR->dhTimesdRMuRM1[1] = F[19];
     FDR->dhTimesdRMuRM1[2] = F[20];

     FDR->dhTimesdRMuR0[0] = F[21];
     FDR->dhTimesdRMuR0[1] = F[22];
     FDR->dhTimesdRMuR0[2] = F[23];

     FDR->dhTimesdRMuR1[0] = F[24];
     FDR->dhTimesdRMuR1[1] = F[25];
     FDR->dhTimesdRMuR1[2] = F[26];

   };

}

/*--------------------------------------------------------------*/
/*- routine for computing frequency-independent panel-panel     */
/*- integrals                                                   */
/*-                                                             */
/*- inputs:                                                     */
/*-                                                             */
/*-  Va[i][j] = jth cartesian coord of ith vertex of panel A    */
/*-  Qa[i]    = ith cartesian coordinate of current source/sink */
/*-             vertex of panel A                               */
/*-  Vb, Qb   = similarly for panel B                           */
/*-  NeedDerivatives = set to 1 or 0                            */
/*-  FDR      = must point to a preallocated FIPPIDataRecord    */
/*-             structure                                       */
/*-                                                             */
/*- the return value is FDR.                                    */
/*--------------------------------------------------------------*/
FIPPIDataRecord *ComputeFIPPIDataRecord(double **Va, double *Qa,
                                        double **Vb, double *Qb,
                                        int NeedDerivatives,
                                        FIPPIDataRecord *FDR)
{ 
  int ncv=AssessPanelPair(Va, Vb);

  /*--------------------------------------------------------------*/
  /*- if there are no common vertices, then use 4-dimensional    -*/
  /*- adaptive cubature over both triangles to compute the FIPPIs-*/
  /*--------------------------------------------------------------*/
  if (ncv==0)
   { ComputeFIPPIDataRecord_Cubature(Va, Qa, Vb, Qb,
                                     NeedDerivatives,
                                     FDR);
     return FDR;
   };

  FDR->HaveDerivatives=0; // no derivatives for common-vertex cases 

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
  
  /*--------------------------------------------------------------*/
  /*- hDot integrals ---------------------------------------------*/
  /*--------------------------------------------------------------*/
  FDR->hDotRM1=real( TaylorMaster(WhichCase, TM_RP, TM_DOT, -1.0,
                                  Va[0], Va[1], Va[2],
                                  Vb[1], Vb[2], Qa, Qb)
                   );

  /* FIXME this one can be computed analytically */ 
  FDR->hDotR0 =real( TaylorMaster(WhichCase, TM_RP, TM_DOT, 0.0, 
                                  Va[0], Va[1], Va[2], 
                                  Vb[1], Vb[2], Qa, Qb)
                   );

  FDR->hDotR1 =real( TaylorMaster(WhichCase, TM_RP, TM_DOT, 1.0, 
                                  Va[0], Va[1], Va[2], 
                                  Vb[1], Vb[2], Qa, Qb)
                   );

  /* FIXME this one can be computed analytically */ 
  FDR->hDotR2 =real( TaylorMaster(WhichCase, TM_RP, TM_DOT, 2.0, 
                                  Va[0], Va[1], Va[2], 
                                  Vb[1], Vb[2], Qa, Qb)
                   );
  
  /*--------------------------------------------------------------*/
  /*- hNabla integrals -------------------------------------------*/
  /*--------------------------------------------------------------*/
  FDR->hNablaRM1 = 4.0*real( TaylorMaster(WhichCase, TM_RP, TM_ONE, -1.0,
                                          Va[0], Va[1], Va[2], 
                                          Vb[1], Vb[2], Qa, Qb)
                           );

  /* FIXME this one can be computed analytically */ 
  FDR->hNablaR0  = 4.0*real( TaylorMaster(WhichCase, TM_RP, TM_ONE, 0.0, 
                                          Va[0], Va[1], Va[2], 
                                          Vb[1], Vb[2], Qa, Qb)
                           );

  FDR->hNablaR1  = 4.0*real( TaylorMaster(WhichCase, TM_RP, TM_ONE, 1.0, 
                                          Va[0], Va[1], Va[2], 
                                          Vb[1], Vb[2], Qa, Qb)
                           );

  /* FIXME this one can be computed analytically */ 
  FDR->hNablaR2  = 4.0*real( TaylorMaster(WhichCase, TM_RP, TM_ONE, 2.0, 
                                          Va[0], Va[1], Va[2], 
                                          Vb[1], Vb[2], Qa, Qb)
                           );
  
  /*--------------------------------------------------------------*/
  /*- hTimes integrals -------------------------------------------*/
  /*--------------------------------------------------------------*/
  FDR->hTimesRM3=real( TaylorMaster(WhichCase, TM_RP, TM_CROSS, -3.0,
                                    Va[0], Va[1], Va[2],
                                    Vb[1], Vb[2], Qa, Qb)
                     );

  FDR->hTimesRM1=real( TaylorMaster(WhichCase, TM_RP, TM_CROSS, -1.0,
                                    Va[0], Va[1], Va[2],
                                    Vb[1], Vb[2], Qa, Qb)
                     );

  /* FIXME this one can be computed analytically */ 
  FDR->hTimesR0=real( TaylorMaster(WhichCase, TM_RP, TM_CROSS, 0.0,
                                   Va[0], Va[1], Va[2],
                                   Vb[1], Vb[2], Qa, Qb)
                    );

  FDR->hTimesR1=real( TaylorMaster(WhichCase, TM_RP, TM_CROSS, 1.0,
                                   Va[0], Va[1], Va[2],
                                   Vb[1], Vb[2], Qa, Qb)
                    );

}
