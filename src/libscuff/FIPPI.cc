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

#define ABSTOL 1.0e-12
#define RELTOL 1.0e-8

/***************************************************************/
/* implementation of the FIPPIDT class *************************/
/***************************************************************/

/*--------------------------------------------------------------*/
/*- class constructor ------------------------------------------*/
/*--------------------------------------------------------------*/
FIPPIDataTable()
{
}

/*--------------------------------------------------------------*/
/*- class destructor  ------------------------------------------*/
/*--------------------------------------------------------------*/
~FIPPIDataTable()
{
} 

/*--------------------------------------------------------------*/
/*- routine for fetching a FIPPI data record from a FIPPIDT:    */
/*- we look through our table to see if we have a record for    */
/*- this panel pair, we return it if we do, and otherwise we    */ 
/*- compute a new FIPPI data record for this panel pair and     */
/*- add it to the table                                         */
/*--------------------------------------------------------------*/
FIPPIDataRecord *FIPPIDataTable::GetFIPPIDataRecord(double **Va, double *Qa,
                                                    double **Vb, double *Qb,
                                                    int NeedDerivatives);

/*--------------------------------------------------------------*/
/*- integrand routine used for evaluating FIPPIs by cubature   -*/
/*- note: CFDR = 'compute FIPPI data record'                   -*/
/*--------------------------------------------------------------*/
typedef struct CFDRData
 {
   double *V0, A[3], B[3], *Q;
   double *V0P, AP[3], BP[3], *QP;
   int NeedDerivatives;
 } CFDRData;

void CFDRIntegrand(unsigned ndim, const double *x, void *params,
                   unsigned fdim, double *fval)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  CFDRData *CFDRD=(CFDRData *)params;

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
  fval[ 9] = hTimes / r1;
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
  int Lower[4]={0.0, 0.0, 0.0, 0.0};
  int Upper[4]={1.0, 1.0, 1.0, 1.0};
  int fdim = NeedDerivatives ? 12 : 27;
  double F[fdim], E[fdim];
  adapt_integrate(fdim, CFDRIntegrand, CFDRD, 4, Lower, Upper,
                  0, ABSTOL, RELTOL, F, E);

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
     
     FDR->dhTimesdRMuRm3[0] = F[15];
     FDR->dhTimesdRMuRm3[1] = F[16];
     FDR->dhTimesdRMuRm3[2] = F[17];

     FDR->dhTimesdRMuRm1[0] = F[18];
     FDR->dhTimesdRMuRm1[1] = F[19];
     FDR->dhTimesdRMuRm1[2] = F[20];

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
  ncv=AssessPanelPair(Va, Vb);

  /*--------------------------------------------------------------*/
  /*- if there are no common vertices, then use 4-dimensional    -*/
  /*- adaptive cubature over both triangles to compute the FIPPIs-*/
  /*--------------------------------------------------------------*/
  if (ncv==0)
   { ComputeFIPPIDataRecord_Cubature(Va, Qa, Vb, Qb,
                                     NeedDerivatives,
                                     FDR);
     return;
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
                                  Vb[1], Vb[2], Qa, Qb);
                   );

  /* FIXME this one can be computed analytically */ 
  FDR->hDotR0 =real( TaylorMaster(WhichCase, TM_RP, TM_DOT, 0.0, 
                                  Va[0], Va[1], Va[2], 
                                  Vb[1], Vb[2], Qa, Qb);
                   );

  FDR->hDotR1 =real( TaylorMaster(WhichCase, TM_RP, TM_DOT, 1.0, 
                                  Va[0], Va[1], Va[2], 
                                  Vb[1], Vb[2], Qa, Qb);
                   );

  /* FIXME this one can be computed analytically */ 
  FDR->hDotR2 =real( TaylorMaster(WhichCase, TM_RP, TM_DOT, 2.0, 
                                  Va[0], Va[1], Va[2], 
                                  Vb[1], Vb[2], Qa, Qb);
                   );
  
  /*--------------------------------------------------------------*/
  /*- hNabla integrals -------------------------------------------*/
  /*--------------------------------------------------------------*/
  FDR->hNablaRM1 = 4.0*real( TaylorMaster(WhichCase, TM_RP, TM_ONE, -1.0,
                                          Va[0], Va[1], Va[2], 
                                          Vb[1], Vb[2], Qa, Qb);
                           );

  /* FIXME this one can be computed analytically */ 
  FDR->hNablaR0  = 4.0*real( TaylorMaster(WhichCase, TM_RP, TM_ONE, 0.0, 
                                          Va[0], Va[1], Va[2], 
                                          Vb[1], Vb[2], Qa, Qb);
                           );

  FDR->hNablaR1  = 4.0*real( TaylorMaster(WhichCase, TM_RP, TM_ONE, 1.0, 
                                          Va[0], Va[1], Va[2], 
                                          Vb[1], Vb[2], Qa, Qb);
                           );

  /* FIXME this one can be computed analytically */ 
  FDR->hNablaR2  = 4.0*real( TaylorMaster(WhichCase, TM_RP, TM_ONE, 2.0, 
                                          Va[0], Va[1], Va[2], 
                                          Vb[1], Vb[2], Qa, Qb);
                           );
  
  /*--------------------------------------------------------------*/
  /*- hTimes integrals -------------------------------------------*/
  /*--------------------------------------------------------------*/
  FDR->hTimesRM3=real( TaylorMaster(WhichCase, TM_RP, TM_CROSS, -3.0,
                                    Va[0], Va[1], Va[2],
                                    Vb[1], Vb[2], Qa, Qb);
                     );

  FDR->hTimesRM1=real( TaylorMaster(WhichCase, TM_RP, TM_CROSS, -1.0,
                                    Va[0], Va[1], Va[2],
                                    Vb[1], Vb[2], Qa, Qb);
                     );

  /* FIXME this one can be computed analytically */ 
  FDR->hTimesR0=real( TaylorMaster(WhichCase, TM_RP, TM_CROSS, 0.0,
                                   Va[0], Va[1], Va[2],
                                   Vb[1], Vb[2], Qa, Qb);
                    );

  FDR->hTImesR1=real( TaylorMaster(WhichCase, TM_RP, TM_CROSS, 1.0,
                                   Va[0], Va[1], Va[2],
                                   Vb[1], Vb[2], Qa, Qb);
                    );

}
