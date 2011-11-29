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
   double V0MC[3], A[3], B[3];
   double V0MCP[3], AP[3], BP[3];
   double R0[3];
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

  double *V0MC=CFDRD->V0MC;
  double *A=CFDRD->A;
  double *B=CFDRD->B;

  double *V0MCP=CFDRD->V0MCP;
  double *AP=CFDRD->AP;
  double *BP=CFDRD->BP;

  double *R0=CFDRD->R0;

  int NeedDerivatives=CFDRD->NeedDerivatives;

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
  double YA[3], YB[3], R[3];
  double r, r2=0.0, r3, r5;
  int Mu;
  for (Mu=0; Mu<3; Mu++)
   { YA[Mu] = V0MC[Mu]  + u*A[Mu]   + v*B[Mu];
     YB[Mu] = V0MCP[Mu] + up*AP[Mu] + vp*BP[Mu];
     R[Mu]  = R0[Mu] + YA[Mu] - YB[Mu];
     r2    += R[Mu]*R[Mu];
   };
  r=sqrt(r2);
  r3=r*r2;
  r5=r3*r2;

  double YAdYB=VecDot(YA, YB);
  double YAxYB[3], YAmYB[3];
  VecCross(YA, YB, YAxYB);
  VecSub(YA, YB, YAmYB);

  /*--------------------------------------------------------------*/
  /*- assemble output vector -------------------------------------*/
  /*--------------------------------------------------------------*/
  int rp, nf;
  double RPowers[6]={uup/r5, uup/r3, uup/r, uup, uup*r, uup*r2};
  nf=0;
  for(rp=2; rp<6; rp++)
   { fval[nf++] = RPowers[rp]*YAdYB;
     fval[nf++] = RPowers[rp]*YA[0];
     fval[nf++] = RPowers[rp]*YA[1];
     fval[nf++] = RPowers[rp]*YA[2];
     fval[nf++] = RPowers[rp]*YB[0];
     fval[nf++] = RPowers[rp]*YB[1];
     fval[nf++] = RPowers[rp]*YB[2];
     fval[nf++] = RPowers[rp];
   };

  for(rp=1; rp<5; rp++)
   { fval[nf++] = RPowers[rp]*YAmYB[0];
     fval[nf++] = RPowers[rp]*YAmYB[1];
     fval[nf++] = RPowers[rp]*YAmYB[2];
     fval[nf++] = RPowers[rp]*YAxYB[0];
     fval[nf++] = RPowers[rp]*YAxYB[1];
     fval[nf++] = RPowers[rp]*YAxYB[2];
   };

  if (NeedDerivatives)
   { 
     for(Mu=0; Mu<3; Mu++)
      for(rp=1; rp<5; rp++)
       { fval[nf++] = R[Mu] * RPowers[rp] * YAdYB;
         fval[nf++] = R[Mu] * RPowers[rp] * YA[0];
         fval[nf++] = R[Mu] * RPowers[rp] * YA[1];
         fval[nf++] = R[Mu] * RPowers[rp] * YA[2];
         fval[nf++] = R[Mu] * RPowers[rp] * YB[0];
         fval[nf++] = R[Mu] * RPowers[rp] * YB[1];
         fval[nf++] = R[Mu] * RPowers[rp] * YB[2];
         fval[nf++] = R[Mu] * RPowers[rp];
       };

     fval[nf++] = RPowers[1]*YA[0];
     fval[nf++] = RPowers[1]*YA[1];
     fval[nf++] = RPowers[1]*YA[2];
     fval[nf++] = RPowers[1]*YB[0];
     fval[nf++] = RPowers[1]*YB[1];
     fval[nf++] = RPowers[1]*YB[2];
     fval[nf++] = RPowers[1];
   
     for(Mu=0; Mu<3; Mu++)
      { fval[nf++] = R[Mu] * RPowers[0] * YAmYB[0];
        fval[nf++] = R[Mu] * RPowers[0] * YAmYB[1];
        fval[nf++] = R[Mu] * RPowers[0] * YAmYB[2];
        fval[nf++] = R[Mu] * RPowers[0] * YAxYB[0];
        fval[nf++] = R[Mu] * RPowers[0] * YAxYB[1];
        fval[nf++] = R[Mu] * RPowers[0] * YAxYB[2];

        fval[nf++] = R[Mu] * RPowers[1] * YAmYB[0];
        fval[nf++] = R[Mu] * RPowers[1] * YAmYB[1];
        fval[nf++] = R[Mu] * RPowers[1] * YAmYB[2];
        fval[nf++] = R[Mu] * RPowers[1] * YAxYB[0];
        fval[nf++] = R[Mu] * RPowers[1] * YAxYB[1];
        fval[nf++] = R[Mu] * RPowers[1] * YAxYB[2];

        fval[nf++] = R[Mu] * RPowers[3] * YAmYB[0];
        fval[nf++] = R[Mu] * RPowers[3] * YAmYB[1];
        fval[nf++] = R[Mu] * RPowers[3] * YAmYB[2];
        fval[nf++] = R[Mu] * RPowers[3] * YAxYB[0];
        fval[nf++] = R[Mu] * RPowers[3] * YAxYB[1];
        fval[nf++] = R[Mu] * RPowers[3] * YAxYB[2];
     };
   };

} 

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
void ComputeFIPPIDataRecord_Cubature(double **Va, double **Vb,
                                     int NeedDerivatives,
                                     FIPPIDataRecord *FDR)
{
  FDR->HaveDerivatives=NeedDerivatives;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double CentroidA[3], CentroidB[3];

  CentroidA[0] = (Va[0][0] + Va[1][0] + Va[2][0])/3.0;
  CentroidA[1] = (Va[0][1] + Va[1][1] + Va[2][1])/3.0;
  CentroidA[2] = (Va[0][2] + Va[1][2] + Va[2][2])/3.0;
  CentroidB[0] = (Vb[0][0] + Vb[1][0] + Vb[2][0])/3.0;
  CentroidB[1] = (Vb[0][1] + Vb[1][1] + Vb[2][1])/3.0;
  CentroidB[2] = (Vb[0][2] + Vb[1][2] + Vb[2][2])/3.0;
  
  /*--------------------------------------------------------------*/
  /* fill in the data structure used to pass parameters to the    */
  /* integrand routine                                            */
  /*--------------------------------------------------------------*/
  CFDRData MyCFDRData, *CFDRD=&MyCFDRData;

  VecSub(Va[0], CentroidA, CFDRD->V0MC);
  VecSub(Va[1], Va[0], CFDRD->A);
  VecSub(Va[2], Va[1], CFDRD->B);

  VecSub(Vb[0], CentroidB, CFDRD->V0MCP);
  VecSub(Vb[1], Vb[0], CFDRD->AP);
  VecSub(Vb[2], Vb[1], CFDRD->BP);

  VecSub(CentroidA, CentroidB, CFDRD->R0);

  CFDRD->NeedDerivatives = NeedDerivatives;

  /*--------------------------------------------------------------*/
  /*- evaluate the adaptive cubature over the pair of triangles  -*/
  /*--------------------------------------------------------------*/
  double Lower[4]={0.0, 0.0, 0.0, 0.0};
  double Upper[4]={1.0, 1.0, 1.0, 1.0};
  int fdim = NeedDerivatives ? 217 : 48;
  double F[fdim], E[fdim];
CFDRD->nCalls=0;
  adapt_integrate(fdim, CFDRIntegrand, CFDRD, 4, Lower, Upper,
                  0, ABSTOL, RELTOL, F, E);
printf("FIPPI cubature: %i calls\n",CFDRD->nCalls);

  /*--------------------------------------------------------------*/
  /*- unpack the results into the output data record -------------*/
  /*--------------------------------------------------------------*/
  int nf=0;

  FDR->YAdYB_RM1    = F[nf++];
  FDR->YA_RM1[0]    = F[nf++];
  FDR->YA_RM1[1]    = F[nf++];
  FDR->YA_RM1[2]    = F[nf++];
  FDR->YB_RM1[0]    = F[nf++];
  FDR->YB_RM1[1]    = F[nf++];
  FDR->YB_RM1[2]    = F[nf++];
  FDR->RM1          = F[nf++];

  FDR->YAdYB_R0     = F[nf++];
  FDR->YA_R0[0]     = F[nf++];
  FDR->YA_R0[1]     = F[nf++];
  FDR->YA_R0[2]     = F[nf++];
  FDR->YB_R0[0]     = F[nf++];
  FDR->YB_R0[1]     = F[nf++];
  FDR->YB_R0[2]     = F[nf++];
  FDR->R0           = F[nf++];

  FDR->YAdYB_R1     = F[nf++];
  FDR->YA_R1[0]     = F[nf++];
  FDR->YA_R1[1]     = F[nf++];
  FDR->YA_R1[2]     = F[nf++];
  FDR->YB_R1[0]     = F[nf++];
  FDR->YB_R1[1]     = F[nf++];
  FDR->YB_R1[2]     = F[nf++];
  FDR->R1           = F[nf++];

  FDR->YAdYB_R2     = F[nf++];
  FDR->YA_R2[0]     = F[nf++];
  FDR->YA_R2[1]     = F[nf++];
  FDR->YA_R2[2]     = F[nf++];
  FDR->YB_R2[0]     = F[nf++];
  FDR->YB_R2[1]     = F[nf++];
  FDR->YB_R2[2]     = F[nf++];
  FDR->R2           = F[nf++];

  FDR->YAmYB_RM3[0] = F[nf++];
  FDR->YAmYB_RM3[1] = F[nf++];
  FDR->YAmYB_RM3[2] = F[nf++];
  FDR->YAxYB_RM3[0] = F[nf++];
  FDR->YAxYB_RM3[1] = F[nf++];
  FDR->YAxYB_RM3[2] = F[nf++];

  FDR->YAmYB_RM1[0] = F[nf++];
  FDR->YAmYB_RM1[1] = F[nf++];
  FDR->YAmYB_RM1[2] = F[nf++];
  FDR->YAxYB_RM1[0] = F[nf++];
  FDR->YAxYB_RM1[1] = F[nf++];
  FDR->YAxYB_RM1[2] = F[nf++];

  FDR->YAmYB_R0[0]  = F[nf++];
  FDR->YAmYB_R0[1]  = F[nf++];
  FDR->YAmYB_R0[2]  = F[nf++];
  FDR->YAxYB_R0[0]  = F[nf++];
  FDR->YAxYB_R0[1]  = F[nf++];
  FDR->YAxYB_R0[2]  = F[nf++];

  FDR->YAmYB_R1[0]  = F[nf++];
  FDR->YAmYB_R1[1]  = F[nf++];
  FDR->YAmYB_R1[2]  = F[nf++];
  FDR->YAxYB_R1[0]  = F[nf++];
  FDR->YAxYB_R1[1]  = F[nf++];
  FDR->YAxYB_R1[2]  = F[nf++];

  if (NeedDerivatives)
   {
     int Mu;
     for(Mu=0; Mu<3; Mu++)
      { 
        FDR->Ri_YAdYB_RM3[Mu]  = F[nf++];
        FDR->Ri_YA_RM3[3*Mu+0] = F[nf++]; 
        FDR->Ri_YA_RM3[3*Mu+1] = F[nf++]; 
        FDR->Ri_YA_RM3[3*Mu+2] = F[nf++]; 
        FDR->Ri_YB_RM3[3*Mu+0] = F[nf++]; 
        FDR->Ri_YB_RM3[3*Mu+1] = F[nf++]; 
        FDR->Ri_YB_RM3[3*Mu+2] = F[nf++]; 
        FDR->Ri_RM3[Mu]        = F[nf++]; 

        FDR->Ri_YAdYB_RM1[Mu]  = F[nf++];
        FDR->Ri_YA_RM1[3*Mu+0] = F[nf++]; 
        FDR->Ri_YA_RM1[3*Mu+1] = F[nf++]; 
        FDR->Ri_YA_RM1[3*Mu+2] = F[nf++]; 
        FDR->Ri_YB_RM1[3*Mu+0] = F[nf++]; 
        FDR->Ri_YB_RM1[3*Mu+1] = F[nf++]; 
        FDR->Ri_YB_RM1[3*Mu+2] = F[nf++]; 
        FDR->Ri_RM1[Mu]        = F[nf++]; 

        FDR->Ri_YAdYB_R0[Mu]   = F[nf++];
        FDR->Ri_YA_R0[3*Mu+0]  = F[nf++]; 
        FDR->Ri_YA_R0[3*Mu+1]  = F[nf++]; 
        FDR->Ri_YA_R0[3*Mu+2]  = F[nf++]; 
        FDR->Ri_YB_R0[3*Mu+0]  = F[nf++]; 
        FDR->Ri_YB_R0[3*Mu+1]  = F[nf++]; 
        FDR->Ri_YB_R0[3*Mu+2]  = F[nf++]; 
        FDR->Ri_R0[Mu]         = F[nf++]; 

        FDR->Ri_YAdYB_R1[Mu]   = F[nf++];
        FDR->Ri_YA_R1[3*Mu+0]  = F[nf++]; 
        FDR->Ri_YA_R1[3*Mu+1]  = F[nf++]; 
        FDR->Ri_YA_R1[3*Mu+2]  = F[nf++]; 
        FDR->Ri_YB_R1[3*Mu+0]  = F[nf++]; 
        FDR->Ri_YB_R1[3*Mu+1]  = F[nf++]; 
        FDR->Ri_YB_R1[3*Mu+2]  = F[nf++]; 
        FDR->Ri_R1[Mu]         = F[nf++]; 
      };

     FDR->YA_RM3[0] = F[nf++];
     FDR->YA_RM3[1] = F[nf++];
     FDR->YA_RM3[2] = F[nf++];
     FDR->YB_RM3[0] = F[nf++];
     FDR->YB_RM3[1] = F[nf++];
     FDR->YB_RM3[2] = F[nf++];
     FDR->RM3       = F[nf++];

     for(Mu=0; Mu<3; Mu++)
      { 
        FDR->Ri_YAmYB_RM5[3*Mu + 0] = F[nf++]; 
        FDR->Ri_YAmYB_RM5[3*Mu + 1] = F[nf++]; 
        FDR->Ri_YAmYB_RM5[3*Mu + 2] = F[nf++]; 
        FDR->Ri_YAxYB_RM5[3*Mu + 0] = F[nf++]; 
        FDR->Ri_YAxYB_RM5[3*Mu + 1] = F[nf++]; 
        FDR->Ri_YAxYB_RM5[3*Mu + 2] = F[nf++]; 

        FDR->Ri_YAmYB_RM3[3*Mu + 0] = F[nf++]; 
        FDR->Ri_YAmYB_RM3[3*Mu + 1] = F[nf++]; 
        FDR->Ri_YAmYB_RM3[3*Mu + 2] = F[nf++]; 
        FDR->Ri_YAxYB_RM3[3*Mu + 0] = F[nf++]; 
        FDR->Ri_YAxYB_RM3[3*Mu + 1] = F[nf++]; 
        FDR->Ri_YAxYB_RM3[3*Mu + 2] = F[nf++]; 

        FDR->Ri_YAmYB_R0[3*Mu + 0]  = F[nf++]; 
        FDR->Ri_YAmYB_R0[3*Mu + 1]  = F[nf++]; 
        FDR->Ri_YAmYB_R0[3*Mu + 2]  = F[nf++]; 
        FDR->Ri_YAxYB_R0[3*Mu + 0]  = F[nf++]; 
        FDR->Ri_YAxYB_R0[3*Mu + 1]  = F[nf++]; 
        FDR->Ri_YAxYB_R0[3*Mu + 2]  = F[nf++]; 

      };
   };

}

/*--------------------------------------------------------------*/
/*- note: on entry it is assumed that a call to AssessPanelPair */
/*- has already been made and that Va, Vb are in the order      */
/*- they were put in by that routine                            */
/*--------------------------------------------------------------*/
FIPPIDataRecord *ComputeFIPPIDataRecord_TaylorDuffy(int ncv, 
                                                    double **Va, double **Vb, 
                                                    FIPPIDataRecord *FDR)
{ 

  FDR->HaveDerivatives=0;

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
   ComputeFIPPIDataRecord_Cubature(Va, Vb, NeedDerivatives, FDR);
  else
   ComputeFIPPIDataRecord_TaylorDuffy(ncv, Va, Vb, FDR);

}
