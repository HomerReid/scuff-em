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
 * PanelPanelInteractions.cc -- libscuff routines for evaluating interactions
 *                           -- between individual panels
 * 
 * homer reid                -- 11/2005 -- 10/2011
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
/* short-wavelength thresholds: we are in the 'short-wavelength'      */
/* regime if |k*PanelRadius| > SWTHRESHOLD, and we are in the         */
/* 'very-short-wavelength' regime if |k*PanelRadius|>VERYSWTHRESHOLD. */
/*                                                                    */
/* note:                                                              */
/* (a) 'short-wavelength' means the low-k taylor-series-expansion     */
/*     method (i.e. the desingularization method) doesn't work, in    */
/*     which case, if there are any common vertices, we have to use   */
/*     the taylor-duffy method for PPIs.                              */
/*                                                                    */
/* (b) 'very-short-wavelength' means that within the taylor-duffy     */
/*     computation we can use the k->infty limit of the Helmholtz     */
/*     kernel, which significantly accelerates that computation.      */
/**********************************************************************/
#define SWTHRESHOLD 1.0*M_PI
#define VERYSWTHRESHOLD 20.0

/**********************************************************************/
// the 'desingularization radius': if the relative distance between   */
// between two panels is < this number, we evaluate the PPIS using    */
// the low-k taylor-series expansion (desingularization).             */
/**********************************************************************/
#define DESINGULARIZATION_RADIUS 4.0

#define AA0 1.0
#define AA1 1.0
#define AA2 (1.0/2.0)
#define AA3 (1.0/6.0)
#define BB0 (-1.0)
#define BB2 (1.0/2.0)
#define BB3 (1.0/3.0)
#define BB4 (1.0/6.0)

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- relative exponential routine. this is not the same routine -*/
/*- as ExpRelV3P0 defined in TaylorDuffy(); the two computations*/
/*- have different normalizations.                              */
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
#define EXPRELTOL  1.0e-8
#define EXPRELTOL2 EXPRELTOL*EXPRELTOL  
cdouble ExpRel(cdouble x, int n)
{
  int m;
  cdouble Term, Sum;
  double mag2Term, mag2Sum;

  for(Term=1.0, m=1; m<n; m++)
   Term*=x/((double)m);

  for(Sum=0.0 ; m<100; m++)
   { Term*=x/((double)m);
     Sum+=Term;
     mag2Term=norm(Term);
     mag2Sum=norm(Sum);
     if ( mag2Term < EXPRELTOL2*mag2Sum )
      break;
   };
  return Sum;
} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AssembleInnerPPIIntegrand(double wp, double *R, double *X,
                               double *F, double *FP, cdouble k,
                               GBarAccelerator *GBA,
                               int DeSingularize,
                               int NumTorqueAxes, double *GammaMatrix,
                               cdouble *HInner, cdouble *GradHInner, cdouble *dHdTInner)
{ 
  /* compute polynomial factors */
  cdouble ik=II*k, ik2=ik*ik;
  cdouble hPlus = VecDot(F,FP) + 4.0/ik2;
  double FxFP[3];
  VecCross(F, FP, FxFP);

  /* compute the Helmholtz Green's function and its derivatives */
  /* either directly (non-PBC case) or via Ewald summation or   */
  /* interpolation (PBC case)                                   */
  if (GBA)
   { 
     cdouble G, dG[3], ddG[9];
     G=GetGBar(R, GBA, dG, (GradHInner ? ddG : 0 ) );

     HInner[0] += wp * hPlus*G;
     HInner[1] += wp * (FxFP[0]*dG[0] + FxFP[1]*dG[1] + FxFP[2]*dG[2]);

     if (GradHInner) // derivatives of the G integral
      { for(int Mu=0; Mu<3; Mu++)
         { GradHInner[ 2*Mu + 0 ] += wp * hPlus * dG[Mu];
           GradHInner[ 2*Mu + 1 ] += wp * (  FxFP[0]*ddG[3*Mu + 0] 
                                           + FxFP[1]*ddG[3*Mu + 1] 
                                           + FxFP[2]*ddG[3*Mu + 2]);
         };
      };
   }
  else
   { 
     double r2=VecNorm2(R), r=sqrt(r2);
     cdouble Phi;
     if (DeSingularize)
      Phi = ExpRel(ik*r,4) / (4.0*M_PI*r);
     else
      Phi = exp(ik*r) / (4.0*M_PI*r);
     if ( !IsFinite(real(Phi)) ) Phi=0.0;

     // put the cubature weight into Phi since Phi is a factor in
     // all integrand components
     Phi*=wp; 

     cdouble Psi  = Phi * (ik - 1.0/r) / r;
     cdouble Zeta = Phi * (ik2 - 3.0*ik/r + 3.0/r2) / r2;

     double hTimes=VecDot(FxFP, R);

     HInner[0] += hPlus * Phi;
     HInner[1] += hTimes * Psi;

     if ( GradHInner )
      for(int Mu=0; Mu<3; Mu++)
       { GradHInner[2*Mu + 0] += R[Mu]*hPlus*Psi;
         GradHInner[2*Mu + 1] += R[Mu]*hTimes*Zeta + FxFP[Mu]*Psi;
       };
   
     /* angular derivatives */
     double dX[3], dF[3], Puv, dFxFP[3];
     if ( NumTorqueAxes>0 && GammaMatrix!=0 )
      for(int nta=0; nta<NumTorqueAxes; nta++)
       { memset(dX,0,3*sizeof(double));
         memset(dF,0,3*sizeof(double));
         for(int Mu=0; Mu<3; Mu++)
          for(int Nu=0; Nu<3; Nu++)
           { dX[Mu]+=GammaMatrix[9*nta + Mu + 3*Nu]*X[Nu];
             dF[Mu]+=GammaMatrix[9*nta + Mu + 3*Nu]*F[Nu];
           };
         Puv=VecDot(R,dX);
         dHdTInner[2*nta + 0] += hPlus*Puv*Psi + VecDot(dF,FP)*Phi;
         dHdTInner[2*nta + 1] += hTimes*Puv*Zeta 
                               + (  VecDot(VecCross(dF,FP,dFxFP),R) 
                                  + VecDot(FxFP,dX) 
                                 )*Psi;
       }; // for(nta= ... )
   };
}

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- PART 1: Routine to compute panel-panel integrals using     -*/
/*-         fixed-order numerical cubature for both panels.    -*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
void GetPPIs_Cubature(GetPPIArgStruct *Args,
                      int DeSingularize, int HighOrder,
                      double **Va, double *Qa,
                      double **Vb, double *Qb)
{ 
  /***************************************************************/
  /* preliminary setup for numerical cubature.                   */
  /* in what follows, X runs over the 'destination triangle' and */
  /* XP runs over the 'source triangle' according to             */
  /*  X  = Va_1 +  u*(Va_2 - Va_1) +  v*(Va_3-Vb_1)              */
  /*  XP = Vb_1 + up*(Va_2 - Vb_1) + vp*(Va_3-Vb_1)              */
  /* where (V_1, V_2, V_3) are the triangle vertices and (u,v)   */
  /* are the cubature points for a 2D numerical cubature rule    */
  /* over the standard triangle with vertices at (0,0)(1,0)(0,1).*/
  /* note that the jacobian of the transformation is 4*A*AP      */
  /* where A and AP are the areas of the triangles; this         */
  /* conveniently cancels the corresponding factor coming from   */
  /* the RWG basis function prefactor.                           */
  /***************************************************************/
  double *V0, A[3], B[3], *Q;
  V0=Va[0];
  VecSub(Va[1], Va[0], A);
  VecSub(Va[2], Va[0], B);
  Q=Qa;

  double *V0P, AP[3], BP[3], *QP;
  V0P=Vb[0];
  VecSub(Vb[1], Vb[0], AP);
  VecSub(Vb[2], Vb[0], BP);
  QP=Qb;

  cdouble H[2]; 
  cdouble *GradH, GradHBuffer[6];
  int NumGradientComponents=Args->NumGradientComponents;
  if (NumGradientComponents>0 && Args->GradH!=0)
   GradH = GradHBuffer;
  else
   GradH = 0;
 
  cdouble *dHdT, dHdTBuffer[6];
  int NumTorqueAxes=Args->NumTorqueAxes;
  double *GammaMatrix=Args->GammaMatrix;
  if (NumTorqueAxes>0 && GammaMatrix!=0 && Args->dHdT!=0)
   dHdT = dHdTBuffer;
  else
   dHdT = 0;

  cdouble HInner[2], GradHInnerBuffer[6], dHdTInnerBuffer[6];
  cdouble *GradHInner = GradH ? GradHInnerBuffer : 0;
  cdouble *dHdTInner  = dHdT  ? dHdTInnerBuffer  : 0;

  /***************************************************************/
  /* choose order of quadrature scheme to use.                   */
  /* TCR ('triangle cubature rule') points to a vector of 3N     */
  /* doubles (for an N-point cubature rule).                     */
  /* TCR[3*n,3*n+1,3*n+2]=(u,v,w), where (u,v)                   */
  /* are the (x,y) coordinates of the nth quadrature point and w */
  /* is its weight.                                              */
  /* note we use the same quadrature rule for both the source    */
  /* and destination triangles.                                  */
  /***************************************************************/
  double *TCR;
  int NumPts;
  if (HighOrder)
   TCR=GetTCR(20, &NumPts);
  else
   TCR=GetTCR(4, &NumPts);

  /***************************************************************/
  /* outer loop **************************************************/
  /***************************************************************/
  cdouble k = Args->k;
  memset(H,0,2*sizeof(cdouble));
  if (GradH) memset(GradH,0,6*sizeof(cdouble));
  if (dHdT) memset(dHdT,0,2*NumTorqueAxes*sizeof(cdouble));
  for(int np=0, ncp=0; np<NumPts; np++)
   { 
     double u=TCR[ncp++];
     double v=TCR[ncp++];
     double w=TCR[ncp++];

     /***************************************************************/
     /* set X and F=X-Q *********************************************/
     /***************************************************************/
     double X[3], F[3];
     for(int Mu=0; Mu<3; Mu++)
      { X[Mu] = V0[Mu] + u*A[Mu] + v*B[Mu];
        F[Mu] = X[Mu] - Q[Mu];
      };

     /***************************************************************/
     /* inner loop to calculate value of inner integrand ************/
     /***************************************************************/
     memset(HInner,0,2*sizeof(cdouble));
     if (GradH) memset(GradHInner,0,6*sizeof(cdouble));
     if (dHdT) memset(dHdTInner,0,2*NumTorqueAxes*sizeof(cdouble));
     for(int npp=0, ncpp=0; npp<NumPts; npp++)
      { 
        double up=TCR[ncpp++];
        double vp=TCR[ncpp++];
        double wp=TCR[ncpp++];

        /***************************************************************/ 
        /* set XP and FP=XP-QP *****************************************/
        /***************************************************************/
        double XP[3], FP[3], R[3];
        for(int Mu=0; Mu<3; Mu++)
         { XP[Mu] = V0P[Mu] + up*AP[Mu] + vp*BP[Mu];
           FP[Mu] = XP[Mu] - QP[Mu];
           R[Mu] = X[Mu] - XP[Mu];
         };

        AssembleInnerPPIIntegrand(wp, R, X, F, FP, k, Args->GBA,
                                  DeSingularize, NumTorqueAxes, GammaMatrix,
                                  HInner, GradHInner, dHdTInner);

      }; /* for(npp=ncpp=0; npp<NumPts; npp++) */

     /*--------------------------------------------------------------*/
     /*- accumulate contributions to outer integral                  */
     /*--------------------------------------------------------------*/
     H[0]+=w*HInner[0];
     H[1]+=w*HInner[1];
     if (GradH)
      for(int Mu=0; Mu<6; Mu++)
       GradH[Mu]+=w*GradHInner[Mu];
     if (dHdT)
      for(int Mu=0; Mu<2*NumTorqueAxes; Mu++)
       dHdT[Mu]+=w*dHdTInner[Mu];

   }; // for(np=ncp=0; np<nPts; np++) 

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  memcpy(Args->H, H, 2*sizeof(cdouble));
  if (GradH) memcpy(Args->GradH, GradH, 6*sizeof(cdouble));
  if (dHdT) memcpy(Args->dHdT, dHdT, 2*NumTorqueAxes*sizeof(cdouble));

}

/***************************************************************/
/* calculate integrals over a single pair of triangles using   */
/* one of several different methods based on how near the two  */
/* triangles are to each other.                                */
/***************************************************************/
void GetPanelPanelInteractions(GetPPIArgStruct *Args)
{ 
  /***************************************************************/
  /* local copies of fields in argument structure ****************/
  /***************************************************************/
  RWGSurface *Sa            = Args->Sa;
  RWGSurface *Sb            = Args->Sb;
  int npa                   = Args->npa; 
  int npb                   = Args->npb; 
  int iQa                   = Args->iQa; 
  int iQb                   = Args->iQb; 
  cdouble k                 = Args->k;
  int NumGradientComponents = Args->NumGradientComponents;
  int NumTorqueAxes         = Args->NumTorqueAxes;
  double *Displacement      = Args->Displacement;
  cdouble *H                = Args->H;
  cdouble *GradH            = Args->GradH;
  cdouble *dHdT             = Args->dHdT;

  /***************************************************************/
  /* extract panel vertices, detect common vertices, measure     */
  /* relative distance                                           */
  /***************************************************************/
  RWGPanel *Pa = Sa->Panels[npa];
  RWGPanel *Pb = Sb->Panels[npb];
  double *Qa   = Sa->Vertices + 3*Pa->VI[iQa];
  double *Qb   = Sb->Vertices + 3*Pb->VI[iQb];
  double *Va[3], *Vb[3];
  double VbDisplaced[3][3];
  double rRel; 
  int ncv;

  if (Displacement==0)
   ncv=AssessPanelPair(Sa,npa,Sb,npb,&rRel,Va,Vb);
  else 
   { 
     Va[0] = Sa->Vertices + 3*Pa->VI[0];
     Va[1] = Sa->Vertices + 3*Pa->VI[1];
     Va[2] = Sa->Vertices + 3*Pa->VI[2];

     VecScaleAdd(Sb->Vertices + 3*Pb->VI[0], 1.0, Displacement, VbDisplaced[0]);
     VecScaleAdd(Sb->Vertices + 3*Pb->VI[1], 1.0, Displacement, VbDisplaced[1]);
     VecScaleAdd(Sb->Vertices + 3*Pb->VI[2], 1.0, Displacement, VbDisplaced[2]);
     Vb[0] = VbDisplaced[0];
     Vb[1] = VbDisplaced[1];
     Vb[2] = VbDisplaced[2];
     Qb    = VbDisplaced[iQb];

     double DC[3]; // 'delta centroid' 
     DC[0] = Pa->Centroid[0] - Pb->Centroid[0] - Displacement[0];
     DC[1] = Pa->Centroid[1] - Pb->Centroid[1] - Displacement[1];
     DC[2] = Pa->Centroid[2] - Pb->Centroid[2] - Displacement[2];

     double rMax = fmax(Pa->Radius, Pb->Radius);
     rRel = VecNorm(DC) / rMax; 

     ncv=AssessPanelPair(Va, Vb, rMax);
   };

  /***************************************************************/
  /* if the panels are far apart, or if we have an interpolator, */
  /* then just use simple low-order non-desingularized cubature  */
  /***************************************************************/
  if ( Args->GBA || (rRel > DESINGULARIZATION_RADIUS) )
   { Args->WhichAlgorithm=PPIALG_LOCUBATURE;
     GetPPIs_Cubature(Args, 0, 0, Va, Qa, Vb, Qb);
     return;
   };

  /***************************************************************/
  /* determine if we are in the short-wavelength regime          */
  /***************************************************************/
  double kR=abs(k*fmax(Pa->Radius, Pb->Radius));
  int InSWRegime = kR > SWTHRESHOLD;
  int InVerySWRegime = kR > VERYSWTHRESHOLD;

  /***************************************************************/
  /* if we are in the short-wavelength regime and there are no   */
  /* common vertices then we use high-order non-adaptive cubature*/
  /***************************************************************/
  if ( InSWRegime && ncv==0 )
   { Args->WhichAlgorithm=PPIALG_HOCUBATURE;
     GetPPIs_Cubature(Args, 0, 1, Va, Qa, Vb, Qb);
     return; 
   };

  /***************************************************************/
  /* if we are in the short-wavelength regime and there are 1 or */
  /* 2 common vertices, or if we are in any regime and there are */
  /* 3 common vertices, or if the called explicitly requested it,*/
  /* we use the Taylor-Duffy method                              */
  /***************************************************************/
  if ( ncv==3 || ( ncv>0 && (InSWRegime || Args->ForceTaylorDuffy) ) )
   { 
     TaylorDuffyArgStruct TDArgStruct, *TDArgs=&TDArgStruct;
     InitTaylorDuffyArgs(TDArgs);

     int PIndex[3]={TD_UNITY, TD_PMCHWG1, TD_PMCHWC};
     int KIndex[3]={TD_HELMHOLTZ, TD_HELMHOLTZ, TD_GRADHELMHOLTZ};
     cdouble KParam[3]={k,k,k};
     cdouble Result[3], Error[3];

     TDArgs->WhichCase=ncv;
     TDArgs->NumPKs = (ncv==3) ? 2 : 3;
     TDArgs->PIndex=PIndex;
     TDArgs->KIndex=KIndex;
     TDArgs->KParam=KParam;
     TDArgs->V1=Va[0];
     TDArgs->V2=Va[1];
     TDArgs->V3=Va[2];
     TDArgs->V2P=Vb[1];
     TDArgs->V3P=Vb[2];
     TDArgs->Q=Qa;
     TDArgs->QP=Qb;
     TDArgs->Result=Result;
     TDArgs->Error=Error;

     if ( InVerySWRegime && RWGGeometry::UseHighKTaylorDuffy )
      { Args->WhichAlgorithm=PPIALG_HKTD;
        KIndex[0]=TD_HIGHK_HELMHOLTZ;
        KIndex[1]=TD_HIGHK_HELMHOLTZ;
        KIndex[2]=TD_HIGHK_GRADHELMHOLTZ;
      }
     else
      Args->WhichAlgorithm=PPIALG_TD;

     TaylorDuffy(TDArgs);

     H[0] = Result[1] - 4.0*Result[0]/(k*k);
     H[1] = (ncv==3) ? 0.0 : Result[2];

     if (GradH) memset(GradH, 0, 2*NumGradientComponents*sizeof(cdouble));
     if (dHdT)  memset(dHdT,  0, 2*NumTorqueAxes*sizeof(cdouble));
     return;
   };

  /*****************************************************************/
  /* ok, we are in the desingularization regime.                   */
  /* if the user requested derivatives, then we make a first call  */
  /* to GetPPIs_Cubature *without* desingularization to get just   */
  /* the derivative integrals, because desingularization of        */
  /* derivative integrals is not implemented and probably never    */
  /* will be. for this call we use high-order cubature but then    */
  /* promptly throw away the H integrals since we proceed to       */
  /* compute those using the more-accurate desingularization       */
  /* method below.                                                 */
  /*****************************************************************/
  cdouble GradHSave[6], dHdTSave[6];
  Args->WhichAlgorithm=PPIALG_DESING;
  if ( NumGradientComponents>0 || NumTorqueAxes>0 )
   { GetPPIs_Cubature(Args, 0, 1, Va, Qa, Vb, Qb);
     memcpy(GradHSave, Args->GradH, 2*NumGradientComponents*sizeof(cdouble));
     memcpy(dHdTSave, Args->dHdT, 2*NumTorqueAxes*sizeof(cdouble));
   };

  /*****************************************************************/
  /* finally, if we made it here, then we use desingularization:   */
  /*                                                               */
  /*  1) use naive cubature to compute panel-panel integral with   */
  /*     singular terms removed                                    */
  /*  2) get the singular terms either by looking them up in the   */
  /*     table (if one was provided) or just computing them on the */
  /*     fly                                                       */
  /*  3) add the singular and non-singular contributions           */
  /*                                                               */
  /*****************************************************************/
  // step 1
  GetPPIs_Cubature(Args, 1, 0, Va, Qa, Vb, Qb);

  // step 2
  QDFIPPIData MyQDFD, *QDFD=&MyQDFD;
  if (Args->opFC)
   GetQDFIPPIData(Va, Qa, Vb, Qb, ncv, Args->opFC, QDFD);
  else
   GetQDFIPPIData(Va, Qa, Vb, Qb, ncv, &GlobalFIPPICache, QDFD);

  // step 3
  // note: PF[n] = (ik)^n / (4\pi)
  cdouble ik=II*k; 
  cdouble OOIK2=1.0/(ik*ik);
  cdouble PF[5];
  PF[0]=1.0/(4.0*M_PI);
  PF[1]=ik*PF[0];
  PF[2]=ik*PF[1];
  PF[3]=ik*PF[2];
  PF[4]=ik*PF[3];

  // add contributions to panel-panel integrals
  H[0] +=  PF[0]*AA0*( QDFD->hDotRM1 + OOIK2*QDFD->hNablaRM1)
          +PF[1]*AA1*( QDFD->hDotR0  + OOIK2*QDFD->hNablaR0 )
          +PF[2]*AA2*( QDFD->hDotR1  + OOIK2*QDFD->hNablaR1 )
          +PF[3]*AA3*( QDFD->hDotR2  + OOIK2*QDFD->hNablaR2 );
  
  H[1] +=  PF[0]*BB0*QDFD->hTimesRM3
          +PF[2]*BB2*QDFD->hTimesRM1
          +PF[3]*BB3*QDFD->hTimesR0 
          +PF[4]*BB4*QDFD->hTimesR1;

  // restore derivative integrals as necessary 
  if (NumGradientComponents>0)
   memcpy(GradH, GradHSave, 2*NumGradientComponents*sizeof(cdouble));
  if (NumTorqueAxes>0)
   memcpy(dHdT, dHdTSave, 2*NumTorqueAxes*sizeof(cdouble));

} 

/***************************************************************/
/* this is an alternate entry point to GetPanelPanelInteractions*/
/* that copies the results out of the structure body into      */
/* user-specified buffers                                      */
/***************************************************************/
void GetPanelPanelInteractions(GetPPIArgStruct *Args,
                               cdouble *H, 
                               cdouble *GradH, 
                               cdouble *dHdT)
{ 
  GetPanelPanelInteractions(Args);

  memcpy(H, Args->H, 2*sizeof(cdouble));
  if(GradH)  
   memcpy(GradH, Args->GradH, 2*Args->NumGradientComponents*sizeof(cdouble));
  if(dHdT)  
   memcpy(dHdT, Args->dHdT, 2*Args->NumTorqueAxes*sizeof(cdouble));
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void InitGetPPIArgs(GetPPIArgStruct *Args)
{
  Args->NumGradientComponents=0;
  Args->NumTorqueAxes=0;
  Args->ForceTaylorDuffy=0;
  Args->GammaMatrix=0;
  Args->opFC=0;
  Args->Displacement=0;
  Args->GBA=0;
}

} // namespace scuff
