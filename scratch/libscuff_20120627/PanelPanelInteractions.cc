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

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
extern int ForceNCV=-1; /* 20120625 */
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

namespace scuff {

#define II cdouble(0.0,1.0)

// the 'short-wavelength threshold:' we are in the short-wavelength
// (high-frequency) regime if |k*PanelRadius| > SWTHRESHOLD
#define SWTHRESHOLD 1.0*M_PI

// the 'desingularization radius': if the relative distance  
// between two panels is < this number, we evaluate the      
// panel-panel integrals using desingularization.
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
/*- relative exponential routine -------------------------------*/
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

  cdouble ik=II*Args->k;
  cdouble ik2=ik*ik;

  cdouble H[2]; 

  cdouble *GradH, GradHBuffer[6];
  int NumGradientComponents=Args->NumGradientComponents;
  if (NumGradientComponents>0 && Args->GradH!=0)
   GradH = GradHBuffer;
  else
   GradH = 0;
 
  cdouble *dHdT, dHdTBuffer[6];
  int nta, NumTorqueAxes=Args->NumTorqueAxes;
  double *GammaMatrix=Args->GammaMatrix;
  if (NumTorqueAxes>0 && GammaMatrix!=0 && Args->dHdT!=0)
   dHdT = dHdTBuffer;
  else
   dHdT = 0;

  cdouble HInner[2], GradHInner[6], dHdTInner[6];

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
  double hDot, hNabla=4.0, hTimes; // note hNabla is constant throughout
  cdouble hPlus;
  int np, ncp, npp, ncpp;
  int Mu, Nu;
  double u, v, w, up, vp, wp;
  double X[3], F[3], XP[3], FP[3], R[3], FxFP[3];
  double dX[3], dF[3], Puv, dFxFP[3];
  double r, r2; 
  cdouble Phi, Psi, Zeta;
  memset(H,0,2*sizeof(cdouble));
  if (GradH) memset(GradH,0,6*sizeof(cdouble));
  if (dHdT) memset(dHdT,0,2*NumTorqueAxes*sizeof(cdouble));
  for(np=ncp=0; np<NumPts; np++) 
   { 
     u=TCR[ncp++]; v=TCR[ncp++]; w=TCR[ncp++];

     /***************************************************************/
     /* set X and F=X-Q *********************************************/
     /***************************************************************/
     for(Mu=0; Mu<3; Mu++)
      { X[Mu] = V0[Mu] + u*A[Mu] + v*B[Mu];
        F[Mu] = X[Mu] - Q[Mu];
      };

     /***************************************************************/
     /* inner loop to calculate value of inner integrand ************/
     /***************************************************************/
     memset(HInner,0,2*sizeof(cdouble));
     if (GradH) memset(GradHInner,0,6*sizeof(cdouble));
     if (dHdT) memset(dHdTInner,0,2*NumTorqueAxes*sizeof(cdouble));
     for(npp=ncpp=0; npp<NumPts; npp++)
      { 
        up=TCR[ncpp++]; vp=TCR[ncpp++]; wp=TCR[ncpp++];

        /***************************************************************/ 
        /* set XP and FP=XP-QP *****************************************/
        /***************************************************************/
        for(Mu=0; Mu<3; Mu++)
         { XP[Mu] = V0P[Mu] + up*AP[Mu] + vp*BP[Mu];
           FP[Mu] = XP[Mu] - QP[Mu];
           R[Mu] = X[Mu] - XP[Mu];
         };
      
        /***************************************************************/
        /* inner integrand  ********************************************/
        /***************************************************************/
        r=VecNorm(R);
        r2=r*r;

        /* compute h factors */
        hDot=VecDot(F, FP);
        hPlus=hDot + hNabla/ik2;
        VecCross(F, FP, FxFP);
        hTimes=VecDot(FxFP, R);
   
        /* compute Phi, Psi, Zeta factors */
        if (DeSingularize)
         Phi = ExpRel(ik*r,4) / (4.0*M_PI*r);
        else
         Phi = exp(ik*r) / (4.0*M_PI*r);
        if ( !isfinite(real(Phi)) ) Phi=0.0;

        // put the cubature weight into Phi since Phi is a factor in
        // all integrand components
        Phi*=wp; 

        Psi = Phi * (ik - 1.0/r) / r;
        Zeta = Phi * (ik2 - 3.0*ik/r + 3.0/r2) / r2;

        // combine h terms with kernel factors as necessary 
        // for the various integrand components
        HInner[0] += hPlus * Phi;
        HInner[1] += hTimes * Psi;

        if ( GradH )
         for(Mu=0; Mu<3; Mu++)
           { GradHInner[2*Mu + 0] += R[Mu]*hPlus*Psi;
             GradHInner[2*Mu + 1] += R[Mu]*hTimes*Zeta + FxFP[Mu]*Psi;
           };

        /* 3. d/dTheta L_{0,1,2} */
        if ( NumTorqueAxes>0 && GammaMatrix!=0 )
         for(nta=0; nta<NumTorqueAxes; nta++)
          { memset(dX,0,3*sizeof(double));
            memset(dF,0,3*sizeof(double));
            for(Mu=0; Mu<3; Mu++)
             for(Nu=0; Nu<3; Nu++)
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

      }; /* for(npp=ncpp=0; npp<NumPts; npp++) */

     /*--------------------------------------------------------------*/
     /*- accumulate contributions to outer integral                  */
     /*--------------------------------------------------------------*/
     H[0]+=w*HInner[0];
     H[1]+=w*HInner[1];
     if (GradH)
      for(Mu=0; Mu<6; Mu++)
       GradH[Mu]+=w*GradHInner[Mu];
     if (dHdT)
      for(Mu=0; Mu<2*NumTorqueAxes; Mu++)
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
  RWGPanel *Pa              = Args->Pa;
  RWGPanel *Pb              = Args->Pb;
  int npa                   = Args->npa; 
  int npb                   = Args->npb; 
  int iQa                   = Args->iQa; 
  int iQb                   = Args->iQb; 
  double *VerticesA         = Args->VerticesA;
  double *VerticesB         = Args->VerticesB;
  cdouble k                 = Args->k;
  int NumGradientComponents = Args->NumGradientComponents;
  int NumTorqueAxes         = Args->NumTorqueAxes;
  cdouble *H                = Args->H;
  cdouble *GradH            = Args->GradH;
  cdouble *dHdT             = Args->dHdT;

  /***************************************************************/
  /* extract panel vertices, detect common vertices, measure     */
  /* relative distance                                           */
  /***************************************************************/
  double *Qa   = VerticesA + 3*Pa->VI[iQa];
  double *Qb   = VerticesB + 3*Pb->VI[iQb];
  double *Va[3], *Vb[3];
  double rRel; 

  int ncv=AssessPanelPair(Pa,VerticesA,Pb,VerticesB,&rRel,Va,Vb);

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
if (ForceNCV>-1 && ncv!=ForceNCV) /*  20120625 */
 { H[0]=H[1]=0.0; return; } /*  20120625 */
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  /***************************************************************/
  /* if the panels are far apart then just use simple low-order  */
  /* non-desingularized cubature                                 */
  /***************************************************************/
  if ( rRel > DESINGULARIZATION_RADIUS )
   { GetPPIs_Cubature(Args, 0, 0, Va, Qa, Vb, Qb);
     return;
   };

  /***************************************************************/
  /* determine if we are in the short-wavelength regime          */
  /***************************************************************/
  int InSWRegime = abs(k*fmax(Pa->Radius, Pb->Radius)) > SWTHRESHOLD;

  /***************************************************************/
  /* if we are in the short-wavelength regime and there are no   */
  /* common vertices then we use high-order non-adaptive cubature*/
  /***************************************************************/
  if ( InSWRegime && ncv==0 )
   { GetPPIs_Cubature(Args, 0, 1, Va, Qa, Vb, Qb);
     return; 
   };

  /***************************************************************/
  /* figure out if we are in one of the situations in which      */
  /* we want to switch off immediately to the taylor-duffy method*/
  /***************************************************************/
  if ( ncv==3 || ( ncv>0 && (InSWRegime || Args->ForceTaylorDuffy) ) )
   { 
     TaylorDuffyArgStruct TDArgStruct, *TDArgs=&TDArgStruct;
     InitTaylorDuffyArgs(TDArgs);

     TDArgs->WhichCase=ncv;
     TDArgs->GParam=k;
     TDArgs->V1=Va[0];
     TDArgs->V2=Va[1];
     TDArgs->V3=Va[2];
     TDArgs->V2P=Vb[1];
     TDArgs->V3P=Vb[2];
     TDArgs->Q=Qa;
     TDArgs->QP=Qb;

     // loosen the default taylor-duffy tolerance in the 
     // short-wavelength regime
     if (InSWRegime)
      TDArgs->RelTol=RWGGeometry::SWPPITol;

     TDArgs->WhichG=TM_EIKR_OVER_R;
     TDArgs->WhichH=TM_DOTPLUS;
     H[0]=TaylorDuffy(TDArgs);

     if (ncv==3)
      H[1]=0.0; /* 'C' kernel integral vanishes for the common-triangle case */
     else
      { TDArgs->WhichG=TM_GRADEIKR_OVER_R;
        TDArgs->WhichH=TM_CROSS;
TDArgs->AbsTol = RWGGeometry::SWPPITol*abs(H[0]); // TODO explain me
        H[1]=TaylorDuffy(TDArgs);
      };

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
}

} // namespace scuff
