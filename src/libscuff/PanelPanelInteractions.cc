/*
 * PanelPanelInteractions.cc -- routines for evaluating interactions
 *                           -- between individual RWG panels
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
#include "TaylorMaster.h"

#define II cdouble(0.0,1.0)

// the 'short-wavelength threshold:' we are in the short-wavelength
// (high-frequency) regime if |k*PanelRadius| > SWTHRESHOLD
#define SWTHRESHOLD 1.0*M_PI

#define AA0 1.0
#define AA1 1.0
#define AA2 (1.0/2.0)
#define AA3 (1.0/6.0)
#define BB0 (-1.0)
#define BB2 (1.0/2.0)
#define BB3 (1.0/3.0)
#define BB4 (1.0/6.0)
#define CC0 3.0
#define CC2 (-1.0/2.0)
#define CC5 (1.0/6.0)

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
  double hDot, hNabla, hTimes; // note hNabla is constant throughout
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
        if ( !finite(real(Phi)) ) Phi=0.0;

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
           { GradHInner[3*Mu + 0] += R[Mu]*hPlus*Psi;
             GradHInner[3*Mu + 1] += R[Mu]*hTimes*Zeta + FxFP[Mu]*Psi; 
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
            dHdTInner[3*nta + 0] += hPlus*Puv*Psi + VecDot(dF,FP)*Phi;
            dHdTInner[3*nta + 1] += hTimes*Puv*Zeta 
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
  RWGObject *Oa             = Args->Oa;
  RWGObject *Ob             = Args->Ob;
  int npa                   = Args->npa; 
  int npb                   = Args->npb; 
  int iQa                   = Args->iQa; 
  int iQb                   = Args->iQb; 
  cdouble k                 = Args->k;
  int NumGradientComponents = Args->NumGradientComponents;
  int NumTorqueAxes         = Args->NumTorqueAxes;
  double *GammaMatrix       = Args->GammaMatrix;
  cdouble *H                = Args->H;
  cdouble *GradH            = Args->GradH;
  cdouble *dHdT             = Args->dHdT;

  /***************************************************************/
  /* extract panel vertices, detect common vertices, measure     */
  /* relative distance                                           */
  /***************************************************************/
  RWGPanel *Pa = Oa->Panels[npa];
  RWGPanel *Pb = Ob->Panels[npb];
  double *Qa   = Oa->Vertices + 3*Pa->VI[iQa];
  double *Qb   = Ob->Vertices + 3*Pb->VI[iQb];
  double *Va[3], *Vb[3];
  double rRel; 
  int ncv=AssessPanelPair(Oa,npa,Ob,npb,&rRel,Va,Vb);

  /***************************************************************/
  /* if the panels are far apart then just use simple low-order  */
  /* non-desingularized cubature.                                */
  /***************************************************************/
  if ( rRel > DESINGULARIZATION_RADIUS )
   { GetPPIs_Cubature(Args, 0, 0, Va, Qa, Vb, Qb);
     return;
   };

  /***************************************************************/
  /* if the panels are identical then we use taylor's scheme for */
  /* the full panel integral no matter what the frequency        */
  /***************************************************************/
  if( ncv==3 )
   { 
     Args->H[0]=TaylorMaster(TM_COMMONTRIANGLE, TM_EIKR_OVER_R, TM_DOTPLUS, k,
                             Va[0], Va[1], Va[2], Vb[1], Vb[2], Qa, Qb);

     Args->H[1]=0.0; /* 'H_\times' vanishes for the common-triangle case */
   
     if (GradH) memset(GradH, 0, 2*NumGradientComponents*sizeof(cdouble));
     if (dHdT)  memset(dHdT,  0, 2*NumTorqueAxes*sizeof(cdouble));

     return;
   };

  /*****************************************************************/
  /* if we are in the high-frequency regime then desingularization */
  /* doesn't work; in this case, if there are any common vertices  */
  /* we use taylor's method for the full panel integral, and       */
  /* otherwise we use high-order naive cubature.                   */
  /*****************************************************************/
  if ( abs(k*fmax(Pa->Radius, Pb->Radius)) > SWTHRESHOLD )
   { 
     if( ncv==2 )
      { 
        /*--------------------------------------------------------------*/
        /* common-edge case                                             */
        /*--------------------------------------------------------------*/
        Args->H[0]=TaylorMaster(TM_COMMONEDGE, TM_EIKR_OVER_R, TM_DOTPLUS, k,
                                Va[0], Va[1], Va[2], Vb[1], Vb[2], Qa, Qb);

        Args->H[1]=TaylorMaster(TM_COMMONEDGE, TM_EIKR_OVER_R, TM_CROSS, k,
                                Va[0], Va[1], Va[2], Vb[1], Vb[2], Qa, Qb);

        if (GradH) memset(GradH, 2, 2*NumGradientComponents*sizeof(cdouble));
        if (dHdT)  memset(dHdT, 2, 2*NumTorqueAxes*sizeof(cdouble));

        return;
      }
     else if( ncv==1 )
      { 
        /*--------------------------------------------------------------*/
        /* common-vertex case                                           */
        /*--------------------------------------------------------------*/
        Args->H[0]=TaylorMaster(TM_COMMONVERTEX, TM_EIKR_OVER_R, TM_DOTPLUS, k,
                                Va[0], Va[1], Va[2], Vb[1], Vb[2], Qa, Qb);

        Args->H[1]=TaylorMaster(TM_COMMONVERTEX, TM_EIKR_OVER_R, TM_CROSS, k,
                                Va[0], Va[1], Va[2], Vb[1], Vb[2], Qa, Qb);

        if (GradH) memset(GradH, 2, 2*NumGradientComponents*sizeof(cdouble));
        if (dHdT)  memset(dHdT, 2, 2*NumTorqueAxes*sizeof(cdouble));

        return;
      }
     else // ncv==0
      { 
        /*--------------------------------------------------------------*/
        /* no common vertices                                           */
        /*--------------------------------------------------------------*/
        GetPPIs_Cubature(Args, 0, 1, Va, Qa, Vb, Qb);
        return;
      };

   };

  /*****************************************************************/
  /* if the user requested angular derivatives, then we make a     */
  /* first call to GetPPIs_Cubature *without* desingularization    */
  /* just to get the angular derivative integrals, because         */
  /* desingularization of angular derivative integrals is not      */
  /* implemented and probably never will be; for this call we      */
  /* use high-order cubature but then promptly throw away the H    */
  /* and GradH integrals since we proceed to compute those using   */
  /* desingularization below.                                      */
  /*****************************************************************/
  cdouble dHdTSave[6];
  if (NumTorqueAxes>0)
   { GetPPIs_Cubature(Args, 0, 1, Va, Qa, Vb, Qb);
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
  int NeedDerivatives=NumGradientComponents>0;
  FIPPIDataTable *FIPPIDT=0;
  FIPPIDataRecord MyFDR, *FDR;
  if (FIPPIDT)
   FDR=FIPPIDT->GetFIPPIDataRecord(Va, Qa, Vb, Qb, NeedDerivatives);
  else
   FDR=ComputeFIPPIDataRecord(Va, Qa, Vb, Qb, NeedDerivatives, &MyFDR);

  // step 3
  // note: PF[n] = (ik)^n / (4\pi)
  cdouble ik=II*k, OOIK2=1.0/(ik*ik);
  cdouble PF[5];
  PF[0]=1.0/(4.0*M_PI);
  PF[1]=ik*PF[0];
  PF[2]=ik*PF[1];
  PF[3]=ik*PF[2];
  PF[4]=ik*PF[3];
  PF[5]=ik*PF[4];

  // add contributions to panel-panel integrals 
  Args->H[0] +=  PF[0]*AA0*(FDR->hDotRM1 + OOIK2*FDR->hNablaRM1)
                +PF[1]*AA1*(FDR->hDotR0  + OOIK2*FDR->hNablaR0 )
                +PF[2]*AA2*(FDR->hDotR1  + OOIK2*FDR->hNablaR1 )
                +PF[3]*AA3*(FDR->hDotR2  + OOIK2*FDR->hNablaR2 );
  
  Args->H[1] +=  PF[0]*BB0*FDR->hTimesRM3
                +PF[2]*BB2*FDR->hTimesRM1
                +PF[3]*BB3*FDR->hTimesR0 
                +PF[4]*BB4*FDR->hTimesR1;

  // restore angular derivatives as necessary 
  if (NumTorqueAxes>0)
   memcpy(Args->dHdT, dHdTSave, 2*NumTorqueAxes*sizeof(cdouble));

  if (NumGradientComponents==0)
   return;

  // add contributions to gradients of panel-panel integrals
  double Rab[3];
  VecSub(Pa->Centroid, Pb->Centroid, Rab);

  cdouble GradH0Scalar, GradH1Scalar;
  GradH0Scalar=  PF[0]*BB0*(FDR->hDotRM3 + OOIK2*FDR->hNablaRM3)
                +PF[2]*BB2*(FDR->hDotRM1 + OOIK2*FDR->hNablaRM1)
                +PF[3]*BB3*(FDR->hDotR0  + OOIK2*FDR->hNablaR0)
                +PF[4]*BB4*(FDR->hDotR1  + OOIK2*FDR->hNablaR1);

  GradH1Scalar=  PF[0]*CC0*FDR->hTimesRM5
                +PF[2]*CC2*FDR->hTimesRM3
                +PF[5]*CC5*FDR->hTimesR0;

  int Mu;
  for(Mu=0; Mu<NumGradientComponents; Mu++)
   { Args->GradH[2*Mu + 0] += Rab[Mu]*GradH0Scalar;
     Args->GradH[2*Mu + 1] += Rab[Mu]*GradH1Scalar
                               +PF[0]*BB0*FDR->dhTimesdRMuRM3[Mu]
                               +PF[2]*BB2*FDR->dhTimesdRMuRM1[Mu]
                               +PF[3]*BB3*FDR->dhTimesdRMuR0[Mu]
                               +PF[4]*BB3*FDR->dhTimesdRMuR1[Mu];
   };

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
  Args->GammaMatrix=0;
}
