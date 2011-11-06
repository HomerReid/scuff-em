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
#include "TaylorMaster.h"

#define II cdouble(0,1)

// the 'common vertex threshold:' two vertices are considered to be 
// the same if their distance is less than CVTHRESHOLD*the panel radius
#define CVTHRESHOLD 1.0e-6

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

/***************************************************************/
/* this routine gathers some information on the pair of panels */
/* (Oa, npa) -- (Ob,npb).                                      */
/*                                                             */
/* on return from this routine,                                */
/*  a) the return value is the # of common vertices (0,1,2,3)  */
/*  b) *rRel is set to the 'relative distance'                 */
/*  c) Va[0..2] and Vb[0..2] are pointers to the vertices of   */
/*     the two panels                                          */
/*  d) if there are any common vertices, then the ordering of  */
/*     the Va and Vb arrays is such that any common vertices   */
/*     come first; for example, if there are 2 common vertices */
/*     then Va[0] = Vb[0] and Va[1] = Vb[1].                   */
/***************************************************************/
int AssessPanelPair(RWGObject *Oa, int npa, RWGObject *Ob, int npb,
                    double *rRel, double **Va, double **Vb)
{
  RWGPanel *Pa=Oa->Panels[npa];
  RWGPanel *Pb=Ob->Panels[npb];

  Va[0] = Oa->Vertices + 3*Pa->VI[0];
  Va[1] = Oa->Vertices + 3*Pa->VI[1];
  Va[2] = Oa->Vertices + 3*Pa->VI[2];

  Vb[0] = Ob->Vertices + 3*Pb->VI[0];
  Vb[1] = Ob->Vertices + 3*Pb->VI[1];
  Vb[2] = Ob->Vertices + 3*Pb->VI[2];

  double rMax=fmax(Pa->Radius, Pb->Radius);

  *rRel = VecDistance(Pa->Centroid, Pb->Centroid) / fmax(Pa->Radius, Pb->Radius);
  if ( *rRel > 2.0 ) // there can be no common vertices in this case 
   return 0;

  /***************************************************************/
  /* look for common vertices.                                   */
  /* NOTE: in an earlier incarnation of this code, i looked for  */
  /*       common vertices by simple integer comparisons         */
  /*       (comparing indices within a table of vertices), but   */
  /*       i specifically DON'T want to do that here for several */
  /*       reasons. ultimately it would be nice to avoid doing   */
  /*       9 separate comparisons here, although in practice it  */
  /*       won't matter much since most cases will be caught by  */
  /*       the rRel>2 check above.                               */
  /***************************************************************/
  int CVIa[3], CVIb[3];
  int ncv=0;
  int ia, ib;
  double ThresholdDistance2=CVTHRESHOLD*CVTHRESHOLD*rMax*rMax;
  for(ia=0; ia<3; ia++)
   for(ib=0; ib<3; ib++)
    { if ( VecDistance2(Va[ia],Vb[ib]) < ThresholdDistance2 )
       { CVIa[ncv]=ia;    
         CVIb[ncv]=ib;
         ncv++;
       };
    };
  
  /***************************************************************/
  /* if there were any common vertices, reorganize the Va and Vb */
  /* arrays so that common vertices appear first                 */
  /***************************************************************/
  if (ncv==0)
   { // vertices are already in acceptable order
   } 
  else if (ncv==1)
   { Va[0] = Oa->Vertices + 3*Pa->VI[  CVIa[0] ];
     Va[1] = Oa->Vertices + 3*Pa->VI[ (CVIa[0]+1)%3 ];
     Va[2] = Oa->Vertices + 3*Pa->VI[ (CVIa[0]+2)%3 ];
     Vb[0] = Ob->Vertices + 3*Pb->VI[  CVIb[0] ];
     Vb[1] = Ob->Vertices + 3*Pb->VI[ (CVIb[0]+1)%3 ];
     Vb[2] = Ob->Vertices + 3*Pb->VI[ (CVIb[0]+2)%3 ];
   }
  else if (ncv==2)
   { Va[0] = Oa->Vertices + 3*Pa->VI[  CVIa[0] ];
     Va[1] = Oa->Vertices + 3*Pa->VI[  CVIa[1] ];
     Va[2] = Oa->Vertices + 3*Pa->VI[  3-CVIa[0]-CVIa[1] ];
     Vb[0] = Ob->Vertices + 3*Pb->VI[  CVIb[0] ];
     Vb[1] = Ob->Vertices + 3*Pb->VI[  CVIb[1] ];
     Vb[2] = Ob->Vertices + 3*Pb->VI[  3-CVIb[0]-CVIb[1] ];
   }
  else if (ncv==3)
   { Va[0] = Oa->Vertices + 3*Pa->VI[  CVIa[0] ];
     Va[1] = Oa->Vertices + 3*Pa->VI[  CVIa[1] ];
     Va[2] = Oa->Vertices + 3*Pa->VI[  CVIa[2] ];
     Vb[0] = Ob->Vertices + 3*Pb->VI[  CVIb[0] ];
     Vb[1] = Ob->Vertices + 3*Pb->VI[  CVIb[1] ];
     Vb[2] = Ob->Vertices + 3*Pb->VI[  CVIb[2] ];
   };

  return ncv;
 
}

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- desingularized exponential routines ------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
#define EXPRELTOL  1.0e-8
#define EXPRELTOL2 EXPRELTOL*EXPRELTOL  
cdouble cExpRel(cdouble x, int n)
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
void GetPPI_Fixed(GPPIArgStruct *Args, int Desingularize, int HighOrder,
                  double **Va, double *Qa, double **Vb, double *Qb)
#if 0
void GetPPIs_Fixed(cdouble Wavenumber, int NeedCross,
                   int NumTorqueAxes, double *GammaMatrix,
                   double *VD[3], int iQD, double *VS[3], int iQS,
                   int Order, int DeSingularize,
                   cdouble *L, cdouble *GradL, cdouble *dLdT)
#endif
{ 
  int np, ncp, npp, ncpp, m, mu, nu, ri, nta;
  double u, v, w, up, vp, wp;
  double AD[3], BD[3], AS[3], BS[3], *QD, *QS;
  double X[3], XmQD[3], XP[3], XPmQS[3], gXh[3], dgXh[3];
  double *TCR;
  int NumPts;
  double r, r2, R[3], dX[3], dG[3];
  double LTerm[3], Puv;
  cdouble Phi, Psi, Zeta, ik;
  cdouble LInner[3], GradLInner[9], dLdTInner[3*NumTorqueAxes];

  /***************************************************************/
  /* preliminary setup for numerical cubature.                   */
  /* in what follows, X runs over the 'destination triangle' and */
  /* XP runs over the 'source triangle' according to             */
  /*  X  = VD_1 +  u*(VD_2 - VD_1) +  v*(VD_3-VD_1)              */
  /*  XP = VS_1 + up*(VS_2 - VS_1) + vp*(VS_3-VS_1)              */
  /* where (V_1, V_2, V_3) are the triangle vertices and (u,v)   */
  /* are the cubature points for a 2D numerical cubature rule    */
  /* over the standard triangle with vertices at (0,0)(1,0)(0,1).*/
  /* note that the jacobian of the transformation is 4*A*AP      */
  /* where A and AP are the areas of the triangles; this         */
  /* conveniently cancels the corresponding factor coming from   */
  /* the RWG basis function prefactor                            */
  /***************************************************************/
  Qa
  VecSub(VD[(iQD+1)%3],VD[iQD],AD);
  VecSub(VD[(iQD+2)%3],VD[iQD],BD);

  QS=VS[iQS];
  VecSub(VS[(iQS+1)%3],VS[iQS],AS);
  VecSub(VS[(iQS+2)%3],VS[iQS],BS);

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
  TCR=GetTCR(Order, &NumPts);

  /***************************************************************/
  /* preliminary setup before entering quadrature loops          */
  /***************************************************************/
  memset(L,0,3*sizeof(cdouble));
  if (GradL) memset(GradL,0,9*sizeof(cdouble));
  if (dLdT) memset(dLdT,0,3*NumTorqueAxes*sizeof(cdouble));

  LTerm[1]=4.0;   // this is constant throughout

  if (!NeedCross)
   LInner[2]=0.0 ;

  ik = II*Wavenumber;

  /***************************************************************/
  /* outer loop **************************************************/
  /***************************************************************/
  for(np=ncp=0; np<NumPts; np++) 
   { 
     u=TCR[ncp++]; v=TCR[ncp++]; w=TCR[ncp++];

     /***************************************************************/
     /* set XmQ and X ***********************************************/
     /***************************************************************/
     for(mu=0; mu<3; mu++)
      { XmQD[mu]=u*AD[mu] + v*BD[mu];
        X[mu]=XmQD[mu] + QD[mu];
      };

     /***************************************************************/
     /* inner loop to calculate value of inner integrand ************/
     /***************************************************************/
     memset(LInner,0,3*sizeof(cdouble));
     if (GradL) memset(GradLInner,0,9*sizeof(cdouble));
     if (dLdT) memset(dLdTInner,0,3*NumTorqueAxes*sizeof(cdouble));
     for(npp=ncpp=0; npp<NumPts; npp++)
      { 
        up=TCR[ncpp++]; vp=TCR[ncpp++]; wp=TCR[ncpp++];

        /***************************************************************/ 
        /* set XPmQp and XP ********************************************/
        /***************************************************************/
        for(mu=0; mu<3; mu++)
         { XPmQS[mu]=up*AS[mu] + vp*BS[mu];
           XP[mu]=XPmQS[mu] + QS[mu];
         };
      
        /***************************************************************/
        /* inner integrand  ********************************************/
        /***************************************************************/
        VecSub(X,XP,R);
        r=VecNorm(R);
        r2=r*r;

        /* compute L factors */
        LTerm[0]=VecDot(XmQD, XPmQS); 
        if (NeedCross)
         { VecCross(XmQD, XPmQS, gXh);
           LTerm[2]=VecDot( gXh, R );
         };
   
        /* compute Phi, Psi, Zeta factors */
        if (DeSingularize)
         Phi = cExpRel(ik*r,4) / (4.0*M_PI*r);
        else
         Phi = exp(ik*r) / (4.0*M_PI*r);

        if ( !finite(real(Phi)) ) Phi=0.0 ;
        Phi*=wp;
        Psi = Phi * (ik - 1.0/r) / r;
        Zeta = Phi * (ik*ik - 3.0*ik/r + 3.0/r2) / r2;

        /* now combine factors as necessary for the various integrands */

        /* 1. L_{0,1,2} */
        LInner[0] += LTerm[0] * Phi;
        LInner[1] += LTerm[1] * Phi;
        if ( NeedCross )
         LInner[2] += LTerm[2] * Psi;

        /* 2. d/dX_mu L_{0,1,2} */
        if ( GradL )
         for(mu=0; mu<3; mu++)
           { GradLInner[3*mu + 0] += LTerm[0]*R[mu]*Psi;
             GradLInner[3*mu + 1] += LTerm[1]*R[mu]*Psi;
             GradLInner[3*mu + 2] += gXh[mu]*Psi + LTerm[2]*R[mu]*Zeta;
           };

        /* 3. d/dTheta L_{0,1,2} */
        if ( NumTorqueAxes>0 && GammaMatrix!=0 && dLdT!=0 )
         for(nta=0; nta<NumTorqueAxes; nta++)
          { memset(dX,0,3*sizeof(double));
            memset(dG,0,3*sizeof(double));
            for(mu=0; mu<3; mu++)
             for(nu=0; nu<3; nu++)
              { dX[mu]+=GammaMatrix[9*nta + mu + 3*nu]*X[nu];
                dG[mu]+=GammaMatrix[9*nta + mu + 3*nu]*XmQD[nu];
              };
            Puv=VecDot(R,dX);
            dLdTInner[3*nta + 0] += LTerm[0]*Puv*Psi + VecDot(dG,XPmQS)*Phi;
            dLdTInner[3*nta + 1] += LTerm[1]*Puv*Psi;
            dLdTInner[3*nta + 2] += LTerm[2]*Puv*Zeta 
                                    + (  VecDot(VecCross(dG,XPmQS,dgXh),R) 
                                       + VecDot(gXh,dX)
                                      )*Psi;
          }; // for(nta= ... )

      }; /* for(npp=ncpp=0; npp<NumPts; npp++) */

     /*--------------------------------------------------------------*/
     /*- accumulate contributions to outer integral                  */
     /*--------------------------------------------------------------*/
     for(m=0; m<3; m++)
      L[m]+=w*LInner[m];
     if (GradL)
      for(m=0; m<9; m++)
       GradL[m]+=w*GradLInner[m];
     if (dLdT)
      for(m=0; m<3*NumTorqueAxes; m++)
       dLdT[m]+=w*dLdTInner[m];

   }; // for(np=ncp=0; np<nPts; np++) 
}

/***************************************************************/
/* calculate integrals over a single pair of triangles using   */
/* one of several different methods based on how near the two  */
/* triangles are to each other.                                */
/***************************************************************/
void GetPanelPanelInteractions(GPPIArgStruct *Args)
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
  cdouble *GC               = Args->GC;
  cdouble *GradGC           = Args->GradGC;
  cdouble *dGCdT            = Args->dGCdT;

  /***************************************************************/
  /* extract panel vertices, detect common vertices, measure     */
  /* relative distance                                           */
  /***************************************************************/
  double rRel; 
  double *Va[3], *Vb[3];
  RWGPanel *Pa = Oa->Panels[npa];
  RWGPanel *Pb = Ob->Panels[npb];
  double *Qa   = Oa->Panels[iQa];
  double *Qb   = Oa->Panels[iQb];
  int ncv=AssessPanelPair(Oa,npa,Ob,npb,&rRel,Va,Vb);
 
  /***************************************************************/
  /* if the panels are far apart then just use naive moderate-   */
  /* order cubature                                              */
  /***************************************************************/
  if ( rRel > DESINGULARIZATION_RADIUS )
   { GetPPI_Fixed(Args, 0, 0, Va, Vb, Qa, Qb);
     return;
   };

  /***************************************************************/
  /* if the panels are identical then we use taylor's scheme for */
  /* the full panel integral no matter what the frequency        */
  /***************************************************************/
  if( ncv==3 )
   { 
     Args->GC[0]=TaylorMaster(TM_COMMONTRIANGLE, TM_EIKR_OVER_R, TM_DOTPLUS, k,
                              Va[0], Va[1], Va[2], Vb[1], Vb[2], Qa, Qb);

     Args->GC[1]=0.0; /* 'C' integral vanishes for the common-triangle case */

     if (GradGC) memset(GradGC, 2*NumGradientComponents, 0*sizeof(cdouble));
     if (dGCdT)  memset(dGCdT, 2*NumTorqueAxes, 0*sizeof(cdouble));

     return;
   };

  /*****************************************************************/
  /* if we are in the high-frequency regime then desingularization */
  /* doesn't work; in this case, if there are any common vertices  */
  /* we use taylor's method for the full panel integral, and       */
  /* otherwise we use high-order naive cubature                    */
  /*****************************************************************/
  cdouble CPreFac=-1.0/(II*k);
  if ( abs(k*fmax(Pa->Radius, Pb->Radius)) > SWTHRESHOLD )
   { 
     if( ncv==2 )
      { 
        /*--------------------------------------------------------------*/
        /* common-edge case                                             */
        /*--------------------------------------------------------------*/
        Args->GC[0]=TaylorMaster(TM_COMMONEDGE, TM_EIKR_OVER_R, TM_DOTPLUS, k,
                                 Va[0], Va[1], Va[2], Vb[1], Vb[2], Qa, Qb);

        Args->GC[1]=CPreFac*TaylorMaster(TM_COMMONEDGE, TM_EIKR_OVER_R, TM_CROSS, k,
                                         Va[0], Va[1], Va[2], Vb[1], Vb[2], Qa, Qb);

        if (GradGC) memset(GradGC, 2*NumGradientComponents, 0*sizeof(cdouble));
        if (dGCdT)  memset(dGCdT, 2*NumTorqueAxes, 0*sizeof(cdouble));

        return;
      }
     else if( ncv==1 )
      { 
        /*--------------------------------------------------------------*/
        /* common-vertex case                                           */
        /*--------------------------------------------------------------*/
        Args->GC[0]=TaylorMaster(TM_COMMONVERTEX, TM_EIKR_OVER_R, TM_DOTPLUS, k,
                                 Va[0], Va[1], Va[2], Vb[1], Vb[2], Qa, Qb);

        Args->GC[1]=CPreFac*TaylorMaster(TM_COMMONVERTEX, TM_EIKR_OVER_R, TM_CROSS, k,
                                         Va[0], Va[1], Va[2], Vb[1], Vb[2], Qa, Qb);

        if (GradGC) memset(GradGC, 2*NumGradientComponents, 0*sizeof(cdouble));
        if (dGCdT)  memset(dGCdT, 2*NumTorqueAxes, 0*sizeof(cdouble));

        return;
      }
     else // ncv==0
      { 
        /*--------------------------------------------------------------*/
        /* no common vertices                                           */
        /*--------------------------------------------------------------*/
        GetPPI_Fixed(Args, 0, 1, Va, Vb, Qa, Qb);
        return;
      };

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
  GetPPI_Fixed(Args, 1, 0, Va, Vb, Qa, Qb);

  // step 2
  FIPPID=FIPPIDT->GetDataRecord(Va, iQa, Vb, iQb);

  // step 3
  // note: PF[n] = (ik)^n / (4\pi)
  PF[0]=1.0/(4.0*M_PI);
  PF[1]=ik*PF[0];
  PF[2]=ik*PF[1];
  PF[3]=ik*PF[2];
  PF[4]=ik*PF[3];
  PF[5]=ik*PF[4];

  // contributions to panel-panel integrals 
  Args->GC[0] += PF[0]*AA0*(FIPPID->hDotRm1 - OOK2*FIPPID->hNablaRm1)
                +PF[1]*AA1*(FIPPID->hDotR0  - OOK2*FIPPID->hNablaR0 )
                +PF[2]*AA2*(FIPPID->hDotR1  - OOK2*FIPPID->hNablaR1 )
                +PF[3]*AA3*(FIPPID->hDotR2  - OOK2*FIPPID->hNablaR2 );
  
  Args->GC[1] += PF[0]*BB0*FIPPID->hTimesRm3
                +PF[2]*BB2*FIPPID->hTimesRm1
                +PF[3]*BB3*FIPPID->hTimesR0 
                +PF[4]*BB4*FIPPID->hTimesR1;

  if (Args->NumGradientComponents==0 && Args->NumTorqueAxes==0) 
   return;

  // contributions to gradients of panel-panel integrals
  double Rab[3];
  VecSub(Pa->Centroid, Pb->Centroid, Rab);

  GradGScalar=  PF[0]*BB0*(FIPPID->hDotRm3 - OOK2*FIPPID->hNablaRm3)
               +PF[2]*BB2*(FIPPID->hDotRm1 - OOK2*FIPPID->hNablaRm1)
               +PF[3]*BB3*(FIPPID->hDotR0  - OOK2*FIPPID->hNablaR0)
               +PF[4]*BB4*(FIPPID->hDotR1  - OOK2*FIPPID->hNablaR1);

  GradCScalar=  PF[0]*CC0*FIPPID->hTimesRm5
               +PF[2]*CC2*FIPPID->hTimesRm3
               +PF[5]*CC5*FIPPID->hTimesR;

  for(Mu=0; Mu<(Args->NumGradientComponents); Mu++)
   { Args->GradGC[2*Mu + 0] += Rab[Mu]*GradGSCalar;
     Args->GradGC[2*Mu + 1] += Rab[Mu]*GradCSCalar; 
   };

  for(Mu=0; Mu<(Args->NumGradientComponents); Mu++)
   { Args->GradGC[2*Mu + 0] +=  PF[0]*BB0*FIPPID->dhTimesdRMuRm3[Mu]
                               +PF[2]*BB2*FIPPID->dhTimesdRMuRm1[Mu]
                               +PF[3]*BB3*FIPPID->dhTimesdRMuR0[Mu]
                               +PF[4]*BB3*FIPPID->dhTimesdRMuR1[Mu];
   };

  // contributions to angular derivatives of panel-panel integrals

} 

/***************************************************************/
/* this is an alternate entry point to GetPanelPanelInteractions*/
/* that copies the results out of the structure body into      */
/* user-specified buffers                                      */
/***************************************************************/
void GetPanelPanelInteractions(GEEIArgStruct *Args,
                               cdouble *GC, 
                               cdouble *GradGC, 
                               cdouble *dGCdT);
{ 
  GetPanelPanelInteractions(Args);

  memcpy(GC, Args->GC, 2*sizeof(cdouble));
  if(GradGC)  
   memcpy(GradGC, Args->GradGC, 2*Args->NumGradientComponents*sizeof(cdouble));
  if(dGCdT)  
   memcpy(dGCdT, Args->dGCdT, 2*Args->NumTorqueAxes*sizeof(cdouble));
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void InitGPPIArgs(GPPIArgStruct *Args)
{
  Args->NumGradientComponents=0;
  Args->NumTorqueAxes=0;
  Args->GammaMatrix=0;
}
