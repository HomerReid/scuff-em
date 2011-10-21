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
  /* the full panel integral                                    */
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
  if ( fabs(k*fmax(Pa->Radius, Pb->Radius)) > SWTHRESHOLD )
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
        GetPPI_Fixed(Args, 0, 0, Va, Vb, Qa, Qb);
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
  PPIData=GetStaticPPIData(SPPIDT, Va, Vb);

  // step 3

  double rCC, rRel;
  double *PV1[3], *PV2[3];
  int iTemp, XXPFlipped;
  int ncv, Index[3], IndexP[3];
  RWGPanel *P1, *P2;

  /***************************************************************/
  /* extract a little more information on the panel pair         */
  /***************************************************************/
  rCC=VecDistance(P1->Centroid, P2->Centroid);
  rRel=rCC/fmax(P1->Radius,P2->Radius);

  if (O1==O2)
   ncv=O1->CountCommonVertices(np1,np2,Index,IndexP);
  else
   ncv=0;

  /***************************************************************/
  /* if there are no common vertices and the imaginary part of   */
  /* kr is >= 25 (i.e. the panel-panel integrals will be down    */
  /* by a factor of e^{-25}) then we just call them zero.        */
  /***************************************************************/
  if ( ncv==0 && (imag(Wavenumber))*rCC>25.0 )
   { memset(L,0,3*sizeof(cdouble));
     if (GradL) memset(GradL,0,9*sizeof(cdouble));
     if (dLdT) memset(dLdT,0,3*NumTorqueAxes*sizeof(cdouble));
     return;
   };

  /***************************************************************/
  /* the next simplest possibility is that the panels are far    */
  /* enough away from each other that we don't need to           */
  /* desingularize the kernel in the integrand.                  */
  /***************************************************************/
  if ( rRel > 2*DESINGULARIZATION_RADIUS )
   PanelPanelInt_Fixed(Wavenumber, NeedCross, NumTorqueAxes, GammaMatrix, 
                       PV1, iQ1, PV2, iQ2, 4, 0, L, GradL, dLdT);
  else if ( rRel > DESINGULARIZATION_RADIUS )
   PanelPanelInt_Fixed(Wavenumber, NeedCross, NumTorqueAxes, GammaMatrix, 
                       PV1, iQ1, PV2, iQ2, 7, 0, L, GradL, dLdT);
  else if ( O1!=O2 && rRel>0.5*DESINGULARIZATION_RADIUS )
   PanelPanelInt_Fixed(Wavenumber, NeedCross, NumTorqueAxes, GammaMatrix, 
                       PV1, iQ1, PV2, iQ2, 14 , 0, L, GradL, dLdT);
  else if ( O1!=O2 ) /* this should happen rarely */
   PanelPanelInt_Fixed(Wavenumber, NeedCross, NumTorqueAxes, GammaMatrix, 
                       PV1, iQ1, PV2, iQ2, 20 , 0, L, GradL, dLdT);
  else
   { 
     /***************************************************************/
     /* panels are relatively near each other on the same object.   */
     /***************************************************************/

     /***************************************************************/
     /* 1. first handle the high-frequency case. for large Kappa,   */
     /*    the desingularized cubature method doesn't work well,    */
     /*    but the 'taylor method' works fine.                      */
     /*    we use the taylor method for panels with 1, 2, or 3      */
     /*    common vertices, and non-desingularized quadrature for   */
     /*    panels with no common vertices. (this latter may not be  */
     /*    very accurate, but in the high-frequency limit the       */
     /*    panel-panel integral over panels with no common vertices */
     /*    decays exponentially with kappa, whereas it decays       */
     /*    algebraically with kappa for panels with common vertices,*/
     /*    so the contributions from panels with no common vertices */
     /*    will be small and do not need to be evaluated to high    */
     /*    accuracy)                                                */
     /***************************************************************/
     if ( (imag(Wavenumber)*fmax(P1->Radius, P2->Radius)) > 2.0 )
      {
        if (ncv==2)
         {
           L[0]=TaylorMaster(TM_COMMONEDGE, TM_EIKR_OVER_R, TM_DOT, Wavenumber, 0.0,
                             PV1[Index[0]],  PV1[Index[1]],  PV1[3-Index[0]-Index[1]],
                             PV2[IndexP[1]], PV2[3-IndexP[0]-IndexP[1]],
                             PV1[iQ1], PV2[iQ2], 1.0);

           L[1]=4.0*TaylorMaster(TM_COMMONEDGE, TM_EIKR_OVER_R,  TM_ONE, Wavenumber, 0.0,
                                 PV1[Index[0]],  PV1[Index[1]],  PV1[3-Index[0]-Index[1]],
                                 PV2[IndexP[1]], PV2[3-IndexP[0]-IndexP[1]],
                                 PV1[iQ1], PV2[iQ2], 1.0);

           L[2]=TaylorMaster(TM_COMMONEDGE, TM_EIKR_OVER_R, TM_CROSS, Wavenumber, 0.0,
                             PV1[Index[0]],  PV1[Index[1]],  PV1[3-Index[0]-Index[1]],
                             PV2[IndexP[1]], PV2[3-IndexP[0]-IndexP[1]],
                             PV1[iQ1], PV2[iQ2], 1.0);
         }
        else if (ncv==1)
         { 
           L[0]=TaylorMaster(TM_COMMONVERTEX, TM_EIKR_OVER_R, TM_DOT, Wavenumber, 0.0,
                             PV1[Index[0]],   PV1[ (Index[0]+1)%3],  PV1[ (Index[0]+2)%3 ],
                             PV2[ (IndexP[0]+1)%3], PV2[ (IndexP[0]+2)%3 ],
                             PV1[iQ1], PV2[iQ2], 1.0);

           L[1]=4.0*TaylorMaster(TM_COMMONVERTEX, TM_EIKR_OVER_R, TM_ONE, Wavenumber, 0.0,
                                 PV1[Index[0]],   PV1[ (Index[0]+1)%3],  PV1[ (Index[0]+2)%3 ],
                                 PV2[ (IndexP[0]+1)%3], PV2[ (IndexP[0]+2)%3 ],
                                 PV1[iQ1], PV2[iQ2], 1.0);

           L[2]=TaylorMaster(TM_COMMONVERTEX, TM_EIKR_OVER_R, TM_CROSS, Wavenumber, 0.0,
                             PV1[Index[0]],   PV1[ (Index[0]+1)%3],  PV1[ (Index[0]+2)%3 ],
                             PV2[ (IndexP[0]+1)%3], PV2[ (IndexP[0]+2)%3 ],
                             PV1[iQ1], PV2[iQ2], 1.0);

         }
        else if (ncv==0)
         PanelPanelInt_Fixed(Wavenumber, NeedCross, 0, 0, 
                             PV1, iQ1, PV2, iQ2, 20 , 0, L, 0, 0); 
        return;
      };

     /***************************************************************/
     /* 2. otherwise, we use the desingularized quadrature method.  */
     /***************************************************************/
     StaticPPIData MySPPID, *SPPID;
     int rp;
     double A[3], AP[3], B[3], BP[3];
     double VmVP[3], VmQ[3], VPmQP[3], QmQP[3], QxQP[3], VScratch[3];
     double gDot[9], gCross[9];
     cdouble ik, ikPowers[RPMAX+2], ikPowers2[RPMAX+1];

     /***************************************************************/
     /* 2a. do some preliminary setup *******************************/
     /***************************************************************/
     double *pV1, *pV2, *pV3, *pV1P, *pV2P, *pV3P, *pQ, *pQP;
     double V1[3], V2[3], V3[3], V1P[3], V2P[3], V3P[3], Q[3], QP[3];
     if (ncv==1)
      { 
        pV1 = O1->Vertices + 3*O1->Panels[np1]->VI[ Index[0] ];
        pV2 = O1->Vertices + 3*O1->Panels[np1]->VI[ (Index[0]+1) %3 ];
        pV3 = O1->Vertices + 3*O1->Panels[np1]->VI[ (Index[0]+2) %3 ];
     
        //V1P = V1; 
        pV1P = O2->Vertices + 3*O2->Panels[np2]->VI[ IndexP[0] ];
        pV2P = O2->Vertices + 3*O2->Panels[np2]->VI[ (IndexP[0]+1) %3 ];
        pV3P = O2->Vertices + 3*O2->Panels[np2]->VI[ (IndexP[0]+2) %3 ];
      }
     else if (ncv==2)
      { 
        pV1 = O1->Vertices  + 3*O1->Panels[np1]->VI[ Index[0] ];
        pV2  = O1->Vertices + 3*O1->Panels[np1]->VI[ Index[1] ];
        pV3  = O1->Vertices + 3*O1->Panels[np1]->VI[ 3 - Index[0] - Index[1] ];
     
        //V1P = V1;
        //V2P = V2;
        pV1P = O2->Vertices + 3*O2->Panels[np2]->VI[ IndexP[0] ];
        pV2P = O2->Vertices + 3*O2->Panels[np2]->VI[ IndexP[1] ];
        pV3P = O2->Vertices + 3*O2->Panels[np2]->VI[ 3 - IndexP[0] - IndexP[1] ];
      }
     else
      { 
        pV1  = O1->Vertices + 3*O1->Panels[np1]->VI[0];
        pV2  = O1->Vertices + 3*O1->Panels[np1]->VI[1];
        pV3  = O1->Vertices + 3*O1->Panels[np1]->VI[2];
        pV1P = O2->Vertices + 3*O2->Panels[np2]->VI[0];
        pV2P = O2->Vertices + 3*O2->Panels[np2]->VI[1];
        pV3P = O2->Vertices + 3*O2->Panels[np2]->VI[2];
      };
     
     pQ=O1->Vertices + 3*O1->Panels[np1]->VI[iQ1];
     pQP=O2->Vertices + 3*O2->Panels[np2]->VI[iQ2];

     memcpy(V1,pV1,3*sizeof(double));
     memcpy(V2,pV2,3*sizeof(double));
     memcpy(V3,pV3,3*sizeof(double));
     memcpy(V1P,pV1P,3*sizeof(double));
     memcpy(V2P,pV2P,3*sizeof(double));
     memcpy(V3P,pV3P,3*sizeof(double));
     memcpy(Q,pQ,3*sizeof(double));
     memcpy(QP,pQP,3*sizeof(double));
     
     /*****************************************************************/
     /* 2b. first compute the panel-panel integral using medium-order */
     /*     cubature with the three most singular terms removed.      */
     /*****************************************************************/
     P1=O1->Panels[np1];
     P2=O2->Panels[np2];

     PV1[0]=O1->Vertices + 3*P1->VI[0];
     PV1[1]=O1->Vertices + 3*P1->VI[1];
     PV1[2]=O1->Vertices + 3*P1->VI[2];

     PV2[0]=O2->Vertices + 3*P2->VI[0];
     PV2[1]=O2->Vertices + 3*P2->VI[1];
     PV2[2]=O2->Vertices + 3*P2->VI[2];

     PanelPanelInt_Fixed(Wavenumber, NeedCross, 0, 0,
                         PV1, iQ1, PV2, iQ2, 4, 1, L, 0, 0);

     /****************************************************************/
     /* 2c. next get static panel-panel integral data, from a lookup */
     /*     table if we have one, or else by computing it on the fly */
     /****************************************************************/
     SPPID=0;
     if( O1->SPPIDTable )
      { SPPID=GetStaticPPIData(O1->SPPIDTable, np1, np2, &MySPPID);
          
        /* this step is important: we calculated the static panel-panel data */
        /* with the object mesh in its original configuration, i.e. as       */
        /* specified in the .msh file.  however, since then we may have      */
        /* transformed (rotated and/or translated) the object. thus, we need */
        /* to do the following computation with the vertices transformed     */
        /* back to their original locations.                                 */

        O1->UnTransformPoint(V1);
        O1->UnTransformPoint(V2);
        O1->UnTransformPoint(V3);
        O1->UnTransformPoint(Q);
        O2->UnTransformPoint(V1P);
        O2->UnTransformPoint(V2P);
        O2->UnTransformPoint(V3P);
        O2->UnTransformPoint(QP);
      }
     else
      { ComputeStaticPPIData(O1, np1, O2, np2, &MySPPID);
        SPPID=&MySPPID;
      };
   
     /***************************************************************/
     /* 2d. finally, augment the desingularized integrals by adding */
     /* the contributions from the integrals of the singular terms. */
     /***************************************************************/
     VecSub(V1, V1P, VmVP);

     VecSub(V1,Q,VmQ);
     VecSub(V2,V1,A);
     VecSub(V3,V2,B);

     VecSub(V1P,QP,VPmQP);
     VecSub(V2P,V1P,AP);
     VecSub(V3P,V2P,BP);

     gDot[0]=VecDot(VmQ,VPmQP);     /* 1   */
     gDot[1]=VecDot(VmQ,AP);        /* up  */
     gDot[2]=VecDot(VmQ,BP);        /* vp  */
     gDot[3]=VecDot(A,VPmQP);       /* u   */
     gDot[4]=VecDot(A,AP);          /* uup */
     gDot[5]=VecDot(A,BP);          /* uvp */
     gDot[6]=VecDot(B,VPmQP);       /* v   */
     gDot[7]=VecDot(B,AP);          /* vup */
     gDot[8]=VecDot(B,BP);          /* vvp */

     if (NeedCross)
      { 
        VecSub(Q,QP,QmQP);
        VecCross(Q,QP,QxQP);

        gCross[0]=  VecDot(VmVP, QxQP) + VecDot(VecCross(V1, V1P, VScratch), QmQP); /* 1   */ 
        gCross[1]= -VecDot(  AP, QxQP) + VecDot(VecCross(V1, AP , VScratch), QmQP); /* up  */
        gCross[2]= -VecDot(  BP, QxQP) + VecDot(VecCross(V1, BP , VScratch), QmQP); /* vp  */
        gCross[3]=  VecDot(   A, QxQP) + VecDot(VecCross(A,  V1P, VScratch), QmQP); /* u   */ 
        gCross[4]=                       VecDot(VecCross(A,  AP,  VScratch), QmQP); /* uup */
        gCross[5]=                       VecDot(VecCross(A,  BP,  VScratch), QmQP); /* uvp */
        gCross[6]=  VecDot(   B, QxQP) + VecDot(VecCross(B,  V1P, VScratch), QmQP); /* v   */ 
        gCross[7]=                       VecDot(VecCross(B,  AP,  VScratch), QmQP); /* vup */
        gCross[8]=                       VecDot(VecCross(B,  BP,  VScratch), QmQP); /* vvp */

        L[2] -= (  VecDot(SPPID->XmXPoR3Int, QxQP) 
                 + VecDot(SPPID->XxXPoR3Int, QmQP) ) / (4.0*M_PI);
      };

     ik = II*Wavenumber;

     ikPowers[0]=1.0;
     ikPowers[1]=ik;
     ikPowers[2]=ik*ikPowers[1] / 2.0;
     ikPowers[3]=ik*ikPowers[2] / 3.0;

     ikPowers2[0]=ikPowers[2];
     ikPowers2[1]=2.0*ikPowers[3];
     ikPowers2[2]=ik*ikPowers[3];

     for(rp=-1; rp<=RPMAX; rp++)
      { 
         L[0] += ikPowers[rp+1]*(    SPPID->hrnInt[rp+1][0]*gDot[0]
                                   + SPPID->hrnInt[rp+1][1]*gDot[1]
                                   + SPPID->hrnInt[rp+1][2]*gDot[2]
                                   + SPPID->hrnInt[rp+1][3]*gDot[3]
                                   + SPPID->hrnInt[rp+1][4]*gDot[4]
                                   + SPPID->hrnInt[rp+1][5]*gDot[5]
                                   + SPPID->hrnInt[rp+1][6]*gDot[6]
                                   + SPPID->hrnInt[rp+1][7]*gDot[7]
                                   + SPPID->hrnInt[rp+1][8]*gDot[8] 
                                ) / (4.0*M_PI);

         L[1] += 4.0*ikPowers[rp+1]*SPPID->hrnInt[rp+1][0] / (4.0*M_PI);

         if(NeedCross && rp<RPMAX)
          L[2] += ikPowers2[rp+1]*(    SPPID->hrnInt[rp+1][0]*gCross[0]
                                     + SPPID->hrnInt[rp+1][1]*gCross[1]
                                     + SPPID->hrnInt[rp+1][2]*gCross[2]
                                     + SPPID->hrnInt[rp+1][3]*gCross[3]
                                     + SPPID->hrnInt[rp+1][4]*gCross[4]
                                     + SPPID->hrnInt[rp+1][5]*gCross[5]
                                     + SPPID->hrnInt[rp+1][6]*gCross[6]
                                     + SPPID->hrnInt[rp+1][7]*gCross[7]
                                     + SPPID->hrnInt[rp+1][8]*gCross[8]
                                  ) / (4.0*M_PI);

      };

   };

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
