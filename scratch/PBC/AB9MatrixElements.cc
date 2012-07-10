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
 * GetAB9MatrixElements()
 *
 * homer reid  -- 4/2011 -- 7/2012
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>
#include <libMDInterp.h>
#include <libTriInt.h>
#include <libscuff.h>
#include <PBCGeometry.h>

namespace scuff{

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetAB9PanelPanelInteraction(double **Va, double *Qa,
                                 double **Vb, double *Qb,
                                 cdouble k,
                                 Interp3D *Interpolator,
                                 cdouble GC[2])
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

  /***************************************************************/
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
  TCR=GetTCR(7, &NumPts);
 // if (HighOrder)
 //  TCR=GetTCR(20, &NumPts);
 // else
 //  TCR=GetTCR(4, &NumPts);

  /***************************************************************/
  /* outer loop **************************************************/
  /***************************************************************/
  memset(GC,0,2*sizeof(cdouble));

  double hDot, hNabla=4.0; // note hNabla is constant throughout
  cdouble hPlus, ik = II*k, ik2=ik*ik;
  int np, ncp, npp, ncpp, i;
  double u, v, w, up, vp, wp;
  double X[3], F[3], XP[3], FP[3], R[3], FxFP[3];
  int ZFlipped;
  double PhiVD[16];
  cdouble GBar, GradGBar[3];
  cdouble GCInner[2];
  for(np=ncp=0; np<NumPts; np++) 
   { 
     u=TCR[ncp++]; v=TCR[ncp++]; w=TCR[ncp++];

     /***************************************************************/
     /* set X and F=X-Q *********************************************/
     /***************************************************************/
     for(i=0; i<3; i++)
      { X[i] = V0[i] + u*A[i] + v*B[i];
        F[i] = X[i] - Q[i];
      };

     /***************************************************************/
     /* inner loop to calculate value of inner integrand ************/
     /***************************************************************/
     memset(GCInner,0,2*sizeof(cdouble));
     for(npp=ncpp=0; npp<NumPts; npp++)
      { 
        up=TCR[ncpp++]; vp=TCR[ncpp++]; wp=TCR[ncpp++];

        /***************************************************************/ 
        /* set XP and FP=XP-QP *****************************************/
        /***************************************************************/
        for(i=0; i<3; i++)
         { XP[i] = V0P[i] + up*AP[i] + vp*BP[i];
           FP[i] = XP[i] - QP[i];
           R[i] = X[i] - XP[i];
         };

        /***************************************************************/
        /***************************************************************/
        /***************************************************************/
        if ( R[2] < 0.0 )
         { R[2] *= -1.0;
           ZFlipped=1;
         }
        else
         ZFlipped=0;

        /***************************************************************/
        /***************************************************************/
	/***************************************************************/    
        Interpolator->EvaluatePlus(R[0], R[1], R[2], PhiVD);
        GBar = cdouble(PhiVD[0],PhiVD[8+0]);
        GradGBar[0] = cdouble(PhiVD[1],PhiVD[8+0]);
        GradGBar[1] = cdouble(PhiVD[2],PhiVD[8+1]);
        GradGBar[2] = cdouble(PhiVD[3],PhiVD[8+2]);

        if (ZFlipped)
         GradGBar[2]*=-1.0; // flip the sign of dG/dz 
      
        /***************************************************************/
        /* inner integrand  ********************************************/
        /***************************************************************/
        /* compute h factors */
        hDot=VecDot(F, FP);
        hPlus=hDot + hNabla/ik2;
        VecCross(F, FP, FxFP);
   
        GCInner[0] += wp * hPlus * GBar;   
        GCInner[1] += wp * (  FxFP[0]*GradGBar[0] 
                             +FxFP[1]*GradGBar[1] 
                             +FxFP[2]*GradGBar[2] );

      }; /* for(npp=ncpp=0; npp<NumPts; npp++) */

     /*--------------------------------------------------------------*/
     /*- accumulate contributions to outer integral                  */
     /*--------------------------------------------------------------*/
     GC[0]+=w*GCInner[0];
     GC[1]+=w*GCInner[1];

   }; // for(np=ncp=0; np<nPts; np++) 


}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetAB9EdgeEdgeInteractions(RWGObject *Oa, int nea, RWGObject *Ob, int neb, 
                                cdouble k, Interp3D *Interpolator, cdouble *GC)
{

  RWGEdge *Ea = Oa->Edges[nea];
  RWGEdge *Eb = Ob->Edges[neb];

  double *Va[3], *Qa, *Vb[3], *Qb;

  Va[1] = Oa->Vertices + 3*(Ea->iV1);
  Va[2] = Oa->Vertices + 3*(Ea->iV2);

  Vb[1] = Oa->Vertices + 3*(Eb->iV1);
  Vb[2] = Ob->Vertices + 3*(Eb->iV2);

  /*- PP ---------------------------------------------------------*/ 
  Va[0] = Qa = Oa->Vertices + 3*(Ea->iQP);
  Vb[0] = Qb = Ob->Vertices + 3*(Eb->iQP);
  GetAB9PanelPanelInteraction(Va, Qa, Vb, Qb, k, Interpolator, GCPP);

  /*- PM ---------------------------------------------------------*/ 
  Va[0] = Qa = Oa->Vertices + 3*(Ea->iQP);
  Vb[0] = Qb = Ob->Vertices + 3*(Eb->iQM);
  GetAB9PanelPanelInteraction(Va, Qa, Vb, Qb, k, Interpolator, GCPM);

  /*- MP ---------------------------------------------------------*/ 
  Va[0] = Qa = Oa->Vertices + 3*(Ea->iQM);
  Vb[0] = Qb = Ob->Vertices + 3*(Eb->iQP);
  GetAB9PanelPanelInteraction(Va, Qa, Vb, Qb, k, Interpolator, GCMP);

  /*- MM ---------------------------------------------------------*/ 
  Va[0] = Qa = Oa->Vertices + 3*(Ea->iQM);
  Vb[0] = Qb = Ob->Vertices + 3*(Eb->iQM);
  GetAB9PanelPanelInteraction(Va, Qa, Vb, Qb, k, Interpolator, GCMM);

  double PreFac = Ea->Length * Eb->Length;
  GC[0] = PreFac * (GCPP[0] - GCPM[0] - GCMP[0] + GCMM[0]);
  GC[1] = PreFac * (GCPP[1] - GCPM[1] - GCMP[1] + GCMM[1]);
  
}

} //namespace scuff
