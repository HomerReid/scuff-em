/*
 * EdgeEdgeInteractions.cc -- routines for evaluating interactions
 *                         -- between RWG edge basis functions
 * 
 * homer reid    -- 11/2005 -- 10/2011
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

/***************************************************************/
/* Evaluate and return the L-functions between two RWG basis   */
/* functions.                                                  */
/*                                                             */
/* Inputs:                                                     */
/*                                                             */
/*     O1,ne1:      object and BF index of first BF            */
/*     O2,ne2:      object and BF index of second BF           */
/*     Wavenumber:  complex wavevector                         */
/*     NeedCross:   set to 1 if you need the cross-product     */
/*                  integrals. if you set this to 0, the       */
/*                  corresponding slot in the return values    */
/*                  is set to 0.                               */
/*                                                             */
/*    NumTorqueAxes: number (0--3) of axes about which to      */
/*                  compute theta derivatives.                 */
/*     GammaMatrix: matrices describing the rotation axes      */
/*                  about which d/dTheta integrals are         */
/*                  computed (for torque computations).        */
/*                  GammaMatrix[ 9*nta + (i+3*j)] = i,j entry  */
/*                  in the Gamma matrix for torque axis #nta.  */
/*                  (Not referenced if NumTorqueAxes==0).      */
/*                                                             */
/* Outputs:                                                    */
/*                                                             */
/*     L[m]:        The inner product of basis functions       */
/*                  (O1,ne1) and (O2,ne2) with operator L_m    */
/*                  (m=0, 1, 2 for bullet, nabla, times).      */
/*                                                             */
/*     GradL[3*mu + m]: d/dx_mu L[m]                           */
/*                  where the derivative is wrt displacement   */
/*                  of the first object (O1, ne1) in the x_mu  */
/*                  direction.                                 */
/*                  If GradL is NULL on input, derivatives     */
/*                  are not computed.                          */
/*                                                             */
/*     dLdT[3*nta + m]: d/dTheta L[m]                          */
/*                  where the derivative is wrt rotation of    */
/*                  the first object (O1, ne1) about an axis   */
/*                  described by GammaMatrix[nta] as           */
/*                  discussed above.                           */
/*                  If dLdT or GammaMatrix are NULL on input   */
/*                  (and/or if NumTorqueAxes==0) then          */
/*                  Theta derivatives are not computed.        */
/***************************************************************/
void GetEEIs(RWGObject *O1, int ne1, RWGObject *O2, int ne2,
             cdouble Wavenumber,
             int NumTorqueAxes, double *GammaMatrix,
             cdouble L[3], cdouble *GradL, cdouble *dLdT)
{ 
  RWGEdge *E1, *E2;
  int m, mu;
  int nta;
  int iPPanel1, iMPanel1, iPPanel2, iMPanel2;  
  int iQP1, iQM1, iQP2, iQM2;
  double PreFac;
  cdouble LPP[3], LPM[3], LMP[3], LMM[3];
  cdouble GradLPP[9], GradLPM[9], GradLMP[9], GradLMM[9];
  cdouble dLPPdT[9], dLPMdT[9], dLMPdT[9], dLMMdT[9];

  /*--------------------------------------------------------------*/
  /*- preliminary setup ------------------------------------------*/
  /*--------------------------------------------------------------*/
  E1=O1->Edges[ne1];
  E2=O2->Edges[ne2];

  iPPanel1=E1->iPPanel;
  iMPanel1=E1->iMPanel;
  iPPanel2=E2->iPPanel;
  iMPanel2=E2->iMPanel;
  iQP1=E1->PIndex;
  iQM1=E1->MIndex;
  iQP2=E2->PIndex;
  iQM2=E2->MIndex;

  PreFac=E1->Length * E2->Length;

  /*--------------------------------------------------------------*/
  /*- compute integrals over each of the four pairs of panels    -*/
  /*--------------------------------------------------------------*/
  PanelPanelInt(Wavenumber, NeedCross, NumTorqueAxes, GammaMatrix,
                O1, iPPanel1, iQP1, O2, iPPanel2, iQP2, 
                LPP, GradLPP, dLPPdT);

  PanelPanelInt(Wavenumber, NeedCross, NumTorqueAxes, GammaMatrix,
                O1, iPPanel1, iQP1, O2, iMPanel2, iQM2, 
                LPM, GradLPM, dLPMdT);

  PanelPanelInt(W, Wavenumber, NeedCross, NumTorqueAxes, GammaMatrix,
                O1, iMPanel1, iQM1, O2, iPPanel2, iQP2, 
                LMP, GradLMP, dLMPdT);

  PanelPanelInt(W, Wavenumber, NeedCross, NumTorqueAxes, GammaMatrix,
                O1, iMPanel1, iQM1, O2, iMPanel2, iQM2, 
                LMM, GradLMM, dLMMdT);

  /*--------------------------------------------------------------*/
  /*- assemble the final quantities ------------------------------*/
  /*--------------------------------------------------------------*/
  for(m=0; m<3; m++)
   L[m]=PreFac*(LPP[m] - LPM[m] - LMP[m] + LMM[m]);

  if (GradL)
   for(m=0; m<9; m++)
    GradL[m]=PreFac*(GradLPP[m] - GradLPM[m] - GradLMP[m] + GradLMM[m]);

  if (dLdT)
   for(m=0; m<3*NumTorqueAxes; m++)
    dLdT[m]=PreFac*(dLPPdT[m] - dLPMdT[m] - dLMPdT[m] + dLMMdT[m]); 

}

/***************************************************************/
/* utility routines used to construct a GammaMatrix that may   */
/* be passed to MatrixElement to compute theta derivatives.    */
/*                                                             */
/* there are three different entry points to this routine      */
/* depending on how the user prefers to specify the axis about */
/* which objects are rotated.                                  */
/*                                                             */
/* in all cases, GammaMatrix must point to a buffer with       */
/* space for 9 doubles, which is filled in with the matrix and */
/* may be subsequently passed to MatrixElement.                */
/*                                                             */
/* Algorithm:                                                  */
/*  1. Construct the matrix Lambda that rotates the Z axis     */
/*     into alignment with the TorqueAxis.                     */
/*  2. Construct the matrix Gamma_0 such that Gamma_0*dTheta is*/
/*     the matrix that rotates through an infinitesimal angle  */
/*     dTheta about the Z axis.                                */
/*  3. Set GammaMatrix = Lambda^{-1} * Gamma_0 * Lambda.       */
/***************************************************************/

/* 1: specify torque axis as a vector of cartesian components */
void CreateGammaMatrix(double *TorqueAxis, double *GammaMatrix)
{ 
  int i, j, k, l;
  double Lambda[3][3], Gamma0[3][3], ct, st, cp, sp;
  double MyTorqueAxis[3];

  /* make a quick copy so we don't modify the user's vector */
  memcpy(MyTorqueAxis,TorqueAxis,3*sizeof(double));

  VecNormalize(MyTorqueAxis);
  ct=MyTorqueAxis[2];
  st=sqrt(1.0-ct*ct);
  cp= ( st < 1.0e-8 ) ? 1.0 : MyTorqueAxis[0] / st;
  sp= ( st < 1.0e-8 ) ? 0.0 : MyTorqueAxis[1] / st;

  Lambda[0][0]=ct*cp;  Lambda[0][1]=ct*sp;  Lambda[0][2]=-st;
  Lambda[1][0]=-sp;    Lambda[1][1]=cp;     Lambda[1][2]=0.0;
  Lambda[2][0]=st*cp;  Lambda[2][1]=st*sp;  Lambda[2][2]=ct;

  Gamma0[0][0]=0.0;    Gamma0[0][1]=-1.0;   Gamma0[0][2]=0.0;
  Gamma0[1][0]=1.0;    Gamma0[1][1]=0.0;    Gamma0[1][2]=0.0;
  Gamma0[2][0]=0.0;    Gamma0[2][1]=0.0;    Gamma0[2][2]=0.0;

  memset(GammaMatrix,0,9*sizeof(double));
  for(i=0; i<3; i++)
   for(j=0; j<3; j++)
    for(k=0; k<3; k++)
     for(l=0; l<3; l++)
      // Gamma[i+3*j] += LambdaInverse[i][k] * Gamma0[k][l] * Lambda[l][j];
      GammaMatrix[i+3*j] += Lambda[k][i] * Gamma0[k][l] * Lambda[l][j];
}


/* 2: specify torque axis as three separate cartesian components */
void CreateGammaMatrix(double TorqueAxisX, double TorqueAxisY, 
                       double TorqueAxisZ, double *GammaMatrix)
{ double TorqueAxis[3];

  TorqueAxis[0]=TorqueAxisX;
  TorqueAxis[1]=TorqueAxisY;
  TorqueAxis[2]=TorqueAxisZ;
  CreateGammaMatrix(TorqueAxis,GammaMatrix);
}

/* 3: specify (Theta,Phi) angles of torque axis */ 
void CreateGammaMatrix(double Theta, double Phi, double *GammaMatrix)
{ 
  double TorqueAxis[3];

  TorqueAxis[0]=sin(Theta)*cos(Phi);
  TorqueAxis[1]=sin(Theta)*sin(Phi);
  TorqueAxis[2]=cos(Theta);
  CreateGammaMatrix(TorqueAxis,GammaMatrix);
} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
void InitGEEIArgs(GEEIArgStruct *Args)
{
  Args->NumGradientComponents=0;
  Args->NumTorqueAxes=0;
  Args->GammaMatrix=0;
  
}
