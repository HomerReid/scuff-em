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

#include "libscuff.h"
#include "libscuffInternals.h"

#define II cdouble(0,1)

// the 'distant basis function' threshold: two basis functions are
// 'distant' if their midpoint--midpoint distance is greater than 
// DBFTHRESHOLD * the larger of the radii of the two basis functions
#define DBFTHRESHOLD 10.0

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetEdgeEdgeInteractions(GEEIArgStruct *Args)
{ 
  /***************************************************************/
  /* local copies of fields in argument structure ****************/
  /***************************************************************/
  RWGObject *Oa             = Args->Oa;
  RWGObject *Ob             = Args->Ob;
  int nea                   = Args->nea; 
  int neb                   = Args->neb; 
  cdouble k                 = Args->k; int NumGradientComponents = Args->NumGradientComponents;
  int NumTorqueAxes         = Args->NumTorqueAxes;
  double *GammaMatrix       = Args->GammaMatrix;

  /***************************************************************/
  /* first look to see if the two basis functions are far enough */
  /* apart that we can use the multipole method                  */
  /***************************************************************/
  RWGEdge *Ea=Oa->Edges[nea];
  RWGEdge *Eb=Oa->Edges[neb];

  #if 0
  double RMax=fmax(Ea->Radius, Eb->Radius);
  if ( VecDistance(Ea->Centroid, Eb->Centroid) > DBFTHRESHOLD*RMax )
   { GetEEIMultipole(Args);
     return;
   };
  #endif

  /***************************************************************/
  /* otherwise, obtain the edge-edge interactions as a sum of    */
  /* four panel-panel interactions                               */
  /***************************************************************/
  cdouble HPP[2], HPM[2], HMP[2], HMM[2];
  cdouble GradHPP[6], GradHPM[6], GradHMP[6], GradHMM[6];
  cdouble dGCdTPP[6], dGCdTPM[6], dGCdTMP[6], dGCdTMM[6];

  /*--------------------------------------------------------------*/
  /*- initialize argument structure for GetPanelPanelInteractions */
  /*--------------------------------------------------------------*/
  GPPIArgStruct MyGPPIArgs, *GPPIArgs=&MyGPPIArgs;
  InitGPPIArgs(&GPPIArgs);

  GPPIArgs->Oa                     = Oa;
  GPPIArgs->Ob                     = Ob;
  GPPIArgs->k                      = k;
  GPPIArgs->NumGradientComponents  = NumGradientComponents;
  GPPIArgs->NumTorqueAxes          = NumTorqueAxes;
  GPPIArgs->GammaMatrix            = GammaMatrix;

  /*--------------------------------------------------------------*/
  /*- positive-positive, positive-negative, etc. -----------------*/
  /*--------------------------------------------------------------*/
  GPPIArgs->npa = Ea->iPPanel;     GPPIArgs->iQa = Ea->PIndex;
  GPPIArgs->npb = Eb->iPPanel;     GPPIArgs->iQb = Eb->PIndex;
  GetPanelPanelInteractions(GPPIArgs, HPP, GradHPP, dHdTPP);

  GPPIArgs->npa = Ea->iPPanel;     GPPIArgs->iQa = Ea->PIndex;
  GPPIArgs->npb = Eb->iMPanel;     GPPIArgs->iQb = Eb->MIndex;
  GetPanelPanelInteractions(GPPIArgs, HPM, GradHPM, dHdTPM);

  GPPIArgs->npa = Ea->iMPanel;     GPPIArgs->iQa = Ea->MIndex;
  GPPIArgs->npb = Eb->iPPanel;     GPPIArgs->iQb = Eb->PIndex;
  GetPanelPanelInteractions(GPPIArgs, HMP, GradHMP, dHdTMP);

  GPPIArgs->npa = Ea->iMPanel;     GPPIArgs->iQa = Ea->MIndex;
  GPPIArgs->npb = Eb->iMPanel;     GPPIArgs->iQb = Eb->MIndex;
  GetPanelPanelInteractions(GPPIArgs, HMM, GradHMM, dHdTMM);

  /*--------------------------------------------------------------*/
  /*- assemble the final quantities ------------------------------*/
  /*--------------------------------------------------------------*/
  double GPreFac = Ea->Length*Eb->Length;
  cdouble CPreFac = Ea->Length*Eb->Length / (II*k);

  Args->GC[0] = PreFac*(HPP[0] - HPM[0] - HMP[0] + HMM[0]);
  Args->GC[1] = PreFac*(HPP[1] - HPM[1] - HMP[1] + HMM[1]);

  for(Mu=0; Mu<NumGradientComponents; Mu++)
   { Args->GradGC[2*Mu+0] = PreFac*( GradHPP[2*Mu+0] - GradHPM[2*Mu+0] - GradHMP[2*Mu+0] + GradHMM[2*Mu+0] );
     Args->GradGC[2*Mu+1] = PreFac*( GradHPP[2*Mu+1] - GradHPM[2*Mu+1] - GradHMP[2*Mu+1] + GradHMM[2*Mu+1] );
   };

  for(Mu=0; Mu<NumTorqueAxes; Mu++)
   { Args->dGCdT[2*Mu+0] = PreFac*( dHdTPP[2*Mu+0] - dHdTPM[2*Mu+0] - dHdTMP[2*Mu+0] + dHdTMM[2*Mu+0]);
     Args->dGCdT[2*Mu+1] = PreFac*( dHdTPP[2*Mu+1] - dHdTPM[2*Mu+1] - dHdTMP[2*Mu+1] + dHdTMM[2*Mu+1]);
   };

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
