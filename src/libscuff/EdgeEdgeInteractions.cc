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
  cdouble k                 = Args->k;
  int NumGradientComponents = Args->NumGradientComponents;
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
  if ( VecDistance(Ea->Centroid, Eb->Centroid) > EEITHRESHOLD*RMax )
   { GetEEIMultipole(Args);
     return;
   };
  #endif

  /***************************************************************/
  /* otherwise, obtain the edge-edge interactions as a sum of    */
  /* four panel-panel interactions                               */
  /***************************************************************/
  cdouble GPP, GPM, GMP, GMM;
  cdouble CPP, CPM, CMP, CMM;
  cdouble GradGPP[3], GradGPM[3], GradGMP[3], GradGMM[3];
  cdouble GradCPP[3], GradCPM[3], GradCMP[3], GradCMM[3];
  cdouble dGdThetaPP[3], dGdThetaPM[3], dGdThetaMP[3], dGdThetaMM[3];
  cdouble dCdThetaPP[3], dCdThetaPM[3], dCdThetaMP[3], dCdThetaMM[3];

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
  GetPanelPanelInteractions(GPPIArgs, GPP, CPP, GradGPP, GradCPP, dGdThetaPP, dCdThetaPP);

  GPPIArgs->npa = Ea->iPPanel;     GPPIArgs->iQa = Ea->PIndex;
  GPPIArgs->npb = Eb->iMPanel;     GPPIArgs->iQb = Eb->MIndex;
  GetPanelPanelInteractions(GPPIArgs, GPM, CPM, GradGPM, GradCPM, dGdThetaPM, dCdThetaPM);

  GPPIArgs->npa = Ea->iMPanel;     GPPIArgs->iQa = Ea->MIndex;
  GPPIArgs->npb = Eb->iPPanel;     GPPIArgs->iQb = Eb->PIndex;
  GetPanelPanelInteractions(GPPIArgs, GMP, CMP, GradGMP, GradCMP, dGdThetaMP, dCdThetaMP);

  GPPIArgs->npa = Ea->iMPanel;     GPPIArgs->iQa = Ea->MIndex;
  GPPIArgs->npb = Eb->iMPanel;     GPPIArgs->iQb = Eb->MIndex;
  GetPanelPanelInteractions(GPPIArgs, GMM, CMM, GradGMM, GradCMM, dGdThetaMM, dCdThetaMM);

  /*--------------------------------------------------------------*/
  /*- assemble the final quantities ------------------------------*/
  /*--------------------------------------------------------------*/
  double PreFac = Ea->Length*Eb->Length;

  Args->GInt = PreFac*(GPP[m] - GPM[m] - GMP[m] + GMM[m]);
  Args->CInt = PreFac*(CPP[m] - CPM[m] - CMP[m] + CMM[m]);

  for(Mu=0; Mu<NumGradientComponents; Mu++)
   { Args->GradGInt[Mu] = PreFac*( GradGPP[Mu] - GradGPM[Mu] - GradGMP[Mu] + GradGMM[Mu]);
     Args->GradCInt[Mu] = PreFac*( GradCPP[Mu] - GradCPM[Mu] - GradCMP[Mu] + GradCMM[Mu]);
   };

  for(Mu=0; Mu<NumTorqueAxes; Mu++)
   { Args->dGIntdTheta[Mu] = PreFac*( dGdThetaPP[Mu] - dGdThetaPM[Mu] - dGdThetaMP[Mu] + dGdThetaMM[Mu]);
     Args->dCIntdTheta[Mu] = PreFac*( dCdThetaPP[Mu] - dCdThetaPM[Mu] - dCdThetaMP[Mu] + dCdThetaMM[Mu]);
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
