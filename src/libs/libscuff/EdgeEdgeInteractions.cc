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

namespace scuff {

#define II cdouble(0,1)

// the 'distant basis function' threshold: two basis functions are
// 'distant' if their midpoint--midpoint distance is greater than 
// DBFTHRESHOLD * the larger of the radii of the two basis functions
#define DBFTHRESHOLD 10.0

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetEdgeEdgeInteractions(GetEEIArgStruct *Args)
{ 
  /***************************************************************/
  /* local copies of fields in argument structure ****************/
  /***************************************************************/
  RWGSurface *Sa            = Args->Sa;
  RWGSurface *Sb            = Args->Sb;
  int nea                   = Args->nea; 
  int neb                   = Args->neb;
  cdouble k                 = Args->k; 
  int NumGradientComponents = Args->NumGradientComponents;
  int NumTorqueAxes         = Args->NumTorqueAxes;

  RWGEdge *Ea=Sa->Edges[nea];
  RWGEdge *Eb=Sb->Edges[neb];

  /***************************************************************/
  /* Since this code doesn't work at DC anyway, we don't bother  */
  /* to compute the edge--edge interactions at k==0, but instead */
  /* just set them to zero. This is actually useful as it gives  */
  /* an easy way to zero out the contributions of individual     */
  /* regions and/or the external medium to the BEM matrix: just  */
  /* set epsilon and/or mu for that region temporarily to 0.     */
  /***************************************************************/
  if ( real(k)==0.0 && imag(k)==0.0 )
   { memset(Args->GC, 0, 2*sizeof(cdouble));
     memset(Args->GradGC, 0, 6*sizeof(cdouble));
     memset(Args->dGCdT, 0, 6*sizeof(cdouble));
     return;
   };

  /***************************************************************/
  /* figure out which method to use, as follows:                 */
  /*  a) if the Args->Force flag was set to force a specific     */
  /*     method, use that method                                 */
  /*  b) otherwise, use spherical multipoles if the relative     */
  /*     distance between the edges exceeds the DBFTHRESHOLD,    */
  /*     and otherwise use the panel-panel integral method       */
  /* // FIXME to take into account the possibility of nonzero    */
  /* displacement when figuring the distance between edge centroids */
  /***************************************************************/
#if 0
  bool UseSMMethod=false;
  if (Args->Force==EEI_FORCE_SM)
   UseSMMethod=true;
  else if (Args->Force==EEI_FORCE_PP)
   UseSMMethod=false;
  else
   { 
     double RMax=fmax(Ea->Radius, Eb->Radius);
     if (VecDistance(Ea->Centroid, Eb->Centroid) > DBFTHRESHOLD*RMax )
      UseSMMethod=true;
   };
#endif

  /***************************************************************/
  /* get edge-edge interactions by spherical multipole method if */
  /* that was the verdict of the above                           */
  /***************************************************************/
  #if 0
  if ( UseSMMethod )
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
  cdouble dHdTPP[6], dHdTPM[6], dHdTMP[6], dHdTMM[6];

  memset(HPM,     0, 2*sizeof(cdouble));
  memset(HMP,     0, 2*sizeof(cdouble));
  memset(HMM,     0, 2*sizeof(cdouble));
  memset(GradHPM, 0, 6*sizeof(cdouble));
  memset(GradHMP, 0, 6*sizeof(cdouble));
  memset(GradHMM, 0, 6*sizeof(cdouble));
  memset(dHdTPM,  0, 6*sizeof(cdouble));
  memset(dHdTMP,  0, 6*sizeof(cdouble));
  memset(dHdTMM,  0, 6*sizeof(cdouble));

  /*--------------------------------------------------------------*/
  /*- initialize argument structure for GetPanelPanelInteractions */
  /*--------------------------------------------------------------*/
  GetPPIArgStruct MyGetPPIArgs, *GetPPIArgs=&MyGetPPIArgs;
  InitGetPPIArgs(GetPPIArgs);

  GetPPIArgs->Sa                     = Sa;
  GetPPIArgs->Sb                     = Sb;
  GetPPIArgs->k                      = k;
  GetPPIArgs->NumGradientComponents  = NumGradientComponents;
  GetPPIArgs->NumTorqueAxes          = NumTorqueAxes;
  GetPPIArgs->GammaMatrix            = Args->GammaMatrix;
  GetPPIArgs->opFC                   = Args->opFC;
  GetPPIArgs->Displacement           = Args->Displacement;
  GetPPIArgs->GBA                    = Args->GBA;    
  GetPPIArgs->ForceFullEwald         = Args->ForceFullEwald;

  /*--------------------------------------------------------------*/
  /*- positive-positive, positive-negative, etc. -----------------*/
  /*--------------------------------------------------------------*/
  GetPPIArgs->npa = Ea->iPPanel;     GetPPIArgs->iQa = Ea->PIndex;
  GetPPIArgs->npb = Eb->iPPanel;     GetPPIArgs->iQb = Eb->PIndex;
  GetPanelPanelInteractions(GetPPIArgs, HPP, GradHPP, dHdTPP);
  Args->PPIAlgorithmCount[GetPPIArgs->WhichAlgorithm]++;

  if ( Eb->iMPanel!=-1 )
   { GetPPIArgs->npa = Ea->iPPanel;     GetPPIArgs->iQa = Ea->PIndex;
     GetPPIArgs->npb = Eb->iMPanel;     GetPPIArgs->iQb = Eb->MIndex;
     GetPanelPanelInteractions(GetPPIArgs, HPM, GradHPM, dHdTPM);
     Args->PPIAlgorithmCount[GetPPIArgs->WhichAlgorithm]++;
   };

  if ( Ea->iMPanel!=-1 )
   { GetPPIArgs->npa = Ea->iMPanel;     GetPPIArgs->iQa = Ea->MIndex;
     GetPPIArgs->npb = Eb->iPPanel;     GetPPIArgs->iQb = Eb->PIndex;
     GetPanelPanelInteractions(GetPPIArgs, HMP, GradHMP, dHdTMP);
     Args->PPIAlgorithmCount[GetPPIArgs->WhichAlgorithm]++;
   };
 
  if ( Ea->iMPanel!=-1 && Eb->iMPanel!=-1 )
   { GetPPIArgs->npa = Ea->iMPanel;     GetPPIArgs->iQa = Ea->MIndex;
     GetPPIArgs->npb = Eb->iMPanel;     GetPPIArgs->iQb = Eb->MIndex;
     GetPanelPanelInteractions(GetPPIArgs, HMM, GradHMM, dHdTMM);
     Args->PPIAlgorithmCount[GetPPIArgs->WhichAlgorithm]++;
   };

  /*--------------------------------------------------------------*/
  /*- assemble the final quantities ------------------------------*/
  /*--------------------------------------------------------------*/
  double GPreFac = Ea->Length*Eb->Length;
  cdouble CPreFac = Ea->Length*Eb->Length / (II*k);
  int Mu;

  Args->GC[0] = GPreFac*(HPP[0] - HPM[0] - HMP[0] + HMM[0]);
  Args->GC[1] = CPreFac*(HPP[1] - HPM[1] - HMP[1] + HMM[1]);

  for(Mu=0; Mu<NumGradientComponents; Mu++)
   { Args->GradGC[2*Mu+0] = GPreFac*( GradHPP[2*Mu+0] - GradHPM[2*Mu+0] - GradHMP[2*Mu+0] + GradHMM[2*Mu+0] );
     Args->GradGC[2*Mu+1] = CPreFac*( GradHPP[2*Mu+1] - GradHPM[2*Mu+1] - GradHMP[2*Mu+1] + GradHMM[2*Mu+1] );
   };

  for(Mu=0; Mu<NumTorqueAxes; Mu++)
   { Args->dGCdT[2*Mu+0] = GPreFac*( dHdTPP[2*Mu+0] - dHdTPM[2*Mu+0] - dHdTMP[2*Mu+0] + dHdTMM[2*Mu+0]);
     Args->dGCdT[2*Mu+1] = CPreFac*( dHdTPP[2*Mu+1] - dHdTPM[2*Mu+1] - dHdTMP[2*Mu+1] + dHdTMM[2*Mu+1]);
   };

  /*--------------------------------------------------------------*/
  /*- 20150224 detect loss of precision in PBC calculations and   */
  /*--------------------------------------------------------------*/
  if (Args->GBA && Args->ForceFullEwald==false)
   { 
     double AbsGC0 = 0.25*GPreFac*( abs(HPP[0]) + abs(HPM[0]) + abs(HMP[0]) + abs(HMM[0]) );
     if ( (AbsGC0 > 1.0e-12) && (abs(Args->GC[0]) < 0.01*AbsGC0 ) )
      { Log("Loss of precision in GetEEIs({%i,%i},{%i,%i}={%e,%e}) (recomputing)",
             Sa->Index, nea, Sb->Index, neb,AbsGC0,abs(Args->GC[0]));
        Args->ForceFullEwald=true;
        GetEdgeEdgeInteractions(Args);
        Args->ForceFullEwald=false;
      };
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
void InitGetEEIArgs(GetEEIArgStruct *Args)
{
  Args->NumGradientComponents=0;
  Args->NumTorqueAxes=0;
  Args->GammaMatrix=0;
  Args->Displacement=0;
  Args->opFC=0;
  Args->Force=EEI_NOFORCE;
  Args->GBA=0;
  Args->ForceFullEwald=false;
  memset(Args->PPIAlgorithmCount, 0, NUMPPIALGORITHMS*sizeof(unsigned));
}

} // namespace scuff
