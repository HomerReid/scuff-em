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
 * libSubstrate/FullWave.cc
 *
 * homer reid             -- 9/2017
 *
 */
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>

#include "libhrutil.h"
#include "libMDInterp.h"
#include "libSGJC.h"
#include "libSpherical.h"
#include "libSubstrate.h"

#define II cdouble(0.0,1.0)

/******************************************************************/
/* on entry                                                       */
/*  G = G_ij in a coordinate system in which (x,y)=(Rho,0)        */
/* on return:                                                     */
/*  G = G_ij in a coordinate system in which (x,y)=(Rho*CP,Rho*SP)*/
/*******************************************************************/
void LayeredSubstrate::RotateG(cdouble G[6][6], int P, int Q, double CP, double SP)
{
  cdouble GP[3][3];

  GP[0][0]=CP*CP*G[P+0][Q+0] - CP*SP*GP[P+0][Q+1] - CP*SP*G[P+1][Q+0] + SP*SP*G[P+1][Q+1];
  GP[0][1]=CP*CP*G[P+0][Q+1] + CP*SP*GP[P+0][Q+0] - CP*SP*G[P+1][Q+1] - SP*SP*G[P+1][Q+0];
  GP[1][0]=CP*CP*G[P+1][Q+0] + CP*SP*GP[P+0][Q+0] - CP*SP*G[P+1][Q+1] - SP*SP*G[P+0][Q+1];
  GP[1][1]=CP*CP*G[P+1][Q+1] + CP*SP*GP[P+0][Q+1] + CP*SP*G[P+1][Q+0] + SP*SP*G[P+0][Q+0];

  GP[0][2]=CP*G[P+0][Q+2] - SP*G[P+1][Q+2];
  GP[2][0]=CP*G[P+2][Q+0] - SP*G[P+2][Q+1];

  GP[1][2]=CP*G[P+1][Q+2] + SP*G[P+0][Q+2];
  GP[2][1]=CP*G[P+2][Q+1] + SP*G[P+2][Q+0];

  GP[2][2]=G[P+2][Q+2];

  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    G[P+Mu][Q+Nu] = GP[Mu][Nu];
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void LayeredSubstrate::RotateG(cdouble Gij[6][6], double Phi)
{
  double CP=cos(Phi), SP=sin(Phi);
  for(int p=0; p<2; p++)
   for(int q=0; q<2; q++)
    RotateG(Gij, 3*p, 3*q, CP, SP);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void LayeredSubstrate::GetSubstrateDGF_StaticLimit(cdouble Omega,
                                                   double *XD, double *XS,
                                                   cdouble Gij[6][6])
{
  /***************************************************************/
  /* EM, ME quadrants determined by E-fields of electrostatic    */
  /* monopole                                                    */
  /***************************************************************/
  double PhiE0[4], *E0=PhiE0+1;
  GetTotalPhiE(XD, XS, PhiE0);
  for(int i=0; i<3; i++)
   { int j=(i+1)%3, k=(i+2)%3;
     Gij[i][3+i] = Gij[3+i][i] = 0.0;
     Gij[i][3+j] = Gij[3+j][i] = +1.0*E0[1+k];
     Gij[j][3+i] = Gij[3+i][j] = -1.0*E0[1+k];
   };

  /***************************************************************/
  /* EE, MM quadrants determined by E-fields of electrostatic    */
  /* dipole = derivatives of monopole fields wrt source point    */
  /***************************************************************/
  for(int j=0; j<3; j++)
   { 
     double XSP[3], XSM[3];
     XSP[0] = XSM[0] = XS[0];
     XSP[1] = XSM[1] = XS[1];
     XSP[2] = XSM[2] = XS[2];

     double DeltaX = 1.0e-4*fabs(XS[j]);
     if (DeltaX==0.0)
      DeltaX=1.0e-4;
     XSP[j] += DeltaX;
     XSM[j] -= DeltaX;

     double PhiEP[4], *EP=PhiEP+1, PhiEM[4], *EM=PhiEM+1, djEi[3];
     GetTotalPhiE(XD, XSP, PhiEP);
     GetTotalPhiE(XD, XSM, PhiEM);
     djEi[0] = (EP[0] - EM[0]) / (2.0*DeltaX);
     djEi[1] = (EP[1] - EM[1]) / (2.0*DeltaX);
     djEi[2] = (EP[2] - EM[2]) / (2.0*DeltaX);

     cdouble EEPrefac = II*ZVAC / Omega;
     cdouble MMPrefac = II/(ZVAC*Omega);
     for(int i=0; i<3; i++)
       { Gij[0+i][0+j] = EEPrefac * djEi[i];
         Gij[3+i][3+j] = MMPrefac * djEi[i];
       };
   };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void LayeredSubstrate::GetSubstrateDGF_StaticLimit(cdouble Omega,
                                                   HMatrix *XMatrix,
                                                   HMatrix *GMatrix)
{
  for(int nx=0; nx<XMatrix->NR; nx++)
   {
     double XD[6], *XS=XD+3;
     XMatrix->GetEntriesD(nx,":",XD);

     cdouble ScriptG[6][6];
     GetSubstrateDGF_StaticLimit(Omega, XD, XS, ScriptG);
     for(int Mu=0; Mu<6; Mu++)
      for(int Nu=0; Nu<6; Nu++)
       GMatrix->SetEntry(nx, 6*Mu+Nu, ScriptG[Mu][Nu] );
   };
}

/***************************************************************/
/* switchboard routine *****************************************/
/***************************************************************/
HMatrix *LayeredSubstrate::GetSubstrateDGF(cdouble Omega,
                                           HMatrix *XMatrix,
                                           HMatrix *GMatrix,
                                           DGFMethod Method)
{
  int NX=XMatrix->NR;
  Log("Computing substrate DGF at %i points...",NX);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if ( GMatrix && (GMatrix->NR!=NX || GMatrix->NC!=36) )
   { Warn("wrong-size GMatrix passed to GetSubstrateDGF: (%i,%i)!=(%i,%i) reallocating...",GMatrix->NR,GMatrix->NC,NX,36);
     delete GMatrix;
     GMatrix=0;
   };
  if (GMatrix==0)
   GMatrix = new HMatrix(NX, 36, LHM_COMPLEX);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (Omega==0.0) Method=STATIC_LIMIT;
  switch(Method)
   { 
     case SURFACE_CURRENT:
       GetSubstrateDGF_SurfaceCurrent(Omega, XMatrix, GMatrix);
       break;

     case STATIC_LIMIT:
     default:
       GetSubstrateDGF_StaticLimit(Omega, XMatrix, GMatrix);
       break;

  //   case PLANE_WAVE:
  //      GetSubstrateDGF_PlaneWave(Omega, XMatrix, GMatrix);
  //     break;
   };
     
  return GMatrix;
}  

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *LayeredSubstrate::GetSubstrateDGF(cdouble Omega,
                                           HMatrix *XMatrix,
                                           DGFMethod Method)
{ return GetSubstrateDGF(Omega, XMatrix, 0, Method); }
