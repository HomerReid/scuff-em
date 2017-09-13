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

/***************************************************************/
/***************************************************************/
/***************************************************************/
void LayeredSubstrate::GetFullWaveDGF_Static(cdouble Omega,
                                             double *XD, double *XS,
                                             cdouble ScriptG[6][6])
{
  int nrDest  = GetRegionIndex(XD[2]);
  cdouble EpsDest, MuDest;
  MPLayer[nrDest]->GetEpsMu(Omega, &EpsDest, &MuDest);
  cdouble k2 = EpsDest*MuDest*Omega*Omega;

  // finite-differencing of electrostatic monopole fields
  // with respect to source location
  cdouble GradG0[3];
  cdouble Gij[3][3];
  for(int Nu=0; Nu<3; Nu++)
   { 
     double XSP[3], XSM[3];
     XSP[0] = XSM[0] = XS[0];
     XSP[1] = XSM[1] = XS[1];
     XSP[2] = XSM[2] = XS[2];

     double DeltaX = 1.0e-4*fabs(XS[Nu]);
     if (DeltaX==0.0)
      DeltaX=1.0e-4;
     XSP[Nu] += DeltaX;
     XSM[Nu] -= DeltaX;

     double PhiEP[4], PhiEM[4];
     GetDeltaPhiE(XD, XSP, PhiEP);
     GetDeltaPhiE(XD, XSM, PhiEM);
     GradG0[Nu] = (PhiEP[0] - PhiEM[0]) / (2.0*DeltaX);

     for(int Mu=0; Mu<3; Mu++)
      Gij[Mu][Nu] = (PhiEP[1+Mu] - PhiEM[1+Mu])/(2.0*DeltaX);
   };

  cdouble ikCij[3][3];
  ikCij[0][1] = GradG0[2];    ikCij[1][0] = -1.0*ikCij[0][1];
  ikCij[1][2] = GradG0[0];    ikCij[2][1] = -1.0*ikCij[1][2];
  ikCij[2][0] = GradG0[1];    ikCij[0][2] = -1.0*ikCij[2][0];
  ikCij[0][0] = ikCij[1][1] = ikCij[2][2] = 0.0;

  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    { ScriptG[0+Mu][0+Nu] = EpsDest*ZVAC*Gij[Mu][Nu] / k2;
      ScriptG[0+Mu][3+Nu] = +ikCij[Mu][Nu];
      ScriptG[3+Mu][0+Nu] = -ikCij[Mu][Nu];
      ScriptG[3+Mu][3+Nu] = MuDest*ZVAC*Gij[Mu][Nu] / k2;
    };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void LayeredSubstrate::GetFullWaveDGF_Static(cdouble Omega,
                                             HMatrix *XMatrix,
                                             HMatrix *GMatrix)
{
  for(int nx=0; nx<XMatrix->NR; nx++)
   {
     double XD[6], *XS=XD+3;
     XMatrix->GetEntriesD(nx,":",XD);
     if (XMatrix->NC==3)
      memcpy(XS, XD, 3*sizeof(double));

     cdouble ScriptG[6][6];
     GetFullWaveDGF_Static(Omega, XD, XS, ScriptG);
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
  if (Omega==0.0) Method=STATIC;
  switch(Method)
   { 
     case PLANE_WAVE:
       GetFullWaveDGF_PlaneWave(Omega, XMatrix, GMatrix);
       break;

     case STATIC:
     default:
       GetFullWaveDGF_Static(Omega, XMatrix, GMatrix);
       break;

  //   case PLANE_WAVE:
  //      GetFullWaveDGF_PlaneWave(Omega, XMatrix, GMatrix);
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
