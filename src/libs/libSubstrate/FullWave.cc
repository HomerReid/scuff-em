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
void LayeredSubstrate::GetSubstrateDGF_StaticLimit(cdouble Omega,
                                                   double *XD, double *XS,
                                                   cdouble Gij[6][6])
{
  /***************************************************************/
  /* explain me **************************************************/
  /***************************************************************/
  bool Total = ( GetRegionIndex(XD[2]) != GetRegionIndex(XS[2]) );

  /***************************************************************/
  /* EM, ME quadrants determined by E-fields of electrostatic    */
  /* monopole                                                    */
  /***************************************************************/
  double PhiE0[4], *E0=PhiE0+1;
  GetPhiE(XD, XS, PhiE0, Total);
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
     GetPhiE(XD, XSP, PhiEP, Total);
     GetPhiE(XD, XSM, PhiEM, Total);
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
     XMatrix->GetEntriesD(nx,"0:5",XD);

     cdouble ScriptG[6][6];
     GetSubstrateDGF_StaticLimit(Omega, XD, XS, ScriptG);
     for(int Mu=0; Mu<6; Mu++)
      for(int Nu=0; Nu<6; Nu++)
       GMatrix->SetEntry(6*Mu+Nu, nx, ScriptG[Mu][Nu] );
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
  if ( GMatrix && (GMatrix->NR!=36 || GMatrix->NC!=NX) )
   { Warn("wrong-size GMatrix passed to GetSubstrateDGF: (%i,%i)!=(%i,%i) reallocating...",GMatrix->NR,GMatrix->NC,36,NX);
     delete GMatrix;
     GMatrix=0;
   };
  if (GMatrix==0)
   GMatrix = new HMatrix(36, NX, LHM_COMPLEX);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (    ForceMethod==AUTO 
       && GetSubstrateDGF_Interp(Omega, XMatrix, GMatrix)
     ) return GMatrix;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (Omega==0.0) Method=STATIC_LIMIT;
  if (ForceMethod!=AUTO) Method=ForceMethod;
  switch(Method)
   { 
     case STATIC_LIMIT:
       GetSubstrateDGF_StaticLimit(Omega, XMatrix, GMatrix);
       break;

     case FULL_SURFACE_CURRENT:
       GetSubstrateDGF_FullSurfaceCurrent(Omega, XMatrix, GMatrix);
       break;

     case FAST_SURFACE_CURRENT:
     case AUTO:
     default:
       GetSubstrateDGF_FastSurfaceCurrent(Omega, XMatrix, GMatrix);
       break;
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

/***************************************************************/
/***************************************************************/
/***************************************************************/
void LayeredSubstrate::GetSubstrateDGF(cdouble Omega,
                                       double XD[3], double XS[3],
                                       cdouble ScriptG[6][6], DGFMethod Method)
{ 
  double XBuffer[6];
  HMatrix XMatrix(1,6,LHM_REAL,XBuffer);
  XMatrix.SetEntriesD(0,"0:2",XD);
  XMatrix.SetEntriesD(0,"3:5",XS);

  cdouble GBuffer[36];
  HMatrix GMatrix(36,1,LHM_COMPLEX,GBuffer);

  GetSubstrateDGF(Omega, &XMatrix, &GMatrix, Method);

  for(int Mu=0; Mu<6; Mu++)
   for(int Nu=0; Nu<6; Nu++)
    ScriptG[Mu][Nu] = GMatrix.GetEntry(6*Mu + Nu,0);  
}
