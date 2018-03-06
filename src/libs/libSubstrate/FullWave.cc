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

namespace scuff{
void CalcGC(double R[3], cdouble Omega,
            cdouble EpsR, cdouble MuR,
            cdouble GMuNu[3][3], cdouble CMuNu[3][3],
            cdouble GMuNuRho[3][3][3], cdouble CMuNuRho[3][3][3]);
}

/***************************************************************/
/* add the 6x6 homogeneous DGF for a medium with given         */
/* relative permittivity and permeability                      */
/* on return, Gamma0 contains the 6x6 matrix as a column-major */
/* array, i.e. Gamma0 could be the ZM field of a 6x6 HMatrix   */
/* (NOT a [6][6] C++ array).                                   */
/***************************************************************/
void AddGamma0(double XD[3], double XS[3], cdouble Omega,
               cdouble EpsRel, cdouble MuRel, cdouble *Gamma0,
               double zGP, bool Image)
{
  double R[3];
  R[0]=XD[0]-XS[0];
  R[1]=XD[1]-XS[1];
  R[2]=XD[2]-XS[2];
  double r2=R[0]*R[0] + R[1]*R[1] + R[2]*R[2], r=sqrt(r2), r3=r*r2;
  cdouble k2=EpsRel*MuRel*Omega*Omega, k=sqrt(k2);
  cdouble ik=II*k, ik2=ik*ik, ikr=ik*r, ikr2=ikr*ikr;
  
  cdouble G[3][3], C[3][3];
  memset( (cdouble *)G, 0, 9*sizeof(cdouble));
  memset( (cdouble *)C, 0, 9*sizeof(cdouble));

  // compute 3x3 G and C matrices
  if (r2==0.0)
   {
     G[0][0] = G[1][1] = G[2][2] = II*real(k)/(6.0*M_PI);
   }
  else
   { cdouble ExpFac=exp(ik*r) / (4.0*M_PI*r3);
     cdouble f2=(1.0-ikr+ikr2)/ik2, f3=(-3.0+3.0*ikr-ikr2)/(ikr2);
     G[0][0] = G[1][1] = G[2][2] = f2*ExpFac;
     for(int i=0; i<3; i++)
      for(int j=0; j<3; j++)
       G[i][j] += f3*ExpFac*R[i]*R[j];
  
     ExpFac*=(-1.0+ikr)/(ik);
     C[0][1]=ExpFac*R[2]; C[1][0]=-C[0][1];
     C[1][2]=ExpFac*R[0]; C[2][1]=-C[1][2];
     C[2][0]=ExpFac*R[1]; C[0][2]=-C[2][0];
   }

  // add quadrants of 6x6 DGF
  cdouble PreFac[2][2];
  PreFac[0][0]=II*Omega*MuRel*ZVAC;
  PreFac[0][1]=+1.0*II*k;
  PreFac[1][0]=-1.0*II*k;
  PreFac[1][1]=II*Omega*EpsRel/ZVAC;
  HMatrix Gamma0Matrix(6,6,LHM_COMPLEX,Gamma0);
  double ImageSign[2][3]={{1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}};
  if (Image) 
   ImageSign[0][0]=ImageSign[0][1]=ImageSign[2][2]=-1.0;
  for(int P=0; P<2; P++)
   for(int Q=0; Q<2; Q++)
    for(int i=0; i<3; i++)
     for(int j=0; j<3; j++)
      Gamma0Matrix.AddEntry(3*P+i, 3*Q+j, ImageSign[Q][j]*PreFac[P][Q]*(P==Q ? G[i][j] : C[i][j]));

  if (!isinf(zGP))
   { R[2] = XD[2] + XS[2] - 2.0*zGP;
     double Zero[3]={0.0, 0.0, 0.0};
     AddGamma0(R, Zero, Omega, EpsRel, MuRel, Gamma0, -1.0*HUGE_VAL, true);
   };

}

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
  bool Total = ( GetLayerIndex(XD[2]) != GetLayerIndex(XS[2]) );

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

     cdouble ScriptGStatic[6][6];
     GetSubstrateDGF_StaticLimit(Omega, XD, XS, ScriptGStatic);

     cdouble *ScriptGBuffer = GMatrix->ZM + nx*36;
     HMatrix ScriptG(6,6,LHM_COMPLEX,ScriptGBuffer);
     for(int Mu=0; Mu<6; Mu++)
      for(int Nu=0; Nu<6; Nu++)
       ScriptG.SetEntry(Mu, Nu, ScriptGStatic[Mu][Nu]);
   };
}

/***************************************************************/
/* given a vector of 22 gFrak values and Theta, fill in the 6x6 tensor ScriptG(Rho) */
/***************************************************************/
void LayeredSubstrate::gFrakToScriptG(cdouble gFrak[NUMGFRAK], double Theta,
                                      cdouble ScriptGBuffer[36])
{
  double CT=cos(Theta), ST=sin(Theta);

  HMatrix ScriptG(6,6,LHM_COMPLEX,ScriptGBuffer);

  ScriptG.SetEntry(0+0, 0+0,      gFrak[_EE2A]*CT*CT + gFrak[_EE0P] + gFrak[_EE2B]);
  ScriptG.SetEntry(0+0, 0+1,      gFrak[_EE2A]*CT*ST);
  ScriptG.SetEntry(0+0, 0+2,      gFrak[_EE1A]*CT);
  ScriptG.SetEntry(0+1, 0+0,      gFrak[_EE2A]*ST*CT);
  ScriptG.SetEntry(0+1, 0+1,      gFrak[_EE2A]*ST*ST + gFrak[_EE0P] + gFrak[_EE2B]);
  ScriptG.SetEntry(0+1, 0+2,      gFrak[_EE1A]*ST);
  ScriptG.SetEntry(0+2, 0+0,      gFrak[_EE1B]*CT);
  ScriptG.SetEntry(0+2, 0+1,      gFrak[_EE1B]*ST);
  ScriptG.SetEntry(0+2, 0+2,      gFrak[_EE0Z]);

  ScriptG.SetEntry(0+0, 3+0,      gFrak[_EM2A]*CT*ST);
  ScriptG.SetEntry(0+0, 3+1,      gFrak[_EM2A]*ST*ST + gFrak[_EM0P]);
  ScriptG.SetEntry(0+0, 3+2, -1.0*gFrak[_EM1A]*ST);
  ScriptG.SetEntry(0+1, 3+0, -1.0*gFrak[_EM2A]*CT*CT - gFrak[_EM0P]);
  ScriptG.SetEntry(0+1, 3+1, -1.0*gFrak[_EM2A]*CT*ST);
  ScriptG.SetEntry(0+1, 3+2,      gFrak[_EM1A]*CT);
  ScriptG.SetEntry(0+2, 3+0, -1.0*gFrak[_EM1B]*ST);
  ScriptG.SetEntry(0+2, 3+1,      gFrak[_EM1B]*CT);
  ScriptG.SetEntry(0+2, 3+2,      gFrak[_EM1A] + gFrak[_EM1B]);

  ScriptG.SetEntry(3+0, 0+0,      gFrak[_ME2A]*CT*ST);
  ScriptG.SetEntry(3+0, 0+1,      gFrak[_ME2A]*ST*ST + gFrak[_ME0P]);
  ScriptG.SetEntry(3+0, 0+2, -1.0*gFrak[_ME1A]*ST);
  ScriptG.SetEntry(3+1, 0+0, -1.0*gFrak[_ME2A]*CT*CT - gFrak[_ME0P]);
  ScriptG.SetEntry(3+1, 0+1, -1.0*gFrak[_ME2A]*CT*ST);
  ScriptG.SetEntry(3+1, 0+2,      gFrak[_ME1A]*CT);
  ScriptG.SetEntry(3+2, 0+0, -1.0*gFrak[_ME1B]*ST);
  ScriptG.SetEntry(3+2, 0+1,      gFrak[_ME1B]*CT);
  ScriptG.SetEntry(3+2, 0+2,      gFrak[_ME1A] + gFrak[_ME1B]);

  ScriptG.SetEntry(3+0, 3+0,      gFrak[_MM2A]*CT*CT + gFrak[_MM0P] + gFrak[_MM2B]);
  ScriptG.SetEntry(3+0, 3+1,      gFrak[_MM2A]*CT*ST);
  ScriptG.SetEntry(3+0, 3+2,      gFrak[_MM1A]*CT);
  ScriptG.SetEntry(3+1, 3+0,      gFrak[_MM2A]*ST*CT);
  ScriptG.SetEntry(3+1, 3+1,      gFrak[_MM2A]*ST*ST + gFrak[_MM0P] + gFrak[_MM2B]);
  ScriptG.SetEntry(3+1, 3+2,      gFrak[_MM1A]*ST);
  ScriptG.SetEntry(3+2, 3+0,      gFrak[_MM1B]*CT);
  ScriptG.SetEntry(3+2, 3+1,      gFrak[_MM1B]*ST);
  ScriptG.SetEntry(3+2, 3+2,      gFrak[_MM0Z]);

}
   
/***************************************************************/
/***************************************************************/
/***************************************************************/
void LayeredSubstrate::GetSubstrateDGF_SurfaceCurrent(cdouble Omega,
                                                      HMatrix *XMatrix,
                                                      HMatrix *ScriptGMatrix,
                                                      bool SubtractQS)
{
  UpdateCachedEpsMu(Omega);

  /***************************************************************/
  /* Get gFrak for all evaluation points (and store it in the    */
  /* ScriptGMatrix buffer for convenience)                       */
  /***************************************************************/
  cdouble *Workspace=0;
  GetgFrak(Omega, XMatrix, ScriptGMatrix->ZM, Workspace, SubtractQS);

  /***************************************************************/
  /*- for each evaluation point, get the full 6x6 substrate DGF  */
  /*- from the vector of gFrak components                        */
  /***************************************************************/
  for(int nx=XMatrix->NR-1; nx>=0; nx--)
   { 
     double Rhox = XMatrix->GetEntryD(nx,0) - XMatrix->GetEntryD(nx,3);
     double Rhoy = XMatrix->GetEntryD(nx,1) - XMatrix->GetEntryD(nx,4);
     double Theta = atan2(Rhoy,Rhox);

     cdouble gFrak[NUMGFRAK];
     memcpy(gFrak, ScriptGMatrix->ZM + nx*NUMGFRAK, NUMGFRAK*sizeof(cdouble));

     gFrakToScriptG(gFrak, Theta, ScriptGMatrix->ZM + nx*36);
   };
}

/***************************************************************/
/* switchboard routine *****************************************/
/***************************************************************/
HMatrix *LayeredSubstrate::GetSubstrateDGF(cdouble Omega,
                                           HMatrix *XMatrix,
                                           HMatrix *GMatrix,
                                           bool AddHomogeneousDGF,
                                           bool SubtractQS)
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
  int Method=SURFACE_CURRENT;
  if (Omega==0.0) Method=STATIC_LIMIT;
  if (ForceMethod!=AUTO) Method=ForceMethod;
  switch(Method)
   { 
     case STATIC_LIMIT:
       GetSubstrateDGF_StaticLimit(Omega, XMatrix, GMatrix);
       break;

     case SURFACE_CURRENT:
     case AUTO:
     default:
       GetSubstrateDGF_SurfaceCurrent(Omega, XMatrix, GMatrix, SubtractQS);
       break;
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
#if 0
  if (AddHomogeneousDGF)
   for(int nx=0; nx<XMatrix->NR; nx++)
    { 
      double XDS[6];
      XMatrix->GetEntriesD(nx,":",XDS);
      int nl=GetLayerIndex(XDS[2]);
      if (nl!=GetLayerIndex(XDS[5])) continue;
      cdouble EpsRel=EpsLayer[nl], MuRel=MuLayer[nl];
      double R[3];
      VecSub(XDS+0, XDS+3, R);
      cdouble G[3][3], C[3][3], dG[3][3][3], dC[3][3][3];
      scuff::CalcGC(R, Omega, EpsRel, MuRel, G, C, dG, dC);

      cdouble k=sqrt(EpsRel*MuRel)*Omega;
      cdouble ZRel=sqrt(MuRel/EpsRel);
      cdouble PreFac[2][2];
      PreFac[0][0] =      II*k*ZVAC*ZRel;   // EE
      PreFac[0][1] =      II*k;             // EM
      PreFac[1][0] = -1.0*II*k;             // ME
      PreFac[1][1] =      II*k/(ZVAC*ZRel); // ME

      HMatrix ScriptG(6,6,LHM_COMPLEX,GMatrix->ZM + nx*36);
      for(int P=0; P<2; P++)
       for(int Q=0; Q<2; Q++)
        for(int i=0; i<3; i++)
         for(int j=0; j<3; j++)
          ScriptG.AddEntry(3*P+i, 3*Q+j, PreFac[P][Q]*(P==Q ? G[i][j] : C[i][j]));
    };
#endif
#if 1
  if (AddHomogeneousDGF && !ForceFreeSpace)
   for(int nx=0; nx<XMatrix->NR; nx++)
    { double XDS[6];
      XMatrix->GetEntriesD(nx,":",XDS);
      if ( XDS[2]<zGP || XDS[5]<zGP )
       continue;
      int nl=GetLayerIndex(XDS[2]);
      if (nl!=GetLayerIndex(XDS[5])) continue;
      cdouble EpsRel=EpsLayer[nl], MuRel=MuLayer[nl];
      double MyzGP = (nl==NumInterfaces ? zGP : -1.0*HUGE_VAL);
      AddGamma0(XDS+0, XDS+3, Omega, EpsRel, MuRel, GMatrix->ZM + nx*36, MyzGP);
    };
#endif
    
    
     
  return GMatrix;
}  

/***************************************************************/
/***************************************************************/
/***************************************************************/
void LayeredSubstrate::GetSubstrateDGF(cdouble Omega, double XD[3], double XS[3],
                                       cdouble ScriptG[6][6],
                                       bool AddHomogeneousDGF, bool SubtractQS)
{ 
  double XBuffer[6];
  HMatrix XMatrix(1,6,LHM_REAL,XBuffer);
  XMatrix.SetEntriesD(0,"0:2",XD);
  XMatrix.SetEntriesD(0,"3:5",XS);

  cdouble GBuffer[36];
  HMatrix GMatrix(36,1,LHM_COMPLEX,GBuffer);

  GetSubstrateDGF(Omega, &XMatrix, &GMatrix, AddHomogeneousDGF, SubtractQS);

  HMatrix ScriptGMatrix(6,6,LHM_COMPLEX,GBuffer);
  for(int Mu=0; Mu<6; Mu++)
   for(int Nu=0; Nu<6; Nu++)
    ScriptG[Mu][Nu] = ScriptGMatrix.GetEntry(Mu,Nu);
}
