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

     cdouble ScriptG[6][6];
     GetSubstrateDGF_StaticLimit(Omega, XD, XS, ScriptG);
     for(int Mu=0; Mu<6; Mu++)
      for(int Nu=0; Nu<6; Nu++)
       GMatrix->SetEntry(6*Mu+Nu, nx, ScriptG[Mu][Nu] );
   };
}

/***************************************************************/
/* Get the substrate Green's function by evaluating the full   */
/* 2D Fourier integral with no fancy accelerations. This is    */
/* too slow for use in practical calculations but offers a     */
/* helpful sanity check for debugging, etc.                    */
/***************************************************************/
void LayeredSubstrate::GetSubstrateDGF_FullSurfaceCurrent(cdouble Omega,
                                                          HMatrix *XMatrix,
                                                          HMatrix *GMatrix)
{
  UpdateCachedEpsMu(Omega);
  (void )XMatrix;
  (void )GMatrix;
#if 0

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  qFunctionData MyData, *Data=&MyData;
  Data->XMatrix    = XMatrix;
  Data->RTwiddle   = new HMatrix(6, 4*NumInterfaces, LHM_COMPLEX);
  Data->WMatrix    = new HMatrix(4*NumInterfaces, 4*NumInterfaces, LHM_COMPLEX);
  Data->STwiddle   = new HMatrix(4*NumInterfaces, 6,               LHM_COMPLEX);
  Data->GTwiddle   = new HMatrix(6,               6,               LHM_COMPLEX);
  Data->byqFile    = WritebyqFiles ? fopen("/tmp/q2D.log","w") : 0;

  int FDim = 36*XMatrix->NR;
  bool ThetaSymmetric=false;
  qIntegrate(Omega, qIntegrandFullSC, (void *)Data, GMatrix->ZM, FDim, ThetaSymmetric);
  
  if (Data->byqFile) fclose(Data->byqFile);

  delete Data->RTwiddle;
  delete Data->WMatrix;
  delete Data->STwiddle;
  delete Data->GTwiddle;
#endif
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetLambdaMatrices(double Theta, double Lambda[7][3][3])
{ 
  double CT=cos(Theta), ST=sin(Theta);
  memset( (double *)Lambda, 0, 7*3*3*sizeof(double));

  Lambda[_LAMBDA0P][0][0]=1.0;
  Lambda[_LAMBDA0P][1][1]=1.0;
  Lambda[_LAMBDA0Z][2][2]=1.0;

  Lambda[_LAMBDA1 ][0][2]=CT;
  Lambda[_LAMBDA1 ][1][2]=ST;

  Lambda[_LAMBDA2 ][0][0]=CT*CT;
  Lambda[_LAMBDA2 ][0][1]=CT*ST;
  Lambda[_LAMBDA2 ][1][0]=ST*CT;
  Lambda[_LAMBDA2 ][1][1]=ST*ST;

  Lambda[_LAMBDA0X][0][1]=1.0;
  Lambda[_LAMBDA0X][1][0]=-1.0;

  Lambda[_LAMBDA1X][0][2]=-ST;
  Lambda[_LAMBDA1X][1][2]=+CT;
  Lambda[_LAMBDA1X][2][2]=1.0;

  Lambda[_LAMBDA2X][0][0]=CT*ST;
  Lambda[_LAMBDA2X][1][1]=-CT*ST;
  Lambda[_LAMBDA2X][0][1]=ST*ST;
  Lambda[_LAMBDA2X][1][0]=-CT*CT;
}

/***************************************************************/
/* compute a, c parameters for the contour-integral portion of */
/* the Sommerfeld integral via the procedure described in      */
/* Golubovic et al, "Efficient Algorithms for Computing        */
/* Sommerfeld Intergral Tails," IEEE Transactions on Antennas  */
/* and Propagation **60** 2409 (2012).                         */
/* DOI: 10.1109/TAP.2012.2189718                               */
/*--------------------------------------------------------------*/
void GetacSommerfeld(LayeredSubstrate *S, cdouble Omega,
                     double Rho, double zDest, double zSource,
                     double *pa, double *pc)
{
  int NumLayers     = S->NumLayers;
  cdouble *EpsLayer = S->EpsLayer;
  cdouble *MuLayer  = S->MuLayer;

  double MaxRealn2 = 0.0;
  for(int nl=0; nl<NumLayers; nl++)
   MaxRealn2 = fmax(MaxRealn2, real(EpsLayer[nl]*MuLayer[nl]));
  double a = 1.5*sqrt(MaxRealn2)*real(Omega);
  double c = real(Omega);
   
  if ( Rho>fabs(zDest-zSource) && real(Omega*Rho)>1.0 )
   c = 1.0/Rho;

  *pa=a;
  *pc=c;
  
}
   

/***************************************************************/
/***************************************************************/
/***************************************************************/
void LayeredSubstrate::GetSubstrateDGF_FastSurfaceCurrent(cdouble Omega,
                                                          HMatrix *XMatrix,
                                                          HMatrix *GMatrix)
{
  UpdateCachedEpsMu(Omega);

  /***************************************************************/
  /* evaluate Sommerfeld integral to get 'gScalar' integrals for  */
  /* all evaluation points                                        */
  /***************************************************************/
  SommerfeldIntegrandData MySID, *SID=&MySID;
  SID->Substrate  = this;
  SID->Omega      = Omega;
  SID->q0         = 0.0;
  SID->uTransform = false;
  SID->XMatrix    = XMatrix;
  SID->RTwiddle   = new HMatrix(6, 4*NumInterfaces, LHM_COMPLEX);
  SID->WMatrix    = new HMatrix(4*NumInterfaces, 4*NumInterfaces, LHM_COMPLEX);
  SID->STwiddle   = new HMatrix(4*NumInterfaces, 6,               LHM_COMPLEX);
  SID->NumPoints  = 0;
  SID->byqFile    = 0;

  double Rhox  = XMatrix->GetEntryD(0,0) - XMatrix->GetEntryD(0,3);
  double Rhoy  = XMatrix->GetEntryD(0,1) - XMatrix->GetEntryD(0,4);
  double Rho   = sqrt(Rhox*Rhox + Rhoy*Rhoy);
  double zDest = XMatrix->GetEntryD(0,2), zSource=XMatrix->GetEntryD(0,5);
  
  //MySID->byqFile=WritebyqFiles ? fopen("/tmp/q2D.log","w") : 0;

  int NX    = XMatrix->NR;
  int zfdim = NUMGSCALAR*NX;
  int fdim  = 2*zfdim;

  cdouble *Error = new cdouble[zfdim];
 
  char *s=getenv("SOMMERFELD_BYQFILE");
  if (s)
   SID->byqFile=fopen(s,"w");

  s=getenv("SOMMERFELD_INTEGRATOR_BYPASS");
  if (s && s[0]=='1')
   { SID->uTransform=true;
     double uMin=0.0, uMax=1.0;
     Log("Computing gScalar integrals via simple quadrature...");
     int ndim=1;
     pcubature(fdim, SommerfeldIntegrand, (void *)SID,
               ndim, &uMin, &uMax, qMaxEval, qAbsTol, qRelTol, ERROR_PAIRED,
               (double *)(GMatrix->ZM), (double *)Error);
     Log("...%i points",SID->NumPoints);
   }
  else
   { bool Verbose=true;
     int xNu=0;
     s=getenv("SCUFF_SOMMERFELD_XNU");
     if (s)
      { xNu = s[0] - '0';
        Log("Setting xNu=%i.\n",xNu);
      };
 
     double a, c;
     GetacSommerfeld(this, Omega, Rho, zDest, zSource, &a, &c);

     Log("Computing gScalar integrals via Sommerfeld integrator...");
     SommerfeldIntegrate(SommerfeldIntegrand, (void *)SID, zfdim,
                         a, c, xNu, Rho, qMaxEval, qMaxEval,
                         qAbsTol, qRelTol, GMatrix->ZM, Error, Verbose);
     Log("...%i points",SID->NumPoints);
   };
  if (SID->byqFile) fclose(SID->byqFile);
  delete[] Error;

  /*--------------------------------------------------------------*/
  /*- for each evaluation point, get the full 6x6 substrate DGF  -*/
  /*- from the gScalar integrals, then add the homogeneous DGF   -*/
  /*- if necessary                                               -*/
  /*--------------------------------------------------------------*/
  for(int nx=NX-1; nx>=0; nx--)
   { 
     Rhox = XMatrix->GetEntryD(nx,0) - XMatrix->GetEntryD(nx,3);
     Rhoy = XMatrix->GetEntryD(nx,1) - XMatrix->GetEntryD(nx,4);
     double Theta = atan2(Rhoy,Rhox);
     double CT=cos(Theta), ST=sin(Theta);

     cdouble g[NUMGSCALAR];
     memcpy(g, GMatrix->ZM + nx*NUMGSCALAR, NUMGSCALAR*sizeof(cdouble));

     cdouble *G = GMatrix->ZM + nx*36;
     G[_EEXX] = g[_EE2A]*CT*CT + g[_EE0P] + g[_EE2B];
     G[_EEXY] = g[_EE2A]*CT*ST;
     G[_EEXZ] = g[_EE1A]*CT;
     G[_EEYX] = g[_EE2A]*CT*ST;
     G[_EEYY] = g[_EE2A]*ST*ST + g[_EE0P] + g[_EE2B];
     G[_EEYZ] = g[_EE1A]*ST;
     G[_EEZX] = g[_EE1B]*CT;
     G[_EEZY] = g[_EE1B]*ST;
     G[_EEZZ] = g[_EE0Z];

     G[_EMXX] =      g[_EM2A]*CT*ST;
     G[_EMXY] =      g[_EM2A]*ST*ST + g[_EM0P];
     G[_EMXZ] = -1.0*g[_EM1A]*ST;
     G[_EMYX] = -1.0*g[_EM2A]*CT*CT - g[_EM0P];
     G[_EMYY] = -1.0*g[_EM2A]*CT*ST;
     G[_EMYZ] =      g[_EM1A]*CT;
     G[_EMZX] = -1.0*g[_EM1B]*ST;
     G[_EMZY] =      g[_EM1B]*CT;
     G[_EMZZ] =      g[_EM1A]+g[_EM1B];

     G[_MEXX] =      g[_ME2A]*CT*ST;
     G[_MEXY] =      g[_ME2A]*ST*ST + g[_ME0P];
     G[_MEXZ] = -1.0*g[_ME1A]*ST;
     G[_MEYX] = -1.0*g[_ME2A]*CT*CT - g[_ME0P];
     G[_MEYY] = -1.0*g[_ME2A]*CT*ST;
     G[_MEYZ] =      g[_ME1A]*CT;
     G[_MEZX] = -1.0*g[_ME1B]*ST;
     G[_MEZY] =      g[_ME1B]*CT;
     G[_MEZZ] =      g[_ME1A]+g[_ME1B];

     G[_MMXX] = g[_MM2A]*CT*CT + g[_MM0P] + g[_MM2B];
     G[_MMXY] = g[_MM2A]*CT*ST;
     G[_MMXZ] = g[_MM1A]*CT;
     G[_MMYX] = g[_MM2A]*CT*ST;
     G[_MMYY] = g[_MM2A]*ST*ST + g[_MM0P] + g[_MM2B];
     G[_MMYZ] = g[_MM1A]*ST;
     G[_MMZX] = g[_MM1B]*CT;
     G[_MMZY] = g[_MM1B]*ST;
     G[_MMZZ] = g[_MM0Z];
   };

  delete SID->RTwiddle;
  delete SID->WMatrix;
  delete SID->STwiddle;
}

/***************************************************************/
/* switchboard routine *****************************************/
/***************************************************************/
HMatrix *LayeredSubstrate::GetSubstrateDGF(cdouble Omega,
                                           HMatrix *XMatrix,
                                           HMatrix *GMatrix,
                                           DGFMethod Method,
                                           bool AddHomogeneousDGF)
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

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (AddHomogeneousDGF)
   for(int nx=0; nx<XMatrix->NR; nx++)
    { double XDS[6];
      XMatrix->GetEntriesD(nx,":",XDS);
      int nl=GetLayerIndex(XDS[2]);
      if (nl!=GetLayerIndex(XDS[5]))
       continue;
      cdouble EpsRel=EpsLayer[nl], MuRel=MuLayer[nl];
      cdouble k=sqrt(EpsRel*MuRel)*Omega;
      cdouble ZRel=sqrt(MuRel/EpsRel);
      double R[3];
      VecSub(XDS+0, XDS+3, R);
      cdouble G[3][3], C[3][3], dG[3][3][3], dC[3][3][3];
      scuff::CalcGC(R, Omega, EpsRel, MuRel, G, C, dG, dC);
      HMatrix ScriptG(6,6,LHM_COMPLEX,GMatrix->ZM + nx*36);
      for(int P=0; P<2; P++)
       for(int Q=0; Q<2; Q++)
        { cdouble PreFac = II*k;
          if (P==0 && Q==0) PreFac*=ZVAC*ZRel;
          if (P==1 && Q==0) PreFac*=-1.0;
          if (P==1 && Q==1) PreFac/=(ZVAC*ZRel);
          for(int Mu=0; Mu<3; Mu++)
           for(int Nu=0; Nu<3; Nu++)
            ScriptG.AddEntry(3*P+Mu,3*Q+Nu,PreFac*(P==Q ? G[Mu][Nu] : C[Mu][Nu]));
        };
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
                                       cdouble ScriptG[6][6], DGFMethod Method, bool AddHomogeneousDGF)
{ 
  double XBuffer[6];
  HMatrix XMatrix(1,6,LHM_REAL,XBuffer);
  XMatrix.SetEntriesD(0,"0:2",XD);
  XMatrix.SetEntriesD(0,"3:5",XS);

  cdouble GBuffer[36];
  HMatrix GMatrix(36,1,LHM_COMPLEX,GBuffer);

  GetSubstrateDGF(Omega, &XMatrix, &GMatrix, Method, AddHomogeneousDGF);

  for(int Mu=0; Mu<6; Mu++)
   for(int Nu=0; Nu<6; Nu++)
    ScriptG[Mu][Nu] = GMatrix.GetEntry(6*Mu + Nu,0);  
}
