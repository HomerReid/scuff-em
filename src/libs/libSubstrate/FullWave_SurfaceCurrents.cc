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
 * FullWave_SurfaceCurrents.cc -- compute the full-wave green's function
 *                             -- above a multi-layer dielectric substrate
 *                             -- using the method of induced surface currents
 *
 * homer reid                  -- 3/2017-9/2017
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
#ifndef ZVAC
#define ZVAC 376.73031346177
#endif

#define SQRT2 1.41421356237309504880

/***************************************************************/
/* compute the 3x3 submatrices that enter the quadrants of the */
/* Fourier transform of the (6x6) homogeneous DGF              */
/***************************************************************/
void GetGC0Twiddle(cdouble k2, double q2D[2], cdouble qz, double Sign,
                   cdouble GT[3][3], cdouble CT[3][3])
{
  double qx  = q2D[0];
  double qy  = q2D[1];

  GT[0][0] = 1.0 - qx*qx/k2;
  GT[1][1] = 1.0 - qy*qy/k2;
  GT[2][2] = 1.0 - qz*qz/k2;

  GT[0][1] = GT[1][0] =      -1.0*qx*qy/k2;
  GT[0][2] = GT[2][0] = -1.0*Sign*qx*qz/k2;
  GT[1][2] = GT[2][1] = -1.0*Sign*qy*qz/k2;

  CT[0][0] = CT[1][1] = CT[2][2] = 0.0;

  if (qz==0.0)
   CT[0][1]=CT[1][0]=CT[0][2]=CT[2][0]=CT[1][2]=CT[2][1]=0.0;
  else
   { CT[0][1] = -1.0*Sign;         CT[1][0]=-1.0*CT[0][1];
     CT[0][2] = +1.0*Sign*qy/qz;   CT[2][0]=-1.0*CT[0][2];
     CT[1][2] = -1.0*Sign*qx/qz;   CT[2][1]=-1.0*CT[1][2];
   };
}

/***************************************************************/
/* add the (Fourier-space) 6x6 homogeneous DGF *****************/
/***************************************************************/
void AddScriptG0Twiddle(cdouble Omega, cdouble EpsRel, cdouble MuRel, double q2D[2],
                        double zDmS, cdouble ScriptG0Twiddle[6][6])
{
  cdouble k2=EpsRel*MuRel*Omega*Omega;
  cdouble qz = sqrt(k2 - q2D[0]*q2D[0] - q2D[1]*q2D[1]);
  double Sign = (zDmS > 0.0) ? 1.0 : -1.0;
  cdouble ExpFac = exp(II*qz*fabs(zDmS));

/*--------------------------------------------------------------*/
  if (Omega==0.0)
   { cdouble q3D[3];
     q3D[0]=q2D[0]; q3D[1]=q2D[1]; q3D[2]=Sign*qz;
     for(int Mu=0; Mu<3; Mu++)
      for(int Nu=0; Nu<3; Nu++)
       { 
         cdouble Factor = (qz==0.0 ? 0.0 : ExpFac*q3D[Mu]*q3D[Nu]/(2.0*qz));
         ScriptG0Twiddle[0+Mu][0+Nu] += Factor*ZVAC/EpsRel;
         ScriptG0Twiddle[3+Mu][3+Nu] += Factor/(ZVAC*MuRel);
       };
     return;
   };
/*--------------------------------------------------------------*/

  cdouble GT[3][3], CT[3][3];
  GetGC0Twiddle(k2, q2D, qz, Sign, GT, CT);
  cdouble Factor = (qz==0.0 ? 0.0 : -0.5*Omega*ExpFac/qz);
  cdouble EEPrefac = Factor*MuRel*ZVAC;
  cdouble EMPrefac = +0.5*ExpFac;
  cdouble MEPrefac = -0.5*ExpFac;
  cdouble MMPrefac = Factor*EpsRel/ZVAC;
  for(int i=0; i<3; i++)
   for(int j=0; j<3; j++)
    { ScriptG0Twiddle[0+i][0+j] += EEPrefac * GT[i][j];
      ScriptG0Twiddle[0+i][3+j] += EMPrefac * CT[i][j];
      ScriptG0Twiddle[3+i][0+j] += MEPrefac * CT[i][j];
      ScriptG0Twiddle[3+i][3+j] += MMPrefac * GT[i][j];
    };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetScriptG0Twiddle(cdouble Omega, cdouble EpsRel, cdouble MuRel, double q2D[2],
                        double zDmS, cdouble ScriptG0Twiddle[6][6])
{ 
  memset(ScriptG0Twiddle, 0, 36*sizeof(cdouble));
  AddScriptG0Twiddle(Omega, EpsRel, MuRel, q2D, zDmS, ScriptG0Twiddle);
}

/**********************************************************************/
/* assemble the matrix that operates on the vector of external-field  */
/* Fourier coefficients in all regions to yield the vector of         */
/* surface-current Fourier coefficients on all material interfaces    */
/*                                                                    */
/* note: UpdateCachedEpsMu() should be called before this routine! We */
/* don't call it here to allow the possibility of passing Omega=0     */
/* to evaluate W in the DC limit (but using Epsilon, Mu values at     */
/* the frequency in question)                                         */
/**********************************************************************/
void LayeredSubstrate::ComputeW(cdouble Omega, double q2D[2], HMatrix *W)
{
  W->Zero();
  for(int a=0; a<NumInterfaces; a++)
   {
     int RowOffset = 4*a;
     for(int b=a-1; b<=a+1; b++)
      { 
        if (b<0 || b>=NumInterfaces) continue;
        int ColOffset = 4*b;

        // contributions of surface currents on interface z_b
        // to tangential-field matching equations at interface z_a
        cdouble ScriptG0Twiddle[6][6];
        double Sign=-1.0;
        if ( b==(a-1) )
         GetScriptG0Twiddle(Omega, EpsLayer[a], MuLayer[a], q2D, zInterface[a]-zInterface[b], ScriptG0Twiddle);
        else if ( b==(a+1) )
         GetScriptG0Twiddle(Omega, EpsLayer[a+1], MuLayer[a+1], q2D, zInterface[a]-zInterface[b], ScriptG0Twiddle);
        else // (b==a)
         { GetScriptG0Twiddle(Omega, EpsLayer[a],   MuLayer[a],   q2D, +1.0e-12, ScriptG0Twiddle);
           AddScriptG0Twiddle(Omega, EpsLayer[a+1], MuLayer[a+1], q2D, -1.0e-12, ScriptG0Twiddle);
           Sign=1.0;
         };

        for(int EH=0; EH<2; EH++)
         for(int KN=0; KN<2; KN++)
          for(int i=0; i<2; i++)
           for(int j=0; j<2; j++)
            W->AddEntry(RowOffset+2*EH+i, ColOffset+2*KN+j, Sign*ScriptG0Twiddle[3*EH+i][3*KN+j]);
      };
   };
{ static int nw=0;
char FileName[100];
sprintf(FileName,"/tmp/W%i.hdf5",++nw);
W->ExportToHDF5(FileName,"W");
}
  W->LUFactorize();

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void LayeredSubstrate::GetScriptGTwiddle(cdouble Omega, double q2D[2],
                                         double zDest, double zSource,
                                         HMatrix *WMatrix, HMatrix *GTwiddle)
{
  UpdateCachedEpsMu(Omega);

  cdouble OmegaArg = (StaticLimit ? 0.0 : Omega);
  cdouble OmegaFac = ((StaticLimit && Omega!=0.0) ? 1.0/Omega : 1.0);

  /**********************************************************************/
  /**********************************************************************/
  /**********************************************************************/
  if (ForceFreeSpace)
   { 
     cdouble ScriptG0Twiddle[6][6];
     cdouble EpsRel=1.0, MuRel=1.0;
     GetScriptG0Twiddle(OmegaArg, EpsRel, MuRel, q2D, zDest-zSource, ScriptG0Twiddle);
     for(int Mu=0; Mu<6; Mu++)
      for(int Nu=0; Nu<6; Nu++)
       GTwiddle->SetEntry(Mu,Nu,OmegaFac*ScriptG0Twiddle[Mu][Nu]);
     return;
   };
 
  /**********************************************************************/
  /* assemble W matrix ***********************************************/
  /**********************************************************************/
  ComputeW(OmegaArg, q2D,  WMatrix);
  WMatrix->LUInvert();
  
  /**********************************************************************/
  /* assemble RHS vector for each (source point, polarization, orientation) */
  /**********************************************************************/
  int DestRegion  = GetRegionIndex(zDest);
  cdouble EpsDest = EpsLayer[DestRegion];
  cdouble MuDest  = MuLayer[DestRegion];
  int pMin        = (DestRegion==0               ? 1 : 0 );
  int pMax        = (DestRegion==NumInterfaces   ? 0 : 1 );
  cdouble ScriptG0TDest[2][6][6];
  for(int p=pMin, nrDest=DestRegion-1+p; p<=pMax; p++, nrDest++)
   GetScriptG0Twiddle(OmegaArg, EpsDest, MuDest, q2D, zDest - zInterface[nrDest], ScriptG0TDest[p]);

  int SourceRegion  = GetRegionIndex(zSource);
  cdouble EpsSource = EpsLayer[SourceRegion];
  cdouble  MuSource = MuLayer[SourceRegion];
  int qMin          = (SourceRegion==0             ? 1 : 0 );
  int qMax          = (SourceRegion==NumInterfaces ? 0 : 1 );
  cdouble ScriptG0TSource[2][6][6];
  for(int q=qMin, nrSource=SourceRegion-1+q; q<=qMax; q++, nrSource++)
   GetScriptG0Twiddle(OmegaArg, EpsSource, MuSource, q2D, zInterface[nrSource]-zSource, ScriptG0TSource[q]);

  /**********************************************************************/
  /**********************************************************************/
  /**********************************************************************/
  GTwiddle->Zero();
  for(int p=pMin, nrDest=DestRegion-1+p; p<=pMax; p++, nrDest++)
   for(int q=qMin, nrSource=SourceRegion-1+q; q<=qMax; q++, nrSource++)
    { 
      double Sign = (p==q) ? 1.0 : -1.0;
      Sign*=-1.0; // EXPLAIN ME 
      int RowOffset = 4*nrDest, ColOffset = 4*nrSource;
      for(int i=0; i<6; i++)
       for(int j=0; j<6; j++)
        for(int k=0; k<4; k++)
         for(int l=0; l<4; l++)
          { int kk = k>=2 ? k+1 : k;
            int ll = l>=2 ? l+1 : l;
            GTwiddle->AddEntry(i,j, Sign*OmegaFac
                                        *ScriptG0TDest[p][i][kk]
                                        *WMatrix->GetEntry(RowOffset+k, ColOffset+l)
                                        *ScriptG0TSource[q][ll][j]
                              );
         };
    };
}

void LayeredSubstrate::GetScriptGTwiddle(cdouble Omega, double qx, double qy,
                                         double zDest, double zSource,
                                         HMatrix *WMatrix, HMatrix *GTwiddle)
{ double q2D[2];
  q2D[0]=qx;
  q2D[1]=qy;
  GetScriptGTwiddle(Omega, q2D, zDest, zSource, WMatrix, GTwiddle);
}


/***************************************************************/
/***************************************************************/
/***************************************************************/
void LayeredSubstrate::Getg0112(cdouble Omega, double qMag,
                                double zDest, double zSource,
                                HMatrix *WMatrix, HMatrix *g0112[4])
{
  // qx=qMag, qy=0
  HMatrix GT10(6,6,LHM_COMPLEX);
  GetScriptGTwiddle(Omega, qMag,       0.0,        zDest, zSource, WMatrix, &GT10);

  // qx=qMag/sqrt[2], qy=qMag/sqrt[2]
  HMatrix GT11(6,6,LHM_COMPLEX);
  GetScriptGTwiddle(Omega, qMag/SQRT2, qMag/SQRT2, zDest, zSource, WMatrix, &GT11);

  // qx=0, qy=qMag
  HMatrix GT01(6,6,LHM_COMPLEX);
  GetScriptGTwiddle(Omega, 0.0,        qMag,       zDest, zSource, WMatrix, &GT01);

  HMatrix *g0 = g0112[0], *g1x=g0112[1], *g1y=g0112[2], *g2=g0112[3];
  g0->Zero();
  g1x->Zero();
  g1y->Zero();
  g2->Zero();
  
  for(int p=0; p<2; p++)
   for(int q=0; q<2; q++)
    { 
      int ROfs = 3*p;
      int COfs = 3*q;

      g0->SetEntry(ROfs+0,COfs+0,GT01.GetEntry(ROfs+0,COfs+0));
      g0->SetEntry(ROfs+0,COfs+1,GT01.GetEntry(ROfs+0,COfs+1));
      g0->SetEntry(ROfs+1,COfs+0,GT01.GetEntry(ROfs+1,COfs+0));
      g0->SetEntry(ROfs+1,COfs+1,GT10.GetEntry(ROfs+1,COfs+1));

      g0->SetEntry(ROfs+0,COfs+2,GT10.GetEntry(ROfs+0,COfs+2) + GT01.GetEntry(ROfs+0, COfs+2) - SQRT2*GT11.GetEntry(ROfs+0, COfs+2));
      g0->SetEntry(ROfs+2,COfs+0,GT10.GetEntry(ROfs+2,COfs+0) + GT01.GetEntry(ROfs+2, COfs+0) - SQRT2*GT11.GetEntry(ROfs+2, COfs+0));
      g0->SetEntry(ROfs+1,COfs+2,GT10.GetEntry(ROfs+1,COfs+2) + GT01.GetEntry(ROfs+1, COfs+2) - SQRT2*GT11.GetEntry(ROfs+1, COfs+2));
      g0->SetEntry(ROfs+2,COfs+1,GT10.GetEntry(ROfs+2,COfs+1) + GT01.GetEntry(ROfs+2, COfs+1) - SQRT2*GT11.GetEntry(ROfs+2, COfs+1));

      g0->SetEntry(ROfs+2,COfs+2,GT01.GetEntry(ROfs+2,COfs+2));

      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      g1x->SetEntry(ROfs+0,COfs+2, GT10.GetEntry(ROfs+0, COfs+2) - g0->GetEntry(ROfs+0, COfs+2) );
      g1x->SetEntry(ROfs+2,COfs+0, GT10.GetEntry(ROfs+2, COfs+0) - g0->GetEntry(ROfs+2, COfs+0) );
      g1x->SetEntry(ROfs+1,COfs+2, GT10.GetEntry(ROfs+1, COfs+2) - g0->GetEntry(ROfs+1, COfs+2) );
      g1x->SetEntry(ROfs+2,COfs+1, GT10.GetEntry(ROfs+2, COfs+1) - g0->GetEntry(ROfs+2, COfs+1) );

      g1y->SetEntry(ROfs+0,COfs+2, GT01.GetEntry(ROfs+0, COfs+2) - g0->GetEntry(ROfs+0, COfs+2) );
      g1y->SetEntry(ROfs+2,COfs+0, GT01.GetEntry(ROfs+2, COfs+0) - g0->GetEntry(ROfs+2, COfs+0) );
      g1y->SetEntry(ROfs+1,COfs+2, GT01.GetEntry(ROfs+1, COfs+2) - g0->GetEntry(ROfs+1, COfs+2) );
      g1y->SetEntry(ROfs+2,COfs+1, GT01.GetEntry(ROfs+2, COfs+1) - g0->GetEntry(ROfs+2, COfs+1) );

      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      g2->SetEntry(ROfs+0,COfs+0,       GT10.GetEntry(ROfs+0, COfs+0) - g0->GetEntry(ROfs+0, COfs+0));
      g2->SetEntry(ROfs+0,COfs+1, 2.0*( GT11.GetEntry(ROfs+0, COfs+1) - g0->GetEntry(ROfs+0, COfs+1)));
      g2->SetEntry(ROfs+1,COfs+0, 2.0*( GT11.GetEntry(ROfs+1, COfs+0) - g0->GetEntry(ROfs+1, COfs+0)));
      g2->SetEntry(ROfs+1,COfs+1,       GT01.GetEntry(ROfs+1, COfs+1) - g0->GetEntry(ROfs+1, COfs+1));

    };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct IntegrandData 
 {
   LayeredSubstrate *Substrate;
   HMatrix *XMatrix;
   HMatrix *WMatrix;
   cdouble Omega;
   double rkLayer2;
   double qScale;
   FILE *LogFile;
   bool Propagating;
   int nCalls;
   int EntryOnly;
   bool OmitZComponents;

 } IntegrandData;

int SubstrateDGFIntegrand_SC(unsigned ndim, const double *uVector,
                             void *UserData, unsigned fdim, double *fval)
{
  (void )ndim; //unused

  IntegrandData *Data         = (IntegrandData *)UserData;
  Data->nCalls++;

  LayeredSubstrate *Substrate = Data->Substrate;
  HMatrix *XMatrix            = Data->XMatrix;
  HMatrix *WMatrix            = Data->WMatrix;
  cdouble Omega               = Data->Omega;
  double rkLayer2             = Data->rkLayer2;
  FILE *LogFile               = Data->LogFile;
  double qScale               = Data->qScale;
  bool OmitZComponents        = Data->OmitZComponents;

  cdouble *Integrand = (cdouble *)fval;
  memset(Integrand, 0, fdim*sizeof(double));

  double Jacobian = 1.0, u=uVector[0], qMag;
  cdouble qz;
#define UMIN 1.0e-3
  if (Data->Propagating)
   { 
     if ( u<=UMIN )
      u=UMIN;
     if ( u>=(1.0-UMIN) )
      u=1.0-UMIN;
     qz        = u*qScale;
     Jacobian  = qScale*real(qz)/(2.0*M_PI);
     qMag      = real(sqrt(rkLayer2 - qz*qz));
   }
  else // evanescent
   { 
     if (u>1.0-UMIN)
      return 0;    // integrand vanishes at qMag=iqz=infinity
     if (u<=UMIN)
      u=UMIN;
     double Denom = 1.0 / (1.0-u);
     double Kappa = qScale*u*Denom;
     qMag         = sqrt(qScale*qScale + Kappa*Kappa);
     qz           = II*Kappa;
     Jacobian     = qScale*Kappa*Denom*Denom/(2.0*M_PI);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int nx=0, ni=0; nx<XMatrix->NR; nx++)
   { 
     double Rho[2];
     Rho[0]          = XMatrix->GetEntryD(nx,0) - XMatrix->GetEntryD(nx,3);
     Rho[1]          = XMatrix->GetEntryD(nx,1) - XMatrix->GetEntryD(nx,4);
     double RhoMag   = sqrt( Rho[0]*Rho[0] + Rho[1]*Rho[1] );
     double CosTheta = (RhoMag==0.0) ? 0.0 : Rho[0] / RhoMag;
     double SinTheta = (RhoMag==0.0) ? 0.0 : Rho[1] / RhoMag;
     double zDest    = XMatrix->GetEntryD(nx,2);
     double zSource  = XMatrix->GetEntryD(nx,5);

     // fetch scalar quantities
     HMatrix g0(6, 6, LHM_COMPLEX);
     HMatrix g1x(6, 6, LHM_COMPLEX);
     HMatrix g1y(6, 6, LHM_COMPLEX);
     HMatrix g2(6, 6, LHM_COMPLEX);
     HMatrix *g0112[4];
     g0112[0]=&g0;
     g0112[1]=&g1x;
     g0112[2]=&g1y;
     g0112[3]=&g2;
     Substrate->Getg0112(Omega, qMag, zDest, zSource, WMatrix, g0112);

     // fetch bessel functions
     cdouble J[3];
     double qRho = qMag*RhoMag;
     double Workspace[12];
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
if(qRho>10000.0)
 memset(J,0,3*sizeof(cdouble));
else
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
     AmosBessel('J', qRho, 0.0, 3, false, J, Workspace);
     cdouble J1oqRho = (qRho<1.0e-3) ? 0.5 : J[1]/qRho;
     cdouble J0      = J[0], J1[2], J2[2][2];
     J1[0]    = II*CosTheta*J[1];
     J1[1]    = II*SinTheta*J[1];
     J2[0][0] = -1.0*CosTheta*CosTheta*J[2] + J1oqRho;
     J2[0][1] = -1.0*CosTheta*SinTheta*J[2];
     J2[1][0] = -1.0*SinTheta*CosTheta*J[2];
     J2[1][1] = -1.0*SinTheta*SinTheta*J[2] + J1oqRho;
     
     // assemble GTwiddle
     cdouble GTwiddle[6][6];
     for(int p=0; p<2; p++)
      for(int q=0; q<2; q++)
       for(int Mu=0; Mu<3; Mu++)
        for(int Nu=0; Nu<3; Nu++)
         { 
           if (OmitZComponents && (Mu==2 || Nu==2))
            { GTwiddle[3*p+Mu][3*q+Nu]=0.0;
              continue;
            };
            
           cdouble g0MuNu  = g0. GetEntry(3*p+Mu, 3*q+Nu);
           cdouble g1xMuNu = g1x.GetEntry(3*p+Mu, 3*q+Nu);
           cdouble g1yMuNu = g1y.GetEntry(3*p+Mu, 3*q+Nu);
           cdouble g2MuNu  = g2. GetEntry(3*p+Mu, 3*q+Nu);

           GTwiddle[3*p+Mu][3*q+Nu]   = g0MuNu * J0;

           if (Mu<2 && Nu<2)
            GTwiddle[3*p+Mu][3*q+Nu] += g2MuNu * J2[Mu][Nu];
           else if (Mu<2 && Nu==2)
            GTwiddle[3*p+Mu][3*q+Nu] += (g1xMuNu*J1[0] + g1yMuNu*J1[1]);
           else if (Mu==2 && Nu<2)
            GTwiddle[3*p+Mu][3*q+Nu] += (g1xMuNu*J1[0] + g1yMuNu*J1[1]);
         };

/*--------------------------------------------------------------*/
if (Data->EntryOnly!=-1)
 { cdouble *GT = (cdouble *)GTwiddle;
   for(int ne=0; ne<36; ne++)
    if (ne!=Data->EntryOnly)
     GT[ne]=0.0;
 };
/*--------------------------------------------------------------*/

     // stamp into integrand vector
     for(int Mu=0; Mu<6; Mu++)
      for(int Nu=0; Nu<6; Nu++)
       Integrand[ni++] = Jacobian * GTwiddle[Mu][Nu];

     if (LogFile)
      { fprintf(LogFile,"%e %e %e ",qMag,real(qz),imag(qz));
        fprintf(LogFile,"%e %e %e %e ",Rho[0],Rho[1],zDest,zSource);
        fprintf(LogFile,"%e %e ",real(J[0]),imag(J[0]));
        fprintf(LogFile,"%e %e ",real(J1oqRho),imag(J1oqRho));
        fprintf(LogFile,"%e %e ",real(J[2]),imag(J[2]));
        fprintVecCR(LogFile,Integrand+ni-36,36);
        fflush(LogFile);
      };

   }; // for(int nx=0, ni=0; nx<XMatrix->NR; nx++)
  return 0;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void LayeredSubstrate::GetSubstrateDGF_SurfaceCurrent(cdouble Omega,
                                                      HMatrix *XMatrix,
                                                      HMatrix *GMatrix)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int NX = XMatrix->NR;
  int IDim = 36*NX;

  HMatrix *WMatrix  = new HMatrix(4*NumInterfaces, 4*NumInterfaces, LHM_COMPLEX);
  HMatrix *kLayer   = new HMatrix(NumLayers, 2, LHM_REAL);
  cdouble *Integral = new cdouble[IDim];
  cdouble *Error    = new cdouble[IDim];

  UpdateCachedEpsMu(Omega);
 
  // compute wavevectors in all layers and sort by real part
  for(int n=0; n<NumLayers; n++)
   { cdouble k=sqrt(EpsLayer[n]*MuLayer[n])*Omega;
     kLayer->SetEntry(n, 0, real(k));
     kLayer->SetEntry(n, 1, imag(k));
   };
  kLayer->Sort(0);

  double Lower=0.0, Upper=1.0;
  IntegrandData MyData, *Data=&MyData;
  Data->Substrate  = this;
  Data->XMatrix    = XMatrix;
  Data->WMatrix    = WMatrix;
  Data->Omega      = Omega;
  Data->LogFile    = 0;

  Data->EntryOnly       = -1;
  Data->OmitZComponents = false;
  char *s;
  if ( (s=getenv("SCUFF_ENTRY_ONLY")) != 0 )
   sscanf(s, "%i", &(Data->EntryOnly));
  if ( (s=getenv("SCUFF_OMIT_ZCOMPONENTS")) != 0)
   Data->OmitZComponents = (s[0]=='1');

  GMatrix->Zero();

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int n=0; n<NumLayers; n++)
   {
     if (LayerOnly!=-1 && LayerOnly!=n)
      continue;

     double knm1       = (n==0) ? 0.0 : kLayer->GetEntryD(n-1,0);
     double kn         = kLayer->GetEntryD(n,0);
     Data->rkLayer2    = kn*kn;
     Data->qScale      = sqrt(kn*kn - knm1*knm1);
     Data->nCalls      = 0;
     Data->Propagating = true;
     Data->LogFile     = vfopen("/tmp/Propagating%i%s.log","w",n+1,StaticLimit ? "_SL" : "");
     pcubature(2*IDim, SubstrateDGFIntegrand_SC, (void *)Data, 1,
               &Lower, &Upper, qMaxEval, qAbsTol, qRelTol, 
               ERROR_PAIRED, (double *)Integral, (double *)Error);
     if (Data->LogFile) fclose(Data->LogFile);
     Data->LogFile=0;
     Log("Integral %i [%g,%g]: %i calls, GEExx,GEEzz,GMMxx,GMMzz=%s,%s",n+1,knm1,kn,Data->nCalls,
          CD2S(Integral[0]),CD2S(Integral[14]), CD2S(Integral[21]),CD2S(Integral[35]));

     for(int nx=0; nx<NX; nx++)
      for(int ng=0; ng<36; ng++)
       GMatrix->AddEntry(nx, ng, Integral[36*nx + ng]);
   };

  if (LayerOnly==-1 || LayerOnly==NumLayers)
   { Data->nCalls      = 0;
     Data->Propagating = false;
     Data->LogFile     = vfopen("/tmp/Evanescent%s.log","w", StaticLimit ? "_SL" : "");
     Data->qScale      = kLayer->GetEntryD(NumLayers-1,0);
     pcubature(2*IDim, SubstrateDGFIntegrand_SC, (void *)Data, 1,
               &Lower, &Upper, qMaxEval, qAbsTol, qRelTol,
               ERROR_PAIRED, (double *)Integral, (double *)Error);
     Log("Integral %i [%g,inf]: %i calls, GEExx,GEEzz,GMMxx,GMMzz=%s,%s",NumLayers,
          Data->qScale, Data->nCalls,CD2S(Integral[0]),CD2S(Integral[14]), CD2S(Integral[21]),CD2S(Integral[35]));
     if (Data->LogFile) fclose(Data->LogFile);
     Data->LogFile=0;

     for(int nx=0; nx<NX; nx++)
      for(int ng=0; ng<36; ng++)
       GMatrix->AddEntry(nx, ng, Integral[36*nx + ng]);
   };
  
  delete[] Error;
  delete[] Integral;
  delete kLayer;
  delete WMatrix;
  
}
