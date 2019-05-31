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
 * libSpherical.cc -- a collection of various utilities useful for 
 *                    working with spherical coordinates 
 *
 * homer reid      -- 4/2005 -- 2/2010
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>

#include <config.h>

#ifdef HAVE_LIBGSL
 #include <gsl/gsl_sf_legendre.h>
#endif

#include "libSpherical.h"

#define II cdouble(0.0,1.0)
#define ROOT2 1.41421356237309504880

/***************************************************************/
/* utility routines for indices of spherical-wave / Maxwell-   */
/* wave functions and radial functions                         */
/***************************************************************/
int GetMWIndex(int L, int M, int T) // maxwell-wave index
{
  int Alpha = LM2ALPHA(L,M) - 1;
  return 2*Alpha + T;
}

int GetSWIndex(int L, int M, int P) // spherical-wave index
{
  return 3*LM2ALPHA(L, M) + P;
}

int GetRFIndex(int L, int P)       // radial-function index
{ return 3*L + P; }



/***************************************************************/
/* AAMatrix[Mu,nsw] = Mu component of angular function #nsw    */
/* sizeof(Workspace) >= 2*NAlpha                               */
/***************************************************************/
void GetAngularFunctionArray(int LMax, double Theta, double Phi,
                             cdouble *Workspace, HMatrix *AAMatrix)
{
  int NAlpha=(LMax+1)*(LMax+1);

  cdouble *YArray         = Workspace;
  cdouble *dYdThetaArray  = Workspace + NAlpha;

  if ( fabs(Theta)<1.0e-6) Theta=1.0e-6;
  if ( fabs(M_PI-Theta)<1.0e-6) Theta=M_PI-1.0e-6;
  GetYlmDerivArray(LMax, Theta, Phi, YArray, dYdThetaArray);

  AAMatrix->Zero();
  for(int L=1, Alpha=1; L<=LMax; L++)
   for(int M=-L; M<=L; M++, Alpha++)
    { 
      cdouble X[3]={0.0, 0.0, 0.0};
      cdouble Z[3]={0.0, 0.0, 0.0};
      cdouble Y[3]={0.0, 0.0, 0.0};

      double RtLLP1=sqrt(L*(L+1.0));

      X[1] = -1.0*M*YArray[Alpha] / ( RtLLP1*sin(Theta) );
      X[2] = -1.0*II*dYdThetaArray[Alpha] / RtLLP1;

      Z[1] = -1.0*II*X[2];
      Z[2] = +1.0*II*X[1];

      Y[0] = YArray[Alpha];

      for(int Mu=0; Mu<3; Mu++)
       { AAMatrix->SetEntry(Mu, 3*Alpha+0, X[Mu]);
         AAMatrix->SetEntry(Mu, 3*Alpha+1, Z[Mu]);
         AAMatrix->SetEntry(Mu, 3*Alpha+2, Y[Mu]);
       };
    };

}

/***************************************************************/
/* like the previous routine but returns cartesian components  */
/***************************************************************/
void GetAngularFunctionArray(int LMax, double X[3],
                             cdouble *Workspace, HMatrix *AAMatrix)
{
  double r, Theta, Phi;
  CoordinateC2S(X, &r, &Theta, &Phi);
  GetAngularFunctionArray(LMax, Theta, Phi, Workspace, AAMatrix);

  for(int nc=0; nc<AAMatrix->NC; AAMatrix++)
   VectorS2C(Theta, Phi, (cdouble *)AAMatrix->GetColumnPointer(nc));
}

/***************************************************************/
/* Compute the radial-function factors in the vector spherical */
/* wavefunctions.                                              */
/* RFArray[3*L + P] = radial-function factor that multiplies   */
/*                    angular function of polarization P for   */
/*                    spherical-wave index L                   */
/*                    where P={0,1,2} for angular functions    */
/*                            {X,Z,Y}                          */
/* RFArray[3*L + 0] = R(kr)                                    */
/* RFArray[3*L + 1] = R(kr)/kr + |dR(z)/dz|_{z=kr}             */
/* RFArray[3*L + 2] = -Sqrt[L*(L+1)]*R(kr)/kr                  */
/*                                                             */
/* where R(kr) = j_\ell(kr)     for WaveType=LS_REGULAR        */
/*               y_\ell(kr)     for WaveType=LS_IRREGULAR      */
/*               h^(1)_\ell(kr) for WaveType=LS_OUTGOING       */
/*               h^(2)_\ell(kr) for WaveType=LS_INCOMING       */
/*                                                             */
/* if TimesrFactor=true, all functions are multiplied by a     */
/* factor of r (if r=0 the limiting value as r->0 is returned).*/
/*                                                             */
/* if Conjugate=true, the radial factor in all functions is    */
/* complex-conjugated.                                         */
/***************************************************************/
void GetVSWRadialFunctions(int LMax, cdouble k, double r,
                           int WaveType, cdouble *RFArray,
                           double *Workspace, bool TimesrFactor,
                           bool Conjugate)
{
  if (r==0.0)
   { memset(RFArray, 0, 3*(LMax+1));
     if (TimesrFactor) return;
     RFArray[3*1+1] = 2.0/3.0;
     RFArray[3*1+2] = -M_SQRT2/3.0; //-1.41421356237309504880/3.0;
     return;
   };

  cdouble z=k*r;

  char Func =   (WaveType==LS_REGULAR)   ? 'j'
              : (WaveType==LS_OUTGOING)  ? 'o'
              : (WaveType==LS_IRREGULAR) ? 'y'
              : (WaveType==LS_INCOMING)  ? 't' : 'x'; 

  if (Func=='x') 
   ErrExit("unknown WaveType=%i in GetVSWRadialFunctions",WaveType);

  AmosBessel(Func,z,0.0,LMax+2,false,RFArray,Workspace);

  for(int L=LMax; L>=0; L--)
   { 
     double RtLLP1 = sqrt(L*(L+1.0));

     cdouble R      = RFArray[L];
     cdouble Rlp1   = RFArray[L+1];
     cdouble Roz    = R/z;
     cdouble dRdz   = ((double )L)*Roz - Rlp1;
     RFArray[3*L+0] = R;
     RFArray[3*L+1] = Roz + dRdz;
     RFArray[3*L+2] = -RtLLP1*Roz;
   };
 
  if (TimesrFactor)
   VecScale(RFArray, r, 3*(LMax+1) );

  if (Conjugate)
   for(int n=0; n<3*(LMax+1); n++)
    RFArray[n]=conj(RFArray[n]);

}

/***************************************************************/
/* Get the full set of spherical waves (MaxwellWaves==false)   */
/* or Maxwell waves (MaxwellWaves==true).                      */
/*                                                             */
/* In the former case,                                         */
/* WaveMatrix[Mu,nsw] = Muth spherical component of nswth      */
/* spherical wave at (r,Theta,Phi).                            */
/*                                                             */
/* In the latter case,                                         */
/* WaveMatrix[Mu,nmw] = Muth spherical component of nmwth      */
/* Maxwell wave at (r,Theta,Phi).                              */
/*                                                             */
/* If WaveMatrix is non-null, it should have dimensions 3xNW   */
/* where NW = 2*( (LMax+1)^2 - 1) (Maxwell waves)              */
/*            3*( (LMax+1)^2 - 1) (spherical waves)            */
/*                                                             */
/* If RConjugate==true, the radial factor is complex-conjugated.*/
/*                                                             */
/* If Workspace is non-null, it should satisfy                 */
/* sizeof(Workspace) >= (7*LMax+7)*(LMax+1) * sizeof(cdouble)  */
/***************************************************************/
HMatrix *GetWaveMatrix(double r, double Theta, double Phi,
                       cdouble k, int LMax, int WaveType,
                       HMatrix *WaveMatrix, bool MaxwellWaves,
                       cdouble *Workspace, bool RConjugate)
{
  int NAlpha = (LMax+1)*(LMax+1);
  int NMW    = 2*(NAlpha-1);
  int NSW    = 3*(NAlpha-1);
  int NW     = MaxwellWaves ? NMW : NSW;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (WaveMatrix && (WaveMatrix->NR!=3 || WaveMatrix->NC!=NW))
   { Warn("wrong-size WaveMatrix in GetWaveMatrix (reallocating...)");
     delete WaveMatrix;
     WaveMatrix=0;
   };
  if (WaveMatrix==0)
   WaveMatrix = new HMatrix(3, NW, LHM_COMPLEX);

  /***************************************************************/
  /* fetch radial-function factors  ******************************/
  /* RFArray[3*L+P] = radial-function factor for polarization P  */
  /***************************************************************/
  bool OwnsWorkspace = false;
  if (Workspace==0)
   { Workspace = (cdouble *)mallocEC( (7*LMax+7)*(LMax+1)*sizeof(cdouble));
     OwnsWorkspace=true;
   };

  cdouble *YArray        = Workspace;
  cdouble *dYdThetaArray = YArray        + NAlpha;
  cdouble *RFArray       = dYdThetaArray + NAlpha;
  cdouble *dWorkspace    = RFArray       + 3*(LMax+1);

  GetVSWRadialFunctions(LMax, k, r, WaveType, RFArray, (double *)dWorkspace, /* TimesrFactor */ false, RConjugate);

  /***************************************************************/
  /* fetch angular functions                                     */
  /***************************************************************/
  cdouble AA[3][3]={ {0.0, 0.0, 0.0},
                     {0.0, 0.0, 0.0},
                     {0.0, 0.0, 0.0} 
                   };
  if ( fabs(Theta)<1.0e-6) Theta=1.0e-6;
  if ( fabs(M_PI-Theta)<1.0e-6) Theta=M_PI-1.0e-6;
  GetYlmDerivArray(LMax, Theta, Phi, YArray, dYdThetaArray);

  WaveMatrix->Zero();
  for(int L=1, Alpha=1; L<=LMax; L++)
   {
     double RtLLP1=sqrt(L*(L+1.0));

     for(int M=-L; M<=L; M++, Alpha++)
      { 
        /* AA[P][Mu] = Mu comp. of angular function for pol P */
        AA[0][1] = -1.0*M*YArray[Alpha] / ( RtLLP1*sin(Theta) );
        AA[0][2] = -1.0*II*dYdThetaArray[Alpha] / RtLLP1;
        AA[1][1] = -1.0*II*AA[0][2];
        AA[1][2] = +1.0*II*AA[0][1];
        AA[2][0] = YArray[Alpha];

        for(int P=0; P<3; P++)
         { int Index =    MaxwellWaves
                       ?  GetMWIndex(L, M, (P==0) ? 0 : 1)
                       :  GetSWIndex(L, M, P);
           for(int Mu=0; Mu<3; Mu++)
             WaveMatrix->AddEntry(Mu,Index,RFArray[3*L+P]*AA[P][Mu]);
         };
      };
   };
  if (OwnsWorkspace) 
   free(Workspace);

  return WaveMatrix;
}

HMatrix *GetWaveMatrix(double X[3],
                       cdouble k, int LMax, int WaveType,
                       HMatrix *WaveMatrix, bool MaxwellWaves,
                       cdouble *Workspace, bool RConjugate)
{
  double r, Theta, Phi;
  CoordinateC2S(X, &r, &Theta, &Phi);
  HMatrix *WM=GetWaveMatrix(r, Theta, Phi, k, LMax, WaveType, 
                            WaveMatrix, MaxwellWaves, Workspace,
                            RConjugate);
  for(int nc=0; nc<WM->NC; nc++)
   VectorS2C(Theta, Phi, (cdouble *)WM->GetColumnPointer(nc));
  return WM;
}

HMatrix *GetMWMatrix(double r, double Theta, double Phi,
                     cdouble k, int LMax, int WaveType,
                     HMatrix *WaveMatrix, cdouble *Workspace,
                     bool RConjugate)
{ 
  return GetWaveMatrix(r, Theta, Phi, k, LMax, WaveType,
                       WaveMatrix, true, Workspace, RConjugate);
}

HMatrix *GetMWMatrix(double X[3], cdouble k, int LMax, int WaveType,
                     HMatrix *WaveMatrix, cdouble *Workspace, 
                     bool RConjugate)
{ 
  return GetWaveMatrix(X, k, LMax, WaveType,
                       WaveMatrix, true, Workspace, RConjugate);
}

/***************************************************************/
/* Get the coefficients in the expansion of the z-derivative   */
/* of a vector Maxwell wave as a linear combination of vector  */
/* Maxwell waves                                               */
/* d/dz F_{LMT} = k*\sum C_{LMT, L'M'T'} F_{L'M'T'}            */
/* where T=0,1 for M,N waves                                   */
/***************************************************************/
static double almCoefficient(int L, int M)
 { return ((double)M)/(L*(L+1.0)); }

static double blmCoefficient(int L, int M)
{ 
  double Num   = L*(L+2.0)*(L-M+1.0)*(L+M+1.0);
  double Denom = (2.0*L+1.0)*(2.0*L+3.0);
  return sqrt(Num/Denom) / (L+1.0);
}

double GetdzVSWCoefficient(int L,  int M,  int T,
                           int LP, int MP, int TP)
{
  if (M!=MP)
   return 0.0;

  if (L==LP && T!=TP)
   return ( T==0 ? 1.0 : -1.0) * almCoefficient(L,M);

  int DL = L-LP;
  if ( T!=TP || abs(DL)!=1 )
   return 0.0;

  if (L>LP)
   return blmCoefficient(LP,M);
  return -1.0*blmCoefficient(L,M);
}
