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
 * libSpherical.h  -- header file for libSpherical, a collection of 
 *                 -- various utilities useful for working with 
 *                 -- spherical coordinates 
 *
 * homer reid      -- 4/2005 -- 2/2010
 */

#ifndef LIBSPHERICAL_H
#define LIBSPHERICAL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex>
#include <cmath>

#include <libhmat.h>

#ifndef cdouble
 typedef std::complex<double> cdouble;
#endif

// values for the 'WaveType' flag passed to functions below
#define LS_REGULAR   0
#define LS_OUTGOING  1
#define LS_INCOMING  2
#define LS_IRREGULAR 3

/***************************************************************/
/* my routines for computing spherical harmonics and other     */
/* quantities store those quantities in the following order:   */
/*  (l,m) = (0,0)                                              */
/*          (1,-1)                                             */
/*          (1,0)                                              */
/*          (1,+1)                                             */
/*          (2,-2)                                             */
/*          ...                                                */
/* and this macro returns the index in the list for a given    */
/* pair (l,m).                                                 */
/***************************************************************/
#define LM2ALPHA(l,m) ( (l)*(l+1) + (m) )

/***************************************************************/
/* conversion routines                                         */
/***************************************************************/
void CoordinateC2S(double X, double Y, double Z, double *r, double *Theta, double *Phi);
void CoordinateC2S(const double X[3], double *r, double *Theta, double *Phi);
void CoordinateC2S(const double X[3], double R[3]);
void CoordinateC2S(const double X[3]);
void CoordinateS2C(double r, double Theta, double Phi, double X[3]);
void CoordinateS2C(double r, double Theta, double Phi, double *X, double *Y, double *Z);
void CoordinateS2C(const double R[3], double X[3]);
void CoordinateS2C(double R[3]);
void VectorC2S(double Theta, double Phi, const cdouble VC[3], cdouble VS[3]);
void VectorC2S(double Theta, double Phi, const double VC[3], double VS[3]);
void VectorC2S(double Theta, double Phi, double V[3]);
void VectorC2S(double Theta, double Phi, cdouble V[3]);
void VectorS2C(double Theta, double Phi, const cdouble VS[3], cdouble VC[3]);
void VectorS2C(double Theta, double Phi, const double VS[3], double VC[3]);
void VectorS2C(double Theta, double Phi, double V[3]);
void VectorS2C(double Theta, double Phi, cdouble V[3]);

/***************************************************************/
/* finite-difference vector-calculus operators                 */
/***************************************************************/
void GetDiv(double r, double Theta, double Phi,
            void (*F)(double r, double Theta, double Phi, cdouble *V),
            cdouble Terms[3], cdouble *DivF);

cdouble GetDiv(double r, double Theta, double Phi,
               void (*F)(double r, double Theta, double Phi, cdouble *V));

void GetCurl(double r, double Theta, double Phi,
             void (*F)(double r, double Theta, double Phi, cdouble *V),
             cdouble *CurlF, cdouble *CurlF2);

/***************************************************************/
/* spherical harmonics                                         */
/***************************************************************/
void GetYlmArray(int lMax, double Theta, double Phi, cdouble *Ylm);
cdouble GetYlm(int l, int m, double Theta, double Phi);
void GetRealYlmArray(int lMax, double Theta, double Phi, double *Ylm);
double GetRealYlm(int l, int m, double Theta, double Phi);
void GetYlmDerivArray(int lMax, double Theta, double Phi, 
                      cdouble *Ylm, cdouble *dYlmdTheta);

/***************************************************************/ 
/* spherical bessel functions needed for scalar and vector     */
/* helmholtz solutions                                         */
/*                                                             */
/* note: R and dRdr must have space for (lMax+2) cdoubles even */
/* though only (lMax+1) functions are calculated.              */
/*                                                             */
/* note: there is no typo here! R and dRdr must have space for */
/* lMax+2 slots, not lMax+1 slots, even though the user is only*/
/* asking for lMax+1 return quantities, because we need one    */
/* extra slot to compute derivatives via recurrence relations. */
/***************************************************************/
void GetRadialFunctions(int lMax, cdouble k, double r, int WaveType,
                        cdouble *R, cdouble *dRdr, double *Workspace=0);

void GetRadialFunction(int l, cdouble k, double r, int WaveType,
                       cdouble *R, cdouble *dRdr=0, cdouble *RlSlash=0);

/***************************************************************/
/* scalar helmholtz solutions **********************************/
/***************************************************************/
void GetScalarHelmholtzSolutions(int lMax, cdouble k,
                                 double r, double Theta, double Phi,
                                 int WaveType, cdouble *Psi);

/***************************************************************/
/* vector helmholtz solutions: spherical waves / Maxwell waves */
/***************************************************************/
int GetMWIndex(int L, int M, int T); // maxwell-wave index
int GetSWIndex(int L, int M, int P); // spherical-wave index
int GetRFIndex(int L, int P);        // radial-function index

void GetAngularFunctionArray(int LMax, double Theta, double Phi,
                             cdouble *Workspace, HMatrix *AAMatrix);

void GetAngularFunctionArray(int LMax, double X[3],
                             cdouble *Workspace, HMatrix *AAMatrix);

void GetVSWRadialFunctions(int LMax, cdouble k, double r,
                           int WaveType, cdouble *RFArray,
                           double *Workspace=0,
                           bool TimesrFactor=false,
                           bool Conjugate=false);

HMatrix *GetWaveMatrix(double r, double Theta, double Phi,
                       cdouble k, int LMax, int WaveType,
                       HMatrix *WaveMatrix=0, bool MaxwellWaves=true,
                       cdouble *Workspace=0, bool RConjugate=false);

HMatrix *GetWaveMatrix(double X[3],
                       cdouble k, int LMax, int WaveType,
                       HMatrix *WaveMatrix=0, bool MaxwellWaves=true,
                       cdouble *Workspace=0, bool RConjugate=0);

HMatrix *GetMWMatrix(double r, double Theta, double Phi,
                     cdouble k, int LMax, int WaveType,
                     HMatrix *WaveMatrix=0, cdouble *Workspace=0,
                     bool RConjugate=false);

HMatrix *GetMWMatrix(double X[3], cdouble k, int LMax, int WaveType,
                     HMatrix *WaveMatrix=0, cdouble *Workspace=0, 
                     bool RConjugate=false);

double GetdzVSWCoefficient(int L,  int M,  int T,
                           int LP, int MP, int TP);

void GetTranslationMatrices(double Xij[3], cdouble k, int lMax, 
                            HMatrix *A, HMatrix *B, HMatrix *C);

/***************************************************************/
/* stuff below this line is legacy, replaced by stuff in new   */
/* file VectorSphericalWaves.cc, and should eventually be      */
/* removed                                                     */
/***************************************************************/

/***************************************************************/
/* vector spherical harmonics                                  */
/***************************************************************/
void GetXlm(int l, int m, double Theta, double Phi, cdouble X[3]);
void GetXlmArray(int lMax, double Theta, double Phi, cdouble *X);

/***************************************************************/
/* vector helmholtz solutions                                  */
/***************************************************************/
void GetMNlm(int l, int m, cdouble k, double r, double Theta, double Phi, 
             int WaveType, cdouble M[3], cdouble N[3]);

void GetMNlmArray(int lMax, cdouble k,
                  double r, double Theta, double Phi,
                  int WaveType, cdouble *M, cdouble *N, 
                  double *Workspace=0, cdouble *LL=0, cdouble *DivLL=0);

/***************************************************************/
/***************************************************************/
/* differential operators **************************************/
/***************************************************************/
/***************************************************************/

/***************************************************************/
/* bessel and airy functions in AmosBessel.cc                  */
/***************************************************************/
int AmosBessel(char WhichFunction, cdouble z,
               double MinOrder, int NumOrders,
               bool Scale, cdouble *f, double *Workspace=0);

int AmosAiry(char WhichFunction, cdouble z, bool Scale, cdouble *f);


#endif
