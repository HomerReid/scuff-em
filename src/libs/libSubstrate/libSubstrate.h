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
 * libSubstrate.h -- header file for C++ library that implements
 *                -- implicit handling of layered material substrates
 *                -- (both electrostatic and full-wave cases)
 *
 * homer reid  -- 3/2017 - 9/2017
 */

#ifndef LIBSUBSTRATE_H
#define LIBSUBSTRATE_H

#include <stdio.h>
#include <math.h>
#include <stdarg.h>

#include <libhrutil.h>
#include <libhmat.h>
#include <libMatProp.h>
#include <libMDInterp.h>
#include <libSGJC.h>

// plane-wave polarizations
#define POL_TE 0
#define POL_TM 1

// methods for full-wave DGF computation
enum DGFMethod {AUTO, FAST_SURFACE_CURRENT, FULL_SURFACE_CURRENT,
                STATIC_LIMIT};

#define LIBSUBSTRATE_NOLOGGING 0
#define LIBSUBSTRATE_TERSE     1
#define LIBSUBSTRATE_VERBOSE   2
#define LIBSUBSTRATE_VERBOSE2  3

#define G0TIME     0
#define BESSELTIME 1
#define WTIME      2
#define SOLVETIME  3
#define STAMPTIME  4
#define NUMTIMES   5
extern const char *TimeNames[];

#ifndef ZVAC
#define ZVAC 376.73031346177
#endif
#define II cdouble(0.0,1.0)
#define SQRT2 1.41421356237309504880

/********************************************************************/
/* labels for gTwiddle functions (of which there are 18) and        */
/* gScalar functions (of which there are 22)                        */
/********************************************************************/
#define _EE0P 0
#define _EE0Z 1
#define _EE1A 2
#define _EE1B 3
#define _EE2A 4

#define _EM0P 5
#define _EM1A 6
#define _EM1B 7
#define _EM2A 8

#define _ME0P 9
#define _ME1A 10
#define _ME1B 11
#define _ME2A 12

#define _MM0P 13
#define _MM0Z 14
#define _MM1A 15
#define _MM1B 16
#define _MM2A 17

#define _EE2B 18
#define _EM2B 19
#define _ME2B 20
#define _MM2B 21

#define NUMGTWIDDLE 18
#define NUMGSCALAR  22

#define _EEXX 0
#define _EEYX 1
#define _EEZX 2
#define _MEXX 3
#define _MEYX 4
#define _MEZX 5

#define _EEXY 6
#define _EEYY 7
#define _EEZY 8
#define _MEXY 9
#define _MEYY 10
#define _MEZY 11

#define _EEXZ 12
#define _EEYZ 13
#define _EEZZ 14
#define _MEXZ 15
#define _MEYZ 16
#define _MEZZ 17

#define _EMXX 18
#define _EMYX 19
#define _EMZX 20
#define _MMXX 21
#define _MMYX 22
#define _MMZX 23

#define _EMXY 24
#define _EMYY 25
#define _EMZY 26
#define _MMXY 27
#define _MMYY 28
#define _MMZY 29

#define _EMXZ 30
#define _EMYZ 31
#define _EMZZ 32
#define _MMXZ 33
#define _MMYZ 34
#define _MMZZ 35

#define _LAMBDA0P 0
#define _LAMBDA0Z 1
#define _LAMBDA1  2
#define _LAMBDA2  3
#define _LAMBDA0X 4
#define _LAMBDA1X 5
#define _LAMBDA2X 6
#define NUMLAMBDA 7
void GetLambdaMatrices(double Theta, double Lambda[7][3][3]);

/***************************************************************/
/* prototype for user-defined Fourier integrand                */
/***************************************************************/
class LayeredSubstrate;
typedef void (*qFunction)(LayeredSubstrate *Substrate,
                          double q2D[2], cdouble Omega,
                          void *UserData, cdouble *Integrand);

/***************************************************************/
/* data structure for layered material substrate               */
/***************************************************************/
class LayeredSubstrate
 {
public:

   // constructor entry point 1: for the case in which we read in the
   // substrate definition from a .Substrate file
   LayeredSubstrate(const char *FileName);

   // constructor entry point 2: for the case in which we read in the
   // substrate definition starting from the second line of a 
   // SUBSTRATE...ENDSUBSTRATE section in an open .scuffgeo file
   LayeredSubstrate(FILE *f, int *pLineNum);

   // actual body of constructor
   // if ErrMsg is nonzero on return, something failed
   void Initialize(FILE *f, const char *FileName, int *pLineNum=0);

   // write a text description of the substrate to f (stdout if f==0)
   void Describe(FILE *f=0);

   // destructor
   ~LayeredSubstrate();

   /*--------------------------------------------------------------*/
   /* electrostatic case: get electrostatic potential and E-field  */
   /* at XD ("Destination") due to monopole at XS ("Source")       */
   /*--------------------------------------------------------------*/
   // substrate contribution
   void GetDeltaPhiE(double XD[3], double XS[3],
                     double PhiE[4], double *pG0Correction=0);

   // total (free-space + substrate)
   void GetTotalPhiE(double XD[3], double XS[3], double PhiE[4]);

   // switchboard between the previous two routines
   void GetPhiE(double XD[3], double XS[3], double PhiE[4],
                bool Total=true);

   /*--------------------------------------------------------------*/
   /* full-wave case: get the 6x6 dyadic Green's function giving   */
   /* the E,H fields at XD due to J, M currents at XS              */
   /*--------------------------------------------------------------*/
   HMatrix *GetSubstrateDGF(cdouble Omega, HMatrix *XMatrix,
                            HMatrix *GMatrix=0, DGFMethod Method=AUTO, bool AddHomogeneousDGF=false);
   HMatrix *GetSubstrateDGF(cdouble Omega, HMatrix *XMatrix, DGFMethod Method);

   void GetSubstrateDGF(cdouble Omega, double XD[3], double XS[3],
                        cdouble ScriptG[6][6], DGFMethod Method=AUTO, bool AddHomogeneousDGF=false);

   // various implementations of the full-wave calculation
   void GetSubstrateDGF_StaticLimit(cdouble Omega,
                                    double *XD, double *XS,
                                    cdouble Gij[6][6]);

   void GetSubstrateDGF_StaticLimit(cdouble Omega, HMatrix *XMatrix,
                                    HMatrix *GMatrix);

   void GetSubstrateDGF_FastSurfaceCurrent(cdouble Omega,
                                           HMatrix *XMatrix,
                                           HMatrix *GMatrix);


   void GetSubstrateDGF_FullSurfaceCurrent(cdouble Omega,
                                           HMatrix *XMatrix,
                                           HMatrix *GMatrix);

   /*--------------------------------------------------------------*/
   /* general-purpose routine for evaluating q (Fourier) integrals */
   /*--------------------------------------------------------------*/
   void GetqIntegral(double RhoMag, double zD, double zS, 
                     double qIntegral[3]);
  
   void qIntegrate(cdouble Omega, qFunction UserFunction,
                   void *UserData, cdouble *Integral, int FDim,
                   bool ThetaIndependent=true);

   /*--------------------------------------------------------------*/
   /*- routines for working with interpolation tables -------------*/
   /*--------------------------------------------------------------*/
   void InitStaticAccelerator1D(double RhoMin, double RhoMax, double z);

   void InitAccelerator1D(cdouble Omega, double RhoMin, double RhoMax, double z);

   void InitAccelerator3D(cdouble Omega,
                          double RhoMin, double RhoMax,
                          double ZDMin, double ZDMax,
                          double ZSMin, double ZSMax);

   bool GetSubstrateDGF_Interp3D(cdouble Omega,
                                 double *XD, double *XS,
                                 cdouble Gij[6][6]);

   bool GetSubstrateDGF_Interp1D(cdouble Omega,
                                 double *XD, double *XS,
                                 cdouble Gij[6][6]);

   bool GetSubstrateDGF_Interp(cdouble Omega, HMatrix *XMatrix,
                               HMatrix *GMatrix);

// private:

// internal ("private") class methods

   void UpdateCachedEpsMu(cdouble Omega);

   // helper functions for static DGF calculation
   double GetStaticG0Correction(double z);
   double GetStaticG0Correction(double zD, double zS);
   void GetSigmaTwiddle(double zS, double q, double *SigmaTwiddle);

   // helper functions for full-wave DGF calculation
   void ComputeW(cdouble Omega, cdouble q[2], HMatrix *W);
   void GetScriptG0Twiddle(cdouble Omega, cdouble q2D[2],
                           double zDest, double zSource,
                           cdouble ScriptG0Twiddle[6][6],
                           int ForceLayer=-1, double ForceSign=0.0,
                           bool Accumulate=false,
                           bool dzDest=false, bool dzSource=false);
   void GetScriptGTwiddle(cdouble Omega, cdouble q2D[2],
                          double zDest, double zSource,
                          HMatrix *RTwiddle, HMatrix *WMatrix,
                          HMatrix *STwiddle, HMatrix *GTwiddle,
                          bool dzDest=false, bool dzSource=false);
   void gTwiddleFromGTwiddle(cdouble Omega, cdouble q,
                             double zDest, double zSource,
                             HMatrix *RTwiddle, HMatrix *WMatrix,
                             HMatrix *STwiddle,
                             cdouble *gTwiddleVD[2][2],
                             bool dzDest=false, bool dzSource=false);
   void gTwiddleHardCoded(cdouble Omega, cdouble q,
                          double zDest, double zSource,
                          cdouble *gTwiddleVD[2][2],
                          bool dzDest=false, bool dzSource=false);

   void RotateG(cdouble Gij[6][6], double Phi);
   void RotateG(cdouble G[6][6], int P, int Q, double CP, double SP);
 
   void GetReflectionCoefficients(double Omega, double *q,
                                  cdouble r[2][2]);
   int GetLayerIndex(double z);

// internal ("private") class data

   char *ErrMsg;

   // info on substrate geometry
   int NumInterfaces;   // number of separating planes
   int NumLayers;       // == NumInterfaces+1
   MatProp **MPLayer;   // MPLayer[n] = properties of layer n
   cdouble  *EpsLayer;  // EpsLayer[n] = permittivity of layer n
   cdouble  *MuLayer;   // MuLayer[n]  = permeability of layer n
   cdouble OmegaCache;  // frequency at which EpsLayer, MuLayer were cached
   double *zInterface;  // z[n] = z-coordinate of layer n lower boundary
   double zGP;          // == z-coordinate of ground plane (or -inf if absent)

   // convergence parameters for q integration
   int qMaxEval;
   double qAbsTol;
   double qRelTol;
   int PPIOrder;
   int PhiEOrder;
   int WhichIntegral;

   // internal storage buffers 
   Interp1D *I1D;
   double I1DRhoMin, I1DRhoMax, I1DZ;
   cdouble I1DOmega;

   Interp3D *I3D;
   double I3DRhoMin, I3DRhoMax;
   double I3DZDMin,  I3DZDMax;
   double I3DZSMin,  I3DZSMax;
   cdouble I3DOmega;

   double Times[NUMTIMES];

   // flags to help in debugging
   int EntryOnly;
   int LayerOnly;
   bool EEOnly;
   bool XYOnly;
   DGFMethod ForceMethod;
   bool ForceFreeSpace;
   bool HardCoded;
   bool StaticLimit;
   int LogLevel;
   bool WritebyqFiles;
   
 };

void AddPhiE0(double XDest[3], double xs, double ys, double zs, double Q, double PhiE[4]);
void AddPhiE0(double XDest[3], double XSource[3], double Q, double PhiE[4]);

LayeredSubstrate *CreateLayeredSubstrate(const char *FileContent);

/***************************************************************/
/* SommerfeldIntegrator.cc                                     */ 
/***************************************************************/
void SommerfeldIntegrate(integrand f, void *fdata, unsigned zfdim,
                         double q0, double qR, int xNu, double Rho,
                         size_t MaxEvalA, size_t MaxEvalB,
                         double AbsTol, double RelTol,
                         cdouble *Integral, cdouble *Error,
                         bool Verbose);

/***************************************************************/
/* SommerfeldIntegrand.cc **************************************/
/***************************************************************/
typedef struct SommerfeldIntegrandData 
 {
   LayeredSubstrate *Substrate;
   cdouble Omega;
   double q0;
   bool uTransform;
   HMatrix *XMatrix;
   HMatrix *RTwiddle;
   HMatrix *WMatrix;
   HMatrix *STwiddle;
   int NumPoints;
   FILE *byqFile;
 } SommerfeldIntegrandData;

int SommerfeldIntegrand(unsigned ndim, const double *x,
                        void *UserData, unsigned fdim, double *fval);

#endif // LIBSUBSTRATE_H
