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
 * libSubstrate/Interpolation.cc -- accelerate calculations of
 *                               -- substrate Green's functions
 *                               -- using interpolation tables
 *
 * homer reid             -- 3/2017
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

/******************************************************************/
/* on entry                                                       */
/*  G = G_ij in a coordinate system in which (x,y)=(Rho,0)        */
/* on return:                                                     */
/*  G = G_ij in a coordinate system in which (x,y)=(Rho*CP,Rho*SP)*/
/* p,q label the quadrant: EE, EM, ME, MM                         */
/*******************************************************************/
void LayeredSubstrate::RotateG(cdouble G[6][6], int p, int q, double CP, double SP)
{
  cdouble GP[3][3];
  
  int P=3*p, Q=3*q;

  GP[0][0]=CP*CP*G[P+0][Q+0] - CP*SP*G[P+0][Q+1] - CP*SP*G[P+1][Q+0] + SP*SP*G[P+1][Q+1];
  GP[0][1]=CP*CP*G[P+0][Q+1] + CP*SP*G[P+0][Q+0] - CP*SP*G[P+1][Q+1] - SP*SP*G[P+1][Q+0];
  GP[1][0]=CP*CP*G[P+1][Q+0] + CP*SP*G[P+0][Q+0] - CP*SP*G[P+1][Q+1] - SP*SP*G[P+0][Q+1];
  GP[1][1]=CP*CP*G[P+1][Q+1] + CP*SP*G[P+0][Q+1] + CP*SP*G[P+1][Q+0] + SP*SP*G[P+0][Q+0];

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
/* rotate full 6x6 tensor **************************************/
/***************************************************************/
void LayeredSubstrate::RotateG(cdouble Gij[6][6], double Phi)
{
  double CP=cos(Phi), SP=sin(Phi);
  for(int p=0; p<2; p++)
   for(int q=0; q<2; q++)
    RotateG(Gij, p, q, CP, SP);
}

/***************************************************************/
/* entry point with the proper prototype for passage to the    */
/* Interp1D() initialization routine.                          */
/***************************************************************/
typedef struct fInterpData
 {
   LayeredSubstrate *Substrate;
   double zD, zS;
   cdouble Omega;
 } fInterpData;

void fInterp1DStatic(double Rho, void *UserData, double *fInterp)
{
  fInterpData *fID   = (fInterpData *)UserData;
  LayeredSubstrate *Substrate = fID->Substrate;
  double zD                   = fID->zD;
  double zS                   = fID->zS;

  double qI[3], dqI[3];
  if (Rho==0.0)
   { double DeltaRho=1.0e-5;
     Substrate->GetqIntegral(Rho,          zD, zS, qI);
     Substrate->GetqIntegral(Rho+DeltaRho, zD, zS, dqI);
     VecPlusEquals(dqI, -1.0, qI, 3);
     VecScale(dqI, 1.0/DeltaRho, 3);
   }
  else
   { double DeltaRho=1.0e-5 * fabs(Rho);
     Substrate->GetqIntegral(Rho+DeltaRho, zD, zS, dqI);
     Substrate->GetqIntegral(Rho-DeltaRho, zD, zS, qI);
     VecPlusEquals(dqI, -1.0, qI, 3);
     VecScale(dqI, 1.0/(2.0*DeltaRho), 3);
     Substrate->GetqIntegral(Rho,          zD, zS, qI);
   };
 
  fInterp[0] =  qI[0];
  fInterp[1] = dqI[0];
  fInterp[2] =  qI[1];
  fInterp[3] = dqI[1];
  fInterp[4] =  qI[2];
  fInterp[5] = dqI[2];

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void LayeredSubstrate::InitStaticAccelerator1D(double RhoMin,
                                               double RhoMax,
                                               double z)
{
  if (      I1D 
       &&  (I1DOmega   ==0.0     )
       &&  (I1DRhoMin  <= RhoMin )
       &&  (I1DRhoMax  >= RhoMax )
       &&  (I1DZ       == z      )
     )
   { Log(" reusing existing substrate interpolation table");
     return;
   };

  I1DOmega  = 0.0;
  I1DRhoMin = fmin(RhoMin, I1DRhoMin);
  I1DRhoMax = fmax(RhoMax, I1DRhoMax);
  
  if (I1D) delete I1D;

  Log(" (re)allocating substrate I1D(%g,%g,%g)...",I1DRhoMin,I1DRhoMax,z);

  int NGrid=100;
  char *s=getenv("SCUFF_SUBSTRATE_NGRID");
  if (s)
  sscanf(s,"%i",&NGrid);
   
  struct fInterpData MyfID, *fID=&MyfID;
  fID->Substrate = this;
  fID->zD        = fID->zS=z;
  I1D= new Interp1D(RhoMin, RhoMax, NGrid, 3, fInterp1DStatic,
                    (void *)fID);
  I1DZ      = z;

  Log(" ...done with substrate I1D");

}

/***************************************************************/
/* GEE_{xx}, GEE_{yy}, GEE_{xz}, GEE_{zx}, GEE_{zz}            */
/* GEM_{xy}, GEM_{xz}, GEM_{zx}, GEM_{zy}                      */
/* GMM_{xx}, GMM_{yy}, GMM_{xz}, GMM_{zx}, GMM_{zz}            */
/***************************************************************/
void LayeredSubstrate::InitAccelerator1D(cdouble Omega,
                                         double RhoMin,
                                         double RhoMax,
                                         double z)
{
  /***************************************************************/
  /* number of nonzero components that must be stored per grid   */
  /* point                                                       */
  /***************************************************************/
  int nzFun;
  if (EEOnly && XYOnly)
   nzFun = 2;
  else if (EEOnly)
   nzFun = 5;
  else
   nzFun = 14;
  int nFun = 2*nzFun; // 2 real-valued functions for each z-valued function

  /***************************************************************/
  /* attempt to reuse existing interpolation table ***************/
  /***************************************************************/
  if (      I1D 
       &&  (I1DOmega   == Omega  )
       &&  (I1DRhoMin  <= RhoMin )
       &&  (I1DRhoMax  >= RhoMax )
       &&  (I1DZ       == z      )
       &&  (I1D->nFun  >= nFun   )
     )
   { Log(" reusing existing substrate interpolation table");
     return;
   };

  /***************************************************************/
  /* create new table ********************************************/
  /***************************************************************/
  if (I1D) delete I1D;
  I1D=0;

  I1DOmega  = Omega;
  I1DRhoMin = RhoMin;
  I1DRhoMax = RhoMax;
  I1DZ      = z;

  int NGrid=100;
  char *s=getenv("SCUFF_SUBSTRATE_NGRID");
  if (s)
   sscanf(s,"%i",&NGrid);
   
  Log(" (re)allocating substrate I1D(%g,%g,%g)",I1DRhoMin,I1DRhoMax,z);

  /***************************************************************/
  /* tabulate function at all data points ************************/
  /***************************************************************/
  HMatrix *XMatrix = new HMatrix(NGrid, 6);
  double DeltaRho = (RhoMax - RhoMin) / (NGrid-1);
  for(int nx=0; nx<XMatrix->NR; nx++)
   { XMatrix->SetEntry(nx, 0, RhoMin + nx*DeltaRho);
     XMatrix->SetEntry(nx, 2, z);
     XMatrix->SetEntry(nx, 5, z);
   };
  HMatrix *GMatrix= new HMatrix(36, NGrid, LHM_COMPLEX);

  Log(" computing DGF at %i points...",NGrid);
  GetSubstrateDGF(Omega, XMatrix, GMatrix);
  Log(" ...done with DGF computation");

  /***************************************************************/
  /* put X, Y points into appropriate arrays for passage to      */
  /* Interp1D constructor                                        */
  /***************************************************************/
  double *XPoints  = (double *)XMatrix->GetColumnPointer(0);
  cdouble *YPoints = (cdouble *)mallocEC(NGrid*nzFun*sizeof(cdouble));
  for(int nx=0, ny=0; nx<NGrid; nx++)
   { 
     cdouble *Gij = (cdouble *)GMatrix->GetColumnPointer(nx);

     YPoints[ny++] = Gij[ (0+0)*6 + (0+0) ];
     YPoints[ny++] = Gij[ (0+1)*6 + (0+1) ];
     if (nzFun==2) continue;

     YPoints[ny++] = Gij[ (0+0)*6 + (0+2) ];
     YPoints[ny++] = Gij[ (0+2)*6 + (0+0) ];
     YPoints[ny++] = Gij[ (0+2)*6 + (0+2) ];
     if (nzFun==5) continue;

     YPoints[ny++] = Gij[ (0+0)*6 + (3+1) ];
     YPoints[ny++] = Gij[ (0+0)*6 + (3+2) ];
     YPoints[ny++] = Gij[ (0+2)*6 + (3+0) ];
     YPoints[ny++] = Gij[ (0+2)*6 + (3+1) ];

     YPoints[ny++] = Gij[ (3+0)*6 + (3+0) ];
     YPoints[ny++] = Gij[ (3+1)*6 + (3+1) ];
     YPoints[ny++] = Gij[ (3+0)*6 + (3+2) ];
     YPoints[ny++] = Gij[ (3+2)*6 + (3+0) ];
     YPoints[ny++] = Gij[ (3+2)*6 + (3+2) ];
   };

  I1D= new Interp1D(XPoints, (double *)YPoints, NGrid, 2*nzFun, LMDI_LOGLEVEL_VERBOSE);

  Log(" ...done initializing substrate");

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
bool LayeredSubstrate::GetSubstrateDGF_Interp1D(cdouble Omega, double *XD, double *XS, cdouble Gij[6][6])
{
  /***************************************************************/
  /* sanity checks to determine whether this will work           */
  /***************************************************************/
  if (!I1D) return false;
  if (Omega!=I1DOmega) return false;

  if ( !EqualFloat(XD[2],XS[2]) || !EqualFloat(XD[2],I1DZ) )
   return false;

  double RhoX = XD[0] - XS[0];
  double RhoY = XD[1] - XS[1];
  double Rho  = sqrt(RhoX*RhoX + RhoY*RhoY);

  if (Rho<I1DRhoMin || Rho>I1DRhoMax)
   return false;

  /***************************************************************/
  /* evaluate interpolator and stamp entries into DGF            */
  /***************************************************************/
  cdouble GEntries[14];
  I1D->Evaluate(Rho, (double *)GEntries);

  // set Gij = DGF in rotated coordinate system in which the
  //         = in-plane separation vector is parallel to the x axis
  memset( (cdouble *)Gij, 0, 36*sizeof(cdouble));
  int ne=0;
  Gij[0+0][0+0]     = GEntries[ne++];   // GEE_{xx}
  Gij[0+1][0+1]     = GEntries[ne++];   // GEE_{yy}
  if(I1D->nFun >= 10)
   { Gij[0+0][0+2]  = GEntries[ne++];   // GEE_{xz}
     Gij[0+2][0+0]  = GEntries[ne++];   // GEE_{zx}
     Gij[0+2][0+2]  = GEntries[ne++];   // GEE_{zz}
   };
  if(I1D->nFun == 28)
   { Gij[0+0][3+1]  = GEntries[ne++];   // GEM_{xy}
     Gij[0+0][3+2]  = GEntries[ne++];   // GEM_{xz}
     Gij[0+2][3+0]  = GEntries[ne++];   // GEM_{zx}
     Gij[0+2][3+1]  = GEntries[ne++];   // GEM_{zy}
     Gij[3+0][3+0]  = GEntries[ne++];   // GMM_{xx}
     Gij[3+1][3+1]  = GEntries[ne++];   // GMM_{yy}
     Gij[3+0][3+2]  = GEntries[ne++];   // GMM_{xz}
     Gij[3+2][3+0]  = GEntries[ne++];   // GMM_{zx}
     Gij[3+2][3+2]  = GEntries[ne++];   // GMM_{zz}
   };

  /***************************************************************/
  /* rotate back to original coordinates *************************/
  /***************************************************************/
  RotateG(Gij, atan2(RhoY, RhoX) );
  return true;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void LayeredSubstrate::InitAccelerator3D(cdouble Omega,
                                         double RhoMin, double RhoMax,
                                         double ZDMin, double ZDMax,
                                         double ZSMin, double ZSMax)
{
  /***************************************************************/
  /* number of nonzero components that must be stored per grid   */
  /* point                                                       */
  /***************************************************************/
  int nzFun;
  if (EEOnly && XYOnly)
   nzFun = 2;
  else if (EEOnly)
   nzFun = 5;
  else
   nzFun = 14;
  //int nFun = 2*nzFun; // 2 real-valued functions for each z-valued function

  /***************************************************************/
  /* attempt to reuse existing interpolation table ***************/
  /***************************************************************/
  if (      I3D 
       &&  (I3DOmega == Omega  )
       &&  ( RhoMin<=I3DRhoMin &&  I3DRhoMax<=RhoMax )
       &&  ( ZDMin <=I3DZDMin  &&  I3DZDMax <=ZDMax  )
       &&  ( ZSMin <=I3DZSMin  &&  I3DZSMax <=ZSMax  )
     )
   { Log(" reusing existing substrate interpolation table");
     return;
   };

  /***************************************************************/
  /* create new table ********************************************/
  /***************************************************************/
  if (I3D) delete I3D;
  I3D=0;

  I3DOmega  = Omega;
  I3DRhoMin = RhoMin;
  I3DRhoMax = RhoMax;
  I3DZDMin  = ZDMin;
  I3DZDMax  = ZDMax;
  I3DZSMin  = ZSMin;
  I3DZSMax  = ZSMax;

  int NGrid=10;
  char *s=getenv("SCUFF_SUBSTRATE_NGRID");
  if (s)
   sscanf(s,"%i",&NGrid);
   
  Log(" (re)allocating substrate I3D[ (%g,%g) (%g,%g) (%g,%g) ]",
        RhoMin,RhoMax,ZDMin,ZDMax,ZSMin,ZSMax);
#if 0
  Log(" computing DGF at %i points...",NGrid);
  GetSubstrateDGF(Omega, XMatrix, GMatrix);
  Log(" ...done with DGF computation");

  /***************************************************************/
  /* put X, Y points into appropriate arrays for passage to      */
  /* Interp1D constructor                                        */
  /***************************************************************/
  double *XPoints  = (double *)XMatrix->GetColumnPointer(0);
  cdouble *YPoints = (cdouble *)mallocEC(NGrid*nzFun*sizeof(cdouble));
  for(int nx=0, ny=0; nx<NGrid; nx++)
   { 
     YPoints[ny++] = GMatrix->GetEntry(nx, (0+0)*6 + (0+0) );
     YPoints[ny++] = GMatrix->GetEntry(nx, (0+1)*6 + (0+1) );
     if (nzFun==2) continue;

     YPoints[ny++] = GMatrix->GetEntry(nx, (0+0)*6 + (0+2) );
     YPoints[ny++] = GMatrix->GetEntry(nx, (0+2)*6 + (0+0) );
     YPoints[ny++] = GMatrix->GetEntry(nx, (0+2)*6 + (0+2) );
     if (nzFun==5) continue;

     YPoints[ny++] = GMatrix->GetEntry(nx, (0+0)*6 + (3+1) );
     YPoints[ny++] = GMatrix->GetEntry(nx, (0+0)*6 + (3+2) );
     YPoints[ny++] = GMatrix->GetEntry(nx, (0+2)*6 + (3+0) );
     YPoints[ny++] = GMatrix->GetEntry(nx, (0+2)*6 + (3+1) );

     YPoints[ny++] = GMatrix->GetEntry(nx, (3+0)*6 + (3+0) );
     YPoints[ny++] = GMatrix->GetEntry(nx, (3+1)*6 + (3+1) );
     YPoints[ny++] = GMatrix->GetEntry(nx, (3+0)*6 + (3+2) );
     YPoints[ny++] = GMatrix->GetEntry(nx, (3+2)*6 + (3+0) );
     YPoints[ny++] = GMatrix->GetEntry(nx, (3+2)*6 + (3+2) );
   };

  I1D= new Interp1D(XPoints, (double *)YPoints, NGrid, 2*nzFun, LMDI_LOGLEVEL_VERBOSE);

  Log(" ...done initializing substrate");
#endif

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
bool LayeredSubstrate::GetSubstrateDGF_Interp3D(cdouble Omega, double *XD, double *XS, cdouble Gij[6][6])
{
  /***************************************************************/
  /* sanity checks to determine whether this will work           */
  /***************************************************************/
  if (!I3D) return false;
  if (Omega!=I3DOmega) return false;

  double RhoX = XD[0] - XS[0];
  double RhoY = XD[1] - XS[1];
  double Rho  = sqrt(RhoX*RhoX + RhoY*RhoY);
  double ZD   = XD[2];
  double ZS   = XS[2];

  if ( Rho<I3DRhoMin || Rho>I3DRhoMax )
   return false;
  if ( ZD <I3DZDMin  || ZD >I3DZDMax  )
   return false;
  if ( ZS <I3DZSMin  || ZS >I3DZSMax  )
   return false;

  /***************************************************************/
  /* evaluate interpolator and stamp entries into DGF            */
  /***************************************************************/
  cdouble GEntries[14];
  I3D->Evaluate(Rho, ZD, ZS, (double *)GEntries);

  // set Gij = DGF in rotated coordinate system in which the
  //         = in-plane separation vector is parallel to the x axis
  memset( (cdouble *)Gij, 0, 36*sizeof(cdouble));
  int ne=0;
  Gij[0+0][0+0]     = GEntries[ne++];   // GEE_{xx}
  Gij[0+1][0+1]     = GEntries[ne++];   // GEE_{yy}
  if(I1D->nFun >= 10)
   { Gij[0+0][0+2]  = GEntries[ne++];   // GEE_{xz}
     Gij[0+2][0+0]  = GEntries[ne++];   // GEE_{zx}
     Gij[0+2][0+2]  = GEntries[ne++];   // GEE_{zz}
   };
  if(I1D->nFun == 28)
   { Gij[0+0][3+1]  = GEntries[ne++];   // GEM_{xy}
     Gij[0+0][3+2]  = GEntries[ne++];   // GEM_{xz}
     Gij[0+2][3+0]  = GEntries[ne++];   // GEM_{zx}
     Gij[0+2][3+1]  = GEntries[ne++];   // GEM_{zy}
     Gij[3+0][3+0]  = GEntries[ne++];   // GMM_{xx}
     Gij[3+1][3+1]  = GEntries[ne++];   // GMM_{yy}
     Gij[3+0][3+2]  = GEntries[ne++];   // GMM_{xz}
     Gij[3+2][3+0]  = GEntries[ne++];   // GMM_{zx}
     Gij[3+2][3+2]  = GEntries[ne++];   // GMM_{zz}
   };

  /***************************************************************/
  /* rotate back to original coordinates *************************/
  /***************************************************************/
  RotateG(Gij, atan2(RhoY, RhoX) );
  return true;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
bool LayeredSubstrate::GetSubstrateDGF_Interp(cdouble Omega, HMatrix *XMatrix, HMatrix *GMatrix)
{ 
  if (!XMatrix || !GMatrix)
   ErrExit("%s:%i: internal error",__FILE__,__LINE__);

  int NX=XMatrix->NR;
  if (GMatrix->NR!=36 || GMatrix->NC!=NX)
   ErrExit("%s:%i: internal error (%i,%i)",__FILE__,__LINE__,
            GMatrix->NR, GMatrix->NC);
 
  Log("Attemping to get %i DGF values via interpolation...",NX);
  for(int nx=0; nx<NX; nx++)
   { double XDS[6]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
     if (XMatrix->NC>=6)
      XMatrix->GetEntriesD(nx,"0:5",XDS);
     else
      XMatrix->GetEntriesD(nx,":",XDS);

     cdouble Gij[6][6];
     if ( GetSubstrateDGF_Interp1D(Omega, XDS+0, XDS+3, Gij) )
      GMatrix->SetEntries(":", nx, (cdouble *)Gij); 
     else if ( GetSubstrateDGF_Interp3D(Omega, XDS+0, XDS+3, Gij) )
      GMatrix->SetEntries(":", nx, (cdouble *)Gij); 
     else
      { Log("...failed! ");
        return false;
      };
   };

  Log("...success!");
  return true;

}
