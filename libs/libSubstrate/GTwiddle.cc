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
 * GTwiddle.cc -- compute the Fourier-space dyadic Green's function
 *             -- function for a layered substrate using the method
 *             -- of induced surface currents
 *
 * homer reid  -- 3/2017-9/2017
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
/* compute the 3x3 submatrices that enter the quadrants of the */
/* Fourier transform of the 6x6 DGF for a homogeneous medium.  */
/***************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
int Gamma0DSTerm=0;
int Gamma0Term=0;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void GetGCTwiddle(cdouble k2, cdouble q2D[2], cdouble qz, double Sign,
                  cdouble GT[3][3], cdouble CT[3][3])
{
  cdouble qx  = q2D[0];
  cdouble qy  = q2D[1];

  if (k2==0.0)
   { cdouble q3D[3];
     q3D[0]=q2D[0]; q3D[1]=q2D[1]; q3D[2]=Sign*qz;
     for(int i=0; i<3; i++)
      for(int j=0; j<3; j++)
       { GT[i][j] = -1.0*q3D[i]*q3D[j];
         CT[i][j] = 0.0;
       }
     return;
   }

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
   { CT[0][1] = -1.0*Sign;    CT[1][0]=-1.0*CT[0][1];
     CT[0][2] = +1.0*qy/qz;   CT[2][0]=-1.0*CT[0][2];
     CT[1][2] = -1.0*qx/qz;   CT[2][1]=-1.0*CT[1][2];
   };

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
if (Gamma0Term==0) return;
memset( (cdouble *)GT, 0, 9*sizeof(cdouble));
memset( (cdouble *)CT, 0, 9*sizeof(cdouble));
if (Gamma0Term==1)
 { 
   GT[0][0]=GT[1][1]=GT[2][2]=1.0;
 }
else if (Gamma0Term==2)
 { GT[0][0] = - qx*qx/k2;
   GT[1][1] = - qy*qy/k2;
   GT[2][2] = - qz*qz/k2;
   GT[0][1] = GT[1][0] =      -1.0*qx*qy/k2;
   GT[0][2] = GT[2][0] = -1.0*Sign*qx*qz/k2;
   GT[1][2] = GT[2][1] = -1.0*Sign*qy*qz/k2;
 }
else if (Gamma0Term==3)
 { CT[0][1] = -1.0*Sign;    CT[1][0]=-1.0*CT[0][1];
   CT[0][2] = +1.0*qy/qz;   CT[2][0]=-1.0*CT[0][2];
   CT[1][2] = -1.0*qx/qz;   CT[2][1]=-1.0*CT[1][2];
 };
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
}

/***************************************************************/
/* Get the (Fourier-space) 6x6 homogeneous DGF for the         */
/* medium occupying substrate layer #nl.                       */
/* If Omega==0, return instead the quantity                    */
/*  \lim_{\omega \to 0 }  -i\omega * Gamma(\omega)             */
/***************************************************************/
void LayeredSubstrate::GetGamma0Twiddle(cdouble Omega, cdouble q2D[2],
                                        double zDest, double zSource,
                                        cdouble Gamma0Twiddle[6][6],
                                        int ForceLayer, double Sign,
                                        bool Accumulate,
                                        bool dzDest, bool dzSource)
{
  if (!Accumulate)
   memset( (cdouble *) Gamma0Twiddle, 0, 36*sizeof(cdouble) );

  // autodetect region unless user forced it
  int nl = ForceLayer;
  if (nl==-1)
   { nl = GetLayerIndex(zDest);
     if (nl!=GetLayerIndex(zSource) ) 
      Warn("%s:%i: internal inconsistency",__FILE__,__LINE__);
   };

  // autodetect sign if not specified 
  if (Sign==0.0)
   Sign = ( (zDest >= zSource) ? 1.0 : -1.0 );

  cdouble EpsRel=EpsLayer[nl], MuRel=MuLayer[nl];
  cdouble k2 = EpsRel*MuRel*Omega*Omega, q2=q2D[0]*q2D[0] + q2D[1]*q2D[1];
  cdouble qz = sqrt(k2 - q2);
  if (imag(qz)<0) qz*=-1.0;
  if ( real(Omega)<0.0 ) // negative Omega means compute quasi-static limit:
   { Omega*=-1.0;        //  retain Omega prefactors, but set Omega->0 when
     qz = II*sqrt(q2);   //  computing qz.
   };

  cdouble GT[3][3], CT[3][3];
  cdouble ExpFac = exp(II*qz*fabs(zDest-zSource));
  if (dzDest)
   ExpFac *= +1.0*Sign*II*qz;
  if (dzSource)
   ExpFac *= -1.0*Sign*II*qz;
  GetGCTwiddle(k2, q2D, qz, Sign, GT, CT);

  cdouble EEPreFac, EMPreFac, MEPreFac, MMPreFac;
  if (qz==0.0)
   { EEPreFac = 0.0;
     EMPreFac = 0.5;
     MEPreFac = -0.5;
     MMPreFac = 0.0;
   }
  else if (Omega==0.0)
   { EEPreFac = 0.5*II*ZVAC/(EpsRel*qz);
     EMPreFac = 0.0;
     MEPreFac = 0.0;
     MMPreFac = 0.5*II/(ZVAC*MuRel*qz);
   }
  else
   { EEPreFac = -0.5*Omega*ZVAC*MuRel/qz;
     EMPreFac = +0.5;
     MEPreFac = -0.5;
     MMPreFac = -0.5*Omega*EpsRel/(ZVAC*qz);
   }

  for(int i=0; i<3; i++)
   for(int j=0; j<3; j++)
    { Gamma0Twiddle[0+i][0+j] += EEPreFac * ExpFac * GT[i][j];
      Gamma0Twiddle[0+i][3+j] += EMPreFac * ExpFac * CT[i][j];
      Gamma0Twiddle[3+i][0+j] += MEPreFac * ExpFac * CT[i][j];
      Gamma0Twiddle[3+i][3+j] += MMPreFac * ExpFac * GT[i][j];
    };
  if ( nl<NumInterfaces || std::isinf(zGP) ) return;

  /***************************************************************/
  /* add image contribution if we are in the bottommost layer    */
  /* and a ground plane is present                               */
  /***************************************************************/
  // only need to refetch if Sign was -1 before
  if (Sign<0.0) GetGCTwiddle(k2, q2D, qz, +1.0, GT, CT);
  ExpFac = exp(II*qz*fabs(zDest + zSource - 2.0*zGP));
  if (dzDest)
   ExpFac *= +1.0*II*qz;
  if (dzSource)
   ExpFac *= +1.0*II*qz;
  const static double ImageSign[3] = {-1.0, -1.0, +1.0};
  for(int i=0; i<3; i++)
   for(int j=0; j<3; j++)
    { Gamma0Twiddle[0+i][0+j] += ImageSign[j] * EEPreFac * ExpFac * GT[i][j];
      Gamma0Twiddle[0+i][3+j] -= ImageSign[j] * EMPreFac * ExpFac * CT[i][j];
      // EXPLAIN ME I don't understand why the following two signs aren't flipped.
      Gamma0Twiddle[3+i][0+j] += ImageSign[j] * MEPreFac * ExpFac * CT[i][j];
      Gamma0Twiddle[3+i][3+j] -= ImageSign[j] * MMPreFac * ExpFac * GT[i][j];
    };
}

/**********************************************************************/
/* assemble the matrix that operates on the vector of external-field  */
/* Fourier coefficients in all regions to yield the vector of         */
/* surface-current Fourier coefficients on all material interfaces    */
/**********************************************************************/
void LayeredSubstrate::ComputeW(cdouble Omega, cdouble q2D[2], HMatrix *W)
{
  UpdateCachedEpsMu(Omega);

  W->Zero();
  for(int a=0; a<NumInterfaces; a++)
   {
     int RowOffset = 4*a;
     double za     = zInterface[a];

     for(int b=a-1; b<=a+1; b++)
      { 
        if (b<0 || b>=NumInterfaces) continue;

        int ColOffset = 4*b;
        double zb     = zInterface[b];

        // contributions of surface currents on interface z_b
        // to tangential-field matching equations at interface z_a
        cdouble Gamma0Twiddle[6][6];
        double Sign=-1.0;
        if ( b==(a-1) )
         GetGamma0Twiddle(Omega, q2D, za, zb, Gamma0Twiddle, a);
        else if ( b==(a+1) )
         GetGamma0Twiddle(Omega, q2D, za, zb, Gamma0Twiddle, b);
        else // (b==a)
         { GetGamma0Twiddle(Omega, q2D, za, za, Gamma0Twiddle, a,   +1.0);
           GetGamma0Twiddle(Omega, q2D, za, za, Gamma0Twiddle, a+1, -1.0, true);
           Sign=1.0;
         };

        for(int EH=0; EH<2; EH++)
         for(int KN=0; KN<2; KN++)
          for(int i=0; i<2; i++)
           for(int j=0; j<2; j++)
            W->AddEntry(RowOffset+2*EH+i, ColOffset+2*KN+j, Sign*Gamma0Twiddle[3*EH+i][3*KN+j]);
      }
   }
  if (LogLevel >= LIBSUBSTRATE_VERBOSE) 
   Log("LU factorizing...");
  W->LUFactorize();
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
cdouble *LayeredSubstrate::CreateScriptGTwiddleWorkspace()
{ int NI = NumInterfaces;
  int RSSize = 6*4*NI, WSize=16*NI*NI, WorkSize=2*RSSize + WSize;
  return (cdouble *)mallocEC(WorkSize*sizeof(cdouble));
}

void LayeredSubstrate::DestroyScriptGTwiddleWorkspace(cdouble *Workspace)
{ if (Workspace) free(Workspace);
}

/***************************************************************/
/* Get the (Fourier-space) 6x6 *inhomogeneous* DGF, i.e. the   */
/* (Fourier components of) the fields due to point sources in  */
/* the presence of the substrate.                              */
/*                                                             */
/* GTwiddle must point to a user-allocated complex-valued 6x6  */
/*  HMatrix.                                                   */
/*                                                             */
/* If Workspace is nonzero, it should be a user-allocated      */
/* cdouble array of length 48*N + 16*N*N                       */
/* where N = number of dielectric interface layers, not        */
/*  counting the ground plane (if present)                     */
/*                                                             */
/* If dzDest and/or dzSource are true, the derivative with     */
/* respect to zDest and/or zSource is returned instead.        */
/*                                                             */
/* If AddGamma0Twiddle is true, the free-space DGF is added    */
/* (only if source and destination points lie in same region). */
/***************************************************************/
void LayeredSubstrate::GetScriptGTwiddle(cdouble Omega, cdouble q2D[2],
                                         double zDest, double zSource,
                                         HMatrix *GTwiddle,
                                         cdouble *Workspace,
                                         bool dzDest, bool dzSource,
                                         bool AddGamma0Twiddle)
{
  if (zDest<zGP || zSource<zGP) // source or destination point below ground plane
   { GTwiddle->Zero();
     return;
   }

  UpdateCachedEpsMu(Omega);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
int Gamma0DestTerm   = Gamma0DSTerm/10;
int Gamma0SourceTerm = Gamma0DSTerm%10;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  /**********************************************************************/
  /**********************************************************************/
  /**********************************************************************/
  int NI=NumInterfaces;
  bool GroundPlaneOnly = (NI==0);
  if (GroundPlaneOnly && !AddGamma0Twiddle) 
   ErrExit("%s:%i:internal error",__FILE__,__LINE__);
  if (ForceFreeSpace || (GroundPlaneOnly && AddGamma0Twiddle))
   { 
     cdouble Gamma0Twiddle[6][6];
     int ForceLayer=0;
     double Sign=0.0;
     bool Accumulate=false;

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
if (Gamma0DSTerm)
 Gamma0Term = (Gamma0DestTerm==Gamma0SourceTerm) ? Gamma0DestTerm : -1;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
      
     GetGamma0Twiddle(Omega, q2D, zDest, zSource, Gamma0Twiddle,
                      ForceLayer, Sign, Accumulate, dzDest, dzSource);

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
 Gamma0Term = 0;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
     for(int Mu=0; Mu<6; Mu++)
      for(int Nu=0; Nu<6; Nu++)
       GTwiddle->SetEntry(Mu,Nu,Gamma0Twiddle[Mu][Nu]);
     return;
   }

  /**********************************************************************/
  /* compute GTwiddle as a sum of matrix-matrix-matrix products via     */
  /* equation (18) in the memo:                                         */
  /*  GTwiddle = \sum RTwiddle * W * STwiddle,                          */
  /* where RTwiddle has dimension 6  x 4N                               */
  /*       W        has dimension 4N x 4N                               */
  /*       STwiddle has dimension 4N x 6                                */
  /**********************************************************************/
  bool OwnsWorkspace = (Workspace==0);
  if (OwnsWorkspace)
   Workspace = CreateScriptGTwiddleWorkspace();
  int RSSize = 6*4*NI, WSize=16*NI*NI;
  HMatrix RTwiddle(6,    4*NI, LHM_COMPLEX, Workspace + 0             );
  HMatrix WMatrix (4*NI, 4*NI, LHM_COMPLEX, Workspace + RSSize        );
  HMatrix STwiddle(4*NI, 6,    LHM_COMPLEX, Workspace + RSSize + WSize);

  /**********************************************************************/
  /* assemble RHS vector for each (source point, polarization, orientation).*/
  /* loops over a and b run over the interfaces above and below the dest*/
  /* or source point.                                                   */
  /**********************************************************************/
  int DestLayer   = GetLayerIndex(zDest);
  int aMin        = (DestLayer==0               ? 1 : 0 );
  int aMax        = (DestLayer==NumInterfaces   ? 0 : 1 );
  cdouble ScriptG0TDest[2][6][6];
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
 Gamma0Term = Gamma0DestTerm;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  for(int a=aMin, nlDest=DestLayer-1+a; a<=aMax; a++, nlDest++)
   GetGamma0Twiddle(Omega, q2D, zDest, zInterface[nlDest],
                    ScriptG0TDest[a], DestLayer, 0.0, false,
                    dzDest, false);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
 Gamma0Term = 0;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  int SourceLayer  = GetLayerIndex(zSource);
  int bMin         = (SourceLayer==0             ? 1 : 0 );
  int bMax         = (SourceLayer==NumInterfaces ? 0 : 1 );
  cdouble ScriptG0TSource[2][6][6];
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
 Gamma0Term = Gamma0SourceTerm;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  for(int b=bMin, nlSource=SourceLayer-1+b; b<=bMax; b++, nlSource++)
   GetGamma0Twiddle(Omega, q2D, zInterface[nlSource], zSource,
                      ScriptG0TSource[b], SourceLayer, 0.0, false,
                      false, dzSource);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
 Gamma0Term = 0;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

 
  /**********************************************************************/
  /* assemble W matrix **************************************************/
  /**********************************************************************/
  ComputeW(Omega, q2D, &WMatrix);

  /**********************************************************************/
  /**********************************************************************/
  /**********************************************************************/
  RTwiddle.Zero();
  STwiddle.Zero();
  double Sign[2]={-1.0, 1.0};
  for(int a=aMin, nlDest=DestLayer-1+a; a<=aMax; a++, nlDest++)
   for(int i=0; i<6; i++)
    { RTwiddle.SetEntry(i, 4*nlDest+0, Sign[a]*ScriptG0TDest[a][i][0]);
      RTwiddle.SetEntry(i, 4*nlDest+1, Sign[a]*ScriptG0TDest[a][i][1]);
      RTwiddle.SetEntry(i, 4*nlDest+2, Sign[a]*ScriptG0TDest[a][i][3]);
      RTwiddle.SetEntry(i, 4*nlDest+3, Sign[a]*ScriptG0TDest[a][i][4]);
    };
  for(int b=bMin, nlSource=SourceLayer-1+b; b<=bMax; b++, nlSource++)
   for(int j=0; j<6; j++)
    { STwiddle.SetEntry(4*nlSource+0, j, -1.0*Sign[b]*ScriptG0TSource[b][0][j]);
      STwiddle.SetEntry(4*nlSource+1, j, -1.0*Sign[b]*ScriptG0TSource[b][1][j]);
      STwiddle.SetEntry(4*nlSource+2, j, -1.0*Sign[b]*ScriptG0TSource[b][3][j]);
      STwiddle.SetEntry(4*nlSource+3, j, -1.0*Sign[b]*ScriptG0TSource[b][4][j]);
    };
  WMatrix.LUSolve(&STwiddle);
  RTwiddle.Multiply(&STwiddle, GTwiddle);

  if (OwnsWorkspace) DestroyScriptGTwiddleWorkspace(Workspace);

  /**********************************************************************/
  /**********************************************************************/
  /**********************************************************************/
  if (AddGamma0Twiddle && (DestLayer==SourceLayer) )
   { 
     cdouble Gamma0Twiddle[6][6];
     int ForceLayer=-1;
     bool ForceSign=0.0;
     bool Accumulate=false;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
if (Gamma0DSTerm)
 Gamma0Term = (Gamma0DestTerm==Gamma0SourceTerm) ? Gamma0DestTerm : -1;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
     GetGamma0Twiddle(Omega, q2D, zDest, zSource, Gamma0Twiddle,
                      ForceLayer, ForceSign, Accumulate, dzDest, dzSource);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
 Gamma0Term = 0;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
     for(int Mu=0; Mu<6; Mu++)
      for(int Nu=0; Nu<6; Nu++)
       GTwiddle->AddEntry(Mu,Nu,Gamma0Twiddle[Mu][Nu]);
     return;
   };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
cdouble *LayeredSubstrate::CreateScriptLTwiddleWorkspace()
{ int NI = NumInterfaces;
  int RSize = 20*4*NI, WSize=16*NI*NI, SSize = 4*NI*8;
  return (cdouble *)mallocEC( (RSize + WSize + SSize)*sizeof(cdouble) );
}

void LayeredSubstrate::DestroyScriptLTwiddleWorkspace(cdouble *Workspace)
{ if (Workspace) free(Workspace); }

void LayeredSubstrate::GetLambda0Twiddle(cdouble Omega, cdouble q2D[2],
                                         double zDest, double zSource, HMatrix *Lambda0Twiddle,
                                         int ForceLayer, double Sign, bool Accumulate,
                                         bool dzDest, bool dzSource)
{
  if (real(Omega)>=0.0)
   UpdateCachedEpsMu(Omega);

  int nl = ForceLayer;
  if (nl==-1)
   { nl = GetLayerIndex(zDest);
     if (nl!=GetLayerIndex(zSource) ) 
      Warn("%s:%i: internal inconsistency",__FILE__,__LINE__);
   };
  cdouble EpsRel = EpsLayer[nl], MuRel=MuLayer[nl];
  cdouble EpsAbs = EpsRel/ZVAC,  MuAbs=MuRel*ZVAC;
  cdouble k2 = EpsRel*MuRel*Omega*Omega, q2=q2D[0]*q2D[0] + q2D[1]*q2D[1];
  cdouble qz = sqrt(k2 - q2);
  if (imag(qz)<0) qz*=-1.0;
  if ( real(Omega)<0.0 ) // negative Omega means compute quasi-static limit:
   { Omega*=-1.0;        //  retain Omega prefactors, but set Omega->0 when
     qz = II*sqrt(q2);   //  computing qz.
   };

  // G0 is the Fourier transform of the scalar Helmholtz Green's function
  cdouble G0 = II*exp(II*qz*fabs(zDest-zSource)) / (4.0*M_PI*qz);
  if (Sign==0.0) 
   Sign = (zDest>zSource ? 1.0 : (zDest<zSource ? -1.0 : 0.0));
  if (dzDest)
   G0 *= Sign*II*qz;
  if (dzSource)
   G0 *= -1.0*Sign*II*qz;

  cdouble dG0[3];
  dG0[0] = II*q2D[0]*G0;
  dG0[1] = II*q2D[1]*G0;
  dG0[2] = Sign*II*qz*G0;

  if (!Accumulate) Lambda0Twiddle->Zero();

  // electric and magnetic vector potentials and their curls
  for(int i=0; i<3; i++)
   {
     int ip1 = (i+1)%3, ip2 = (i+2)%3;
     Lambda0Twiddle->AddEntry(_AEX + i, _KX + i, MuAbs*G0);
     Lambda0Twiddle->AddEntry(_AMX + i, _NX + i, EpsAbs*G0);
     Lambda0Twiddle->AddEntry(_CURLAEX + i, _KX + ip2,      MuAbs* dG0[ip1]);
     Lambda0Twiddle->AddEntry(_CURLAEX + i, _KX + ip1, -1.0*MuAbs* dG0[ip2]);
     Lambda0Twiddle->AddEntry(_CURLAMX + i, _NX + ip2,      EpsAbs*dG0[ip1]);
     Lambda0Twiddle->AddEntry(_CURLAMX + i, _NX + ip1, -1.0*EpsAbs*dG0[ip2]);
   };

  // electric and magnetic scalar potentials and their gradients
  Lambda0Twiddle->AddEntry(_PHIE, _DIVK, G0 / EpsAbs);
  Lambda0Twiddle->AddEntry(_PHIM, _DIVN, G0 / MuAbs);
  for(int i=0; i<3; i++)
   { Lambda0Twiddle->AddEntry(_GRADPHIEX + i, _DIVK, dG0[i]/EpsAbs);
     Lambda0Twiddle->AddEntry(_GRADPHIMX + i, _DIVN, dG0[i]/MuAbs);
   };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void LayeredSubstrate::GetScriptLTwiddle(cdouble Omega, cdouble q2D[2],
                                         double zDest, double zSource,
                                         HMatrix *LTwiddle, cdouble *Workspace,
                                         bool dzDest, bool dzSource,
                                         bool AddLambda0Twiddle)
{
  if (zDest<zGP || zSource<zGP) // source or destination point below ground plane
   { LTwiddle->Zero();
     return;
   }

  UpdateCachedEpsMu(Omega);

  /**********************************************************************/
  /**********************************************************************/
  /**********************************************************************/
  int NI=NumInterfaces;
  bool GroundPlaneOnly = (NI==0);
  if (GroundPlaneOnly && !AddLambda0Twiddle)
   ErrExit("%s:%i:internal error",__FILE__,__LINE__);
  if (ForceFreeSpace || (GroundPlaneOnly && AddLambda0Twiddle))
   { 
     int ForceLayer=0;
     double Sign=0.0;
     bool Accumulate=false;
     GetLambda0Twiddle(Omega, q2D, zDest, zSource, LTwiddle,
                       ForceLayer, Sign, Accumulate, dzDest, dzSource);
     return;
   }

  /**********************************************************************/
  /* compute ScriptLTwiddle as a sum of matrix-matrix-matrix products   */
  /* equation (18) in the memo:                                         */
  /*  LTwiddle = \sum RTwiddle * W * STwiddle,                          */
  /* where RTwiddle has dimension 20 x 4N                               */
  /*       W        has dimension 4N x 4N                               */
  /*       STwiddle has dimension 4N x 8                                */
  /*                                                                    */
  /* RTwiddle[Mu,4*n+I] = contribution of Ith surface-current component */
  /*                      (I=0,1,2,3 for Kx, Ky, Nx, Ny) on interface #n*/
  /*                      to Muth component of extended potential       */
  /*                      vector at destination point                   */
  /*                                                                    */
  /* STwiddle[4*n+J,Nu] = Jth tangential field component                */
  /*                      (J=0,1,2,3 for Ex, Ey, Hx, Hy) on interface #n*/
  /*                      due to Nuth source component at source point  */
  /**********************************************************************/
  bool OwnsWorkspace = (Workspace==0);
  if (OwnsWorkspace)
   Workspace = CreateScriptLTwiddleWorkspace();
  int RSize = 20*4*NI, WSize=16*NI*NI; // SSize = 4*NI*8;
  HMatrix RTwiddle(20,   4*NI, LHM_COMPLEX, Workspace + 0            );
  HMatrix WMatrix (4*NI, 4*NI, LHM_COMPLEX, Workspace + RSize        );
  HMatrix STwiddle(4*NI, 8,    LHM_COMPLEX, Workspace + RSize + WSize);

  /**********************************************************************/
  /* assemble RHS vector for each (source point, polarization, orientation).*/
  /* loops over a and b run over the interfaces above and below the dest*/
  /* or source point.                                                   */
  /**********************************************************************/
  int DestLayer   = GetLayerIndex(zDest);
  int aMin        = (DestLayer==0               ? 1 : 0 );
  int aMax        = (DestLayer==NumInterfaces   ? 0 : 1 );
  cdouble L0TDestBuffer[2*20*8];
  HMatrix Lambda0TDest0(20,8,LHM_COMPLEX,L0TDestBuffer+0);
  HMatrix Lambda0TDest1(20,8,LHM_COMPLEX,L0TDestBuffer+20*8);
  HMatrix *Lambda0TDest[2];
  Lambda0TDest[0] = &Lambda0TDest0;
  Lambda0TDest[1] = &Lambda0TDest1;
  for(int a=aMin, nlDest=DestLayer-1+a; a<=aMax; a++, nlDest++)
   GetLambda0Twiddle(Omega, q2D, zDest, zInterface[nlDest],
                     Lambda0TDest[a], DestLayer, 0.0, false,
                     dzDest, false);

  int SourceLayer  = GetLayerIndex(zSource);
  int bMin         = (SourceLayer==0             ? 1 : 0 );
  int bMax         = (SourceLayer==NumInterfaces ? 0 : 1 );
  cdouble L0TSourceBuffer[2*20*8];
  HMatrix Lambda0TSource0(20,8,LHM_COMPLEX,L0TSourceBuffer+0);
  HMatrix Lambda0TSource1(20,8,LHM_COMPLEX,L0TSourceBuffer+20*8);
  HMatrix *Lambda0TSource[2];
  Lambda0TSource[0] = &Lambda0TSource0;
  Lambda0TSource[1] = &Lambda0TSource1;
  for(int b=bMin, nlSource=SourceLayer-1+b; b<=bMax; b++, nlSource++)
   GetLambda0Twiddle(Omega, q2D, zInterface[nlSource], zSource,
                     Lambda0TSource[b], SourceLayer, 0.0, false,
                     false, dzSource);
 
  /**********************************************************************/
  /* assemble W matrix **************************************************/
  /**********************************************************************/
  ComputeW(Omega, q2D, &WMatrix);

  /**********************************************************************/
  /**********************************************************************/
  /**********************************************************************/
  RTwiddle.Zero();
  STwiddle.Zero();
  double Sign[2]={-1.0, 1.0};
  cdouble iw = II*Omega;
  for(int a=aMin, nlDest=DestLayer-1+a; a<=aMax; a++, nlDest++)
   for(int Mu=0; Mu<20; Mu++)
    for(int i=0; i<2; i++)
     { RTwiddle.SetEntry(Mu, 4*nlDest+0+i, Sign[a]*(            Lambda0TDest[a]->GetEntry(Mu,_KX+i)
                                                    + II*q2D[i]*Lambda0TDest[a]->GetEntry(Mu,_DIVK)
                                                   )
                        );
       RTwiddle.SetEntry(Mu, 4*nlDest+2+i, Sign[a]*(            Lambda0TDest[a]->GetEntry(Mu,_NX+i)
                                                    + II*q2D[i]*Lambda0TDest[a]->GetEntry(Mu,_DIVN)
                                                   )
                        );
    };
  for(int b=bMin, nlSource=SourceLayer-1+b; b<=bMax; b++, nlSource++)
   for(int i=0; i<2; i++)
    { STwiddle.SetEntry(4*nlSource+0+i, _KX+i, -1.0*Sign[b]*iw*Lambda0TSource[b]->GetEntry(_AEX+i,_KX+i));
      STwiddle.SetEntry(4*nlSource+0+i, _DIVK, -1.0*Sign[b]*(-1.0/iw)*Lambda0TSource[b]->GetEntry(_GRADPHIEX+i,_DIVK));
      STwiddle.SetEntry(4*nlSource+2+i, _NX+i, -1.0*Sign[b]*iw*Lambda0TSource[b]->GetEntry(_AMX+i,_NX+i));
      STwiddle.SetEntry(4*nlSource+2+i, _DIVN, -1.0*Sign[b]*(-1.0/iw)*Lambda0TSource[b]->GetEntry(_GRADPHIMX+i,_DIVN));
    };
  WMatrix.LUSolve(&STwiddle);
  RTwiddle.Multiply(&STwiddle, LTwiddle);

  if (OwnsWorkspace) DestroyScriptLTwiddleWorkspace(Workspace);

  /**********************************************************************/
  /**********************************************************************/
  /**********************************************************************/
  if (AddLambda0Twiddle && (DestLayer==SourceLayer) )
   { int ForceLayer=-1;
     bool ForceSign=0.0;
     bool Accumulate=true;
     GetLambda0Twiddle(Omega, q2D, zDest, zSource, LTwiddle,
                       ForceLayer, ForceSign, Accumulate, dzDest, dzSource);
     return;
   };
}
