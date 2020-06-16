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
 * GetDyadicGF.cc  -- compute the scattering parts of the electric and
 *                    magnetic dyadic Green's functions
 *
 * homer reid      -- 5/2012
 *                 -- 20151212 overhaul and major acceleration
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <ctype.h>

#include <libhrutil.h>

#include "libscuff.h"
#include "libscuffInternals.h"
#include "PanelCubature.h"

#define II cdouble(0.0,1.0)

namespace scuff { 

/***************************************************************/
/* 20151212 new routine for computing dyadic GFs that uses a   */
/* much faster strategy and computes DGFs at many points       */
/* simultaneously.                                             */
/* There are two options:                                      */
/*  (a) XMatrix has dimension NX x 3. In this case the source  */
/*      and destination points are the same, with cartesian    */
/*      components X[nx,0:2].                                  */
/*  (b) XMatrix has dimension NX x 6. In this case the source  */
/*      and destination points can be different, and we have   */
/*       X[nx,0:2] = destination point                         */
/*       X[nx,3:5] = source points                             */
/***************************************************************/
HMatrix *RWGGeometry::GetDyadicGFs(cdouble Omega, double *kBloch,
                                   HMatrix *XMatrix, HMatrix *M,
                                   HMatrix *GMatrix,
                                   bool ScatteringOnly)
{ 
  int NBF = TotalBFs;
  int NX  = XMatrix->NR;
  Log("Getting DGFs at %i eval points...",NX);
  
  /* Decide whether real RFSource and real RFDest suffice: */
  /* Condition: for compact geo at imag freq, despite  k=0 */
  bool RealRF = (!LDim && (real(Omega)==0)); 
  bool HavekBloch = false;
  if (kBloch)
   for(int d=0; d<LDim; d++)
    if (kBloch[d]!=0.0) HavekBloch=true;  

  /*--------------------------------------------------------------*/
  /* allocate storage for RFSource, RFDest matrices. I keep these */
  /* on hand as statically-allocated buffers on the assumption    */
  /* that the routine will be called many times with the same     */
  /* number of evaluation points, for example in Brillouin-zone   */
  /* integrations                                                 */
  /*--------------------------------------------------------------*/
  static HMatrix *RFSource=0, *RFDest=0;
  if ( RFSource==0 || RFSource->NR!=NBF || RFSource->NC!=(6*NX) )
   { 
     if (RFSource) delete RFSource;
     if (RFDest)   delete RFDest;
     RFSource = new HMatrix(NBF, 6*NX, RealRF?LHM_REAL:LHM_COMPLEX);
     RFDest   = new HMatrix(NBF, 6*NX, RealRF?LHM_REAL:LHM_COMPLEX);
   };

  /*--------------------------------------------------------------*/
  /*- allocate an output matrix of the right size if necessary   -*/
  /*--------------------------------------------------------------*/
  if (GMatrix==0 || GMatrix->NR!=NX || GMatrix->NC!=18)
   { 
     if (GMatrix) 
      { Warn("wrong-size GMatrix passed to GetDyadicGFs (reallocating)");
        delete GMatrix;
      };
     GMatrix=new HMatrix(NX, 18, LHM_COMPLEX);
   };

  /*--------------------------------------------------------------*/
  /*- if we are computing two-point DGFs, initialize a PointSource*/
  /*- structure to compute the direct (non-scattering) contributions */
  /*--------------------------------------------------------------*/
  bool TwoPointDGF = (XMatrix->NC >= 6);
  bool AddDirectContribution = TwoPointDGF && !ScatteringOnly;
  PointSource PSBuffer, *PS = (AddDirectContribution ? &PSBuffer : 0);

  /*--------------------------------------------------------------*/
  /*- precompute 'reduced-field' vectors -------------------------*/
  /*--------------------------------------------------------------*/
  GetRFMatrix(Omega, kBloch, XMatrix, RFDest, true);
  if ( TwoPointDGF || HavekBloch )
   GetRFMatrix(Omega, kBloch, XMatrix, RFSource, false, TwoPointDGF ? 3 : 0);
  else
   RFSource->Copy(RFDest);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Log(" LUSolving...");
  
  M->LUSolve(RFSource); //Set a break point here  

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Log(" Computing VMVPs...");
  for(int nx=0; nx<NX; nx++)
   { 
     double XDest[3];
     XMatrix->GetEntriesD(nx,"0:2",XDest);
     int nr=GetRegionIndex(XDest);
     if (nr==-1) continue;

     cdouble Eps, Mu;
     RegionMPs[nr]->GetEpsMu(Omega, &Eps, &Mu);
     cdouble k    = Omega * sqrt(Eps*Mu);
     cdouble ZRel = sqrt(Mu/Eps);
     cdouble GEScatNormFac = -1.0/(II*k*ZVAC*ZVAC*ZRel);
     cdouble GMScatNormFac = +ZRel/(II*k);

     cdouble GEDirect[3][3], GMDirect[3][3];
     double LastXSource[3];
     if (AddDirectContribution)
      { 
        cdouble GEDirectNormFac = k*k/Eps;
        cdouble GMDirectNormFac = k*k/Mu;
        double XSource[3];
        XMatrix->GetEntriesD(nx,"3:5",XSource);
        PS->SetX0(XSource);
        if (nx==0 || !VecEqualFloat(XSource,LastXSource) )
         UpdateIncFields(PS, Omega, kBloch);
        memcpy(LastXSource,XSource,3*sizeof(double));
        for(int i=0; i<3; i++)
         { cdouble EH[6];
           cdouble P[3]={0.0, 0.0, 0.0};
           P[i]=1.0;
           PS->SetP(P);
           PS->SetType(LIF_ELECTRIC_DIPOLE);
           PS->GetFields(XDest, EH);
           for(int j=0; j<3; j++)
            GEDirect[j][i] = EH[0+j] / GEDirectNormFac;
           PS->SetType(LIF_MAGNETIC_DIPOLE);
           PS->GetFields(XDest, EH);
           for(int j=0; j<3; j++)
            GMDirect[j][i] = EH[3+j] / GMDirectNormFac;
         };
      };

     for(int i=0; i<3; i++)
      for(int j=0; j<3; j++)
       { 
         cdouble GEScat=0.0, GMScat=0.0;
         for(int nbf=0; nbf<NBF; nbf++)
          { GEScat += RFDest->GetEntry(nbf, 6*nx+0+i) * RFSource->GetEntry(nbf, 6*nx+0+j);
            GMScat += RFDest->GetEntry(nbf, 6*nx+3+i) * RFSource->GetEntry(nbf, 6*nx+3+j);
          };
         GEScat *= GEScatNormFac;
         GMScat *= GMScatNormFac;

         GMatrix->SetEntry(nx, 0 + 3*i + j, GEScat);
         GMatrix->SetEntry(nx, 9 + 3*i + j, GMScat);

         if (AddDirectContribution)
          { GMatrix->AddEntry(nx, 0 + 3*i + j, GEDirect[i][j]);
            GMatrix->AddEntry(nx, 9 + 3*i + j, GMDirect[i][j]);
          };
       };
   };

  return GMatrix;

}

/***************************************************************/
/* wrapper to support legacy interface to GetDyadicGFs         */
/***************************************************************/
void RWGGeometry::GetDyadicGFs(double XEval[3], double XSource[3],
                               cdouble Omega, double *kBloch,
                               HMatrix *M, HVector *KN,
                               cdouble GEScat[3][3], cdouble GMScat[3][3],
                               cdouble GETot[3][3], cdouble GMTot[3][3])
{ 
  (void )KN; // not used, retained for backward compatibility

  double XBuffer[6];
  memcpy(XBuffer+0, XEval,   3*sizeof(double));
  HMatrix *XMatrix = new HMatrix(1,(XEval == XSource)?3:6, LHM_REAL,LHM_NORMAL,XBuffer);
  if (XEval != XSource) memcpy(XBuffer+3, XSource, 3*sizeof(double));
  
  if ( (LBasis && !kBloch) || (!LBasis && kBloch) )
   ErrExit("%s:%i: incorrect kBloch specification",__FILE__,__LINE__);

  cdouble GBuffer[18];
  HMatrix GMatrix(1,18,LHM_COMPLEX,LHM_NORMAL,GBuffer);
  
  GetDyadicGFs(Omega, kBloch, XMatrix, M, &GMatrix, true);

  for(int i=0; i<3; i++)
   for(int j=0; j<3; j++)
    { GEScat[i][j] = GMatrix.GetEntry(0, 0 + 3*i + j );
      GMScat[i][j] = GMatrix.GetEntry(0, 9 + 3*i + j );
    };
  
  /***************************************************************/
  /* add direction contributions of point source                 */
  /***************************************************************/
  cdouble EpsRel, MuRel;
  int nr=GetRegionIndex(XSource);
  RegionMPs[nr]->GetEpsMu(Omega, &EpsRel, &MuRel);
  cdouble k=Omega*sqrt(EpsRel*MuRel);
  cdouble EFactor = k*k/EpsRel;
  cdouble MFactor = k*k/MuRel;

  cdouble P[3]={1.0, 0.0, 0.0};
  PointSource PS(XSource, P);
  if (LDim>0)
   { PS.SetLattice(LBasis);
     PS.SetkBloch(kBloch);
   };

  for(int j=0; j<3; j++)
   { 
     cdouble EH[6];

     // set point source to point in the jth direction
     memset(P, 0, 3*sizeof(cdouble));
     P[j]=1.0;
     PS.SetP(P);

     PS.SetType(LIF_ELECTRIC_DIPOLE);
     PS.GetFields(XEval, EH);
     GETot[0][j] = GEScat[0][j] + EH[0] / EFactor;
     GETot[1][j] = GEScat[1][j] + EH[1] / EFactor;
     GETot[2][j] = GEScat[2][j] + EH[2] / EFactor;

     PS.SetType(LIF_MAGNETIC_DIPOLE);
     PS.GetFields(XEval, EH);
     GMTot[0][j] = GMScat[0][j] + EH[3] / MFactor;
     GMTot[1][j] = GMScat[1][j] + EH[4] / MFactor;
     GMTot[2][j] = GMScat[2][j] + EH[5] / MFactor;
   };

   delete XMatrix;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RWGGeometry::GetDyadicGFs(double X[3], cdouble Omega, double *kBloch,
                               HMatrix *M, HVector *KN,
                               cdouble GEScat[3][3], cdouble GMScat[3][3])
{
  cdouble Dummy[3][3];
  GetDyadicGFs(X, X, Omega, kBloch, M, KN, GEScat, GMScat, Dummy, Dummy);
}

} // namespace scuff
