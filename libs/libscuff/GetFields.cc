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
 * You should have received a copy of the GNU General Public License * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * GetFields.cc  -- libscuff class methods for computing scattered
 *               -- electric and magnetic fields
 *
 * homer reid    -- 3/2007  -- 11/2009
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <libhrutil.h>
#include <libTriInt.h>
#include <libMDInterp.h>

#include "libscuff.h"
#include "libscuffInternals.h"
#include "PanelCubature.h"

#ifdef USE_PTHREAD
#  include <pthread.h>
#endif

#define NUMFIELDS 6 // Ex, Ey, Ez, Hx, Hy, Hz
#define II cdouble(0,1)

namespace scuff {

cdouble GetG(double R[3], cdouble k, cdouble *dG, cdouble *ddG=0);

void GetReducedFields_Nearby(RWGSurface *S, const int ne,
                             const double X0[3], const cdouble k,
                             cdouble e[3], cdouble h[3]);

/***************************************************************/
/* integrand for computation of reduced fields                 */
/***************************************************************/
void GetGCBar2D_Fourier(cdouble k, double *kBloch,
                        HMatrix *RLBasis, double RLVolume,
                        double *XDest, double *XSource,
                        cdouble G[3][3], cdouble C[3][3]);

typedef struct RFIData
 { cdouble k;
   double *X0;
   GBarAccelerator *GBA;
   HMatrix *RLBasis;
   double RLVolume;
   bool NewMethod;
 } RFIData;

void RFIntegrand(double X[3], double b[3], double Divb,
                 void *UserData, double W, double *Integral)
{
  RFIData *Data       = (RFIData *)UserData;
  cdouble k            = Data->k;
  double *X0           = Data->X0;
  GBarAccelerator *GBA = Data->GBA;

  cdouble *GC = (cdouble *)Integral;

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
if (Data->NewMethod)
{
cdouble GG[3][3], CC[3][3];
GetGCBar2D_Fourier(k, GBA->kBloch, Data->RLBasis, Data->RLVolume, X, X0, GG, CC);
for(int Mu=0; Mu<3; Mu++)
 for(int Nu=0; Nu<3; Nu++)
  { GC[0 + Mu] += W*GG[Mu][Nu]*b[Nu];
    GC[3 + Mu] -= W*CC[Mu][Nu]*b[Nu];
  };
return;
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  double XmX0[3];
  XmX0[0] = X[0] - X0[0];
  XmX0[1] = X[1] - X0[1];
  XmX0[2] = X[2] - X0[2];

  cdouble G0, dG[3];
  if (GBA)
   G0=GetGBar(XmX0, GBA, dG);
  else
   G0=GetG(XmX0, k, dG);

  cdouble k2=k*k, ik=II*k;
  for(int i=0; i<3; i++)
   GC[i] += W*(G0*b[i] - Divb*dG[i]/k2);

  GC[3+0] += W*(b[1]*dG[2] - b[2]*dG[1]) / (-1.0*ik);
  GC[3+1] += W*(b[2]*dG[0] - b[0]*dG[2]) / (-1.0*ik);
  GC[3+2] += W*(b[0]*dG[1] - b[1]*dG[0]) / (-1.0*ik);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
#if 0
void GetDistanceToPanel(double X0[3], RWGSurface *S, int np)
{
  RWGPanel *P=S->Panels[np];
  double XdN=VecDot(X0, P->ZHat);

  double XPerp[3], XPar[3];
  l
}

void GetDistanceToBF(double X0[3], RWGSurface *S, int ne)
{
  RWGEdge *E=S->Edges[ne];
  double Distance = GetDistanceToPanel(X0, S, E->iPPanel);
  if (E->iMPanel!=-1)
   Distance = fmin(Distance, GetDistanceToPanel(X0, S, E->iMPanel);
   
}
#endif

/***************************************************************/
/* RFMatrix is a matrix of "reduced fields", i.e. a matrix     */
/* whose columns may be dot-producted with the KN vector (BEM  */
/* system solution vector) to yield components of the          */
/* scattered E and H fields.                                   */
/* More specifically, for Mu=0...5, the (6*nx + Mu)th column   */
/* of RFMatrix is dotted into KN to yield the Muth component   */
/* of the field six-vector F=\{ E \choose H \}.                */
/***************************************************************/
HMatrix *RWGGeometry::GetRFMatrix(cdouble Omega, double *kBloch0,
                                  HMatrix *XMatrix, HMatrix *RFMatrix,
                                  bool MinuskBloch, int ColumnOffset)
{
  /* Decide whether real RFSource and real RFDest suffice */
  bool RealRF = (!LDim && (real(Omega)==0)); 

  /* Flip kBloch if necessary */
  double *kBloch=kBloch0;
  double kBlochBuffer[3];
  if (kBloch && MinuskBloch)
   { kBloch = kBlochBuffer;
     kBloch[0] = kBloch[1] = kBloch[2] = 0.0;
     for(int d=0; d<LDim; d++)
      kBloch[d] = -1.0*kBloch0[d];
   };
  
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NE  = TotalEdges;
  int NBF = TotalBFs;
  int NX  = XMatrix->NR;
  if (     RFMatrix==0 
       || (RFMatrix->NR != NBF) 
       || (RFMatrix->NC != 6*NX) 
     )
   { if (RFMatrix) 
      { Warn("wrong-size RFMatrix passed to GetRFMatrix; reallocating");
        delete RFMatrix;
      };
     RFMatrix = new HMatrix(NBF, 6*NX, RealRF?LHM_REAL:LHM_COMPLEX);
   };
  RFMatrix->Zero();

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double rRelOuterThreshold=4.0;
  double rRelInnerThreshold=1.0;
  int LowOrder=7;
  int HighOrder=20;
  char *s1=getenv("SCUFF_RREL_OUTER_THRESHOLD");
  char *s2=getenv("SCUFF_RREL_INNER_THRESHOLD");
  char *s3=getenv("SCUFF_LOWORDER");
  char *s4=getenv("SCUFF_HIGHORDER");
  if (s1) sscanf(s1,"%le",&rRelOuterThreshold);
  if (s2) sscanf(s2,"%le",&rRelInnerThreshold);
  if (s3) sscanf(s3,"%i",&LowOrder);
  if (s4) sscanf(s4,"%i",&HighOrder);
  if (s1||s2||s3||s4)
   Log("({O,I}rRelThreshold | LowOrder | HighOrder)=(%e,%e,%i,%i)",
       rRelOuterThreshold,rRelInnerThreshold,LowOrder,HighOrder);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble *ZRels   = new cdouble[NumRegions];
  cdouble *ks      = new cdouble[NumRegions];
  for(int nr=0; nr<NumRegions; nr++)
   { cdouble EpsRel, MuRel;
     RegionMPs[nr]->GetEpsMu(Omega, &EpsRel, &MuRel);
     ZRels[nr] = sqrt(MuRel/EpsRel);
     ks[nr]    = sqrt(MuRel*EpsRel) * Omega;
   };

  /***************************************************************/
  /* For the periodic-boundary-condition case, we need to        */
  /* initialize accelerator objects to accelerate computation of */
  /* the periodic Green's function in each extended region of    */
  /* the geometry.                                               */
  /***************************************************************/
  GBarAccelerator **RegionGBAs=0;
  if (LBasis)
   { RegionGBAs=
      (GBarAccelerator **)mallocEC(NumRegions*sizeof(RegionGBAs[0]));
     for(int nr=0; nr<NumRegions; nr++)
      if ( ! ( RegionMPs[nr]->IsPEC() ) )
       RegionGBAs[nr]=CreateRegionGBA(nr, Omega, kBloch, XMatrix);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  static bool Initialized=false, UseNewMethod=false;
  if (Initialized==false)
   { Initialized=true;
     char *s = getenv("SCUFF_NEW_RFMETHOD");
     if(s && s[0]=='1')
      { printf("Using new RF method.\n");
        UseNewMethod=true;
      };
   };
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NENX=NE*NX;
#ifndef USE_OPENMP
  if (LogLevel>SCUFF_VERBOSELOGGING)
   Log("Computing RFMatrix entries at %i points",NX);
#else
  int NumThreads=GetNumThreads();
  if (LogLevel>SCUFF_VERBOSELOGGING)
   Log("Computing RFMatrix entries (%i threads) at %i points",NumThreads,NX);
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int nenx=0; nenx<NENX; nenx++)
   { 
     int nx     = nenx / NE;
     int neFull = nenx % NE;

     int ns, ne, nbf;
     RWGSurface *S = ResolveEdge(neFull, &ns, &ne, &nbf);
     RWGEdge *E    = S->Edges[ne];

     double X[3];
     X[0]=XMatrix->GetEntryD(nx,ColumnOffset+0);
     X[1]=XMatrix->GetEntryD(nx,ColumnOffset+1);
     X[2]=XMatrix->GetEntryD(nx,ColumnOffset+2);
     int RegionIndex = GetRegionIndex(X);
     if (RegionIndex==-1) continue; // inside a closed PEC surface
   
     double Sign=0.0;
     if      (Surfaces[ns]->RegionIndices[0]==RegionIndex) 
      Sign=+1.0;
     else if (Surfaces[ns]->RegionIndices[1]==RegionIndex)
      Sign=-1.0;
     else 
      continue;
   
     cdouble k    = ks[RegionIndex];
     cdouble ZRel = ZRels[RegionIndex];

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     cdouble GC[6]; 
     RFIData MyData, *Data=&MyData;
     Data->X0  = X;
     Data->k   = ks[RegionIndex];
     Data->GBA = RegionGBAs ? RegionGBAs[RegionIndex] : 0;
     Data->RLBasis = RLBasis;
     Data->RLVolume= RLVolume;
     Data->NewMethod = UseNewMethod;

     double rRel = VecDistance(X, E->Centroid) / E->Radius;
     const int IDim=12;
     if (rRel >= rRelOuterThreshold)
      { 
        GetBFCubature2(this, ns, ne, RFIntegrand, (void *)Data,
                       IDim, LowOrder, (double *)GC);
      }
     else if (rRel>=rRelInnerThreshold)
      { 
        GetBFCubature2(this, ns, ne, RFIntegrand, (void *)Data,
                       IDim, HighOrder, (double *)GC);
      }
     else
      { 
        GetReducedFields_Nearby(S, ne, X, k, GC+0, GC+3);
        GC[3] /= (-II*k);
        GC[4] /= (-II*k);
        GC[5] /= (-II*k);

        if (RegionGBAs)
         { cdouble GC1[6], GC2[6];
           int Order=4;
           GetBFCubature2(this, ns, ne, RFIntegrand, (void *)Data,
                          IDim, Order, (double *)GC1);
           Data->GBA = 0;
           GetBFCubature2(this, ns, ne, RFIntegrand, (void *)Data,
                          IDim, Order, (double *)GC2);
           for(int Mu=0; Mu<6; Mu++) 
            GC[Mu] += (GC1[Mu] - GC2[Mu]);
         };
      };

     /***************************************************/
     /*                                                 */
     /* E = ik*Z0 * Zr * k*g + ik*n*c                   */
     /*   = ik*Z0 * Zr * k*g - ik*Z0*nScuff*c           */
     /* H =        -ik * k*c + (ik/(Z0*Zr)) * n*c       */
     /*   =        -ik * k*c - (ik/Zr) *nScuff*c        */
     /***************************************************/
     cdouble *GG=GC+0, *CC=GC+3;
     cdouble EKFactor =      Sign*II*k*ZRel*ZVAC;
     cdouble HKFactor = -1.0*Sign*II*k;
     cdouble ENFactor = -1.0*Sign*II*k*ZVAC;
     cdouble HNFactor = -1.0*Sign*II*k/ZRel;

     RFMatrix->SetEntry(nbf, 6*nx + 0, EKFactor * GG[0] );
     RFMatrix->SetEntry(nbf, 6*nx + 1, EKFactor * GG[1] );
     RFMatrix->SetEntry(nbf, 6*nx + 2, EKFactor * GG[2] );
     RFMatrix->SetEntry(nbf, 6*nx + 3, HKFactor * CC[0] );
     RFMatrix->SetEntry(nbf, 6*nx + 4, HKFactor * CC[1] );
     RFMatrix->SetEntry(nbf, 6*nx + 5, HKFactor * CC[2] );

     if ( !(S->IsPEC) )
      { RFMatrix->SetEntry(nbf+1, 6*nx + 0, ENFactor * CC[0] );
        RFMatrix->SetEntry(nbf+1, 6*nx + 1, ENFactor * CC[1] );
        RFMatrix->SetEntry(nbf+1, 6*nx + 2, ENFactor * CC[2] );
        RFMatrix->SetEntry(nbf+1, 6*nx + 3, HNFactor * GG[0] );
        RFMatrix->SetEntry(nbf+1, 6*nx + 4, HNFactor * GG[1] );
        RFMatrix->SetEntry(nbf+1, 6*nx + 5, HNFactor * GG[2] );
      };

   }; // for(int nenx=0; nenx<NENX; nenx++)

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (RegionGBAs)
   { for(int nr=0; nr<NumRegions; nr++)
      if (RegionGBAs[nr])
       DestroyGBarAccelerator(RegionGBAs[nr]);
     free(RegionGBAs);
   };

  delete[] ZRels;
  delete[] ks;

  return RFMatrix;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *RWGGeometry::GetFields(IncField *IFList, HVector *KN,
                                cdouble Omega, double *kBloch,
                                HMatrix *XMatrix, HMatrix *FMatrix)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if ( XMatrix==0 || XMatrix->NC<3 || XMatrix->NR==0 )
   ErrExit("wrong-size XMatrix (%ix%i) passed to GetFields",
            XMatrix->NR,XMatrix->NC);

  int NX=XMatrix->NR;
  if (LogLevel >= SCUFF_VERBOSELOGGING)
   Log("Computing fields at %i evaluation points...",NX);

  /***************************************************************/
  /* (re)allocate output matrix as necessary *********************/
  /***************************************************************/
  if (FMatrix==0 || FMatrix->NR!=NX || FMatrix->NC!=NUMFIELDS)
   { if (FMatrix)
      { Warn(" ** warning: wrong-size FMatrix passed to GetFields(); reallocating");
        delete FMatrix;
      };
     FMatrix=new HMatrix(NX, NUMFIELDS, LHM_COMPLEX);
   };
  FMatrix->Zero();

  /***************************************************************/
  /* the incident fields will most likely have been updated at   */
  /* the current frequency already by an earlier call to         */
  /* AssembleRHSVector(), but someone might call GetFields()     */
  /* to get information on just the incident fields before       */
  /* before setting up and solving the BEM problem, so we should */
  /* do this just to make sure.                                  */
  /***************************************************************/
  if (IFList)
   UpdateIncFields(IFList, Omega, kBloch);

  /***************************************************************/
  /* get contributions of surface currents if present ************/
  /***************************************************************/
  if (KN)
   {
     HMatrix *RFMatrix = GetRFMatrix(Omega, kBloch, XMatrix, 0, true);
     HMatrix KNMatrix(1, TotalBFs, LHM_COMPLEX, LHM_NORMAL, (void *)KN->ZV);
     HMatrix *FMatrixT = new HMatrix(1, 6*NX, LHM_COMPLEX);
     KNMatrix.Multiply(RFMatrix, FMatrixT);
     for(int nx=0; nx<NX; nx++)
      for(int Mu=0; Mu<6; Mu++)
       FMatrix->SetEntry(nx, Mu, FMatrixT->GetEntry(0, 6*nx + Mu));

     delete RFMatrix;
     delete FMatrixT;
   };

  /***************************************************************/
  /* add contributions of incident fields if present *************/
  /***************************************************************/
  if (IFList)
   for(int nx=0; nx<NX; nx++)
    { 
      double X[3];
      XMatrix->GetEntriesD(nx,"0:2",X);
      int RegionIndex = GetRegionIndex(X);
      if (RegionIndex==-1) continue; // inside a closed PEC surface

      for(IncField *IF=IFList; IF; IF=IF->Next)
       if ( IF->RegionIndex == RegionIndex )
        { cdouble EH[6];
          IF->GetFields(X, EH);
          for(int Mu=0; Mu<6; Mu++)
           FMatrix->AddEntry(nx, Mu, EH[Mu]);
        };
    };

  return FMatrix;
         
}

/***************************************************************/
/* alternative entry points to GetFields                       */
/***************************************************************/
HMatrix *RWGGeometry::GetFields(IncField *IF, HVector *KN, cdouble Omega,
                                HMatrix *XMatrix, HMatrix *FMatrix)
{ return GetFields(IF, KN, Omega, 0, XMatrix, FMatrix); }


void RWGGeometry::GetFields(IncField *IF, HVector *KN, cdouble Omega, 
                            double *kBloch, double *X, cdouble *EH)
{
  HMatrix XMatrix(1, 3, LHM_REAL, LHM_NORMAL, (void *)X);
  HMatrix FMatrix(1, 6, LHM_COMPLEX, LHM_NORMAL, (void *)EH);
  GetFields(IF, KN, Omega, kBloch, &XMatrix, &FMatrix);
} 

void RWGGeometry::GetFields(IncField *IF, HVector *KN, cdouble Omega, 
                            double *X, cdouble *EH)
{ GetFields(IF, KN, Omega, 0, X, EH); }

} // namespace scuff
