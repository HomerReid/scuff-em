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
 * GetPFT.cc   -- master (switchboard) routines for SCUFF algorithms
 *             -- for computing power, force, and torque
 *
 * homer reid  -- 1/2015
 */

#include "libscuff.h"
#include "PFTOptions.h"

namespace scuff { 

/***************************************************************/
/* function prototypes for the various PFT algorithms.         */
/*                                                             */
/* Note: In an effort to improve modularity and readability,   */
/*       the specific PFT algorithms are implemented as        */
/*       standalone (non-class-method) functions.              */
/*       Only the master GetPFT() routine is a class method in */
/*       RWGGeometry; it is just a switchboard routine that    */
/*       hands off to the various non-class-method functions   */
/*       to do the computation.                                */
/***************************************************************/

// PFT by overlap method
void GetOPFT(RWGGeometry *G, int SurfaceIndex, cdouble Omega,
             HVector *KNVector, HVector *RHS, HMatrix *DRMatrix,
             double PFT[NUMPFT], double **ByEdge=0);

// PFT by cartesian multipole-moment method
HMatrix *GetMomentPFTMatrix(RWGGeometry *G, cdouble Omega,
                            IncField *IF,
                            HVector *KNVector, HMatrix *DRMatrix=0,
                            HMatrix *PFTMatrix=0, bool Itemize=false);

// PFT by displaced-surface-integral method
void GetDSIPFT(RWGGeometry *G, cdouble Omega, double *kBloch,
               HVector *KN, IncField *IF, double PFT[NUMPFT],
               char *BSMesh, double R, int NumPoints,
               bool FarField, char *PlotFileName, 
               GTransformation *GT1, GTransformation *GT2);

HMatrix *GetEMTPFTMatrix(RWGGeometry *G, cdouble Omega, IncField *IF,
                         HVector *KNVector, HMatrix *DRMatrix=0,
                         HMatrix *PFTMatrix=0, bool Interior=false,
                         int EMTPFTIMethod=SCUFF_EMTPFTI_DEFAULT,
                         bool Itemize=false);

// matrix-trace version of DSIPFT for non-equilibrium calculations
void GetDSIPFTTrace(RWGGeometry *G, cdouble Omega, HMatrix *DRMatrix,
                    double PFT[NUMPFT],
                    char *DSIMesh=0, double DSIRadius=0.0, int DSIPoints=302,
                    bool FarField=false, char *PlotFileName=0,
                    GTransformation *GT1=0, GTransformation *GT2=0,
                    HMatrix *RFMatrix=0, bool RFMatrixDirty=true);

/***************************************************************/
/***************************************************************/
/***************************************************************/
double AutodetectDSIRadius(RWGSurface *S, GTransformation *GT1, GTransformation *GT2)
{
  double X0[3]={0.0, 0.0, 0.0};
  if (GT1) GT1->Apply(X0);
  if (GT2) GT2->Apply(X0);
  double Radius=0.0;
  for(int np=0; np<S->NumPanels; np++)
    { double *XP=S->Panels[np]->Centroid;
      Radius = fmax(Radius, VecDistance(X0, XP));
    };
  return 1.5*Radius;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetKNBilinears(HVector *KNVector, HMatrix *DRMatrix,
                    bool IsPECA, int KNIndexA,
                    bool IsPECB, int KNIndexB,
                    cdouble Bilinears[4])
{
  cdouble KK, KN, NK, NN;
  if (KNVector && (IsPECA || IsPECB) )
   { 
     cdouble kAlpha =  KNVector->GetEntry(KNIndexA);
     cdouble kBeta  =  KNVector->GetEntry(KNIndexB);
     KK = conj(kAlpha) * kBeta;
     KN = NK = NN = 0.0;
   }
  else if (KNVector)
   { 
     cdouble kAlpha =       KNVector->GetEntry(KNIndexA+0);
     cdouble nAlpha = -ZVAC*KNVector->GetEntry(KNIndexA+1);
     cdouble kBeta  =       KNVector->GetEntry(KNIndexB+0);
     cdouble nBeta  = -ZVAC*KNVector->GetEntry(KNIndexB+1);

     KK = conj(kAlpha) * kBeta;
     KN = conj(kAlpha) * nBeta;
     NK = conj(nAlpha) * kBeta;
     NN = conj(nAlpha) * nBeta;
   }
  else
   {
     KK = DRMatrix->GetEntry(KNIndexB+0, KNIndexA+0);
     KN = DRMatrix->GetEntry(KNIndexB+1, KNIndexA+0);
     NK = DRMatrix->GetEntry(KNIndexB+0, KNIndexA+1);
     NN = DRMatrix->GetEntry(KNIndexB+1, KNIndexA+1);
   }

  Bilinears[0]=KK;
  Bilinears[1]=KN;
  Bilinears[2]=NK;
  Bilinears[3]=NN;
}

void GetKNBilinears(RWGGeometry *G, int neA, int neB,
                    HVector *KNVector, HMatrix *DRMatrix,
                    cdouble Bilinears[4])
{
  int nbfA, nbfB;
  RWGSurface *SA=G->ResolveEdge(neA, 0, 0, &nbfA);
  RWGSurface *SB=G->ResolveEdge(neB, 0, 0, &nbfB);

  GetKNBilinears(KNVector, DRMatrix, 
                 SA->IsPEC, nbfA, SB->IsPEC, nbfB,
                 Bilinears);
}

/***************************************************************/
/* Get the power, force, and torque on surface #SurfaceIndex.  */
/***************************************************************/
void RWGGeometry::GetPFT(int SurfaceIndex, HVector *KN,
                         cdouble Omega, double PFT[NUMPFT],
                         PFTOptions *Options)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  PFTOptions DefaultOptions;
  if (Options==0)
   { Options=&DefaultOptions;
     InitPFTOptions(Options);
   };
  int PFTMethod      = Options->PFTMethod;
  char *FluxFileName = Options->FluxFileName;
  RWGSurface *S      = Surfaces[SurfaceIndex];
  IncField *IF       = Options->IF;
  HMatrix *DRMatrix  = Options->DRMatrix;

  /***************************************************************/
  /* allocate arrays for the edge-by-edge contributions if the   */
  /* user requested flux plots                                   */
  /***************************************************************/
  double **ByEdge=0;
  if (FluxFileName && PFTMethod!=SCUFF_PFT_DSI )
   { 
     int NE = S->NumEdges;
     ByEdge=(double **)mallocEC(NUMPFT*sizeof(double *));
     ByEdge[0]=(double *)mallocEC(NUMPFT*NE*sizeof(double));
     for(int nq=1; nq<NUMPFT; nq++)
      ByEdge[nq] = ByEdge[nq-1] + NE;
   };

  /***************************************************************/
  /* hand off to the individual PFT algorithms to do the         */
  /* computation                                                 */
  /***************************************************************/
  if ( PFTMethod==SCUFF_PFT_OVERLAP )
   { 
     HVector *RHSVector = Options->RHSVector;
     GetOPFT(this, SurfaceIndex, Omega, KN,
             RHSVector, DRMatrix, PFT, ByEdge);
   }
  else if ( PFTMethod==SCUFF_PFT_DSI )
   { 
     char *DSIMesh        = Options->DSIMesh;
     double DSIRadius     = Options->DSIRadius;
     int DSIPoints        = Options->DSIPoints;
     bool DSIFarField     = Options->DSIFarField;
     double *kBloch       = Options->kBloch;
     GTransformation *GT1 = S->OTGT;
     GTransformation *GT2 = S->GT;

     if (DSIRadius==0)
      DSIRadius=AutodetectDSIRadius(S, GT1, GT2);

     if (DRMatrix==0)
      GetDSIPFT(this, Omega, kBloch, KN, IF, PFT,
                DSIMesh, DSIRadius, DSIPoints,
                DSIFarField, FluxFileName, GT1, GT2);
     else 
      GetDSIPFTTrace(this, Omega, DRMatrix, PFT,
                     DSIMesh, DSIRadius, DSIPoints,
                     DSIFarField, FluxFileName, GT1, GT2);
   }
  else if (PFTMethod==SCUFF_PFT_MOMENTS)
   { 
     HMatrix *PFTMatrix 
      = GetMomentPFTMatrix(this, Omega, IF, KN, DRMatrix);
     PFTMatrix->GetEntriesD(SurfaceIndex, ":", PFT);
     delete PFTMatrix;
   }
  else if (PFTMethod==SCUFF_PFT_EMT)
   { 
     bool Interior      = Options->Interior;
     int  Method        = Options->EMTPFTIMethod;  
     HMatrix *PFTMatrix 
      = GetEMTPFTMatrix(this, Omega, IF, KN, DRMatrix, 0, Interior, Method);
     PFTMatrix->GetEntriesD(SurfaceIndex, ":", PFT);
     delete PFTMatrix;
   };


  /***************************************************************/
  /* produce flux plots if that was requested ********************/
  /***************************************************************/
  if (ByEdge)
   { 
     static const char *PFTNames[NUMPFT]
      ={"PAbs","PScat","FX","FY","FZ","TX","TY","TZ"};

     for(int nq=0; nq<NUMPFT; nq++)
      { char Tag[20];
        snprintf(Tag,20,"%s(%s)",PFTNames[nq],z2s(Omega));
        S->PlotScalarDensity(ByEdge[nq],2,FluxFileName,Tag);
      };

     free(ByEdge[0]);
     free(ByEdge);

   };
 
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *RWGGeometry::GetPFTMatrix(HVector *KN, cdouble Omega,
                                   PFTOptions *Options, 
                                   HMatrix *PFTMatrix)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  PFTOptions DefaultOptions;
  if (Options==0)
   { Options=&DefaultOptions;
     InitPFTOptions(Options);
   };
  IncField *IF=Options->IF;
  HMatrix *DRMatrix=Options->DRMatrix;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (    PFTMatrix
       && ( PFTMatrix->NR != NumSurfaces || PFTMatrix->NC != NUMPFT)
     ) 
   { Warn("wrong-size PFTMatrix in GetPFTMatrix (reallocating...)");
     delete PFTMatrix;
     PFTMatrix=0;
   };
  if (PFTMatrix==0)
   PFTMatrix=new HMatrix(NumSurfaces, NUMPFT);
  
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (Options->PFTMethod==SCUFF_PFT_EMT)
   { 
     GetEMTPFTMatrix(this, Omega, IF, KN, DRMatrix,
                     PFTMatrix, Options->Interior, Options->EMTPFTIMethod);
   }
  else if (Options->PFTMethod==SCUFF_PFT_MOMENTS)
   { 
     GetMomentPFTMatrix(this, Omega, IF, KN, DRMatrix, PFTMatrix);
   }
  else
   { 
     for(int ns=0; ns<NumSurfaces; ns++)
      { double PFT[NUMPFT];
        GetPFT(ns, KN, Omega, PFT, Options);
        PFTMatrix->SetEntriesD(ns,":",PFT);
      };
   };

  return PFTMatrix;

}

/***************************************************************/
/* PFTBySurface is an NSxNQ input matrix whose (ns,nq) entry is*/
/* the #nqth PFT quantity for surface #ns.                     */
/* PFTByRegion is an NRxNQ output matrix whose (nr,nq) entry is*/
/* the #nqth PFT quantity for region #nr, calculated by        */
/* summing the surface-resolved PFTs in PFTBySurface, with     */
/* proper signs, for all surfaces bounding region #nr.         */
/***************************************************************/
HMatrix *GetPFTByRegion(RWGGeometry *G, HMatrix *PFTBySurface,
                        HMatrix *PFTByRegion)
{
  int NS = G->NumSurfaces;
  int NR = G->NumRegions;
  if (PFTByRegion && (PFTByRegion->NR!=NR || PFTByRegion->NC!=NUMPFT))
   { Warn("wrong-size matrix in GetPFTByRegion (reallocating)");
     delete PFTByRegion;
     PFTByRegion=0;
   };
  if (PFTByRegion==0)
   PFTByRegion = new HMatrix(NR, NUMPFT);

  PFTByRegion->Zero();
  for(int nr=0; nr<NR; nr++)
   for(int ns=0; ns<NS; ns++)
    { 
      RWGSurface *S=G->Surfaces[ns];
      if (S->IsPEC) continue;
      double Sign = 0.0;
      if (S->RegionIndices[0]==nr) 
       Sign=1.0;
      else if (S->RegionIndices[1]==nr) 
       Sign=-1.0;
      else
       continue;
      for(int nq=0; nq<NUMPFT; nq++)
       PFTByRegion->AddEntry(nr, nq, Sign*PFTBySurface->GetEntry(ns, nq));
    };
  return PFTByRegion;
} 

/***************************************************************/
/* routine for initializing a PFTOptions structure to default  */
/* values; creates and returns a new default structure if      */
/* called with Options=NULL or with no argument                */
/***************************************************************/
PFTOptions *InitPFTOptions(PFTOptions *Options)
{
  if (Options==0)
   Options = (PFTOptions *)mallocEC( sizeof(PFTOptions) );

  // general options
  Options->PFTMethod = SCUFF_PFT_DEFAULT;
  Options->FluxFileName=0;
  Options->DRMatrix=0;
  Options->IF=0;
  Options->kBloch=0;

  // options affecting overlap PFT computation
  Options->RHSVector = 0;

  // options affecting DSI PFT computation
  // (note: DSIRadius=0 means autodetect the radius)
  Options->DSIMesh=0;
  Options->DSIRadius=0.0;
  Options->DSIPoints=302;
  Options->DSIFarField=false;
 
  // options affecting EMTPFT power computation
  Options->Itemize=false;
  Options->Interior=false;
  Options->EMTPFTIMethod=SCUFF_EMTPFTI_DEFAULT;
  Options->TInterior=Options->TExterior=0;

  // options affecting moment PFT computation
  Options->KeepQpTerm=false;
  Options->MomentFileBase=0;

  Options->GetRegionPFTs=false;

  return Options;
}

} // namespace scuff
