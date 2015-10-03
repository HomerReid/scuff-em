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
 * SurfaceSurfaceInteractions.cc -- libscuff routine for computing the 
 *                               -- matrix describing the interactions 
 *                               -- of two objects in the geometry)
 *                               --
 *                               -- (cf. 'libscuff Implementation and Technical
 *                               --  Details', section 8.3, 'Structure of the BEM
 *                               --  Matrix.')
 *                               --
 * homer reid                    -- 10/2006 -- 11/2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libhmat.h>
#include <libhrutil.h>

#include "libscuff.h"
#include "libscuffInternals.h"

#include "cmatheval.h"

#ifdef USE_PTHREAD
#  include <pthread.h>
#endif
#ifdef USE_OPENMP
#  include <omp.h>
#endif

namespace scuff {

#define II cdouble(0,1)

/***************************************************************/
/* Given two surfaces, identify whether they bound zero, one,  */
/* or two common regions. If there are any common regions,     */
/* identify their indices and the relative sign between the    */
/* contributions of surface currents on the two surfaces to    */
/* fields in those regions.                                    */
/***************************************************************/
int CountCommonRegions(RWGSurface *Sa, RWGSurface *Sb, 
                       int CommonRegionIndices[2], double Signs[2])
{
  int NumCommonRegions=0;

  if ( Sa->RegionIndices[0] == Sb->RegionIndices[0] )
   { CommonRegionIndices[NumCommonRegions] = Sa->RegionIndices[0];
     Signs[NumCommonRegions]=+1.0;
     NumCommonRegions++;
   }
  else if ( Sa->RegionIndices[0] == Sb->RegionIndices[1] )
   { CommonRegionIndices[NumCommonRegions] = Sa->RegionIndices[0];
     Signs[NumCommonRegions]=-1.0;
     NumCommonRegions++;
   }
  if ( Sa->RegionIndices[1] == Sb->RegionIndices[0] )
   { CommonRegionIndices[NumCommonRegions] = Sa->RegionIndices[1];
     Signs[NumCommonRegions]=-1.0;
     NumCommonRegions++;
   }
  else if ( !Sa->IsPEC && !Sb->IsPEC && Sa->RegionIndices[1] == Sb->RegionIndices[1] )
   { CommonRegionIndices[NumCommonRegions] = Sa->RegionIndices[1];
     Signs[NumCommonRegions]=+1.0;
     NumCommonRegions++;
   };

  return NumCommonRegions;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct ThreadData
 { 
   GetSSIArgStruct *Args;
   unsigned PPIAlgorithmCount[NUMPPIALGORITHMS];
   int nt, NumTasks;

 } ThreadData;

/***************************************************************/
/* 'GetSurfaceSurfaceInteractionThread'                        */
/***************************************************************/
void *GSSIThread(void *data)
{ 
  /***************************************************************/
  /* extract local copies of fields in argument structure */
  /***************************************************************/
  ThreadData *TD=(ThreadData *)data;
  GetSSIArgStruct *Args= TD->Args;
  RWGGeometry *G       = Args->G;
  RWGSurface *Sa       = Args->Sa;
  RWGSurface *Sb       = Args->Sb;
  cdouble Omega        = Args->Omega;
  int NumTorqueAxes    = Args->NumTorqueAxes;
  double *GammaMatrix  = Args->GammaMatrix;
  int RowOffset        = Args->RowOffset;
  int ColOffset        = Args->ColOffset;
  bool Symmetric       = Args->Symmetric;
  double *Displacement = Args->Displacement;
  HMatrix *B           = Args->B;
  HMatrix **GradB      = Args->GradB;
  HMatrix **dBdTheta   = Args->dBdTheta;
  cdouble EpsA         = Args->EpsA;
  cdouble EpsB         = Args->EpsB;
  cdouble MuA          = Args->MuA;
  cdouble MuB          = Args->MuB;
  double SignA         = Args->SignA;
  double SignB         = Args->SignB;
  bool SaIsPEC         = Args->SaIsPEC;
  bool SbIsPEC         = Args->SbIsPEC;

#ifdef USE_PTHREAD
  SetCPUAffinity(TD->nt);
#endif

  /***************************************************************/
  /* initialize an argument structure to be passed to            */
  /* GetEdgeEdgeInteractions() below                             */
  /***************************************************************/
  GetEEIArgStruct MyGetEEIArgs, *GetEEIArgs=&MyGetEEIArgs;
  InitGetEEIArgs(GetEEIArgs);

  GetEEIArgs->Sa=Sa;
  GetEEIArgs->Sb=Sb;
  GetEEIArgs->NumGradientComponents = GradB ? 3 : 0;
  GetEEIArgs->NumTorqueAxes=NumTorqueAxes;
  GetEEIArgs->GammaMatrix=GammaMatrix;
  GetEEIArgs->Displacement=Displacement;

  /* pointers to arrays inside the structure */
  cdouble *GC=GetEEIArgs->GC;
  cdouble *GradGC=GetEEIArgs->GradGC;
  cdouble *dGCdT=GetEEIArgs->dGCdT;

  /***************************************************************/
  /* precompute the constant prefactors that multiply the        */
  /* integrals returned by GetEdgeEdgeInteractions()             */
  /***************************************************************/
  cdouble kA, PreFac1A, PreFac2A, PreFac3A;
  cdouble kB, PreFac1B, PreFac2B, PreFac3B;

  kA=csqrt2(EpsA*MuA)*Omega;
  PreFac1A =  SignA*II*MuA*Omega;
  PreFac2A = -SignA*II*kA;
  PreFac3A = -SignA*II*EpsA*Omega;

  if (EpsB!=0.0)
   { 
     kB=csqrt2(EpsB*MuB)*Omega;
     PreFac1B =  SignB*II*MuB*Omega;
     PreFac2B = -SignB*II*kB;
     PreFac3B = -SignB*II*EpsB*Omega;
   };

  /***************************************************************/
  /* loop over all internal edges on both objects.               */
  /***************************************************************/
  int nea, NEa=Sa->NumEdges;
  int neb, NEb=Sb->NumEdges;
  int X, Y, Mu, nt=0;
  int NumGradientComponents = GradB ? 3 : 0;
  int nebStart = Symmetric ? 1 : 0;
  for(nea=0; nea<NEa; nea++)
   for(neb=nebStart*nea; neb<NEb; neb++)
    { 
      nt++;
      if (nt==TD->NumTasks) nt=0;
      if (nt!=TD->nt) continue;

      if (G->LogLevel>=SCUFF_VERBOSE2 && (neb==nebStart*nea) )
       LogPercent(nea, NEa);

      /*--------------------------------------------------------------*/
      /*- contributions of first medium (EpsA, MuA)  -----------------*/
      /*--------------------------------------------------------------*/
      GetEEIArgs->nea  = nea;
      GetEEIArgs->neb  = neb;
      GetEEIArgs->k    = kA;
      GetEEIArgs->GBA  = Args->GBA1;
      GetEdgeEdgeInteractions(GetEEIArgs);

      if ( SaIsPEC && SbIsPEC )
       { 
         X=RowOffset + nea;
         Y=ColOffset + neb;  

         B->AddEntry( X, Y, PreFac1A*GC[0] );

         for(Mu=0; Mu<NumGradientComponents; Mu++)
          if (GradB[Mu]) GradB[Mu]->AddEntry( X, Y, PreFac1A*GradGC[2*Mu+0]);

         for(Mu=0; Mu<NumTorqueAxes; Mu++)
          dBdTheta[Mu]->AddEntry( X, Y, PreFac1A*dGCdT[2*Mu+0]);

       }
      else if ( SaIsPEC && !SbIsPEC )
       { 
         X=RowOffset + nea;
         Y=ColOffset + 2*neb;  

         B->AddEntry( X, Y,   PreFac1A*GC[0] );
         B->AddEntry( X, Y+1, PreFac2A*GC[1] );

         for(Mu=0; Mu<NumGradientComponents; Mu++)
          { if (!GradB[Mu]) continue;
            GradB[Mu]->AddEntry( X, Y,   PreFac1A*GradGC[2*Mu+0]);
            GradB[Mu]->AddEntry( X, Y+1, PreFac2A*GradGC[2*Mu+1]);
          };

         for(Mu=0; Mu<NumTorqueAxes; Mu++)
          { dBdTheta[Mu]->AddEntry( X, Y, PreFac1A*dGCdT[2*Mu+0]);
            dBdTheta[Mu]->AddEntry( X, Y+1, PreFac2A*dGCdT[2*Mu+1]);
          };
       }
      else if ( !SaIsPEC && SbIsPEC )
       {
         X=RowOffset + 2*nea;
         Y=ColOffset + neb;  

         B->AddEntry( X,   Y, PreFac1A*GC[0] );
         B->AddEntry( X+1, Y, PreFac2A*GC[1] );

         for(Mu=0; Mu<NumGradientComponents; Mu++)
          { if (!GradB[Mu]) continue;
            GradB[Mu]->AddEntry( X, Y,   PreFac1A*GradGC[2*Mu+0]);
            GradB[Mu]->AddEntry( X+1, Y, PreFac2A*GradGC[2*Mu+1]);
          };

         for(Mu=0; Mu<NumTorqueAxes; Mu++)
          { dBdTheta[Mu]->AddEntry( X, Y,   PreFac1A*dGCdT[2*Mu+0]);
            dBdTheta[Mu]->AddEntry( X+1, Y, PreFac2A*dGCdT[2*Mu+1]);
          };
       }
      else if ( !SaIsPEC && !SbIsPEC )
       { 
         X=RowOffset + 2*nea;
         Y=ColOffset + 2*neb;

         B->AddEntry( X, Y,   PreFac1A*GC[0]);
         B->AddEntry( X, Y+1, PreFac2A*GC[1]);
         if ( !Symmetric || (nea!=neb) )
          B->AddEntry( X+1, Y, PreFac2A*GC[1]);
         B->AddEntry( X+1, Y+1, PreFac3A*GC[0]);

         for(Mu=0; Mu<NumGradientComponents; Mu++)
          { 
            if (!GradB[Mu]) continue;
            GradB[Mu]->AddEntry( X, Y,   PreFac1A*GradGC[2*Mu+0]);
            GradB[Mu]->AddEntry( X, Y+1, PreFac2A*GradGC[2*Mu+1]);
            if ( !Symmetric || (nea!=neb) )
             GradB[Mu]->AddEntry( X+1, Y, PreFac2A*GradGC[2*Mu+1]);
            GradB[Mu]->AddEntry( X+1, Y+1, PreFac3A*GradGC[2*Mu+0]);
          };

         for(Mu=0; Mu<NumTorqueAxes; Mu++)
          { 
            dBdTheta[Mu]->AddEntry( X, Y,   PreFac1A*dGCdT[2*Mu+0]);
            dBdTheta[Mu]->AddEntry( X, Y+1, PreFac2A*dGCdT[2*Mu+1]);
            if ( !Symmetric || (nea!=neb) )
             dBdTheta[Mu]->AddEntry( X+1, Y, PreFac2A*dGCdT[2*Mu+1]);
            dBdTheta[Mu]->AddEntry( X+1, Y+1, PreFac3A*dGCdT[2*Mu+0]);
          };

       }; // if ( OaIsPEC && ObIsPEC ) ... else ... 

      /*--------------------------------------------------------------*/
      /*- contributions of second medium if present.                  */
      /*- note this case we already know we are in the fourth case    */
      /*- of the above if...else statement.                           */
      /*--------------------------------------------------------------*/
      if (EpsB!=0.0)
       { 
         GetEEIArgs->k   = kB;
         GetEEIArgs->GBA = Args->GBA2;
         GetEdgeEdgeInteractions(GetEEIArgs);

         X=RowOffset + 2*nea;
         Y=ColOffset + 2*neb;

         B->AddEntry( X, Y,   PreFac1B*GC[0]);
         B->AddEntry( X, Y+1, PreFac2B*GC[1]);
         if ( !Symmetric || (nea!=neb) )
          B->AddEntry( X+1, Y, PreFac2B*GC[1]);
         B->AddEntry( X+1, Y+1, PreFac3B*GC[0]);

         for(Mu=0; Mu<NumGradientComponents; Mu++)
          { 
            if (!GradB[Mu]) continue;
            GradB[Mu]->AddEntry( X, Y,   PreFac1B*GradGC[2*Mu+0]);
            GradB[Mu]->AddEntry( X, Y+1, PreFac2B*GradGC[2*Mu+1]);
            if ( !Symmetric || (nea!=neb) )
             GradB[Mu]->AddEntry( X+1, Y, PreFac2B*GradGC[2*Mu+1]);
            GradB[Mu]->AddEntry( X+1, Y+1, PreFac3B*GradGC[2*Mu+0]);
          };

         for(Mu=0; Mu<NumTorqueAxes; Mu++)
          { 
            dBdTheta[Mu]->AddEntry( X, Y,   PreFac1B*dGCdT[2*Mu+0]);
            dBdTheta[Mu]->AddEntry( X, Y+1, PreFac2B*dGCdT[2*Mu+1]);
            if ( !Symmetric || (nea!=neb) )
             dBdTheta[Mu]->AddEntry( X+1, Y, PreFac2B*dGCdT[2*Mu+1]);
            dBdTheta[Mu]->AddEntry( X+1, Y+1, PreFac3B*dGCdT[2*Mu+0]);
          };
       }; // if (EpsB!=0.0)

    }; // for(nea=0; nea<NEa; nea++), for(neb=nebStart*nea; neb<NEb; neb++) ... 

  memcpy(TD->PPIAlgorithmCount, GetEEIArgs->PPIAlgorithmCount, NUMPPIALGORITHMS*sizeof(unsigned));
  return 0;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddSurfaceSigmaContributionToBEMMatrix(GetSSIArgStruct *Args)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if ( (Args->Sa != Args->Sb)    ) return;
  if ( !(Args->Sa->SurfaceSigma) && !(Args->Sa->SurfaceSigmaMP) ) return;
  if ( Args->Displacement        ) return;
  if ( Args->GBA1 || Args->GBA2  ) return;

  RWGSurface *S=Args->Sa;
 
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  HMatrix *B    = Args->B;
  int Offset    = Args->RowOffset;

  if (Offset!=Args->ColOffset)
   ErrExit("%s:%i: internal error",__FILE__,__LINE__);

  /*--------------------------------------------------------------*/
  /*- evaluate the surface conductivity for this surface at this -*/
  /*- frequency.                                                 -*/
  /*- The scaling by ZVac here is accounting for the SCUFF-EM    -*/
  /*- convention in which the upper half of the PMCHWT system    -*/
  /*- (the equations specifying tangential E-field continuity)   -*/
  /*- is divided by ZVac as part of the procedure for obtaining  -*/
  /*- a symmetric system of equations.                           -*/
  /*--------------------------------------------------------------*/
  cdouble GZ=Args->Sa->SurfaceSigmaMP->GetEps(Args->Omega);
  Log("Surface conductivity for surface %s is (%e,%e) at Omega=%e",
       Args->Sa->Label,real(GZ),imag(GZ),real(Args->Omega));
  GZ*=ZVAC;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int neAlpha, neBeta;
  double Overlap;
  for(neAlpha=0; neAlpha<S->NumEdges; neAlpha++)
   for(neBeta=neAlpha; neBeta<S->NumEdges; neBeta++)
    { 
      Overlap=S->GetOverlap(neAlpha, neBeta);
      if (Overlap==0.0) continue;

#if 0
#20130426 this code figures out the centroid of the panel in 
          question so we can evaluate the surface conductivity there;
          not needed for the time being as we are using spatially-
          constant surface conductivity. 

      char *SSParmNames[4]={ const_cast<char *>("w"), const_cast<char *>("x"), 
                         const_cast<char *>("y"), const_cast<char *>("z") };
      cdouble SSParmValues[4];
      SSParmValues[0]=Args->Omega*MatProp::FreqUnit;

      EAlpha=S->Edges[neAlpha];
      EBeta=S->Edges[neBeta];

      // if there was a nonzero overlap, get the value 
      // of the surface conductivity at the centroid
      // of the common panel (if there was only one common panel)
      // or of the common edge if there were two common panels.
      if (neAlpha==neBeta)
       { SSParmValues[1] = EAlpha->Centroid[0];
         SSParmValues[2] = EAlpha->Centroid[1];
         SSParmValues[3] = EAlpha->Centroid[2];
       }
      else if ( EAlpha->iPPanel==EBeta->iPPanel || EAlpha->iPPanel==EBeta->iMPanel )
       { P=S->Panels[EAlpha->iPPanel];
         SSParmValues[1] = P->Centroid[0];
         SSParmValues[2] = P->Centroid[1];
         SSParmValues[3] = P->Centroid[2];
       }
      else if ( (EAlpha->iMPanel!=-1) && ( (EAlpha->iMPanel==EBeta->iPPanel) || (EAlpha->iMPanel==EBeta->iMPanel) ) ) 
       { P=S->Panels[EAlpha->iMPanel];
         SSParmValues[1] = P->Centroid[0];
         SSParmValues[2] = P->Centroid[1];
         SSParmValues[3] = P->Centroid[2];
       };

      GZ=cevaluator_evaluate(S->SurfaceSigma, 4, SSParmNames, SSParmValues);
#endif

      if ( S->IsPEC )
       { B->AddEntry(Offset+neAlpha, Offset+neBeta, -1.0*Overlap/GZ);
         if (neAlpha!=neBeta)
          B->AddEntry(Offset+neBeta, Offset+neAlpha, -1.0*Overlap/GZ);
       }
      else
       { B->AddEntry(Offset + 2*neAlpha+1, Offset + 2*neBeta+1, +GZ*Overlap);
         if (neAlpha!=neBeta)
          B->AddEntry(Offset + 2*neBeta+1, Offset + 2*neAlpha+1, +GZ*Overlap);
       }
      
    };

}

/***************************************************************/  
/***************************************************************/  
/***************************************************************/
void GetSurfaceSurfaceInteractions(GetSSIArgStruct *Args)
{ 
  RWGGeometry *G = Args->G;
  cdouble Omega = Args->Omega;
  RWGSurface  *Sa = Args->Sa;
  RWGSurface  *Sb = Args->Sb;

  /*--------------------------------------------------------------*/
  /*- make sure the cached epsilon and mu values for all regions  */
  /*- are up-to-date for the present frequency                    */
  /*--------------------------------------------------------------*/
  G->UpdateCachedEpsMuValues(Omega);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if ( Args->Accumulate==false )
   { Args->B->ZeroBlock(Args->RowOffset, Sa->NumBFs, Args->ColOffset, Sb->NumBFs);
     if (Args->GradB && Args->GradB[0])
      Args->GradB[0]->ZeroBlock(Args->RowOffset, Sa->NumBFs, Args->ColOffset, Sb->NumBFs);
     if (Args->GradB && Args->GradB[1])
      Args->GradB[1]->ZeroBlock(Args->RowOffset, Sa->NumBFs, Args->ColOffset, Sb->NumBFs);
     if (Args->GradB && Args->GradB[2])
      Args->GradB[2]->ZeroBlock(Args->RowOffset, Sa->NumBFs, Args->ColOffset, Sb->NumBFs);
     if (Args->dBdTheta && Args->dBdTheta[0])
      Args->dBdTheta[0]->ZeroBlock(Args->RowOffset, Sa->NumBFs, Args->ColOffset, Sb->NumBFs);
     if (Args->dBdTheta && Args->dBdTheta[1])
      Args->dBdTheta[1]->ZeroBlock(Args->RowOffset, Sa->NumBFs, Args->ColOffset, Sb->NumBFs);
     if (Args->dBdTheta && Args->dBdTheta[2])
      Args->dBdTheta[2]->ZeroBlock(Args->RowOffset, Sa->NumBFs, Args->ColOffset, Sb->NumBFs);
   };

  /***************************************************************/
  /* figure out if the two surfaces have 0, 1, or 2 regions in   */
  /* common and fill in some relevant fields in the Args struct  */
  /***************************************************************/
  double Signs[2];
  int CommonRegions[2]; 
  int NumCommonRegions=CountCommonRegions(Sa, Sb, CommonRegions, Signs);
  if (NumCommonRegions==0)
   return;

  Args->EpsA  = G->EpsTF[ CommonRegions[0] ];
  Args->MuA   = G->MuTF[  CommonRegions[0] ];
  Args->SignA = Signs[0];
  if ( NumCommonRegions==2 )
   { Args->EpsB = G->EpsTF[ CommonRegions[1] ];
     Args->MuB  = G->MuTF[  CommonRegions[1] ];
     Args->SignB = Signs[1];
   }
  else
   Args->EpsB = Args->MuB = Args->SignB = 0.0;

  Args->SaIsPEC = (Sa->IsPEC==1);
  Args->SbIsPEC = (Sb->IsPEC==1);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (Args->OmitRegion1)
   Args->EpsA=Args->MuA=Args->SignA=0.0;

  if (Args->OmitRegion2)
   Args->EpsB=Args->MuB=Args->SignB=0.0;

  if ( Args->EpsA==0.0 && Args->EpsB==0.0 )
   return;

  /***************************************************************/
  /* fire off threads ********************************************/
  /***************************************************************/
  GlobalFIPPICache.Hits=GlobalFIPPICache.Misses=0;

  int nt, NumTasks, NumThreads = GetNumThreads();
  unsigned PPIAlgorithmCount[NUMPPIALGORITHMS];  
  memset(PPIAlgorithmCount, 0, NUMPPIALGORITHMS*sizeof(unsigned));

#ifdef USE_PTHREAD
  ThreadData *TDs = new ThreadData[NumThreads], *TD;
  pthread_t *Threads = new pthread_t[NumThreads];
  for(nt=0; nt<NumThreads; nt++)
   { 
     TD=&(TDs[nt]);
     TD->nt=nt;
     TD->NumTasks=NumThreads;
     TD->Args=Args;
     if (nt+1 == NumThreads)
       GSSIThread((void *)TD);
     else
       pthread_create( &(Threads[nt]), 0, GSSIThread, (void *)TD);
   }
  for(nt=0; nt<NumThreads-1; nt++)
   { pthread_join(Threads[nt],0);
     for(int n=0; n<NUMPPIALGORITHMS; n++)
      PPIAlgorithmCount[n] += TD->PPIAlgorithmCount[n];
   };
  delete[] Threads;
  delete[] TDs;

#else 
#ifndef USE_OPENMP
  NumTasks=NumThreads=1;
  if (G->LogLevel>=SCUFF_VERBOSE2)
   Log(" no multithreading...");
#else
  NumTasks=NumThreads*100;
  if (NumTasks>Sa->NumEdges) 
   NumTasks=Sa->NumEdges;
  if (G->LogLevel>=SCUFF_VERBOSE2)
   Log(" OpenMP multithreading (%i threads,%i tasks)...",NumThreads,NumTasks);
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(nt=0; nt<NumTasks; nt++)
   { 
     ThreadData TD1;
     TD1.nt=nt;
     TD1.NumTasks=NumTasks;
     TD1.Args=Args;
     GSSIThread((void *)&TD1);
     for(int n=0; n<NUMPPIALGORITHMS; n++)
      PPIAlgorithmCount[n] += TD1.PPIAlgorithmCount[n];
   };
#endif

  if (G->LogLevel>=SCUFF_VERBOSE2)
   { Log("  %i/%i cache hits/misses",GlobalFIPPICache.Hits,GlobalFIPPICache.Misses);
     Log("  PPIs: LOC(%u), HOC(%u), TD(%u), HK(%u), D(%u)",
            PPIAlgorithmCount[PPIALG_LOCUBATURE],
            PPIAlgorithmCount[PPIALG_HOCUBATURE],
            PPIAlgorithmCount[PPIALG_TD],
            PPIAlgorithmCount[PPIALG_HKTD],
            PPIAlgorithmCount[PPIALG_DESING]);
   };

  /***************************************************************/
  /* 20120526 handle objects with finite surface conductivity    */
  /***************************************************************/
  if ( (Args->Sa == Args->Sb) && (Args->Sa->SurfaceSigma!=0 || Args->Sa->SurfaceSigmaMP!=0) )
   AddSurfaceSigmaContributionToBEMMatrix(Args);

  /***************************************************************/
  /* if the caller specified the matrix as symmetric, then so far*/
  /* we have only computed the upper triangle, so now we need to */
  /* fill in the lower triangle. The exception is if the matrix  */
  /* uses packed storage, in which case only the upper triangle  */
  /* is stored anyway.                                           */
  /***************************************************************/
  if ( Args->Symmetric && (Args->B->StorageType==LHM_NORMAL) )
   { int N=Args->B->NR;
     for(int nr=1; nr<N; nr++)
      for(int nc=0; nc<nr; nc++)
       Args->B->SetEntry(nr,nc,Args->B->GetEntry(nc,nr));
   };

}

/***************************************************************/
/* initialize an argument structure for                        */
/* GetSurfaceSurfaceInteractions.                              */
/*                                                             */
/* note: what this routine does is to fill in default values   */
/* for some of the lesser-used fields in the structure, while  */
/* leaving several other fields uninitialized; any caller of   */
/* AssembleBEMMatrixBlock must fill in those remaining fields  */
/* before the call.                                            */
/*                                                             */
/* Note: If Accumulate is set to true before the call to       */
/*       AssembleBEMMatrixBlock(), then the entries of the     */
/*       HMatrix will be augmented, not replaced, by the new   */
/*       matrix entries. As far as I know this is only ever    */
/*       useful for accumulating the contributions of          */
/*       neighboring lattice cells when assembling the BEM     */
/*       matrix with periodic boundary conditions.             */
/***************************************************************/
void InitGetSSIArgs(GetSSIArgStruct *Args)
{
  Args->NumTorqueAxes=0;
  Args->GammaMatrix=0;
  
  Args->RowOffset=0;
  Args->ColOffset=0;

  Args->Symmetric=false;

  Args->GBA1=Args->GBA2=0;

  Args->OmitRegion1=false;
  Args->OmitRegion2=false;

  Args->Displacement = 0;

  Args->GradB=0;
  Args->dBdTheta=0;

  Args->Accumulate=false;

}

} // namespace scuff
