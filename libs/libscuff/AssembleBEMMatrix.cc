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
 * AssembleBEMMatrix.cc -- libscuff routines for assembling the BEM matrix.
 *                      --
 *                      -- (cf. 'libscuff Implementation and Technical
 *                      --  Details', section 8.3, 'Structure of the BEM
 *                      --  Matrix.')
 *                      --
 * homer reid           -- 10/2006 -- 10/2011
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libhmat.h>
#include <libhrutil.h>
#include <sys/stat.h>
#include <unistd.h>

#include "libscuff.h"
#include "libscuffInternals.h"

#include "cmatheval.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#ifdef USE_OPENMP
#  include <omp.h>
#endif

#define II cdouble(0,1)

namespace scuff {

HMatrix *AssembleBEMMatrix2018(RWGGeometry *G, cdouble Omega, double *kBloch, HMatrix *M);

/***************************************************************/
/* This routine adds the contents of matrix block B, which has */
/* dimensions NRxNC, to the block of M whose upper-left corner */
/* has indices (RowOffset, ColOffset).                         */
/*                                                             */
/* Each entry of B is scaled by BPF ('bloch phase factor')     */
/* before being added to the corresponding entry of M.         */ 
/*                                                             */
/* If UseSymmetry is true, then the routine additionally adds  */
/* the contents of B' times the complex conjugate of BPF to    */
/* the destination block of M.                                 */
/*                                                             */
/* If GradB[Mu] (Mu=0,1,2) is non-null, the same stamping      */
/* operation is used to stamp GradB[Mu] into GradM[Mu].        */
/***************************************************************/
void StampInNeighborBlock(HMatrix *B, HMatrix **GradB,
                          int NR, int NC,
                          HMatrix *M, HMatrix **GradM,
                          int RowOffset, int ColOffset,
                          double L[2], double kBloch[2], bool UseSymmetry)
{ 
  HMatrix *BList[4];
  HMatrix *MList[4];

  // bloch phase factor 
  cdouble BPF=exp( II*(kBloch[0]*L[0] + kBloch[1]*L[1]) );

  BList[0] = B;
  BList[1] = GradB ? GradB[0] : 0;
  BList[2] = GradB ? GradB[1] : 0;
  BList[3] = GradB ? GradB[2] : 0;

  MList[0] = M;
  MList[1] = GradM ? GradM[0] : 0;
  MList[2] = GradM ? GradM[1] : 0;
  MList[3] = GradM ? GradM[2] : 0;
  
  for(int n=0; n<4; n++)
   { 
     HMatrix *BB = BList[n]; 
     HMatrix *MM = MList[n]; 
     if ( !BList[n] || !MList[n] )
      continue;
     
     if (UseSymmetry)
      { for(int nr=0; nr<NR; nr++)
         for(int nc=0; nc<NC; nc++)
          MM->AddEntry(RowOffset + nr, ColOffset + nc,
                       BPF*BB->GetEntry(nr,nc) + conj(BPF)*BB->GetEntry(nc,nr));
      }
     else
      { for(int nr=0; nr<NR; nr++)
         for(int nc=0; nc<NC; nc++)
          MM->AddEntry(RowOffset + nr, ColOffset + nc, 
                       BPF*BB->GetEntry(nr,nc));
      };
   };
}

/***************************************************************/
/* this is a simple quick-and-dirty mechanism for saving       */
/* and retrieving diagonal T-matrix blocks to disk.            */
/***************************************************************/
#define TBCOP_READ  0
#define TBCOP_WRITE 1
bool TBlockCacheOp(int Op, RWGGeometry *G, int ns,
                   cdouble Omega, double *kBloch,
                   HMatrix *M, int RowOffset, int ColOffset)
{
  char *Dir = getenv("SCUFF_TBLOCK_PATH");
  if (Dir==0 && Op==TBCOP_READ)
   Dir = getenv("SCUFF_TBLOCK_READPATH");
  if (!Dir) return false;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  char K4VStr[50];
  if (imag(Omega)==0.0)
   sprintf(K4VStr,"%.6e",real(Omega));
  else if (real(Omega)==0.0)
   sprintf(K4VStr,"%.6eI",imag(Omega));
  else 
   sprintf(K4VStr,"%.6e+%.6eI",real(Omega),imag(Omega));
  if (G->LDim>=1)
   vstrappend(K4VStr,"_%.6e",kBloch[0]);
  if (G->LDim>=2)
   vstrappend(K4VStr,"_%.6e",kBloch[1]);

  RWGSurface *S=G->Surfaces[ns];
  int nr1=S->RegionIndices[0];
  int nr2=S->RegionIndices[1];
  if (G->RegionMPs[nr1]->Zeroed)
   strcat(K4VStr,"_Interior");
  else if (nr2!=-1 && G->RegionMPs[nr2]->Zeroed)
   strcat(K4VStr,"_Exterior");
    
  char *FileBase = GetFileBase(S->MeshFileName); 
  int MeshTag    = S->MeshTag;
  char FileName[200];
  if (MeshTag != -1)
   snprintf(FileName,200,"%s/%s_%i_%s.hdf5",Dir,FileBase,MeshTag,K4VStr);
  else
   snprintf(FileName,200,"%s/%s_%s.hdf5",Dir,FileBase,K4VStr);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NBF        = G->Surfaces[ns]->NumBFs;
  bool IsComplex = (M->RealComplex==LHM_COMPLEX);
  off_t DataSize = NBF*NBF*(IsComplex ? sizeof(double) : sizeof(cdouble) );

  /***************************************************************/
  /* if the block we are reading/writing is a proper subblock of */
  /* M, then we will need to allocate an HMatrix for just that   */
  /* block, but we hold off on doing that until necessary        */
  /***************************************************************/
  HMatrix *B=0;
  bool OwnsB=false;
  if ( RowOffset==0 && ColOffset==0 && (M->NR==NBF) && (M->NC==NBF) )
   B=M;

  bool Success=false;
  FILE *f=0;
  if (Op==TBCOP_READ)	
   {
     Log("Attempting to read T-block (%s,%s) from file %s...",FileBase,z2s(Omega),FileName);
     struct stat st;
     int Status=stat(FileName, &st);
     f=fopen(FileName,"r");
     if (Status!=0 || f==0)
      Log("...could not access file");
     else if (st.st_size!=DataSize)
      Log("...file had incorrect size %lu (should be %lu)",st.st_size,DataSize);
     else 
      { 
        if (B==0)
         { OwnsB=true;
           B=new HMatrix(NBF,NBF,M->RealComplex);
           if (B->ErrMsg)
            { Log("...could not allocate matrix (%s)",B->ErrMsg);
              delete B;
              fclose(f);
              return false;
            };
         };

        void *DataBuffer = (IsComplex ? ((void *)(B->ZM)) : ((void *)(B->DM)));
        ssize_t SizeToRead = (ssize_t)DataSize;
        ssize_t SizeRead = read(fileno(f), DataBuffer, SizeToRead);
        Success = (SizeRead==SizeToRead);
        if (Success)
         { Log("...success!");
           if (OwnsB)
            M->InsertBlock(B, RowOffset, ColOffset);
         }
        else
         Log("...failed to read full size (%lu)",SizeToRead);
      };
   }
  else //(Op==TBCOP_WRITE)
   { 
     Log("Attempting to write T-block (%s,%s) to file %s...",FileBase,z2s(Omega),FileName);
     if (B==0)
      { OwnsB=true;
        B=new HMatrix(NBF,NBF,M->RealComplex);
        if (B->ErrMsg)
         { Log("...could not allocate matrix (%s)",B->ErrMsg);
           delete B;
           return false;
         };
        M->ExtractBlock(RowOffset, ColOffset, B);
      };
     void *DataBuffer = (IsComplex ? ((void *)(B->ZM)) : ((void *)(B->DM)));
     f=fopen(FileName,"w");
     if (!f) 
      Log("...could not access file");
     else if (DataSize != write(fileno(f), DataBuffer, DataSize))
      Log("...failed to write full file");
     else
      { Log("...success!");
        Success=true;
      };
   };

  if (B && OwnsB) delete B;
  if (f) fclose(f);
  return Success; 
}

/***************************************************************/
/* KBIMBCache = 'kBloch-independent matrix-block cache.'       */
/***************************************************************/
typedef struct KBMIMBCache
 {
   cdouble Omega;
   int NumMatrices;
   bool NeedZDerivative;
   void *Storage;
   HMatrix *B[9], *dBdZ[9];

 } KBIMBCache;

void *RWGGeometry::CreateABMBAccelerator(int nsa, int nsb,
                                         bool PureImagFreq,
                                         bool NeedZDerivative)
{
  if (LBasis==0)
   return 0;

  /*--------------------------------------------------------------*/
  /*- gather some information about the problem ------------------*/
  /*--------------------------------------------------------------*/
  bool OneDLattice = (LDim==1);
  bool TwoDLattice = !OneDLattice;
  bool SameSurface = (nsa==nsb);
  int NR = Surfaces[nsa]->NumBFs;
  int NC = Surfaces[nsb]->NumBFs;

  size_t MatrixSize = PureImagFreq ? NR*NC*sizeof(double) : NR*NC*sizeof(cdouble);
  int RC = PureImagFreq ? LHM_REAL : LHM_COMPLEX;

  int NumMatrices=0;
  if      (  OneDLattice &&  SameSurface ) NumMatrices=2;
  else if (  OneDLattice && !SameSurface ) NumMatrices=3;
  else if (  TwoDLattice &&  SameSurface ) NumMatrices=5;
  else if (  TwoDLattice && !SameSurface ) NumMatrices=9;

  /*--------------------------------------------------------------*/
  /*- attempt to allocate enough storage for the full cache       */
  /*--------------------------------------------------------------*/
  Log("Trying to allocate accelerator for (%i,%i) matrix-block assembly...",nsa,nsb);
  int MatricesNeeded = NeedZDerivative ? 2*NumMatrices : NumMatrices;
  void *Storage = malloc( MatricesNeeded * MatrixSize );
  if (Storage==0) 
   { Log("...failed! not enough memory.");
     return 0;
   };
  Log("...success!");

  void *Buffers[18];
  Buffers[0] = Storage;
  for (int nb=1; nb<18; nb++)
   Buffers[nb] = (void *)( (char *)Buffers[nb-1] + MatrixSize);

  /*--------------------------------------------------------------*/
  /*- only if that succeeded, proceed to allocate the actual cache*/
  /*--------------------------------------------------------------*/
  KBIMBCache *Cache = (KBIMBCache *)mallocEC( sizeof *Cache );
  Cache->Omega           = -1.0;            // initially dirty
  Cache->NumMatrices     = NumMatrices;
  Cache->NeedZDerivative = NeedZDerivative;
  Cache->Storage         = Storage;
  for(int nm=0, nb=0; nm<NumMatrices; nm++)
   { Cache->B[nm] = new HMatrix(NR, NC, RC, LHM_NORMAL, Buffers[nb++]);
     if (NeedZDerivative)
      Cache->dBdZ[nm] = new HMatrix(NR, NC, RC, LHM_NORMAL, Buffers[nb++]);
   };

  return (void *)Cache;
 
}

void RWGGeometry::DestroyABMBAccelerator(void *pCache)
{
  if (pCache==0) return;

  KBIMBCache *Cache    = (KBIMBCache *)pCache;
  int NumMatrices      = Cache->NumMatrices;
  bool NeedZDerivative = Cache->NeedZDerivative;
  for(int nm=0; nm<NumMatrices; nm++)
   { delete Cache->B[nm];
     if (NeedZDerivative) delete Cache->dBdZ[nm];
   };
  free(Cache->Storage);
  free(Cache);

}

/***************************************************************/
/* This routine computes the block of the BEM matrix that      */
/* describes the interaction between surfaces nsa and nsb.     */
/* This block is stamped into M in such a way that the upper-  */ 
/* left element of the block is at the (RowOffset, ColOffset)  */ 
/* entry of M. If GradM is non-null and GradM[Mu] is non-null  */ 
/* (Mu=0,1,2) then the X_{Mu} derivative of BEM matrix is      */ 
/* similarly stamped into GradM[Mu]. If NumTorqueAxes>0 and    */ 
/* dMdT and GammaMatrix are non-null, then the derivative of   */ 
/* M with respect to rotation angle Theta about the Muth torque*/ 
/* axis described by GammaMatrix (Mu=0,...,NumTorqueAxes-1) is */ 
/* similarly stamped into dMdT[Mu].                             */ 
/***************************************************************/
void RWGGeometry::AssembleBEMMatrixBlock(int nsa, int nsb,
                                         cdouble Omega, double *kBloch,
                                         HMatrix *M, HMatrix **GradM,
                                         int RowOffset, int ColOffset,
                                         void *Accelerator, bool TransposeAccelerator,
                                         int NumTorqueAxes, HMatrix **dMdT,
                                         double *GammaMatrix)
{
  if (TransposeAccelerator)
   ErrExit("%s:%i: TransposeAccelerator not implemented");

  if (    nsa==nsb
       && GradM==0
       && TBlockCacheOp(TBCOP_READ, this, nsa, Omega, kBloch, M, RowOffset, ColOffset)
     ) return;

  if (LogLevel>=SCUFF_VERBOSELOGGING)
   Log("Assembling BEM matrix block (%i,%i)",nsa,nsb);

  /***************************************************************/
  /* handle the compact-object case first since it is so simple  */
  /***************************************************************/
  if (LBasis==0)
   {  
     GetSSIArgStruct GetSSIArgs, *Args=&GetSSIArgs;
     InitGetSSIArgs(Args);
     Args->G=this;
     Args->Sa=Surfaces[nsa];
     Args->Sb=Surfaces[nsb];
     Args->Omega=Omega;
     Args->NumTorqueAxes=NumTorqueAxes;
     Args->GammaMatrix=GammaMatrix;
     Args->Symmetric = (nsa==nsb);
     Args->B=M;
     Args->GradB=GradM;
     Args->dBdTheta=dMdT;
     Args->RowOffset=RowOffset;
     Args->ColOffset=ColOffset;
     GetSurfaceSurfaceInteractions(Args);
     if (nsa==nsb)
      TBlockCacheOp(TBCOP_WRITE, this, nsa, Omega, kBloch, M, RowOffset, ColOffset);
     return;
   }

  /***************************************************************/
  /* The remainder of this routine is now for the PBC case only, */
  /* and it consists of two main steps: (a) assemble the kBloch- */
  /* independent contributions of the innermost grid cells and   */
  /* stamp them appropriately into the matrix; then (b) add the  */
  /* contributions of outer grid cells.                          */
  /***************************************************************/
  if ( NumTorqueAxes>0 )
   ErrExit("angular derivatives of BEM matrix not supported for periodic geometries");
  if ( GradM && (GradM[0] || GradM[1]) )
   ErrExit("x,y derivatives of BEM matrix not supported for periodic geometries");

  KBIMBCache *Cache = (KBIMBCache *)Accelerator;
  bool HaveCache = (Cache!=0);
  bool HaveCleanCache = HaveCache && EqualFloat(Cache->Omega, Omega);
  if (HaveCache) Cache->Omega=Omega;

  int NumCommonRegions, CRIndices[2];
  double Signs[2];
  NumCommonRegions=CountCommonRegions(Surfaces[nsa], Surfaces[nsb], CRIndices, Signs);
  if (NumCommonRegions==0) 
   return;
  int nr1=CRIndices[0];
  int nr2=NumCommonRegions==2 ? CRIndices[1] : -1;

  int NBFA=Surfaces[nsa]->NumBFs;
  int NBFB=Surfaces[nsb]->NumBFs;

  bool UseSymmetry = (nsa==nsb);

  double L[3]={0.0, 0.0, 0.0};

  bool OneDLattice = (LDim==1);
  double LBV[3][3];
  LBV[0][0] = LBasis->GetEntryD(0,0);
  LBV[0][1] = LBasis->GetEntryD(1,0);
  LBV[0][2] = LBasis->GetEntryD(2,0);
  LBV[1][0] = LDim > 1 ? LBasis->GetEntryD(0,1) : 0.0;
  LBV[1][1] = LDim > 1 ? LBasis->GetEntryD(1,1) : 0.0;
  LBV[1][2] = LDim > 1 ? LBasis->GetEntryD(2,1) : 0.0;

  /***************************************************************/
  /* pre-initialize arguments for GetSurfaceSurfaceInteractions **/
  /***************************************************************/
  GetSSIArgStruct GetSSIArgs, *Args=&GetSSIArgs;
  InitGetSSIArgs(Args);

  Args->G            = this;
  Args->Sa           = Surfaces[nsa];
  Args->Sb           = Surfaces[nsb];
  Args->Omega        = Omega;
  Args->GBA1         = 0;
  Args->GBA2         = 0;
  Args->Accumulate   = false;
  Args->Displacement = L;

  /*--------------------------------------------------------------*/
  /*- If the caller didn't provide a cache, we need temporary     */
  /*- storage for the kBloch-independent matrix blocks.           */
  /*--------------------------------------------------------------*/
  HMatrix *GradBBuffer[3]={0, 0, 0};
  Args->GradB = (GradM && GradM[2]) ? GradBBuffer : 0;
  if (!HaveCache)
   { 
     Args->B = new HMatrix(NBFA, NBFB, LHM_COMPLEX);
     if (Args->GradB)
      Args->GradB[2] = new HMatrix(NBFA, NBFB, LHM_COMPLEX);
   };

  /***************************************************************/
  /* Assemble and stamp in contributions of innermost grid cells.*/
  /***************************************************************/
  M->ZeroBlock(RowOffset, NBFA, ColOffset, NBFB);
  if (GradM && GradM[2]) GradM[2]->Zero();
  if (LogLevel>=SCUFF_VERBOSELOGGING)
   Log(" Step 1: Contributions of innermost grid cells...");
  for(int n1=+1, nb=0; n1>=-1; n1--)
   for(int n2=+1; n2>=-1; n2--)
    { 
      if ( OneDLattice && n2!=0 ) continue;

      L[0] = n1*LBV[0][0] + n2*LBV[1][0];
      L[1] = n1*LBV[0][1] + n2*LBV[1][1];
 
      if (HaveCache)
       { Args->B = Cache->B[ nb ];
         if (Args->GradB) Args->GradB[2] = Cache->dBdZ[ nb ];
         nb++;
       };

      if ( !HaveCleanCache )
       { 
         // detect extendedness of regions and omit contributions of
         // any regions that are not extended
         if ( n1==0 && n2==0 )
          { Args->OmitRegion1 = false;
            Args->OmitRegion2 = (nr2==-1);
          }
         else 
          { Args->OmitRegion1 = Args->OmitRegion2 = true;
            if ( n1!=0 && RegionIsExtended[0][nr1] ) Args->OmitRegion1=false;
            if ( n1!=0 && nr2!=-1 && RegionIsExtended[0][nr2] ) Args->OmitRegion2=false;
            if ( n2!=0 && RegionIsExtended[1][nr1] ) Args->OmitRegion1=false;
            if ( n2!=0 && nr2!=-1 && RegionIsExtended[1][nr2] ) Args->OmitRegion2=false;
          };
 
         if (LogLevel>=SCUFF_VERBOSE2)
          Log("  ...(%i,%i) block...",n1,n2);
         Args->Symmetric = (nsa==nsb && n1==0 && n2==0); 
         GetSurfaceSurfaceInteractions(Args);
       };

      StampInNeighborBlock(Args->B, Args->GradB, NBFA, NBFB,
                           M, GradM, RowOffset, ColOffset, L, kBloch, 
                           UseSymmetry && !(n1==0 && n2==0) );

      if ( UseSymmetry && (n1==0 && n2==0) )
       goto done; // want break, but need to break out of both loops
    };

done: 

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (!HaveCache)
   { delete Args->B;
     if (Args->GradB) delete Args->GradB[2];
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (LogLevel>=SCUFF_VERBOSELOGGING)
   Log(" Step 2: Contributions of outer grid cells...");
  Args->Displacement = 0;
  Args->Symmetric    = false;
  Args->Accumulate   = true;
  Args->B            = M;
  Args->GradB        = GradM;
  Args->RowOffset    = RowOffset;
  Args->ColOffset    = ColOffset;
  Args->OmitRegion1  = false;
  Args->OmitRegion2  = (nr2==-1);

  Args->GBA1=CreateRegionGBA(nr1, Omega, kBloch, nsa, nsb);
  if (Args->GBA1==0)
   { if (LogLevel>=SCUFF_VERBOSELOGGING)
      Log("Skipping interpolation table for region %i",nr1);
     Args->OmitRegion1  = true;
   }

  if ( !(Args->OmitRegion2) )
   { Args->GBA2=CreateRegionGBA(nr2, Omega, kBloch, nsa, nsb);
     if (Args->GBA2==0)
      { if (LogLevel>=SCUFF_VERBOSELOGGING)
         Log("Skipping interpolation table for region %i",nr2);
        Args->OmitRegion2  = true;
      }
   }
  else
   Args->GBA2=0;

  if ( Args->OmitRegion1==false || Args->OmitRegion2==false )
   GetSurfaceSurfaceInteractions(Args);

  if (Args->GBA1) DestroyGBarAccelerator(Args->GBA1);
  if (Args->GBA2) DestroyGBarAccelerator(Args->GBA2);

  if (nsa==nsb)
   TBlockCacheOp(TBCOP_WRITE, this, nsa, Omega, kBloch, M, RowOffset, ColOffset);

}

/***************************************************************/
/* determine whether or not edge #ne on surface #ns is part of */
/* a multi-material junction.                                  */
/* if not, return false.                                       */
/* if so, return true and set                                  */
/*  *pnMMJ           = index of MMJ                            */
/*  *pnEdgeWithinMMJ = index of edge within MMJ                */
/***************************************************************/
bool EdgeInMMJ(RWGGeometry *G, int ns, int ne, int *pnMMJ, int *pnEdgeWithinMMJ)
{
  for(int nMMJ=0; nMMJ<G->NumMMJs; nMMJ++)
   { MMJData *MMJ = G->MultiMaterialJunctions[nMMJ];
     int MMJSize         = MMJ->NumEdges;
     int *SurfaceIndices = MMJ->SurfaceIndices;
     int *EdgeIndices    = MMJ->EdgeIndices;   

     for(int n=0; n < MMJSize; n++)
      if ( SurfaceIndices[n]==ns && EdgeIndices[n]==ne )
       { if (pnMMJ) *pnMMJ=nMMJ; 
         if (pnEdgeWithinMMJ) *pnEdgeWithinMMJ=n;
         return true;
       };
   };
  return false;
   
}
/***************************************************************/
/***************************************************************/
/***************************************************************/
#define MMJALG_NONE  0
#define MMJALG_MINUS 1
#define MMJALG_PLUS  2
void RWGGeometry::ApplyMMJTransformation(HMatrix *M, HVector *RHS)
{
  int Method=MMJALG_NONE;

  char *s=getenv("SCUFF_MMJ_ALGORITHM");
  if (s)
   Method = s[0] - '0';
  if (Method<0 || Method>2)
   { Warn("invalid MMJ algorithm %s (skipping)",s);
     return;
   };

  if (Method==MMJALG_NONE) return;

  Log("Applying MMJ transformation via algorithm %i.",Method);

  double Sign = (Method==MMJALG_MINUS) ? -1.0 : 1.0;
  for(int nMMJ=0; nMMJ<NumMMJs; nMMJ++)
   { 
     MMJData *Data       = MultiMaterialJunctions[nMMJ];
     int MMJSize         = Data->NumEdges;
     int *SurfaceIndices = Data->SurfaceIndices;
     int *EdgeIndices    = Data->EdgeIndices;

     int ns0             = SurfaceIndices[0];
     int ne0             = EdgeIndices[0];
     bool IsPEC          = Surfaces[ns0]->IsPEC;
     int iBF0            = BFIndexOffset[ns0] + (IsPEC ? 1 : 2)*ne0;
     
     // subtract/add BEM system row #iBF0 from rows #iBF1...#iBFN
     for(int n=1; n<MMJSize; n++)
      { 
         int nsN     = SurfaceIndices[n];
         int neN     = EdgeIndices[n];
         int iBFN    = BFIndexOffset[nsN] + (IsPEC ? 1 : 2)*neN;

         if (M)
          for(int nc=0; nc<M->NC; nc++)
           { M->AddEntry(iBFN, nc, Sign*M->GetEntry(iBF0, nc));
             if (!IsPEC)
              M->AddEntry(iBFN+1, nc, Sign*M->GetEntry(iBF0+1, nc));
           };
         if (RHS)
          { RHS->AddEntry(iBFN, Sign*RHS->GetEntry(iBF0));
            if ( !IsPEC )
             RHS->AddEntry(iBFN+1, Sign*RHS->GetEntry(iBF0+1));
          };
      }; // for(int n=1; n<MMJSize; n++)

     // replace BEM system row #iBF0 with condition \sum K_n = 0
     if (M)
      { for(int nc=0; nc<M->NC; nc++)
         { M->SetEntry(iBF0, nc, 0.0);
           if (!IsPEC) M->SetEntry(iBF0+1, nc, 0.0);
         };
        for(int n=0; n<MMJSize; n++)
         { 
           int ns  = SurfaceIndices[n];
           int ne  = EdgeIndices[n];
           int iBF = BFIndexOffset[ns] + (IsPEC ? 1 : 2)*ne;
           M->SetEntry(iBF0, iBF, 1.0);
           if (!IsPEC) 
            M->SetEntry(iBF0+1, iBF+1, 1.0);
         };
      }; 

     if (RHS)
      { RHS->SetEntry(iBF0, 0.0);
        if (!IsPEC) 
         RHS->SetEntry(iBF0+1,0.0);
      };

   }; // for(int nMMJ=0; nMMJ<NumMMJs; nMMJ++)

}

/***************************************************************/
/* this is the actual API-exposed routine for assembling the   */
/* BEM matrix, which is pretty simple and really just calls    */
/* AssembleBEMMatrixBlock() to do all the dirty work.          */
/*                                                             */
/* If the M matrix is NULL on entry, a new HMatrix of the      */
/* appropriate size is allocated and returned. Otherwise, the  */
/* return value is M.                                          */
/***************************************************************/
HMatrix *RWGGeometry::AssembleBEMMatrix(cdouble Omega, double *kBloch, HMatrix *M)
{ 
  if (CheckEnv("SCUFF_MATRIX_2018") && LDim==0 )
   return AssembleBEMMatrix2018(this, Omega, kBloch, M);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if ( LBasis==0 && kBloch!=0 && (kBloch[0]!=0.0 || kBloch[1]!=0.0) )
   ErrExit("%s:%i: Bloch wavevector is undefined for compact geometries");
  if ( LBasis!=0 && kBloch==0 )
   ErrExit("%s:%i: Bloch wavevector must be specified for PBC geometries");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (M==NULL)
   M=AllocateBEMMatrix();
  else if ( M->NR != TotalBFs || M->NC != TotalBFs )
   { Warn("wrong-size matrix passed to AssembleBEMMatrix; reallocating...");
     M=AllocateBEMMatrix();
   }

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (LDim==0)
   Log("Assembling BEM matrix at Omega=%s",z2s(Omega));
  else if (LDim==1)
   Log("Assembling BEM matrix at {Omega,kx}={%s,%g}",z2s(Omega),kBloch[0]);
  else if (LDim==2)
   Log("Assembling BEM matrix at {Omega,kx,ky}={%s,%g,%g}",z2s(Omega),kBloch[0],kBloch[1]);

  // the overall BEM matrix is symmetric as long as we 
  // don't have a nonzero bloch wavevector.
  bool MatrixIsSymmetric = ( !kBloch || (kBloch[0]==0.0 && kBloch[1]==0.0) );

  /***************************************************************/
  /* loop over all pairs of objects to assemble the diagonal and */
  /* above-diagonal blocks of the matrix                         */
  /***************************************************************/
  int nsm; // 'number of surface mate'
  int nspStart = MatrixIsSymmetric ? 1 : 0;
  for(int ns=0; ns<NumSurfaces; ns++)
   for(int nsp=nspStart*ns; nsp<NumSurfaces; nsp++)
    { 
      // attempt to reuse the diagonal block of an identical previous object
      if (ns==nsp && (nsm=Mate[ns])!=-1)
       { int ThisOffset = BFIndexOffset[ns];
         int MateOffset = BFIndexOffset[nsm];
         int Dim = Surfaces[ns]->NumBFs;
         Log("Block(%i,%i) is identical to block (%i,%i) (reusing)",ns,ns,nsm,nsm);
         M->InsertBlock(M, ThisOffset, ThisOffset, Dim, Dim, MateOffset, MateOffset);
       }
      else
       AssembleBEMMatrixBlock(ns, nsp, Omega, kBloch, M, 0,
                              BFIndexOffset[ns], BFIndexOffset[nsp]);
    }

  /***************************************************************/
  /* if the matrix is symmetric, then the computations above have*/
  /* only filled in its upper triangle, so we need to go back and*/
  /* fill in the lower triangle. (The exception is if the matrix */
  /* is defined to use packed storage, in which case only the    */
  /* upper triangle is needed anyway.)                           */
  /* Note: Technically the lower-triangular parts of the diagonal*/
  /* blocks should already have been filled in, so this code is  */
  /* slightly redundant because it re-fills-in those entries.    */
  /***************************************************************/
  if (MatrixIsSymmetric && M->StorageType==LHM_NORMAL)
   { 
     for(int nr=1; nr<TotalBFs; nr++)
      for(int nc=0; nc<nr; nc++)
       M->SetEntry(nr, nc, M->GetEntry(nc, nr) );
   };

  if (UseHRWGFunctions && NumMMJs>0 )
   ApplyMMJTransformation(M, 0);

  return M;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *RWGGeometry::AssembleBEMMatrix(cdouble Omega, HMatrix *M)
{
  return AssembleBEMMatrix(Omega, 0, M); 
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *RWGGeometry::AllocateBEMMatrix(bool PureImagFreq, bool Packed)
{
  int Storage = Packed ? LHM_SYMMETRIC : LHM_NORMAL;
  int DataType = (!LBasis && PureImagFreq) ? LHM_REAL : LHM_COMPLEX;
  return new HMatrix(TotalBFs, TotalBFs, DataType, Storage);
}

} // namespace scuff
