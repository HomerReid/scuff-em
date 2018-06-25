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
 * AssembleBEMMatrix2018.cc -- the `modern' approach to assembling the system matrix
 *
 * homer reid
 */
#include "libscuff.h"
#include "libscuffInternals.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#ifdef USE_OPENMP
#  include <omp.h>
#endif

#define II cdouble(0,1)

namespace scuff {

/***************************************************************/
/***************************************************************/
/***************************************************************/
char *ToStr(cdouble Omega, double *kBloch, int LDim, char *s=0)
{ static char sBuffer[50];
  if (s==0) s=sBuffer;
  if (LDim==0) snprintf(s,50,"Omega=%s",z2s(Omega));
  if (LDim==1) snprintf(s,50,"Omega,kx={%s,%g}",z2s(Omega),kBloch[0]);
  if (LDim==2) snprintf(s,50,"Omega,kx,ky={%s,%g,%g}",z2s(Omega),kBloch[0],kBloch[1]);
  return s;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
#if 0
void AddSubstrateContributionToSIEMatrixElements(RWGGeometry *G, GetGCMEArgStruct *Args,
                                                 int nea, int neb,
                                                 cdouble Gab, cdouble ikCab, cdouble GPhiTerm,
                                                 cdouble iwEps, cdouble iwMu, HMatrix *MEs)
{
}
#endif

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetSIEMatrixElements(RWGGeometry *G, GetGCMEArgStruct *Args, int nea, int neb,
                          cdouble iwEps[2], cdouble iwMu[2], dVec Signs, HMatrix *MEs)
{
  cdouble Gab[2][NUMGCMES], ikCab[2][NUMGCMES];
  GetGCMatrixElements(G, Args, nea, neb, Gab, ikCab);
  Gab[0][0]*=Signs[0];   ikCab[0][0]*=Signs[0];
  Gab[1][0]*=Signs[1];   ikCab[1][0]*=Signs[1];
  MEs->SetEntry(0, 0, iwMu[0]*Gab[0][0] + iwMu[1]*Gab[1][0]);
  if (MEs->NC>1)
   MEs->SetEntry(0, 1, -1.0*(ikCab[0][0] + ikCab[1][0]));
  if (MEs->NR>1)
   MEs->SetEntry(1, 0, -1.0*(ikCab[0][0] + ikCab[1][0]));
  if (MEs->NR>1 && MEs->NC>1)
   MEs->SetEntry(1, 1, -1.0*(iwEps[0]*Gab[0][0] + iwEps[1]*Gab[1][0]));
}

/***************************************************************/
/* Given two surfaces, identify whether they bound zero, one,  */
/* or two common regions. If there are any common regions,     */
/* identify their indices and the relative sign between the    */
/* contributions of surface currents on the two surfaces to    */
/* fields in those regions.                                    */
/***************************************************************/
iVec GetCommonRegions(RWGSurface *Sa, RWGSurface *Sb, dVec &Signs)
{
  iVec CommonRegions;
  int NRA = (Sa->RegionIndices[1]==-1 ? 2 : 1), NRB = (Sb->RegionIndices[1]==-1 ? 2 : 1);
  for(int nra=0; nra<NRA; nra++)
   for(int nrb=0; nrb<NRB; nrb++)
    if ( Sa->RegionIndices[nra] == Sb->RegionIndices[nrb] )
     { CommonRegions.push_back(Sa->RegionIndices[nra]);
       Signs.push_back( (nra==0)==(nrb==0) ? +1.0 : -1.0 );
     }
  return CommonRegions;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AssembleBEMMatrixBlock2018(RWGGeometry *G, int nsa, int nsb,
                                cdouble Omega, HMatrix *Block, 
                                int OffsetA, int OffsetB)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  EquivalentEdgePairTable *EEPTable=0;
  if ( G->EEPTables.size()>0 && G->EEPTables[nsa][nsb]==0 )
   EEPTable = G->EEPTables[nsa][nsb] = new EquivalentEdgePairTable(G,nsa,nsb);
  
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  RWGSurface *Sa = G->Surfaces[nsa], *Sb=G->Surfaces[nsb];
  int NEA=Sa->NumEdges, NBFPEA = Sa->IsPEC ? 1 : 2;
  int NEB=Sb->NumEdges, NBFPEB = Sb->IsPEC ? 1 : 2;
  int NumPairs = (Sa==Sb ? NEA*(NEA+1)/2 : NEA*NEB);

  cdouble iwEps[2]={0.0, 0.0}, iwMu[2]={0.0, 0.0};
  dVec Signs;
  iVec CommonRegions   = GetCommonRegions(Sa, Sb, Signs);
  int NumCommonRegions = CommonRegions.size();
  for(int nr=0; nr<NumCommonRegions; nr++)
   G->RegionMPs[CommonRegions[nr]]->GetEpsMu(Omega, iwEps+nr, iwMu+nr);

  GetGCMEArgStruct Args;
  InitGetGCMEArgs(&Args);
  Args.nsa        = nsa;
  Args.nsb        = nsb;
  Args.NumRegions = CommonRegions.size();
  Args.k[0]       = sqrt(iwEps[0]*iwMu[0]*Omega);
  if (NumCommonRegions==2)
   Args.k[1]       = sqrt(iwEps[1]*iwMu[1]*Omega);
  Args.NeedGC = true;
  Args.FIBBICache = (nsa==nsb) ? G->FIBBICaches[nsa] : 0;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NumThreads=GetNumThreads();
  Log("Computing SIE matrix block (%i,%i) (%i pairs) (%i threads)",nsa,nsb,NumPairs,NumThreads);
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int nPair=0; nPair<NEA*NEB; nPair++)
   { 
     if (G->LogLevel>=SCUFF_VERBOSE2) LogPercent(nPair, NEA*NEB);

     int nea = nPair/NEB;
     int neb = nPair%NEB;
     if (nsa==nsb && neb<nea) continue;

     if (EEPTable && EEPTable->HasParent(nea,neb)) continue;

     cdouble MEBuffer[4];
     HMatrix MEs(NBFPEA, NBFPEB, MEBuffer);
     GetSIEMatrixElements(G, &Args, nea, neb, iwEps, iwMu, Signs, &MEs);
     Block->InsertBlock(&MEs, OffsetA + NBFPEA*nea, OffsetB + NBFPEB*neb);
   }

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (EEPTable)
   { if (G->LogLevel>=SCUFF_VERBOSE2) Log("Filling in matrix entries for child edge pairs ...");
     ParentPairList Parents = EEPTable->GetParents();
     for(size_t np=0; np<Parents.size(); np++)
      { cdouble ParentMEBuffer[4]={0.0, 0.0, 0.0, 0.0};
        HMatrix ParentMEs(NBFPEA, NBFPEB, ParentMEBuffer);
        Block->ExtractBlock(OffsetA + NBFPEA*Parents[np].nea, OffsetB + NBFPEB*Parents[np].neb, &ParentMEs);
        ChildPairList Children = EEPTable->GetChildren(Parents[np].nea, Parents[np].neb);
        for(size_t nc=0; nc<Children.size(); nc++)
         { int neaChild=Children[nc].nea, nebChild=Children[nc].neb;
           double GFlipped = Children[nc].Signs.Flipped[0], CFlipped = Children[nc].Signs.Flipped[1];
           cdouble ChildMEBuffer[4];
           HMatrix ChildMEs(NBFPEA, NBFPEB, ChildMEBuffer);
           memcpy(ChildMEBuffer, ParentMEBuffer, 4*sizeof(cdouble));
           if (GFlipped) { ChildMEBuffer[0]*=-1.0; ChildMEBuffer[3]*=-1.0; }
           if (CFlipped) { ChildMEBuffer[1]*=-1.0; ChildMEBuffer[2]*=-1.0; }
           Block->InsertBlock(&ChildMEs, OffsetA+NBFPEA*neaChild, OffsetB+NBFPEB*nebChild);
        }
      }
   }

  /***************************************************************/
  /* if this is a diagonal block (nsa==nsb) we only filled in the*/
  /* upper triangle, so go back and fill in the lower part now   */
  /***************************************************************/
  if (nsa==nsb)
   for(int nr=1, Offset=OffsetA; nr<NEA; nr++)
    for(int nc=0; nc<nr; nc++)
     Block->SetEntry(Offset+nr, Offset+nc, Block->GetEntry(Offset+nc, Offset+nr) );
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *AssembleBEMMatrix2018(RWGGeometry *G, cdouble Omega, double *kBloch, HMatrix *M)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if ( G->LBasis==0 && kBloch!=0 && (kBloch[0]!=0.0 || kBloch[1]!=0.0) )
   ErrExit("%s:%i: Bloch wavevector is undefined for compact geometries");
  if ( G->LBasis!=0 && kBloch==0 )
   ErrExit("%s:%i: Bloch wavevector must be specified for PBC geometries");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (M==NULL)
   M=G->AllocateBEMMatrix();
  else if ( M->NR != G->TotalBFs || M->NC != G->TotalBFs )
   { Warn("wrong-size matrix passed to AssembleBEMMatrix; reallocating...");
     M=G->AllocateBEMMatrix();
   }

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  Log("Assembling BEM matrix at %s",ToStr(Omega, kBloch, G->LDim));

  // the full system matrix is symmetric as long as we don't have a nonzero bloch wavevector.
  bool MatrixSymmetric = ( !kBloch || (kBloch[0]==0.0 && kBloch[1]==0.0) );

  /***************************************************************/
  /* loop over all pairs of objects to assemble the diagonal and */
  /* above-diagonal blocks of the matrix                         */
  /***************************************************************/
  for(int nsa=0; nsa<G->NumSurfaces; nsa++)
   for(int nsb = (MatrixSymmetric ? nsa : 0); nsb<G->NumSurfaces; nsb++)
    { 
      // attempt to reuse the diagonal block of an identical previous object
      if (nsa==nsb && (G->Mate[nsa])!=-1)
       { int nsMate = G->Mate[nsa];
         int ThisOffset = G->BFIndexOffset[nsa];
         int MateOffset = G->BFIndexOffset[nsMate];
         int Dim = G->Surfaces[nsa]->NumBFs;
         Log("Block(%i,%i) is identical to block (%i,%i) (reusing)",nsa,nsa,nsMate,nsMate);
         M->InsertBlock(M, ThisOffset, ThisOffset, Dim, Dim, MateOffset, MateOffset);
       }
      else
       AssembleBEMMatrixBlock2018(G, nsa, nsb, Omega, M, G->BFIndexOffset[nsa], G->BFIndexOffset[nsb]);
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
  if (MatrixSymmetric && M->StorageType==LHM_NORMAL)
   for(int nr=1; nr<G->TotalBFs; nr++)
    for(int nc=0; nc<nr; nc++)
    M->SetEntry(nr, nc, M->GetEntry(nc, nr) );

  return M;
}

} // namespace scuff
