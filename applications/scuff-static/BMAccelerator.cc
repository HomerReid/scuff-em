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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

#include <scuff-static.h>

using namespace scuff;

/***************************************************************/
/* return true if GT is the identity transformation            */
/***************************************************************/
bool IsIdentity(GTransformation *GT,
                double LengthScale=1.0, double Tol=1.0e-8)
{ 
  if ( fabs(GT->DX[0]) > Tol*LengthScale ) return false;
  if ( fabs(GT->DX[1]) > Tol*LengthScale ) return false;
  if ( fabs(GT->DX[2]) > Tol*LengthScale ) return false;
  if ( fabs(GT->M[0][0] - 1.0) > Tol     ) return false;
  if ( fabs(GT->M[1][1] - 1.0) > Tol     ) return false;
  if ( fabs(GT->M[2][2] - 1.0) > Tol     ) return false;
  return true;
}

// index of off-diagonal matrix block (ns,nsp) where nsp>ns
int OffDiagonalBlockIndex(int NS, int ns, int nsp)
{ return ns*NS - ns*(ns+1)/2 + nsp - ns - 1; }

/***************************************************************/
/* Detect equivalent surface-surface pairs to allow reuse of   */
/*  off-diagonal BEM matrix blocks.                            */
/*                                                             */
/* (SAlpha, SBeta) is equivalent to (SA, SB) if the following  */
/*  three conditions are all satified.                         */
/*                                                             */
/*  1. SAlpha and SA are both the result of geometrical        */
/*     transformations applied to a common surface S1:         */
/*       SAlpha = TAlpha(S1)                                   */
/*       SA     = TA    (S1)                                   */
/*    where TAlpha, TA are geometrical transforms (which may   */
/*    be the identity transform)                               */
/*                                                             */
/*  2. SBeta and SB are both the result of geometrical         */
/*     transformations applied to a common surface S2:         */
/*       SBeta  = TBeta(S2)                                    */
/*       SB     = TB    S2)                                    */
/*                                                             */
/*  3. We have TAlpha^{-1} * TBeta = TA^{-1} TB                */
/*      or                                                     */
/*             TAlpha^{-1} * TBeta * TB^{-1} * TA = identity   */
/*                                                             */
/* If (SAlpha, SBeta) <==> (SA, SB) then the returned matrix   */
/*  M satifies                                                 */
/*                                                             */
/*  M[SAlpha, SBeta] = index of off-diagonal block (SA,SB)     */
/*                                                             */
/* (with NS==number of surfaces)                               */
/***************************************************************/
HMatrix *GetEquivalentPairMatrix(RWGGeometry *G, HMatrix *M=0)
{
  int NS    = G->NumSurfaces;
  int *Mate = G->Mate;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (M && ((M->NR!=NS) || (M->NC!=NS)) )
   { delete M;
     M=0;
   };
  if (!M)
   M=new HMatrix(NS, NS);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int nsAlpha=1; nsAlpha<NS; nsAlpha++)
   for(int nsBeta=nsAlpha+1; nsBeta<NS; nsBeta++)
    { 
      int nsAlphaMate = Mate[nsAlpha];
      int nsBetaMate  = Mate[nsBeta];
      if (nsAlphaMate==-1 && nsBetaMate==-1)
       continue;

      if (nsAlphaMate==-1) nsAlphaMate = nsAlpha;
      if (nsBetaMate==-1)  nsBetaMate  = nsBeta;

      for(int nsA=0; nsA<=nsAlpha; nsA++)
       for(int nsB=1; nsB<=nsBeta; nsB++)
        {
          if (nsA==nsAlpha && nsB==nsBeta) continue;

          int nsAMate = (Mate[nsA]==-1) ? nsA: Mate[nsA];
          int nsBMate = (Mate[nsB]==-1) ? nsB: Mate[nsB];

          if (nsAlphaMate!=nsAMate || nsBetaMate!=nsBMate)
           continue;

          RWGSurface *SAlpha = G->Surfaces[nsAlpha];
          RWGSurface *SBeta  = G->Surfaces[nsBeta];
          RWGSurface *SA     = G->Surfaces[nsA];
          RWGSurface *SB     = G->Surfaces[nsB];
          
          GTransformation T; // identity transformation
    
          if (SAlpha->GT)   T = T - *(SAlpha->GT);
          if (SAlpha->OTGT) T = T - *(SAlpha->OTGT);
    
          if (SBeta->OTGT)  T = T + *(SBeta->OTGT);
          if (SBeta->GT)    T = T + *(SBeta->GT);
    
          if (SB->GT)       T = T - *(SB->GT);
          if (SB->OTGT)     T = T - *(SB->OTGT);
    
          if (SA->OTGT)     T = T + *(SA->OTGT);
          if (SA->GT)       T = T + *(SA->GT);
    
          double Lengthscale = fmin( VecDistance(SA->RMax, SA->RMin),
                                     VecDistance(SB->RMax, SB->RMin)
                                   );
          if ( IsIdentity(&T, Lengthscale) )
           { double Sign=1.0;
             if (nsB<nsA)
              { int temp = nsA; nsA=nsB; nsB=temp;
                SA=G->Surfaces[nsA];
                SB=G->Surfaces[nsB];
                Sign=-1.0;
              };
             M->SetEntry(nsAlpha, nsBeta, Sign * (double)OffDiagonalBlockIndex(NS, nsA, nsB));
             Log(" %10s<-->%-10s === %10s<-->%-10s",SAlpha->Label,SBeta->Label,SA->Label,SB->Label);
             nsA=nsAlpha; nsB=nsBeta; // to break out of inner loop
           };

        };
    }; 
  return M;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
BMAccelerator *CreateBMAccelerator(StaticSolver *SS, GTCList GTCs,
                                   bool UsePairEquivalence)
{
  RWGGeometry *G     = SS->G;
  int NS             = G->NumSurfaces;
  int NADB           = NS*(NS-1)/2; // number of above-diagonal blocks
  int NT             = GTCs.size();

  BMAccelerator *BMA = (BMAccelerator *)mallocEC(sizeof *BMA);
  BMA->TBlocks       = (HMatrix **)mallocEC(NS*sizeof(HMatrix *));
  BMA->UBlocks       = (HMatrix **)mallocEC(NADB*sizeof(HMatrix *));
  
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int ns=0, nb=0; ns<NS; ns++)
   { 
     int NBF=G->Surfaces[ns]->NumPanels;

     // allocate diagonal block for surface #ns, unless
     // it has a mate in which case we reuse the mate's block
     int nsMate = G->Mate[ns];
     if ( nsMate!=-1 )
      BMA->TBlocks[ns] = BMA->TBlocks[nsMate];
     else
      BMA->TBlocks[ns] = new HMatrix(NBF, NBF, LHM_REAL);

     // allocate off-diagonal blocks for surfaces ns,nsp>ns
     // TODO exploit pair equivalence to reduce number of 
     // blocks allocated
     for(int nsp=ns+1; nsp<NS; nsp++, nb++)
      { int NBFp = G->Surfaces[nsp]->NumPanels;
        BMA->UBlocks[nb] = new HMatrix(NBF, NBFp, LHM_REAL);
      };
   };
  
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if( !UsePairEquivalence )
   { 
     BMA->PairMatrices=0;
   }
  else
   { 
     BMA->PairMatrices = (HMatrix **)mallocEC(NT*sizeof(HMatrix *));
     for(int nt=0; nt<NT; nt++)
      { Log("Looking for equivalent surface pairs (transform %s)...",GTCs[nt]->Tag);
        G->Transform(GTCs[nt]);
        BMA->PairMatrices[nt]=GetEquivalentPairMatrix(G);
        G->UnTransform();
      }
   }

  return BMA;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ReassembleBEMMatrix(StaticSolver *SS, HMatrix **pM,
                         BMAccelerator *BMA, int nt)
{
  HMatrix *M = *pM;

  if (M==0)
   M=*pM=SS->AllocateBEMMatrix();

  /*******************************************************************/
  /* assemble BEM matrix, either (a) all at once if we have only one */
  /* geometric transformation, or (b) with the diagonal and off-     */
  /* diagonal blocks computed separately so that the former can be   */
  /* reused for multiple geometric transformations                   */
  /*******************************************************************/
  if (BMA==0)
   { 
     SS->AssembleBEMMatrix(M);
   }
  else
   { 
     Log("Assembling BEM matrix blocks...");
     RWGGeometry *G      = SS->G;
     int NS              = G->NumSurfaces;
     HMatrix **TBlocks   = BMA->TBlocks;
     HMatrix **UBlocks   = BMA->UBlocks;
     HMatrix *PairMatrix = (BMA->PairMatrices ? BMA->PairMatrices[nt] : 0);
 
     /*******************************************************************/
     /* recompute T blocks only on first transformation *****************/
     /*******************************************************************/
     if (nt==0)
      for(int ns=0; ns<NS; ns++)
       if (G->Mate[ns]!=-1)
        Log("reusing diagonal block %i == %i",ns,G->Mate[ns]);
       else
        SS->AssembleBEMMatrixBlock(ns, ns, TBlocks[ns]);

     /*******************************************************************/
     /* recompute U blocks only for surfaces that have moved and for    */
     /* which there is no equivalent pair                               */
     /*******************************************************************/
     for(int ns=0, nb=0; ns<NS; ns++)
      for(int nsp=ns+1; nsp<NS; nsp++, nb++)
       { if (PairMatrix && (PairMatrix->GetEntryD(ns, nsp)!=0.0) )
          { Log("reusing off-diagonal block %i == %i ",nb,(int)PairMatrix->GetEntryD(ns,nsp));
            continue;
          };
         if (nt==0 || G->SurfaceMoved[ns] || G->SurfaceMoved[nsp] )
          SS->AssembleBEMMatrixBlock(ns, nsp, UBlocks[nb]);
       };

     /*******************************************************************/
     /* stamp blocks into BEM matrix ************************************/
     /*******************************************************************/
     for(int ns=0, nb=0; ns<NS; ns++)
      { 
        int RowOffset=G->PanelIndexOffset[ns];
        M->InsertBlock(TBlocks[ns], RowOffset, RowOffset);
        for(int nsp=ns+1; nsp<NS; nsp++, nb++)
         { int ColOffset=G->PanelIndexOffset[nsp];

           int nbEffective=nb;
           bool Flip=false;
           if (PairMatrix && (PairMatrix->GetEntryD(ns,nsp)!=0.0))
            { nbEffective=(int)(fabs(PairMatrix->GetEntryD(ns,nsp)));
              Flip = (PairMatrix->GetEntryD(ns,nsp) < 0.0 );
            };

           if (Flip)
            { M->InsertBlockTranspose(UBlocks[nbEffective], RowOffset, ColOffset);
              M->InsertBlock(UBlocks[nbEffective], ColOffset, RowOffset);
            }
           else
            { M->InsertBlock(UBlocks[nbEffective], RowOffset, ColOffset);
              M->InsertBlockTranspose(UBlocks[nbEffective], ColOffset, RowOffset);
            };
         };
      };
   };

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  M->LUFactorize();

}
