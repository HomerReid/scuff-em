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
 * FrequencyIntegrand.cc -- evaluate the spectral density of 
 *                       -- power transfer at a single frequency
 *
 * homer reid            -- 2/2012
 *
 */

#include "scuff-heat.h"
#include "libscuffInternals.h"

/***************************************************************/
/***************************************************************/
/***************************************************************/
void CreateFluxPlot(SHData *SHD, cdouble Omega, char *Tag)
{ 
  RWGGeometry *G = SHD->G;

  /***************************************************************/
  /* sanity check ************************************************/
  /***************************************************************/
  int ns;
  for(ns=0; ns<G->NumSurfaces; ns++)
   if ( G->Surfaces[ns]->IsPEC )
    ErrExit("flux plot not available for geometries containing PEC surfaces");

  /***************************************************************/
  /* allocate a vector with enough slots to store one double     */
  /* value per panel in the mesh                                 */
  /***************************************************************/
  double *PFV=(double *)mallocEC( (G->TotalPanels)*sizeof(double) ); // panel flux vector
  memset(PFV, 0, (G->TotalPanels)*sizeof(double));

  /***************************************************************/
  /* fill in the panel flux vector with the flux on each panel   */
  /***************************************************************/
  RWGSurface *S;
  RWGEdge *E;
  int ne, BFIndex, PanelIndex;
  double Value;
  int Offset = (G->NumSurfaces == 1) ? 0 : G->Surfaces[0]->NumBFs;
  for(ns=(G->NumSurfaces==1) ? 0 : 1; ns<G->NumSurfaces; ns++)
   for(S=G->Surfaces[ns], ne=0; ne<S->NumEdges; ne++)
    { 
      E=S->Edges[ne];
      
      BFIndex = G->BFIndexOffset[ns] - Offset + 2*ne;
      Value = 0.5 * ( SHD->DV->GetEntryD(BFIndex + 0) + SHD->DV->GetEntryD(BFIndex + 1) ); 

      PanelIndex = G->PanelIndexOffset[ns] + E->iPPanel;
      PFV[ PanelIndex ] += 0.5*Value / (S->Panels[E->iPPanel]->Area);

      PanelIndex = G->PanelIndexOffset[ns] + E->iMPanel;
      PFV[ PanelIndex ] += 0.5*Value / (S->Panels[E->iMPanel]->Area);

    };
  
  /***************************************************************/
  /* create a GMSH postprocessing file containing a single 'view'*/
  /* that plots the flux on each panel                           */
  /***************************************************************/
  FILE *f=vfopen("%s.%g.%s.flux.pp","w",GetFileBase(G->GeoFileName),real(Omega),Tag);
  fprintf(f,"View \"Flux\"{\n");
  int np;
  RWGPanel *P;
  double *PV[3];
  for(PanelIndex=ns=0; ns<G->NumSurfaces; ns++)
   for(S=G->Surfaces[ns], np=0, P=S->Panels[0]; np<S->NumPanels; np++, PanelIndex++)
    {
      P=S->Panels[np];
      PV[0]=S->Vertices + 3*P->VI[0];
      PV[1]=S->Vertices + 3*P->VI[1];
      PV[2]=S->Vertices + 3*P->VI[2];

      if ( PFV[PanelIndex]!=0.0 )
       fprintf(f,"ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n",
                  PV[0][0], PV[0][1], PV[0][2],
                  PV[1][0], PV[1][1], PV[1][2],
                  PV[2][0], PV[2][1], PV[2][2],
                  PFV[PanelIndex],PFV[PanelIndex],PFV[PanelIndex]);
    };
  fprintf(f,"};\n\n");
  fclose(f);


}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void FlipSignOfMagneticColumns(HMatrix *B)
{ 
  int nr, nc;
  for (nr=0; nr<B->NR; nr++)
   for (nc=1; nc<B->NC; nc+=2)
    B->SetEntry(nr, nc, -1.0*B->GetEntry(nr, nc) );
}

void FlipSignOfMagneticRows(HMatrix *B)
{
  int nr, nc;
  for (nr=1; nr<B->NR; nr+=2)
   for (nc=0; nc<B->NC; nc++)
    B->SetEntry(nr, nc, -1.0*B->GetEntry(nr, nc) );
}

void FillInLowerTriangle(HMatrix *B)
{
  int nr, nc;
  for (nr=1; nr<B->NR; nr++)
   for (nc=0; nc<nr; nc++)
    B->SetEntry(nr, nc, B->GetEntry(nc, nr) );
}

/***************************************************************/
/* this routine mimics the effect of saying                    */
/*  A->InsertBlock(B, RowOffset, ColOffset)                    */
/* with the difference that the symmetrized version of the B   */
/* matrix, i.e. (B+B^\dagger) / 2, is inserted instead.        */
/***************************************************************/
void InsertSymmetrizedBlock(HMatrix *A, HMatrix *B, int RowOffset, int ColOffset)
{ 
  int nr, nc;
  cdouble Z;

  for(nr=0; nr<B->NR; nr++)
   { 
     A->SetEntry(RowOffset+nr, ColOffset+nr, real(B->GetEntry(nr,nr)) );

     for(nc=nr+1; nc<B->NC; nc++)
      { Z = 0.5*( B->GetEntry(nr,nc) + conj(B->GetEntry(nc,nr)) );
        A->SetEntry( RowOffset+nr, ColOffset+nc, Z       );
        A->SetEntry( RowOffset+nc, ColOffset+nr, conj(Z) );
      };
   };

} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetFrequencyIntegrand(SHData *SHD, cdouble Omega, double *FI)
{

  /***************************************************************/
  /* extract fields from SHData structure ************************/
  /***************************************************************/
  RWGGeometry *G     = SHD->G;
  int N1             = SHD->N1;
  HMatrix **TSelf    = SHD->TSelf;
  HMatrix **TMedium  = SHD->TMedium;
  HMatrix **UMedium  = SHD->UMedium;
  HMatrix *SymG1     = SHD->SymG1;
  HMatrix *SymG2     = SHD->SymG2;
  HMatrix *W         = SHD->W;
  HMatrix *W21       = SHD->W21;
  HMatrix *W21SymG1  = SHD->W21SymG1;
  HMatrix *W21DSymG2 = SHD->W21DSymG2;
  HMatrix *Scratch   = SHD->Scratch;
  HVector *DV        = SHD->DV;
  int PlotFlux       = SHD->PlotFlux;

  Log("Computing heat radiation/transfer at omega=%s...",z2s(Omega));

  /***************************************************************/
  /* before entering the loop over transformations, we first     */
  /* assemble the (transformation-independent) T matrix blocks.  */
  /***************************************************************/
  int ns, nsp, nb, nr, NS=G->NumSurfaces;
  for(ns=0; ns<G->NumSurfaces; ns++)
   { 
     Log(" Assembling self contributions to T(%i)...",ns);
     G->RegionMPs[0]->Zero();
     G->AssembleBEMMatrixBlock(ns, ns, Omega, 0, TSelf[ns]);
     FillInLowerTriangle(TSelf[ns]);
     FlipSignOfMagneticColumns(TSelf[ns]);
     G->RegionMPs[0]->UnZero();

     Log(" Assembling medium contributions to T(%i)...",ns);
     G->RegionMPs[ns+1]->Zero();
     G->AssembleBEMMatrixBlock(ns, ns, Omega, 0, TMedium[ns]);
     FillInLowerTriangle(TMedium[ns]);
     FlipSignOfMagneticColumns(TMedium[ns]);
     G->RegionMPs[ns+1]->UnZero();

   };

  /***************************************************************/
  /* pause to dump out the scuff cache to disk. this will be     */
  /* overwritten later if everything goes well, but if something */
  /* goes wrong and we never make it to the later step it's nice */
  /* to dump it out here so that we can at least have a partial  */
  /* cache to accelerate a subsequent calculation involving      */
  /* any of the same surfaces                                    */
  /***************************************************************/
  if ( SHD->WriteCache ) 
   StoreCache( SHD->WriteCache );

  /***************************************************************/
  /* the SymG1 matrix can be formed and stored ahead of time     */
  /* (unlike SymG2)                                              */
  /***************************************************************/
  InsertSymmetrizedBlock(SymG1, TSelf[0], 0, 0 );

  /***************************************************************/
  /* now loop over transformations. ******************************/
  /* note: 'gtc' stands for 'geometrical transformation complex' */
  /***************************************************************/
  int nt;
  char *Tag;
  int RowOffset, ColOffset;
  for(nt=0; nt<SHD->NumTransformations; nt++)
   { 
     /*--------------------------------------------------------------*/
     /*- transform the geometry -------------------------------------*/
     /*--------------------------------------------------------------*/
     Tag=SHD->GTCList[nt]->Tag;
     G->Transform(SHD->GTCList[nt]);
     Log(" Computing quantities at geometrical transform %s",Tag);

     /*--------------------------------------------------------------*/
     /* assemble off-diagonal matrix blocks.                         */
     /* note that not all off-diagonal blocks necessarily need to    */
     /* be recomputed for all transformations; this is what the 'if' */
     /* statement here is checking for.                              */
     /*--------------------------------------------------------------*/
     for(nb=0, ns=0; ns<NS; ns++)
      for(nsp=ns+1; nsp<NS; nsp++, nb++)
       if ( nt==0 || G->SurfaceMoved[ns] || G->SurfaceMoved[nsp] )
        { 
          Log("  Assembling U(%i,%i)...",ns,nsp);
          G->AssembleBEMMatrixBlock(ns, nsp, Omega, 0, UMedium[nb]);
          FlipSignOfMagneticColumns(UMedium[nb]);
        };

     /*--------------------------------------------------------------*/
     /*- put together the full BEM matrix by stamping the T0, TN,    */
     /*- and U blocks in their appropriate places, then LU-factorize */
     /*- and invert it to get the W matrix.                          */
     /*--------------------------------------------------------------*/
     for(nb=0, ns=0; ns<NS; ns++)
      { 
        RowOffset=G->BFIndexOffset[ns];
        W->InsertBlock(TSelf[ns], RowOffset, RowOffset);
        W->AddBlock(TMedium[ns], RowOffset, RowOffset);

        for(nsp=ns+1; nsp<NS; nsp++, nb++)
         { ColOffset=G->BFIndexOffset[nsp];
           W->InsertBlock(UMedium[nb], RowOffset, ColOffset);
           
           FlipSignOfMagneticColumns(UMedium[nb]);
           FlipSignOfMagneticRows(UMedium[nb]);
           W->InsertBlockTranspose(UMedium[nb], ColOffset, RowOffset);
           FlipSignOfMagneticRows(UMedium[nb]);
           FlipSignOfMagneticColumns(UMedium[nb]);
         };
      };
     Log("  LU factorizing M...");
     W->LUFactorize();

     /*--------------------------------------------------------------*/
     /*- invert the W matrix and extract the lower-left subblock W21.*/
     /*--------------------------------------------------------------*/
#if 0 // old (20120306) slower method
     Log("  LU inverting M...");
     W->LUInvert();
     W->ExtractBlock(N1, 0, W21);
#else // new (20120307) hopefully faster method: instead of LUSolving
      // with the full identity matrix to get the full matrix inverse,
      // we LUSolve with just the first N1 columns of the identity matrix
      // since this gives us the only chunk of the inverse that we 
      // need.
      // note: we could achieve a further speedup by truncating 
      // the back-substitution so that we only carry it out far
      // enough to extract the bottommost N2 entries in each
      // row, but this would involve tweaking the lapack routines,
      // so leave it TODO.
     Log("  Partially LU-inverting M...");
     Scratch->Zero();
     for(nr=0; nr<N1; nr++)
      Scratch->SetEntry(nr, nr, 1.0);
     W->LUSolve(Scratch);
     if (NS==1)
      Scratch->ExtractBlock(0, 0, W21);
     else
      Scratch->ExtractBlock(N1, 0, W21);
#endif

     /*--------------------------------------------------------------*/
     /*- fill in the SymG2 matrix. this is just what we did for the  */
     /*- SymG1 matrix above, except that there may be more than one  */
     /*- block to insert, and also if PlotFlux==1 then we don't want */
     /*- to symmetrize.                                              */
     /*--------------------------------------------------------------*/
     SymG2->Zero();
     if (NS==1)
      { if (PlotFlux)
         SymG2->InsertBlock(TMedium[0], 0, 0);
        else 
         InsertSymmetrizedBlock(SymG2, TMedium[0], 0, 0);
      }
     else
      { for(ns=1; ns<NS; ns++)
         { RowOffset=G->BFIndexOffset[ns] - N1;
           if (PlotFlux)
            SymG2->InsertBlock(TSelf[ns], RowOffset, RowOffset );
           else 
            InsertSymmetrizedBlock(SymG2, TSelf[ns], RowOffset, RowOffset );
         };
      };

     /*--------------------------------------------------------------*/
     /*- compute the products W21*sym(G1) and W21^{\dagger} * sym(G2)*/
     /*--------------------------------------------------------------*/
     Log("  Multiplication 1...");
     W21->Multiply(SymG1, W21SymG1);

     Log("  Multiplication 2...");
     W21->Adjoint();
     W21->Multiply(SymG2, W21DSymG2);

     // we have to do this again so that W21 will be the correct
     // size on the next go-round
     W21->Adjoint();

     /*--------------------------------------------------------------*/
     /*- compute the diagonal elements of the matrix--matrix product */ 
     /*- W21*sym(G1)*W21^{\dagger}*sym(G2)                           */
     /*--------------------------------------------------------------*/
     Log("  Multiplication 3...");
     //W21SymG1->Multiply(W21DSymG2, SymG2);
     W21SymG1->GetMatrixProductDiagonal(W21DSymG2, DV);

     /*--------------------------------------------------------------*/
     /*- if we are plotting the flux, then extract the diagonal of   */
     /*- the matrix we just computed and plot it; otherwise, compute */
     /*- the trace and put it in the next slot in the output vector. */
     /*--------------------------------------------------------------*/
     if (PlotFlux)
      { 
        CreateFluxPlot(SHD, Omega, Tag);
      }
     else
      { 
        for(FI[nt]=0.0, nr=0; nr<DV->N; nr++)
         FI[nt] += DV->GetEntryD(nr);
        FI[nt] /= 8.0;

        /***************************************************************/
        /* write the result to the frequency-resolved output file ******/
        /***************************************************************/
        FILE *f=fopen(SHD->ByOmegaFile, "a");
        fprintf(f,"%s %s %e\n",Tag,z2s(Omega),FI[nt]);
        fclose(f);
      };

     /*--------------------------------------------------------------*/
     /* and untransform the geometry                                 */
     /*--------------------------------------------------------------*/
     G->UnTransform();
     Log(" ...done!");

   }; // for (nt=0; nt<SHD->NumTransformations... )

  /*--------------------------------------------------------------*/
  /*- at the end of the first successful frequency calculation,  -*/
  /*- we dump out the cache to disk, and then tell ourselves not -*/
  /*- to dump the cache to disk again (explain me)               -*/
  /*--------------------------------------------------------------*/
  if ( SHD->WriteCache ) 
   { StoreCache( SHD->WriteCache );
     SHD->WriteCache=0;
   };

}
