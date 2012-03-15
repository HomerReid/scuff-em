/*
 * FrequencyIntegrand.cc -- evaluate the spectral density of 
 *                       -- power transfer at a single frequency
 *
 * homer reid            -- 2/2012
 *
 */

#include "scuff-cas3D.h"
#include "libscuffInternals.h"

/***************************************************************/
/***************************************************************/
/***************************************************************/
void CreateFluxPlot(SC3Data *SHD, cdouble Omega, char *Tag)
{ 
  RWGGeometry *G = SC3D->G;

  /***************************************************************/
  /* sanity check ************************************************/
  /***************************************************************/
  int no;
  for(no=0; no<G->NumObjects; no++)
   if ( G->Objects[no]->MP->IsPEC() )
    ErrExit("flux plot not available for geometries containing PEC objects");

  /***************************************************************/
  /* allocate a vector with enough slots to store one double     */
  /* value per panel in the mesh                                 */
  /***************************************************************/
  double *PFV=(double *)mallocEC( (G->TotalPanels)*sizeof(double) ); // panel flux vector
  memset(PFV, 0, (G->TotalPanels)*sizeof(double));

  /***************************************************************/
  /* fill in the panel flux vector with the flux on each panel   */
  /***************************************************************/
  RWGObject *O;
  RWGEdge *E;
  int ne, BFIndex, PanelIndex;
  double Value;
  for(no=0; no<G->NumObjects; no++)
   for(O=G->Objects[no], ne=0; ne<O->NumEdges; ne++)
    { 
      E=O->Edges[ne];
      
      BFIndex = G->BFIndexOffset[no] + 2*ne;
      Value = 0.5 * ( SC3D->DV->GetEntryD(BFIndex + 0) + SHD->DV->GetEntryD(BFIndex + 1) ); 

      PanelIndex = G->PanelIndexOffset[no] + E->iPPanel;
      PFV[ PanelIndex ] += 0.5*Value / (O->Panels[E->iPPanel]->Area);

      PanelIndex = G->PanelIndexOffset[no] + E->iMPanel;
      PFV[ PanelIndex ] += 0.5*Value / (O->Panels[E->iMPanel]->Area);

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
  for(PanelIndex=no=0; no<G->NumObjects; no++)
   for(O=G->Objects[no], np=0, P=O->Panels[0]; np<O->NumPanels; np++, PanelIndex++)
   {
      P=O->Panels[np];
      PV[0]=O->Vertices + 3*P->VI[0];
      PV[1]=O->Vertices + 3*P->VI[1];
      PV[2]=O->Vertices + 3*P->VI[2];

      fprintf(f,"ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n",
                 PV[0][0], PV[0][1], PV[0][2],
                 PV[1][0], PV[1][1], PV[1][2],
                 PV[2][0], PV[2][1], PV[2][2],
                 PFV[PanelIndex],PFV[PanelIndex],PFV[PanelIndex]);
    };
  fprintf(f,"};\n\n");
  fclose(f);


}

#if 0
/***************************************************************/
/***************************************************************/
/***************************************************************/

/***************************************************************/
/* how this works: consider a 2x2 block decomposition of the   */
/* matrix, where                                               */
/*                                                             */
/*  block 1 = basis functions on object 1                      */
/*  block 2 = basis functions on all other objects, 2..N       */
/*                                                             */
/* let W be the inverse of the full BEM matrix and let its     */
/* block structure be W = (A B ; C D).                         */
/*                                                             */
/* the matrices G1 and G2 have the block structure             */
/*                                                             */
/*  G1 = [X1 0; 0 0]                                           */
/*  G2 = [ 0 0; 0 X2]                                          */
/*                                                             */
/* (here G1, G2 include the effect of the symmetrization       */
/* operation, G1=sym(G1), G2=sym(G2))                          */
/*                                                             */
/* using the block decompositions, the trace in question       */
/* becomes:                                                    */
/*                                                             */
/*  tr ( W * G1 * W^T * G2 )                                   */
/*   = tr ( C * X1 * B^T * X2 )                                */
/*                                                             */
/*                                                             */
/***************************************************************/

  /***************************************************************/
  /* assemble the three separate blocks that contribute to the   */
  /* BEM matrix:                                                 */
  /*  0: contributions of environment only                       */
  /*  1: contributions of object 1 only                          */
  /*  2: contributions of objects 2...N only                     */
  /***************************************************************/

  /*--------------------------------------------------------------*/
  /*- 0: contributions of environment only -----------------------*/
  /*--------------------------------------------------------------*/
  int no;
  for(no=0; no<G->NumObjects; no++)
   G->Objects[no]->MP->Zero();

  Log("Assembling M0 matrix...");
  G->AssembleBEMMatrix(Omega, nThread, M0);

  // multiply by the 'S' matrix 
  int nr, nc;
  for(nr=1; nr<M0->NR; nr+=2)
   for(nc=0; nc<M0->NC; nc++)
    M0->SetEntry(nr, nc, -1.0*M0->GetEntry(nr,nc));

  /*--------------------------------------------------------------*/
  /*- 1: contributions of object 1 only    -----------------------*/
  /*--------------------------------------------------------------*/
  G->ExteriorMP->Zero();
  G->Objects[0]->MP->UnZero();

  Log("Assembling M1 matrix...");
  G->AssembleBEMMatrix(Omega, nThread, M1);

  for(nr=1; nr<M1->NR; nr+=2)
   for(nc=0; nc<M1->NC; nc++)
    M1->SetEntry(nr, nc, -1.0*M1->GetEntry(nr,nc));

  /*--------------------------------------------------------------*/
  /*- 2: contributions of objects 2--N only.                      */
  /*-                                                             */
  /*-  NOTE: if we only have a single object then we set M2 = M0; */
  /*-        the calculation them amounts to computing the heat   */
  /*-        transfer the single body to the environment.         */
  /*--------------------------------------------------------------*/
  if (G->NumObjects==1)
   M2->Copy(M0);
  else // (G->NumObjects>1)
   { 
     G->Objects[0]->MP->Zero();
     for(no=1; no<G->NumObjects; no++)
      G->Objects[no]->MP->UnZero();

     Log("Assembling M2 matrix...");
     G->AssembleBEMMatrix(Omega, nThread, M2);

     for(nr=1; nr<M2->NR; nr+=2)
      for(nc=0; nc<M2->NC; nc++)
       M2->SetEntry(nr, nc, -1.0*M2->GetEntry(nr,nc));
   };
  
  // undo the zeroing out of the environment and object 1
  G->ExteriorMP->UnZero();
  G->Objects[0]->MP->UnZero();

  /***************************************************************/
  /* assemble and LU-factorize the full BEM matrix (for now we   */
  /* do this in-place using the M0 matrix)                       */
  /***************************************************************/
  if (G->NumObjects==1)
   { for(nr=0; nr<M0->NR; nr++)
      for(nc=0; nc<M0->NC; nc++)
       M0->AddEntry(nr, nc, M1->GetEntry(nr,nc) );
   }
  else
   { for(nr=0; nr<M0->NR; nr++)
      for(nc=0; nc<M0->NC; nc++)
       M0->AddEntry(nr, nc, M1->GetEntry(nr,nc) + M2->GetEntry(nr,nc) );
   };

  Log("LU-factorizing MTotal...");
  M0->LUFactorize();

  /***************************************************************/
  /* set M1 = sym(M1) and M2 = sym(M2)                           */
  /*  (note that the loop runs over only the upper triangle of   */
  /*  the matrices)                                              */
  /***************************************************************/
  cdouble Sym, SymT;
  for(nr=0; nr<M1->NR; nr++)
   for(nc=nr; nc<M1->NC; nc++)
    { 
      Sym = 0.5*(M1->GetEntry(nr, nc) + conj(M1->GetEntry(nc, nr)));
      M1->SetEntry(nr, nc, Sym );
      if(nc>nr) M1->SetEntry(nc, nr, conj(Sym) );
    }; 

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (SC3D->PlotFlux)
   {
     for(nr=0; nr<M2->NR; nr++)
      for(nc=nr; nc<M2->NC; nc++)
       { Sym  = conj(M2->GetEntry(nr, nc));
         SymT = conj(M2->GetEntry(nc, nr));
         M2->SetEntry(nc, nr, Sym );
         if (nc>nr) M2->SetEntry(nr, nc, SymT );
       };

     Log("LU-solving M1...");
     M0->LUSolve(M1,'N');
     Log("LU-solving M2...");
     M0->LUSolve(M2,'C');
     Log("Multipliying...");
     M1->Multiply(M2, M0);
     for(nr=0; nr<M1->NR; nr++)
      SC3D->DV->SetEntry(nr, M1->GetEntry(nr,nr));
     Log("Plotting flux vector...");
     PlotFlux(SC3D, Omega);

     *FI=0.0;
   }
  else
   {
     for(nr=0; nr<M2->NR; nr++)
      for(nc=nr; nc<M2->NC; nc++)
       { Sym = 0.5*(M2->GetEntry(nr, nc) + conj(M2->GetEntry(nc, nr)));
         M2->SetEntry(nr, nc, Sym );
         if(nc>nr) M2->SetEntry(nc, nr, conj(Sym) );
       };

     /***************************************************************/
     /* set M1 <= M^{-1'} * M1 **************************************/
     /* set M2 <= M^{-1}  * M2 **************************************/
     /***************************************************************/
     Log("LU-solving M1...");
     M0->LUSolve(M1,'C');
     Log("LU-solving M2...");
     M0->LUSolve(M2,'N');
   
     /***************************************************************/
     /* set M0 = M1*M2                                              */
     /***************************************************************/
     Log("Multiplying M1*M2...");
     M1->Multiply(M2, M0);
   
     /***************************************************************/
     /* the value of the frequency integrand is now the trace of M0 */
     /***************************************************************/
     *FI = real( M0->GetTrace() ) / 8.0 ;
     Log("...done!");

     /***************************************************************/
     /* write the result to the frequency-resolved output file ******/
     /***************************************************************/
     FILE *f=fopen(SC3D->ByOmegaFile, "a");
     fprintf(f,"%s %e\n",z2s(Omega),*FI);
     fclose(f);
   };
   
}
#endif

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
void GetFrequencyIntegrand(SC3Data *SHD, cdouble Omega, double *FI)
{
  /***************************************************************/
  /* extract fields from SC3Data structure ************************/
  /***************************************************************/
  RWGGeometry *G     = SC3D->G;
  int N1             = SC3D->N1;
  int N2             = SC3D->N2;
  HMatrix **TSelf    = SC3D->TSelf;
  HMatrix **TMedium  = SC3D->TMedium;
  HMatrix **UMedium  = SC3D->UMedium;
  HMatrix *SymG1     = SC3D->SymG1;
  HMatrix *SymG2     = SC3D->SymG2;
  HMatrix *W         = SC3D->W;
  HMatrix *W21       = SC3D->W21;
  HMatrix *W21SymG1  = SC3D->W21SymG1;
  HMatrix *W21DSymG2 = SC3D->W21DSymG2;
  HMatrix *Scratch   = SC3D->Scratch;
  HVector *DV        = SC3D->DV;
  int PlotFlux       = SC3D->PlotFlux;

  /***************************************************************/
  /* preinitialize an argument structure for the BEM matrix      */
  /* block assembly routine                                      */
  /***************************************************************/
  ABMBArgStruct MyABMBArgStruct, *Args=&MyABMBArgStruct;
  InitABMBArgs(Args);
  Args->G         = SC3D->G;
  Args->Omega     = Omega;
  Args->nThread   = SC3D->nThread;

  Log("Computing heat radiation/transfer at omega=%s...",z2s(Omega));

  /***************************************************************/
  /* before entering the loop over transformations, we first     */
  /* assemble the (transformation-independent) T matrix blocks.  */
  /***************************************************************/
  Args->Symmetric=1;
  int no, nop, nb, nr, NO=G->NumObjects;
  for(no=0; no<G->NumObjects; no++)
   { 
     Args->Oa = Args->Ob = G->Objects[no];

     Log(" Assembling self contributions to T(%i)...",no);
     G->ExteriorMP->Zero();
     Args->B = TSelf[no];
     AssembleBEMMatrixBlock(Args);
     G->ExteriorMP->UnZero();

     Log(" Assembling medium contributions to T(%i)...",no);
     G->Objects[no]->MP->Zero();
     Args->B = TMedium[no];
     AssembleBEMMatrixBlock(Args);
     G->Objects[no]->MP->UnZero();

   };

  /***************************************************************/
  /* pause to dump out the scuff cache to disk. this will be     */
  /* overwritten later if everything goes well, but if something */
  /* goes wrong and we never make it to the later step it's nice */
  /* to dump it out here so that we can at least have a partial  */
  /* cache to accelerate a subsequent calculation involving      */
  /* any of the same objects                                     */
  /***************************************************************/
  if ( SC3D->WriteCache ) 
   StoreGlobalFIPPICache( SC3D->WriteCache );

  /***************************************************************/
  /* the SymG1 matrix can be formed and stored ahead of time     */
  /* (unlike SymG2)                                              */
  /***************************************************************/
  InsertSymmetrizedBlock(SymG1, TSelf[0], 0, 0 );
  FlipSignOfMagneticColumns(SymG1);

  /***************************************************************/
  /* now loop over transformations. ******************************/
  /* note: 'gtc' stands for 'geometrical transformation complex' */
  /***************************************************************/
  int nt;
  char *Tag;
  int RowOffset, ColOffset;
  for(nt=0; nt<SC3D->NumTransformations; nt++)
   { 
     /*--------------------------------------------------------------*/
     /*- transform the geometry -------------------------------------*/
     /*--------------------------------------------------------------*/
     Tag=SC3D->GTCList[nt]->Tag;
     G->Transform(SC3D->GTCList[nt]);
     Log(" Computing quantities at geometrical transform %s",Tag);

     /*--------------------------------------------------------------*/
     /* assemble off-diagonal matrix blocks.                         */
     /* note that not all off-diagonal blocks necessarily need to    */
     /* be recomputed for all transformations; this is what the 'if' */
     /* statement here is checking for.                              */
     /*--------------------------------------------------------------*/
     Args->Symmetric=0;
     for(nb=0, no=0; no<NO; no++)
      for(nop=no+1; nop<NO; nop++, nb++)
       if ( nt==0 || G->ObjectMoved[no] || G->ObjectMoved[nop] )
        { 
          Log("  Assembling U(%i,%i)...",no,nop);
          Args->Oa = G->Objects[no];
          Args->Ob = G->Objects[nop];
          Args->B  = UMedium[nb];
          AssembleBEMMatrixBlock(Args);
        };

     /*--------------------------------------------------------------*/
     /*- put together the full BEM matrix by stamping the T0, TN,    */
     /*- and U blocks in their appropriate places, then LU-factorize */
     /*- and invert it to get the W matrix.                          */
     /*--------------------------------------------------------------*/
     for(nb=0, no=0; no<NO; no++)
      { 
        RowOffset=G->BFIndexOffset[no];
        W->InsertBlock(TSelf[no], RowOffset, RowOffset);
        W->AddBlock(TMedium[no], RowOffset, RowOffset);

        for(nop=no+1; nop<NO; nop++, nb++)
         { ColOffset=G->BFIndexOffset[nop];
           W->InsertBlock(UMedium[nb], RowOffset, ColOffset);
           W->InsertBlockTranspose(UMedium[nb], ColOffset, RowOffset);
         };
      };
     FlipSignOfMagneticColumns(W);
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
     Scratch->ExtractBlock(N1, 0, W21);
#endif

     /*--------------------------------------------------------------*/
     /*- fill in the SymG2 matrix. this is just what we did for the  */
     /*- SymG1 matrix above, except that there may be more than one  */
     /*- block to insert, and also if PlotFlux==1 then we don't want */
     /*- to symmetrize.                                              */
     /*--------------------------------------------------------------*/
     SymG2->Zero();
     for(no=1; no<NO; no++)
      { RowOffset=G->BFIndexOffset[no] - N1;
        if (PlotFlux)
         SymG2->InsertBlock(TSelf[no], RowOffset, RowOffset );
        else 
         InsertSymmetrizedBlock(SymG2, TSelf[no], RowOffset, RowOffset );
      };
     FlipSignOfMagneticColumns(SymG2);

     /*--------------------------------------------------------------*/
     /*- compute the products W21*sym(G1) and W21^{\dagger} * sym(G2)*/
     /*--------------------------------------------------------------*/
     Log("  Multiplication 1...");
     W21->Multiply(SymG1, W21SymG1);

     Log("  Multiplication 2...");
     W21->Adjoint();
     W21->Multiply(SymG2, W21DSymG2);

     // we have to do this again so that W21 will be the correct
     // size again on the next go-round
     W21->Adjoint();

     /*--------------------------------------------------------------*/
     /*- compute the product W21*sym(G1)*W21^{\dagger}*sym(G2)      -*/
     /*--------------------------------------------------------------*/
     Log("  Multiplication 3...");
     W21SymG1->Multiply(W21DSymG2, SymG2);

     /*--------------------------------------------------------------*/
     /*- if we are plotting the flux, then extract the diagonal of   */
     /*- the matrix we just computed and plot it; otherwise, compute */
     /*- the trace and put it in the next slot in the output vector. */
     /*--------------------------------------------------------------*/
     if (PlotFlux)
      { 
        DV->Zero();
        for(nr=0; nr<N2; nr++)
         DV->SetEntry(N1+nr, SymG2->GetEntryD(nr, nr) );
        CreateFluxPlot(SC3D, Omega, Tag);
      }
     else
      { 
        FI[nt] = real(SymG2->GetTrace());

        /***************************************************************/
        /* write the result to the frequency-resolved output file ******/
        /***************************************************************/
        FILE *f=fopen(SC3D->ByOmegaFile, "a");
        fprintf(f,"%s %s %e\n",Tag,z2s(Omega),FI[nt]);
        fclose(f);
      };

     /*--------------------------------------------------------------*/
     /* and untransform the geometry                                 */
     /*--------------------------------------------------------------*/
     G->UnTransform();
     Log(" ...done!");

   }; // for (nt=0; nt<SC3D->NumTransformations... )

  /*--------------------------------------------------------------*/
  /*- at the end of the first successful frequency calculation,  -*/
  /*- we dump out the cache to disk, and then tell ourselves not -*/
  /*- to dump the cache to disk again (explain me)               -*/
  /*--------------------------------------------------------------*/
  if ( SC3D->WriteCache ) 
   { StoreGlobalFIPPICache( SC3D->WriteCache );
     SC3D->WriteCache=0;
   };

}