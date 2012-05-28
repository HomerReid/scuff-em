/*
 * FrequencyIntegrand.cc -- evaluate spectral density of power/momentum 
 *                       -- transfer at a single frequency
 *
 * homer reid            -- 2/2012
 *
 */

#include "scuff-neq.h"
#include "libscuffInternals.h"

/***************************************************************/
/***************************************************************/
/***************************************************************/
void UndoMagneticRenormalization(HMatrix *M)
{ 
#if 0
  int nr, nc;
  for (nr=0; nr<B->NR; nr++)
   for (nc=1; nc<B->NC; nc+=2)
    M->SetEntry(nr, nc, -1.0*M->GetEntry(nr, nc) );
#endif
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AssembleOverlapMatrices(SNEQData *SNEQD)
{
  RWGGeometry *G=SNEQD->G;

  RWGObject *O;
  int no, neAlpha, neBeta;
  double Overlaps[5];
  int Offset;

  for(no=0; no<G->NumObjects; no++)
   { O=G->Objects[no];
     for(neAlpha=0; neAlpha<O->NumEdges; neAlpha++)
      for(neBeta=0; neBeta<O->NumEdges; neBeta++)
       { O->GetOverlaps(ne, 
   };

}

/***************************************************************/
/* evaluate the four-matrix trace formula for the given        */
/* quantity on the given object                                */
/***************************************************************/
double GetTrace(SNEQData *SNEQD, int WhichQuantity, int WhichObject)
{
  int qIndex 

  RWGGeometry *G=SNEQD->G;
  int NO=G->NumObjects;

  /***************************************************************/
  /***************************************************************/
  int nr1, Offset1, no2, nr2, Offset2;

  O1 = G->Objects[WhichObject];
  Offset1 = G->BFIndexOffset[WhichObject];
  for(nr1=0; nr1<O1->NumEdges; nr1++)
   for(no2=0; no2<NO; no2++)
    { 
      if (no2==WhichObject) continue;

      O2=G->Objects[no2];
      Offset2 = G->BFIndexOffset[no2];
      for(nr2=0; nr2<O2->NumEdges; nr2++)

    };

} 

/***************************************************************/
/* the computed quantities are ordered in the output vector    */
/* like this:                                                  */
/*  FI[ nt*NONQ + nq*NQ + no ]                                 */
/*   = nqth quantity for noth object under ntth transform      */
/*  where   NQ = number of quantities computed (between 1-4)   */
/*  where NONQ = number of objects * number of quantities      */
/***************************************************************/
void GetFrequencyIntegrand(SNEQData *SNEQD, cdouble Omega, double *FI)
{
  Log("Computing neq quantities at omega=%s...",z2s(Omega));

  /***************************************************************/
  /* extract fields from SNEQData structure ************************/
  /***************************************************************/
  RWGGeometry *G      = SNEQD->G;
  HMatrix *W          = SNEQD->W;
  HMatrix **T         = SNEQD->T;
  HMatrix **U         = SNEQD->U;
  HVector *OPF;       = SNEQD->OPF;
  HVector *OiMF;      = SNEQD->OiMF;
  int PlotFlux        = SNEQD->PlotFlux;
  int WhichQuantities = SNEQD->WhichQuantities;

  /***************************************************************/
  /* initialize overlap matrices *********************************/
  /***************************************************************/
  AssembleOverlapMatrices(SNEQD);

  /***************************************************************/
  /* preinitialize an argument structure for the BEM matrix      */
  /* block assembly routine                                      */
  /***************************************************************/
  ABMBArgStruct MyABMBArgStruct, *Args=&MyABMBArgStruct;
  InitABMBArgs(Args);
  Args->G         = SNEQD->G;
  Args->Omega     = Omega;

  /***************************************************************/
  /* before entering the loop over transformations, we first     */
  /* assemble the (transformation-independent) T matrix blocks.  */
  /***************************************************************/
  int no, nop, nb, nr, NO=G->NumObjects;
  for(no=0; no<G->NumObjects; no++)
   { 
     Log(" Assembling self contributions to T(%i)...",no);

     Args->Oa = Args->Ob = G->Objects[no];
     Args->B = T[no];
     Args->Symmetric=1;
     AssembleBEMMatrixBlock(Args);
   };

  /***************************************************************/
  /* now loop over transformations. ******************************/
  /* note: 'gtc' stands for 'geometrical transformation complex' */
  /***************************************************************/
  int nt;
  char *Tag;
  int RowOffset, ColOffset;
  for(nt=0; nt<SNEQD->NumTransformations; nt++)
   { 
     /*--------------------------------------------------------------*/
     /*- transform the geometry -------------------------------------*/
     /*--------------------------------------------------------------*/
     Tag=SNEQD->GTCList[nt]->Tag;
     G->Transform(SNEQD->GTCList[nt]);
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
          Args->B  = U[nb];
          Args->Symmetric=0;
          AssembleBEMMatrixBlock(Args);
          FlipSignOfMagneticColumns(U[nb]);
        };

     /*--------------------------------------------------------------*/
     /*- stamp all blocks into the BEM matrix and invert it         -*/
     /*--------------------------------------------------------------*/
     for(nb=0, no=0; no<NO; no++)
      { 
        RowOffset=G->BFIndexOffset[no];
        W->InsertBlock(T[no], RowOffset, RowOffset);

        for(nop=no+1; nop<NO; nop++, nb++)
         { ColOffset=G->BFIndexOffset[nop];
           W->InsertBlock(U[nb], RowOffset, ColOffset);
           W->InsertBlockTranspose(U[nb], ColOffset, RowOffset);
         };
      };
     UndoMagneticRenormalization(W);
     W->LUFactorize();
     W->LUInvert();

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     FILE *f=fopen(SNEQD->ByOmegaFile, "a");
     if ( WhichQuantities & QUANTITY_ENERGY )
      FI[ntnq++] = GetTrace(W, 

     /*--------------------------------------------------------------*/
     /* write the result to the frequency-resolved output file ------*/
     /*--------------------------------------------------------------*/
     fprintf(f,"%s %s %e\n",Tag,z2s(Omega),FI[nt]);
     fclose(f);

     /*--------------------------------------------------------------*/
     /* and untransform the geometry                                 */
     /*--------------------------------------------------------------*/
     G->UnTransform();
     Log(" ...done!");

   }; // for (nt=0; nt<SNEQD->NumTransformations... )

  /*--------------------------------------------------------------*/
  /*- at the end of the first successful frequency calculation,  -*/
  /*- we dump out the cache to disk, and then tell ourselves not -*/
  /*- to dump the cache to disk again (explain me)               -*/
  /*--------------------------------------------------------------*/
  if ( SNEQD->WriteCache ) 
   { StoreCache( SNEQD->WriteCache );
     SNEQD->WriteCache=0;
   };

}
