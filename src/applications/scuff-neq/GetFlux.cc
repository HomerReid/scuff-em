/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This fileis part of SCUFF-EM.
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
 * GetFlux.cc  -- evaluate the four-matrix-trace formulas that give
 *             -- the spectral density of power/momentum flux at a 
 *             -- single frequency
 *
 * homer reid  -- 2/2012
 *
 */

#include "scuff-neq.h"
#include "libscuffInternals.h"

#define II cdouble(0.0,1.0)

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SpyPlot(HMatrix *M, char *Title, char *FileName)
{
  if (!FileName)
   FileName=const_cast<char *>("/tmp/.spyplot.dat");
   
  FILE *f=fopen(FileName,"w");
  if (!f) return;

  cdouble ME;
  for(int nr=0; nr<M->NR; nr++)
   for(int nc=0; nc<M->NC; nc++)
    { ME=M->GetEntry(nr,nc);
      if ( abs(ME) != 0.0 )
       fprintf(f,"%i %i %e %e \n",nr,nc,real(ME),imag(ME));
    };

  fclose(f);

  f=popen("gnuplot -persist","w");
  if (!f) return;
  if (Title)
   fprintf(f,"set title '%s'\n",Title);
  fprintf(f,"set xlabel 'Rows'\n");
  fprintf(f,"set ylabel 'Columns'\n");
  fprintf(f,"plot '%s' w p pt 7 ps 1\n",FileName);
 // pclose(f);
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void UndoSCUFFMatrixTransformation(HMatrix *M)
{ 
  for (int nr=0; nr<M->NR; nr+=2)
   for (int nc=0; nc<M->NC; nc+=2)
    { M->SetEntry(nr,   nc,   ZVAC*M->GetEntry(nr,   nc)   );
      M->SetEntry(nr,   nc+1, -1.0*M->GetEntry(nr,   nc+1) );
      M->SetEntry(nr+1, nc+1, -1.0*M->GetEntry(nr+1, nc+1)/ZVAC );
    };
}

/***************************************************************/
/* the GetFlux() routine computes a large number of quantities.*/
/* each quantity is characterized by values of four indices:   */
/*                                                             */
/*  (a) nt, the geometrical transform                          */
/*  (b) nss the source surface (object)                        */
/*  (c) nsd, the destination surface (object)                  */
/*  (d) nq, the physical quantity (i.e. power, xforce, etc.)   */
/*                                                             */
/* Given values for quantities (a)--(d), this routine computes */
/* a unique index into a big vector storing all the quantities.*/
/***************************************************************/
int GetIndex(SNEQData *SNEQD, int nt, int nss, int nsd, int nq)
{
  int NS = SNEQD->G->NumSurfaces;
  int NQ = SNEQD->NQ;
  int NSNQ = NS*NQ;
  int NS2NQ = NS*NS*NQ;
  return nt*NS2NQ + nss*NSNQ + nsd*NQ + nq; 
}

/***************************************************************/
/* compute the four-matrix-trace formula for the contribution  */
/* of sources inside SourceSurface to the fluxes of power      */
/* and/or momentum through DestSurface.                        */
/* If SourceSurface==DestSurface and SelfTerm==true, then the  */
/* calculation is performed as if the surface were in          */
/* isolation (no other objects in the geometry).               */
/***************************************************************/
void GetTrace(SNEQData *SNEQD, int SourceSurface, int DestSurface,
              double *Results, bool SelfTerm=false)
{
  RWGGeometry *G      = SNEQD->G;

  int DimS            = G->Surfaces[SourceSurface]->NumBFs;
  int OffsetS         = G->BFIndexOffset[SourceSurface];

  int DimD            = G->Surfaces[DestSurface]->NumBFs;
  int OffsetD         = G->BFIndexOffset[DestSurface];
   
  /*--------------------------------------------------------------*/
  /*- Set W_{SD} = S,D subblock of W matrix (SelfTerm==false)     */
  /*-            = inverse T matrix of object S (SelfTerm==true)  */
  /*--------------------------------------------------------------*/
  HMatrix *WSD = new HMatrix(DimS, DimD, LHM_COMPLEX, LHM_NORMAL, SNEQD->Buffer[0]);
  if (SelfTerm)
   { WSD->Copy(SNEQD->T[SourceSurface]);
     UndoSCUFFMatrixTransformation(WSD);
     WSD->LUFactorize();
     WSD->LUInvert();
   }
  else
   SNEQD->W->ExtractBlock(OffsetS, OffsetD, WSD);

  /*--------------------------------------------------------------*/
  /*- set WDOW = W_{SD}^\dagger  * O_S * W_{SD}                  -*/
  /*--------------------------------------------------------------*/
  HMatrix *OW = new HMatrix(DimS, DimD, LHM_COMPLEX, LHM_NORMAL, SNEQD->Buffer[1]);
  if (SNEQD->SymGSource)
   { 
     HMatrix *SymG = new HMatrix(DimS, DimS, LHM_COMPLEX, LHM_NORMAL, SNEQD->Buffer[2]);
     for(int nr=0; nr<DimS; nr++)
      for(int nc=0; nc<DimS; nc++)
       SymG->SetEntry(nr, nc, 0.5*(       SNEQD->TSelf[SourceSurface]->GetEntry(nr,nc)
                                    +conj(SNEQD->TSelf[SourceSurface]->GetEntry(nc,nr))
                                  ));
     SymG->Multiply(WSD, OW);
     delete SymG; 
   }
  else
   { 
     SMatrix *OMatrixS   = SNEQD->SArray[SourceSurface][ 1 + QINDEX_POWER ];
     OMatrixS->Apply(WSD, OW);
   };
  HMatrix *WDOW = new HMatrix(DimD, DimD, LHM_COMPLEX, LHM_NORMAL, SNEQD->Buffer[2]);
  WSD->Multiply(OW, WDOW, "--transA C");

  delete WSD;
  delete OW;

  /*--------------------------------------------------------------*/
  /*- for each quantity requested, compute trace(O*WDOW) where    */
  /*- O is the overlap matrix for the quantity in question        */
  /*--------------------------------------------------------------*/
  int *CIndices;       // column indices
  cdouble *Entries;    // column entries   
  for(int nq=0, QIndex=0; QIndex<MAXQUANTITIES; QIndex++)
   { 
      int QFlag = 1<<QIndex;
      if ( !(SNEQD->QuantityFlags & QFlag) )
       continue;

      // we want to set the self-term for the power to 0
      if ( QFlag==QFLAG_POWER && SelfTerm==true )
       { Results[nq++]=0.0;
         continue;
       };

      //
      if ( QFlag==QFLAG_POWER && SNEQD->SymGDest )
       {
         HMatrix *T=SNEQD->TSelf[DestSurface];
         double FMPTrace=0.0; //'four-matrix-product trace'
         for(int nr=0; nr<DimD; nr++)
          for(int nc=0; nc<DimD; nc++)
           { cdouble SymG = 0.5 * (T->GetEntry(nr,nc) + conj(T->GetEntry(nc,nr)) );
             FMPTrace += real(SymG * WDOW->GetEntry(nc,nr));
           };
         Results[nq++] = (1.0/8.0) * FMPTrace;
         continue;
       };
 
      SMatrix *OMatrixD = SNEQD->SArray[DestSurface][ 1 + QIndex ];
      double FMPTrace=0.0; //'four-matrix-product trace'
      for(int ri=0; ri<DimD; ri++) // 'row index'
       { 
         int nnz=OMatrixD->GetRow(ri, &CIndices, (void **)&Entries);
         for(int nci=0; nci<nnz; nci++)
          FMPTrace += real( Entries[nci] * WDOW->GetEntry(CIndices[nci], ri) );
       };
      Results[nq++] = (-1.0/16.0) * FMPTrace;

   };

 delete WDOW;

} 

/***************************************************************/
/* return false on failure *************************************/
/***************************************************************/
bool CacheRead(SNEQData *SNEQD, cdouble Omega, double *Flux)
{
  if (SNEQD->UseExistingData==false)
   return false;

  FILE *f=vfopen("%s.flux","r",SNEQD->FileBase);
  if (!f) return false;
  Log("Attempting to cache-read flux data for Omega=%e...",real(Omega));

  int NT=SNEQD->NumTransformations;
  int NS=SNEQD->G->NumSurfaces;
  int NQ=SNEQD->NQ;
  int nt, nss, nsd, nq;
  GTComplex **GTCList=SNEQD->GTCList;
  char *FirstTag = GTCList[0]->Tag;
  int ErrorCode, LineNum=0;
  double FileOmega, rOmega=real(Omega);

  char Line[1000];
  char *Tokens[50];
  int NumTokens, MaxTokens=50;
  bool FoundFirstTag=false;
  while( FoundFirstTag==false && fgets(Line,1000,f) )
   { 
     LineNum++;
     NumTokens=Tokenize(Line, Tokens, MaxTokens, " ");
     if (NumTokens<4) continue;
     if (strcmp(Tokens[1],FirstTag)) continue;
     sscanf(Tokens[0],"%lf",&FileOmega);
     if ( fabs(FileOmega - rOmega) > 1.0e-6*rOmega ) continue;
     FoundFirstTag=true;
   
   };
  if(FoundFirstTag==false)
   { ErrorCode=1; goto fail;}

  for(nt=0; nt<NT; nt++)
   for(nss=0; nss<NS; nss++)
    for(nsd=0; nsd<NS; nsd++)
     { 
       if ( !(nt==0 && nss==0 && nsd==0) )
        { if ( !fgets(Line,1000,f) )
           { ErrorCode=2; goto fail; }
          LineNum++;
          NumTokens=Tokenize(Line, Tokens, MaxTokens, " ");
          if ( strcmp(Tokens[1],GTCList[nt]->Tag) )
           { ErrorCode=3; goto fail; }
          sscanf(Tokens[0],"%lf",&FileOmega);
          if ( fabs(FileOmega - rOmega) > 1.0e-6*rOmega )
           { ErrorCode=4; goto fail; }
        };

       if ( NumTokens < 3+NQ ) 
        { ErrorCode=5; goto fail; }
       for(nq=0; nq<NQ; nq++)
        sscanf(Tokens[3+nq],"%le", Flux+GetIndex(SNEQD, nt, nss, nsd, nq) );
     };

  // success:
   Log("...success!");
   fclose(f); 
   return true;

  fail:
   switch(ErrorCode)
    { case 1: Log("could not find first tag (fail)"); break;
      case 2: Log("line %i: unexpected end of file (fail)",LineNum); break;
      case 3: Log("line %i: wrong tag (fail)",LineNum); break;
      case 4: Log("line %i: wrong frequency (fail)",LineNum); break;
      case 5: Log("line %i: too few quantities (fail)",LineNum); break;
    };
   fclose(f); 
   return false;

}

/***************************************************************/
/* the computed quantities are ordered in the output vector    */
/* like this:                                                  */
/*                                                             */
/*  Flux[ nt*NS2NQ + ns*NSNQ + nsp*NQ + nq ]                   */
/*   = contribution of sources inside surface #nsp to flux of  */
/*     quantity #nq into surface #ns, all at transformation #nt*/
/*                                                             */
/*  where    NQ = number of quantities (1--4)                  */
/*  where  NSNQ = number of surface * NQ                       */
/*  where NS2NQ = (number of surfaces)^2* NQ                   */
/***************************************************************/
void GetFlux(SNEQData *SNEQD, cdouble Omega, double *Flux)
{
  if ( CacheRead(SNEQD, Omega, Flux) )
   return;

  Log("Computing neq quantities at omega=%s...",z2s(Omega));

  /***************************************************************/
  /* extract fields from SNEQData structure ************************/
  /***************************************************************/
  RWGGeometry *G      = SNEQD->G;
  HMatrix *W          = SNEQD->W;
  HMatrix **T         = SNEQD->T;
  HMatrix **TSelf     = SNEQD->TSelf;
  HMatrix **U         = SNEQD->U;
  int NQ              = SNEQD->NQ;
  bool *NeedMatrix    = SNEQD->NeedMatrix;
  SMatrix ***SArray   = SNEQD->SArray;

  /***************************************************************/
  /* preinitialize an argument structure for the BEM matrix      */
  /* block assembly routine                                      */
  /***************************************************************/
  GetSSIArgStruct MyGSSIArgStruct, *Args=&MyGSSIArgStruct;
  InitGetSSIArgs(Args);
  Args->G         = G;
  Args->Omega     = Omega;

  /***************************************************************/
  /* before entering the loop over transformations, we first     */
  /* assemble the (transformation-independent) T matrix blocks.  */
  /***************************************************************/
  int NS=G->NumSurfaces;
  for(int ns=0; ns<NS; ns++)
   { 
     if (G->Mate[ns]!=-1)
      { Log(" Block %i is identical to %i (reusing T matrix)",ns,G->Mate[ns]);
        continue;
      }
     else
      Log(" Assembling self contributions to T(%i)...",ns);

     Args->Sa = Args->Sb = G->Surfaces[ns];
     Args->B = T[ns];
     Args->Symmetric=1;
     GetSurfaceSurfaceInteractions(Args);
   };

  /***************************************************************/
  /* also before entering the loop over transformations, we      */
  /* pause to assemble the overlap matrices.                     */
  /***************************************************************/
  for(int ns=0; ns<NS; ns++)
   G->Surfaces[ns]->GetOverlapMatrices(NeedMatrix, SArray[ns], Omega, G->RegionMPs[0]);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (SNEQD->SymGSource || SNEQD->SymGDest)
   { for(int nr=0; nr<G->NumRegions; nr++)
      G->RegionMPs[nr]->Zero();
     for(int ns=0; ns<NS; ns++)
      {
        if (G->Mate[ns]!=-1)
         { Log(" Block %i is identical to %i (reusing TSelf matrix)",ns,G->Mate[ns]);
           continue;
         }
        else
         Log(" Assembling self contributions to TSelf(%i)...",ns);
   
        Args->Sa = Args->Sb = G->Surfaces[ns];
        Args->B = TSelf[ns];
        Args->Symmetric=1;
        G->RegionMPs[ G->Surfaces[ns]->RegionIndices[1] ]->UnZero();
        GetSurfaceSurfaceInteractions(Args);
        G->RegionMPs[ G->Surfaces[ns]->RegionIndices[1] ]->Zero();
        UndoSCUFFMatrixTransformation(TSelf[ns]);
      };
     for(int nr=0; nr<G->NumRegions; nr++)
      G->RegionMPs[nr]->UnZero();
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (NS>50) ErrExit("%s:%i: internal error",__FILE__,__LINE__);
  double SelfContributions[50][MAXQUANTITIES];
  for(int ns=0; ns<NS; ns++)
   { GetTrace(SNEQD, ns, ns, SelfContributions[ns], true);
     Log("ns=%i, Omega=%e: self contributions = (%e,%e,%e,%e)",ns,real(Omega),
          SelfContributions[0], SelfContributions[1], SelfContributions[2], SelfContributions[3]);
   };
       
  /***************************************************************/
  /* now loop over transformations.                              */
  /* note: 'gtc' stands for 'geometrical transformation complex' */
  /***************************************************************/
  char *Tag;
  int RowOffset, ColOffset;
  for(int nt=0; nt<SNEQD->NumTransformations; nt++)
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
     for(int nb=0, ns=0; ns<NS; ns++)
      for(int nsp=ns+1; nsp<NS; nsp++, nb++)
       if ( nt==0 || G->SurfaceMoved[ns] || G->SurfaceMoved[nsp] )
        { 
          Log("  Assembling U(%i,%i)...",ns,nsp);
          Args->Sa = G->Surfaces[ns];
          Args->Sb = G->Surfaces[nsp];
          Args->B  = U[nb];
          Args->Symmetric=0;
          GetSurfaceSurfaceInteractions(Args);
        };

     /*--------------------------------------------------------------*/
     /*- stamp all blocks into the BEM matrix and invert it         -*/
     /*--------------------------------------------------------------*/
     for(int nb=0, ns=0; ns<NS; ns++)
      { 
        RowOffset=G->BFIndexOffset[ns];
        W->InsertBlock(T[ns], RowOffset, RowOffset);

        for(int nsp=ns+1; nsp<NS; nsp++, nb++)
         { ColOffset=G->BFIndexOffset[nsp];
           W->InsertBlock(U[nb], RowOffset, ColOffset);
           W->InsertBlockTranspose(U[nb], ColOffset, RowOffset);
         };
      };
     UndoSCUFFMatrixTransformation(W);
     Log("LU factorizing/inverting...");
     W->LUFactorize();
     W->LUInvert();
     Log("Done with linear algebra...");

     /*--------------------------------------------------------------*/
     /*- compute the requested quantities for all objects -----------*/
     /*- note: nss = 'num surface, source'                          -*/
     /*-       nsd = 'num surface, destination'                     -*/
     /*--------------------------------------------------------------*/
     FILE *f=vfopen("%s.flux","a",SNEQD->FileBase);
     double Quantities[7];
     for(int nss=0; nss<NS; nss++)
      for(int nsd=0; nsd<NS; nsd++)
       { 
         fprintf(f,"%e %s %i%i ",real(Omega),Tag,nss+1,nsd+1);
         GetTrace(SNEQD, nss, nsd, Quantities, false);
         for(int nq=0; nq<NQ; nq++)
          { int Index=GetIndex(SNEQD, nt, nss, nsd, nq);
            Flux[Index] = Quantities[nq]; 
            if ( nss==nsd )
             Flux[Index] -= SelfContributions[nsd][nq]; 
            fprintf(f,"%.8e ",Flux[Index]);
          };
         fprintf(f,"\n");
    
       };
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
