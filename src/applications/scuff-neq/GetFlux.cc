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
/* evaluate the four-matrix trace formula for the contribution */
/* of fluctuations within SourceSurface to the flux of quantity*/
/* QIndex into DestSurface.                                    */
/*                                                             */
/* QIndex      = 0,1,2,3 for power, {x,y,z} momentum           */
/*                                                             */
/* Note: the four-matrix trace is (WD=W^\dagger)               */
/*                                                             */
/*  FMT = tr (O1 * W * O2 * WD )                               */
/*      = \sum_{pqrs} O1_{pq} W_{qr} O2_{rs} WD_{sp}           */
/*      = \sum_{pqrs} O1_{pq} W_{qr} O2_{rs} (W_{ps})^*        */
/*                                                             */
/* because O1 and O2 are sparse, we evaluate the sum by        */
/* looping over the rows of O1 and O2 (which we index by       */
/* p and r respectively). for each row, we extract the         */
/* list of nonzero column indices in each row (that is,        */
/* the nonzero values of q and s) and then evaluate the        */
/* contribution to the sum of (p,q,r,s).                       */
/***************************************************************/
#if 0
double GetTrace(SNEQData *SNEQD, int QIndex, int SourceSurface, int DestSurface)
{
  RWGGeometry *G      = SNEQD->G;
  HMatrix *W          = SNEQD->W;
  int NN              = W->NR;
  SMatrix ***SArray   = SNEQD->SArray;

  int Offset1       = G->BFIndexOffset[DestSurface];
  SMatrix *OMatrix1 = SArray[DestSurface][ 1 + QIndex ];

  int Offset2       = G->BFIndexOffset[SourceSurface];
  SMatrix *OMatrix2 = SArray[SourceSurface][ 1 + QINDEX_POWER];

  int p, q, r, s; 
  int nnzq, nq, nnzs, ns;
  int *qValues, *sValues;
  cdouble *O1Entries, *O2Entries;

  cdouble FMPTrace=0.0; //'four-matrix-product trace'

  for(p=0; p<OMatrix1->NR; p++)
   { 
     nnzq=OMatrix1->GetRow(p, &qValues, (void **)&O1Entries);

     for(r=0; r<OMatrix2->NR; r++)
      { 
        nnzs=OMatrix2->GetRow(r, &sValues, (void **)&O2Entries);

        for(nq=0, q=qValues[0]; nq<nnzq; q=qValues[++nq] )
         for(ns=0, s=sValues[0]; ns<nnzs; s=sValues[++ns] )
          FMPTrace +=  O1Entries[nq]
                      *W->ZM[ (Offset1+q) + NN*(Offset2+r) ]
                      *O2Entries[ns]
                      *conj( W->ZM[ (Offset1+p) + NN*(Offset2+s) ] );
#if 0
          FMPTrace +=  O1Entries[nq]
                      *W->GetEntry( Offset1+q, Offset2+r )
                      *O2Entries[ns]
                      *conj( W->GetEntry( Offset1+p, Offset2+s ));
#endif

      }; // for (r=0...
   }; // for (p=0... 
  FMPTrace *= (-1.0/16.0);

 return real(FMPTrace);

} 
#endif

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

  int DimD            = G->Surfaces[DestSurface]->NumBFs;
  int OffsetD         = G->BFIndexOffset[DestSurface];

  int DimS            = G->Surfaces[SourceSurface]->NumBFs;
  int OffsetS         = G->BFIndexOffset[SourceSurface];
   
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
#if 0
  SMatrix *OMatrixS   = SNEQD->SArray[SourceSurface][ 1 + QINDEX_POWER ];
  HMatrix *OW = new HMatrix(DimS, DimD, LHM_COMPLEX, LHM_NORMAL, SNEQD->Buffer[1]);
  HMatrix *WDOW = new HMatrix(DimS, DimD, LHM_COMPLEX, LHM_NORMAL, SNEQD->Buffer[2]);
  OMatrixS->Apply(WSD, OW);
  WSD->Multiply(OW, WDOW, "--transA C");
#else
  HMatrix *WDOW = new HMatrix(DimS, DimS, LHM_COMPLEX, LHM_NORMAL, SNEQD->Buffer[2]);
  for(int nr=0; nr<DimS; nr++)
   for(int nc=0; nc<DimS; nc++)
    WDOW->SetEntry(nr, nc, 0.5*(      SNEQD->SymG[SourceSurface]->GetEntry(nr,nc)
                                +conj(SNEQD->SymG[SourceSurface]->GetEntry(nc,nr))
                               ));

  HMatrix *OW = new HMatrix(DimS, DimD, LHM_COMPLEX, LHM_NORMAL, SNEQD->Buffer[1]);
  WDOW->Multiply(WSD, OW);
  delete WDOW; 
  WDOW = new HMatrix(DimD, DimD, LHM_COMPLEX, LHM_NORMAL, SNEQD->Buffer[2]);
  WSD->Multiply(OW, WDOW, "--transA C");
#endif
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
  HMatrix **SymG      = SNEQD->SymG;
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
  for(int nr=0; nr<G->NumRegions; nr++)
   G->RegionMPs[nr]->Zero();
  for(int ns=0; ns<NS; ns++)
   {
     if (G->Mate[ns]!=-1)
      { Log(" Block %i is identical to %i (reusing SymG matrix)",ns,G->Mate[ns]);
        continue;
      }
     else
      Log(" Assembling self contributions to SymG(%i)...",ns);

     Args->Sa = Args->Sb = G->Surfaces[ns];
     Args->B = SymG[ns];
     Args->Symmetric=0;
     G->RegionMPs[ G->Surfaces[ns]->RegionIndices[1] ]->UnZero();
     GetSurfaceSurfaceInteractions(Args);
     G->RegionMPs[ G->Surfaces[ns]->RegionIndices[1] ]->Zero();
     UndoSCUFFMatrixTransformation(SymG[ns]);
   };
  for(int nr=0; nr<G->NumRegions; nr++)
   G->RegionMPs[nr]->UnZero();

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (NS>50) ErrExit("%s:%i: internal error",__FILE__,__LINE__);
  double SelfContributions[50][MAXQUANTITIES];
  for(int ns=0; ns<NS; ns++)
   if ( G->Mate[ns]==-1 )
{
    GetTrace(SNEQD, ns, ns, SelfContributions[ns], true);
Log("ns=%i, Omega=%e: self contributions = (%e,%e,%e,%e)",ns,real(Omega),
SelfContributions[0], SelfContributions[1], 
SelfContributions[2], SelfContributions[3]);
}
   else
    memcpy(SelfContributions[ns], SelfContributions[G->Mate[ns]], 
           MAXQUANTITIES*sizeof(double));
       
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
     double Quantities[4];
     int nfc=0;
     for(int nss=0; nss<NS; nss++)
      for(int nsd=0; nsd<NS; nsd++)
       { 
         fprintf(f,"%e %s %i%i ",real(Omega),Tag,nss+1,nsd+1);
         GetTrace(SNEQD, nss, nsd, Quantities, false);
         for(int nq=0; nq<NQ; nq++)
          { Flux[nfc] = Quantities[nq]; 
            if (nss==nsd) Flux[nfc] -= SelfContributions[nsd][nq]; 
            fprintf(f,"%.8e ",Flux[nfc]);
            nfc++;
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
