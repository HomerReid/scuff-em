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
/* the GetFlux() routine computes a large number of flux       */
/* quantities, including both spatially-integrated (SI) and    */
/* spatially-resolved (SR) quantities. these fluxes are stored */
/* in big one-dimensional arrays for passage to numerical      */
/* cubature codes. The following two routines compute unique   */
/* indices into these vectors for individual quantities.       */
/*                                                             */ 
/* SI quantities are characterized by values of four indices:  */
/*                                                             */
/*  (a) nt, the geometrical transform                          */
/*  (b) nss, the source surface (object)                       */
/*  (c) nsd, the destination surface (object)                  */
/*  (d) nq, the physical quantity (i.e. power, xforce, etc.)   */
/*                                                             */
/* SR quantities are characterized by values of four indices:  */
/*                                                             */
/*  (a) nt, the geometrical transform                          */
/*  (b) nss, the source surface (object)                       */
/*  (c) nx, the index of the evaluation point in the EP file   */
/*  (d) nq, the physical quantity (i.e. power, xforce, etc.)   */
/***************************************************************/
int GetSIQIndex(SNEQData *SNEQD, int nt, int nss, int nsd, int nq)
{
  int NS = SNEQD->G->NumSurfaces;
  int NQ = SNEQD->NQ;
  return nt*(NS*NS*NQ) + nss*(NS*NQ) + nsd*NQ + nq;
}

int GetSRQIndex(SNEQData *SNEQD, int nt, int nss, int nx, int nq)
{
  int NumSIQs = SNEQD->NumSIQs;

  int NS = SNEQD->G->NumSurfaces;
  int NX = SNEQD->NX;
  int NQ = SNEQD->NQ;

  return NumSIQs + nt*(NS*NX*NQ) + nss*(NX*NQ) + nx*NQ + nq;
}

/***************************************************************/
/* Compute the Sigma matrix for sources contained in body      */
/* SourceSurface. The matrix is stored in the Sigma field of   */
/* the SNEQD structure.                                        */
/***************************************************************/
void ComputeSigmaMatrix(SNEQData *SNEQD, int SourceSurface)
{
  RWGGeometry *G = SNEQD->G;
  HMatrix *W     = SNEQD->W;
  HMatrix *Sigma = SNEQD->Sigma;
  void **Buffer  = SNEQD->Buffer;

  int s=SourceSurface;
  int NBFS = G->Surfaces[s]->NumBFs;
  int OffsetS = G->BFIndexOffset[s];

  int NS=SNEQD->G->NumSurfaces;
  for(int a=0; a<NS; a++)
   for(int b=0; b<NS; b++)
    {
      // extract (a,s) and (b,s) subblocks of W
      int NBFA    = G->Surfaces[a]->NumBFs;
      int OffsetA = G->BFIndexOffset[a];
  
      int NBFB    = G->Surfaces[b]->NumBFs;
      int OffsetB = G->BFIndexOffset[b];

      // get a,s subblock of W 
      HMatrix Was(NBFA, NBFS, LHM_COMPLEX, LHM_NORMAL, Buffer[0]);
      W->ExtractBlock(OffsetA, OffsetS, &Was);

      // get Sym G_{s} / 4.
      // (The factor of 0.5 here is present in the definition
      //  of SymG, while the 0.25 prefactor appears in the 
      //  definition of Sigma; we roll it into SymG for convenience.
      HMatrix SymGs(NBFS, NBFS, LHM_COMPLEX, LHM_NORMAL, Buffer[1]);
      HMatrix *Ts=SNEQD->TSelf[SourceSurface];
      for(int nr=0; nr<NBFS; nr++)
       for(int nc=0; nc<NBFS; nc++)
        SymGs.SetEntry(nr, nc, 0.25*0.5*(        Ts->GetEntry(nr,nc)
                                          +conj( Ts->GetEntry(nc,nr) )
                                        ) 
                      ); 

      // Temp1 <- W_{a,s}*SymG_{s}
      HMatrix Temp1(NBFA, NBFS, LHM_COMPLEX, LHM_NORMAL, Buffer[2]);
      Was.Multiply(&SymGs, &Temp1);
  //    Was.Multiply(&SymGs, &Temp1, "--TransA C");

      // get b,s subblock of W 
      HMatrix Wbs(NBFB, NBFS, LHM_COMPLEX, LHM_NORMAL, Buffer[0]);
      W->ExtractBlock(OffsetB, OffsetS, &Wbs);

      // Temp2 <- W_{a,s}*SymG_{s}*(W_{b,s}^\dagger)
      HMatrix Temp2(NBFA, NBFB, LHM_COMPLEX, LHM_NORMAL, Buffer[1]);
      Temp1.Multiply(&Wbs, &Temp2, "--TransB C");
  //    Temp1.Multiply(&Wbs, &Temp2);

      Sigma->InsertBlock(&Temp2, OffsetA, OffsetB);
    };

}

/***************************************************************/
/* evaluate trace formulas for the contribution of sources     */
/* inside SourceSurface to the fluxes of all spatially-        */
/* integrated quantities.                                      */
/*                                                             */
/* somewhat confusing: the PFT routines in libscuff compute    */
/* all 7 PFT quantities, but the SIFlux[] output vector needs  */
/* to be filled in with just the number of *requested*         */
/* quantities.                                                 */
/***************************************************************/
#define METHOD_DSIPFT 0
#define METHOD_OPFT   1
#define METHOD_EPPFT  2
void GetSIFlux(SNEQData *SNEQD, int SourceSurface, int DestSurface,
               cdouble Omega, double SIFlux[NUMPFT])
{
  if ( (SourceSurface==DestSurface) && SNEQD->OmitSelfTerms )
   { memset(SIFlux, 0, NUMPFT*sizeof(double));
     return;
   };

  RWGGeometry *G  = SNEQD->G;
  HMatrix *Sigma  = SNEQD->Sigma;
  bool DSISelf    = SNEQD->DSISelf;
  bool DSIOther   = SNEQD->DSIOther;
  bool EPOther    = SNEQD->EPOther;

  double **ByEdge=0;
  if (SNEQD->ByEdge)
   ByEdge=SNEQD->ByEdge[DestSurface];

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int Method;
  if (SourceSurface==DestSurface)
   Method = DSISelf ? METHOD_DSIPFT : METHOD_EPPFT;
  else
   Method = DSIOther ? METHOD_DSIPFT : ( EPOther ? METHOD_EPPFT : METHOD_OPFT );

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double AllFlux[NUMPFT];
  switch(Method)
   { 
     case METHOD_DSIPFT:
      G->GetDSIPFTTrace(DestSurface, Omega, 0, Sigma, AllFlux, ByEdge,
                        SNEQD->DSIMesh, SNEQD->DSIRadius, SNEQD->DSIPoints,
                        SNEQD->Lebedev, SNEQD->FarField);
      break;

     case METHOD_EPPFT:
      G->GetEPPFTTrace(DestSurface, Omega, 0, Sigma, AllFlux, ByEdge,
                       (SourceSurface==DestSurface ? true : false));
      break;

     case METHOD_OPFT:
      G->GetOPFTTrace(DestSurface, Omega, 0, Sigma, AllFlux, ByEdge);
      break;

   };

  /*--------------------------------------------------------------*/
  /*- collapse the full vector of 7 PFTs to just the entries the -*/
  /*- user requested                                             -*/
  /*--------------------------------------------------------------*/
  for(int nq=0, nqq=0; nq<NUMPFT; nq++)
   if (SNEQD->NeedQuantity[nq])
    SIFlux[nqq++] = -4.0*AllFlux[nq];

  /*--------------------------------------------------------------*/
  /*- generate panel-resolved flux plots if that was requested   -*/
  /*--------------------------------------------------------------*/
  if (ByEdge)
   { 
     char FileName[100];
     snprintf(FileName,100,"%s.PFTFlux.pp",G->Surfaces[DestSurface]->Label);

     for(int nq=0; nq<NUMPFT; nq++)
      if (ByEdge[nq])
       G->Surfaces[DestSurface]->PlotScalarDensity(ByEdge[nq],
                                                   FileName,
                                                   "%s_%g",
                                                   QuantityNames[nq],
                                                   real(Omega));
   };

} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetSRFluxes(SNEQData *SNEQD, int SourceSurface, HMatrix *XMatrix,
                 cdouble Omega, double *Results)
{ 
 // for(int nx=0; nx<
}

/***************************************************************/
/* return false on failure *************************************/
/***************************************************************/
bool CacheRead(SNEQData *SNEQD, cdouble Omega, double *kBloch, double *Flux)
{
  if (SNEQD->UseExistingData==false)
   return false;

  // i haven't yet implemented caching of spatially-resolved data
  if (SNEQD->NumSRQs>0)
   return false; 

  FILE *f=vfopen("%s.SIFlux","r",SNEQD->FileBase);
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
        sscanf(Tokens[3+nq],"%le", Flux+GetSIQIndex(SNEQD, nt, nss, nsd, nq) );
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
/* for spatially-integrated fluxes, the computed quantities    */
/* are ordered in the output vector like this:                 */ 
/*                                                             */
/*  SIFlux[ nt*NS2NQ + ns*NSNQ + nsp*NQ + nq ]                 */
/*   = contribution of sources inside surface #nsp to flux of  */
/*     quantity #nq into surface #ns, all at transformation #nt*/
/*                                                             */
/*  where    NQ = number of quantities (1--7)                  */
/*  where  NSNQ = number of surface * NQ                       */
/*  where NS2NQ = (number of surfaces)^2* NQ                   */
/*                                                             */
/* for spatially-resolved fluxes, the computed quantities are  */
/* ordered in the output vector like this:                     */
/*                                                             */
/*  SRFlux[ nt*NPNSNQ + np*NSNQ +  ns*NQ + nq ]                */
/*   = contribution of sources inside surface #ns to flux of   */
/*     quantity #nq at evaluation point np, all at             */
/*     transformation NT                                       */
/*                                                             */
/*  where    NQ  = number of quantities (1--7)                 */
/*  where  NSNQ  = number of surfaces * NQ                     */
/*  where NPNSNQ = number of evaluation points * NPNSNQ        */
/***************************************************************/
void GetFlux(SNEQData *SNEQD, cdouble Omega, double *kBloch, double *Flux)
{
  SetLogFileName("scuff-neq.log");

  if ( CacheRead(SNEQD, Omega, kBloch, Flux) )
   return;

  /***************************************************************/
  /* extract fields from SNEQData structure **********************/
  /***************************************************************/
  RWGGeometry *G      = SNEQD->G;
  HMatrix *W          = SNEQD->W;
  HMatrix **T         = SNEQD->T;
  HMatrix **TSelf     = SNEQD->TSelf;
  HMatrix **U         = SNEQD->U;
  int NS              = SNEQD->G->NumSurfaces;
  int NQ              = SNEQD->NQ;
  int NX              = SNEQD->NX;
  int NumSIQs         = SNEQD->NumSIQs;
  char *FileBase      = SNEQD->FileBase;
  bool *NeedQuantity  = SNEQD->NeedQuantity;

  Log("Computing neq quantities at omega=%s...",z2s(Omega));

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
  for(int ns=0; ns<NS; ns++)
   { 
     if (G->Mate[ns]!=-1)
      { Log(" Block %i is identical to %i (reusing T matrix)",ns,G->Mate[ns]);
        continue;
      }
     else
      Log(" Assembling self contributions to T(%i)...",ns);

     G->AssembleBEMMatrixBlock(ns, ns, Omega, kBloch, T[ns]);
   };

  /***************************************************************/
  /* assemble the TSelf matrices *********************************/
  /***************************************************************/
  for(int nr=0; nr<G->NumRegions; nr++)
   G->RegionMPs[nr]->Zero();
  for(int ns=0; ns<NS; ns++)
   {
     if (G->Mate[ns]!=-1)
      { Log(" Block %i is identical to %i (reusing TSelf matrix)",ns,G->Mate[ns]);
        continue;
      }
     else
      Log(" Assembling self contributions to TSelf(%i)...",ns);

     G->RegionMPs[ G->Surfaces[ns]->RegionIndices[1] ]->UnZero();
     G->AssembleBEMMatrixBlock(ns, ns, Omega, kBloch, TSelf[ns]);
     UndoSCUFFMatrixTransformation(TSelf[ns]);
     G->RegionMPs[ G->Surfaces[ns]->RegionIndices[1] ]->Zero();

   }; // for(int ns=0; ns<NS; ns++)
  for(int nr=0; nr<G->NumRegions; nr++)
   G->RegionMPs[nr]->UnZero();

  /***************************************************************/
  /* now loop over transformations.                              */
  /* note: 'gtc' stands for 'geometrical transformation complex' */
  /***************************************************************/
  for(int nt=0; nt<SNEQD->NumTransformations; nt++)
   { 
     /*--------------------------------------------------------------*/
     /*- insert here code to skip this transformation if the integral*/
     /*- has already converged                                       */
     /*--------------------------------------------------------------*/

     /*--------------------------------------------------------------*/
     /*- transform the geometry -------------------------------------*/
     /*--------------------------------------------------------------*/
     char *Tag=SNEQD->GTCList[nt]->Tag;
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
        G->AssembleBEMMatrixBlock(ns, nsp, Omega, kBloch, U[nb]);

     /*--------------------------------------------------------------*/
     /*- stamp all blocks into the BEM matrix and invert it         -*/
     /*--------------------------------------------------------------*/
     for(int nb=0, ns=0; ns<NS; ns++)
      { 
        int RowOffset=G->BFIndexOffset[ns];
        W->InsertBlock(T[ns], RowOffset, RowOffset);

        for(int nsp=ns+1; nsp<NS; nsp++, nb++)
         { int ColOffset=G->BFIndexOffset[nsp];
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
     /*- compute the requested quantities for all objects           -*/
     /*- note: nss = 'num surface, source'                          -*/
     /*-       nsd = 'num surface, destination'                     -*/
     /*--------------------------------------------------------------*/
     for(int nss=0; nss<NS; nss++)
      {
        // compute the Sigma matrix for this source
        ComputeSigmaMatrix(SNEQD, nss);

        // compute spatially-integrated flux quantities for
        // all destination objects
        if (NumSIQs > 0)
         { 
           FILE *f=vfopen("%s.SIFlux","a",FileBase);
           for(int nsd=0; nsd<NS; nsd++)
            { 
              double SIFlux[7];
              GetSIFlux(SNEQD, nss, nsd, Omega, SIFlux);

              fprintf(f,"%e %s ",real(Omega),Tag);
              if (kBloch) fprintf(f,"%e %e ",kBloch[0],kBloch[1]);
              fprintf(f,"%i%i ",nss+1,nsd+1);
              for(int nq=0; nq<NQ; nq++)
               { int Index= GetSIQIndex(SNEQD, nt, nss, nsd, nq);
                 Flux[Index]=SIFlux[nq];
                 fprintf(f,"%.8e ",Flux[Index]);
               };
              fprintf(f,"\n");
            };
           fclose(f);
         };

        // compute spatially-resolved flux quantities for
        // all evaluation points
        //if (NumSRQs > 0)
        //  { GetSRFlux(
        // };
      };

     /*--------------------------------------------------------------*/
     /* untransform the geometry                                     */
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

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetFlux(SNEQData *SNEQD, cdouble Omega, double *Flux)
 { GetFlux(SNEQD, Omega, 0, Flux); }
