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
void ProcessDRMatrix(SNEQData *SNEQD,
                        cdouble Omega,
                        int SourceSurface)
{ 
  RWGGeometry *G    = SNEQD->G;
  HMatrix *DRMatrix = SNEQD->DRMatrix;
  void **Buffer     = SNEQD->Buffer;

  int N = G->Surfaces[SourceSurface]->NumBFs;

  HMatrix _RCopy(N, N, LHM_COMPLEX, LHM_NORMAL, Buffer[0]);
  HMatrix _U(N, N, LHM_COMPLEX, LHM_NORMAL, Buffer[1]);
  HVector _Lambda(N, LHM_REAL, Buffer[2]);
  HMatrix *RCopy  = &_RCopy;
  HMatrix *U      = &_U;
  HVector *Lambda = &_Lambda;
  RCopy->Copy(DRMatrix);
  Log("Computing eigenvectors of dressed Rytov matrix...");
  RCopy->Eig(Lambda, U);

  char FileName[100];
  snprintf(FileName,100,"%s.Lambda",GetFileBase(G->GeoFileName));
  Lambda->ExportToText(FileName);

  HVector _KN(N, LHM_COMPLEX, Buffer[0]);
  HVector *KN=&_KN;
  for(int nv=0; nv<SNEQD->PlotRytovVectors; nv++)
   { 
     Log("Plotting dressed-Rytov surface currents (%i)...",nv);
     double RtLambda=sqrt(fabs(Lambda->GetEntryD(nv)));
     for(int m=0; m<N; m++)
      KN->SetEntry(m, RtLambda*U->GetEntry(m,nv));
     G->PlotSurfaceCurrents(KN, Omega, "%s.%s.Rytov%i.pp",
                            GetFileBase(G->GeoFileName),z2s(Omega),nv);
   };

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

void ApplySCUFFMatrixTransformation(HMatrix *M)
{ 
  for (int nr=0; nr<M->NR; nr+=2)
   for (int nc=0; nc<M->NC; nc+=2)
    { M->SetEntry(nr,   nc, (1.0/ZVAC)*M->GetEntry(nr,   nc)    );
      M->SetEntry(nr,   nc+1, -1.0*M->GetEntry(nr,   nc+1)      );
      M->SetEntry(nr+1, nc+1, -1.0*ZVAC*M->GetEntry(nr+1, nc+1) );
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
/*  (d) nfc, the index of the flux component                   */
/*      =0,1,2,3,4,...,11 for Px, Py, Pz, T_xx, T_xy, ... T_zz */
/***************************************************************/
int GetSIQIndex(SNEQData *SNEQD, int nt, int nss, int nsd, int nq)
{
  int NS = SNEQD->G->NumSurfaces;
  int NQ = SNEQD->NPFT;
  return nt*(NS*NS*NQ) + nss*(NS*NQ) + nsd*NQ + nq;
}

int GetSRQIndex(SNEQData *SNEQD, int nt, int nss, int nx, int nfc)
{
  int NumSIQs = SNEQD->NumSIQs;

  int NS  = SNEQD->G->NumSurfaces;
  int NX  = SNEQD->NX;
  int NFC = NUMSRFLUX;

  return NumSIQs + nt*(NS*NX*NFC) + nss*(NX*NFC) + nx*NFC + nfc;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define AVAC (1.0/ZVAC) // admittance of vacuum
void TSelfToSymG(HMatrix *TSelf, HMatrix *SymG)
{
  int NE = TSelf->NR / 2;
  for(int nea=0; nea<NE; nea++)
   for(int neb=0; neb<NE; neb++)
    { 
       cdouble TEE_NRNC =  ZVAC*TSelf->GetEntry(2*nea + 0 ,2*neb + 0 );
       cdouble TEH_NRNC =  -1.0*TSelf->GetEntry(2*nea + 0 ,2*neb + 1 );
       cdouble THE_NRNC =       TSelf->GetEntry(2*nea + 1 ,2*neb + 0 );
       cdouble THH_NRNC = -AVAC*TSelf->GetEntry(2*nea + 1 ,2*neb + 1 );

       cdouble TEE_NCNR =  ZVAC*TSelf->GetEntry(2*nea + 0 ,2*neb + 0 );
       cdouble TEH_NCNR =  -1.0*TSelf->GetEntry(2*nea + 0 ,2*neb + 1 );
       cdouble THE_NCNR =       TSelf->GetEntry(2*nea + 1 ,2*neb + 0 );
       cdouble THH_NCNR = -AVAC*TSelf->GetEntry(2*nea + 1 ,2*neb + 1 );

       SymG->SetEntry(2*nea+0, 2*neb+0, 0.25*0.5*(TEE_NRNC + conj(TEE_NCNR)));
       SymG->SetEntry(2*nea+0, 2*neb+1, 0.25*0.5*(TEH_NRNC + conj(TEH_NCNR)));
       SymG->SetEntry(2*nea+1, 2*neb+0, 0.25*0.5*(THE_NRNC + conj(THE_NCNR)));
       SymG->SetEntry(2*nea+1, 2*neb+1, 0.25*0.5*(THH_NRNC + conj(THH_NCNR)));
    };
}

/***************************************************************/
/* Compute the dressed Rytov matrix for sources contained in   */
/* SourceSurface. The matrix is stored in the DRMatrix         */
/* field of the SNEQD structure.                               */
/***************************************************************/
void ComputeDRMatrix(SNEQData *SNEQD, int SourceSurface)
{
  RWGGeometry *G    = SNEQD->G;
  HMatrix *W        = SNEQD->W;
  HMatrix *DRMatrix = SNEQD->DRMatrix;
  void **Buffer     = SNEQD->Buffer;

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
      //  definition of Rytov; we roll it into SymG for convenience.
      HMatrix SymGs(NBFS, NBFS, LHM_COMPLEX, LHM_NORMAL, Buffer[1]);
#if 0
      HMatrix *Ts=SNEQD->TInt[SourceSurface];
#endif
//FIXME
#if 1
HMatrix TSelf(NBFS, NBFS, LHM_COMPLEX, LHM_NORMAL, Buffer[2]);
HMatrix *Ts=&TSelf;
Ts->Copy(SNEQD->TInt[SourceSurface]);
UndoSCUFFMatrixTransformation(Ts);
      for(int nr=0; nr<NBFS; nr++)
       for(int nc=0; nc<NBFS; nc++)
        SymGs.SetEntry(nr, nc, 0.25*0.5*(        Ts->GetEntry(nr,nc)
                                          +conj( Ts->GetEntry(nc,nr) )
                                        ) 
                      );
#else
TSelfToSymG(SNEQD->TInt[SourceSurface], SymGs);
#endif


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

      DRMatrix->InsertBlock(&Temp2, OffsetA, OffsetB);
    };

}

/***************************************************************/
/* evaluate trace formulas for the contribution of sources     */
/* inside SourceSurface to the fluxes of all spatially-        */
/* integrated quantities on all surfaces.                      */
/***************************************************************/
void GetSIFlux(SNEQData *SNEQD, int SourceSurface, cdouble Omega,
               int PFTMethod, HMatrix *PFTMatrix)
{
  RWGGeometry *G      = SNEQD->G;
  int NumSurfaces     = G->NumSurfaces;
  HMatrix *DRMatrix   = SNEQD->DRMatrix;
  bool OmitSelfTerms  = SNEQD->OmitSelfTerms;

  /*--------------------------------------------------------------*/
  /*- initialize PFT options structure ---------------------------*/
  /*--------------------------------------------------------------*/
  PFTOptions *PFTOpts = &(SNEQD->PFTOpts);
  PFTOpts->DRMatrix   = DRMatrix;

  if (PFTMethod>10)
   { PFTOpts->PFTMethod = SCUFF_PFT_DSI;
     PFTOpts->DSIPoints = PFTMethod;
   }
  else
   PFTOpts->PFTMethod  = PFTMethod;

  char FFNBuffer[200];
  if (SNEQD->PlotFlux)
   PFTOpts->FluxFileName=FFNBuffer;

  /*--------------------------------------------------------------*/
  /*- do the PFT calculation           ---------------------------*/
  /*--------------------------------------------------------------*/
  if (PFTMethod==SCUFF_PFT_EMT)
   { 
     G->GetPFTMatrix(0, Omega, PFTOpts, PFTMatrix);
   }
  else
   {
     for(int DestSurface=0; DestSurface<NumSurfaces; DestSurface++)
      {
        if ( (SourceSurface==DestSurface) && OmitSelfTerms )
         continue;

        if(SNEQD->PlotFlux)
         snprintf(FFNBuffer,200,"%s.%sTo%s.PFTFlux.pp",
                                 GetFileBase(G->GeoFileName),
                                 G->Surfaces[SourceSurface]->Label,
                                 G->Surfaces[DestSurface]->Label);

        double Flux[NUMPFT];
        PFTOpts->TInterior=SNEQD->TInt[DestSurface];
        PFTOpts->TExterior=SNEQD->TExt[DestSurface];
        G->GetPFT(DestSurface, 0, Omega, Flux, PFTOpts);
        if (PFTMethod==SCUFF_PFT_DSI)
         Flux[PFT_PSCAT]=-Flux[PFT_PABS];
        PFTMatrix->SetEntriesD(DestSurface, ":", Flux);
      };
   };

  PFTMatrix->Scale(-4.0); // where does this factor come from?

} 

/***************************************************************/
/* return false on failure *************************************/
/***************************************************************/
bool CacheRead(SNEQData *SNEQD, cdouble Omega, double *kBloch, double *Flux)
{
  (void)kBloch;

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
  int NQ=SNEQD->NPFT;
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
/***************************************************************/
/***************************************************************/
bool DoDSIAtThisFrequency(SNEQData *SNEQD, cdouble Omega)
{
  HVector *OmegaPoints=SNEQD->DSIOmegaPoints;
  if (OmegaPoints==0)
   return true;
  for(int n=0; n<OmegaPoints->N; n++)
   if (EqualFloat(Omega, OmegaPoints->GetEntry(n)))
    return true;
  return false;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetFlux(SNEQData *SNEQD, cdouble Omega, double *kBloch, double *Flux)
{
  if ( CacheRead(SNEQD, Omega, kBloch, Flux) )
   return;

  /***************************************************************/
  /* extract fields from SNEQData structure **********************/
  /***************************************************************/
  RWGGeometry *G      = SNEQD->G;
  HMatrix *W          = SNEQD->W;
  HMatrix **TExt      = SNEQD->TExt;
  HMatrix **TInt      = SNEQD->TInt; 
  HMatrix **U         = SNEQD->U;
  int NS              = SNEQD->G->NumSurfaces;
  int NX              = SNEQD->NX;
  int NumSRQs         = SNEQD->NumSRQs;
  char *FileBase      = SNEQD->FileBase;

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
  for(int nr=0; nr<G->NumRegions; nr++)
   G->RegionMPs[nr]->Zero();
  for(int ns=0; ns<NS; ns++)
   { 
     if (G->Mate[ns]!=-1)
      { Log(" Block %i is identical to %i (reusing T matrices)",ns,G->Mate[ns]);
        continue;
      }
     else
      Log(" Assembling self contributions to T(%i)...",ns);

     G->RegionMPs[ G->Surfaces[ns]->RegionIndices[1] ]->UnZero();
     G->AssembleBEMMatrixBlock(ns, ns, Omega, kBloch, TInt[ns]);
     G->RegionMPs[ G->Surfaces[ns]->RegionIndices[1] ]->Zero();

     G->RegionMPs[ G->Surfaces[ns]->RegionIndices[0] ]->UnZero();
     G->AssembleBEMMatrixBlock(ns, ns, Omega, kBloch, TExt[ns]);
     G->RegionMPs[ G->Surfaces[ns]->RegionIndices[0] ]->Zero();

   };
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
Log("...SN done with ABMB");

     /*--------------------------------------------------------------*/
     /*- stamp all blocks into the BEM matrix and invert it         -*/
     /*--------------------------------------------------------------*/
     for(int nb=0, ns=0; ns<NS; ns++)
      { 
        int RowOffset=G->BFIndexOffset[ns];
        W->InsertBlock(TInt[ns], RowOffset, RowOffset);
        W->AddBlock(TExt[ns], RowOffset, RowOffset);

        for(int nsp=ns+1; nsp<NS; nsp++, nb++)
         { int ColOffset=G->BFIndexOffset[nsp];
           W->InsertBlock(U[nb], RowOffset, ColOffset);
           W->InsertBlockTranspose(U[nb], ColOffset, RowOffset);
         };
      };
     UndoSCUFFMatrixTransformation(W);
     Log("LU factorizing...");
     W->LUFactorize();
     Log("LU inverting...");
     W->LUInvert();
     Log("Done with linear algebra...");

     /*--------------------------------------------------------------*/
     /*- compute the requested quantities for all objects           -*/
     /*- note: nss = 'num surface, source'                          -*/
     /*-       nsd = 'num surface, destination'                     -*/
     /*--------------------------------------------------------------*/
     int NumPFTMethods  = SNEQD->NumPFTMethods;
     int *PFTMethods    = SNEQD->PFTMethods;
     HMatrix *PFTMatrix = SNEQD->PFTMatrix;
     for(int nss=0; nss<NS; nss++)
      {
        if (SNEQD->OmitZeroTemperatureFlux && SNEQD->TSurfaces[nss]==0.0)
         continue;

        // compute the "dressed Rytov" matrix for this source
Log("SN computing DR matrix");
        ComputeDRMatrix(SNEQD, nss);
Log("...done with DR matrix");
        if (SNEQD->PlotRytovVectors)
         ProcessDRMatrix(SNEQD, Omega, nss);

        // compute spatially-integrated flux quantities for
        // all destination objects using all requested 
        // calculation methods
        for(int npm=0; npm<NumPFTMethods; npm++)
         { 
           if (    PFTMethods[npm]==SCUFF_PFT_DSI 
                && !DoDSIAtThisFrequency(SNEQD, Omega) 
              ) continue;

           GetSIFlux(SNEQD, nss, Omega, PFTMethods[npm], PFTMatrix);

           FILE *f=vfopen(SNEQD->SIFluxFileNames[npm],"a");
           for(int nsd=0; nsd<NS; nsd++)
            { fprintf(f,"%s %e ",Tag,real(Omega));
              if (kBloch) fprintf(f,"%e %e ",kBloch[0],kBloch[1]);
              fprintf(f,"%i%i ",nss+1,nsd+1);
              for(int nq=0; nq<NUMPFT; nq++)
               fprintf(f,"%+.8e ",PFTMatrix->GetEntryD(nsd,nq));
              fprintf(f,"\n");
            };
           fclose(f);

           if (npm==0)
            for(int nsd=0; nsd<NS; nsd++)
             for(int nq=0; nq<NUMPFT; nq++)
              if (SNEQD->NeedQuantity[nq])
               Flux[ GetSIQIndex(SNEQD, nt, nss, nsd, nq) ]
                = PFTMatrix->GetEntryD(nsd, nq);
         };

        // compute spatially-resolved flux quantities for
        // all evaluation points
        if (NumSRQs > 0)
          { 
            HMatrix *SRXMatrix = SNEQD->SRXMatrix;
            HMatrix *SRFMatrix = SNEQD->SRFMatrix;
            HMatrix *DRMatrix  = SNEQD->DRMatrix;
            GetSRFluxTrace(G, SRXMatrix, Omega, DRMatrix, SRFMatrix);

            FILE *f=vfopen("%s.SRFlux","a",FileBase);
            for(int nx=0; nx<NX; nx++)
             {
               double X[3], SRFlux[NUMSRFLUX];
               SRXMatrix->GetEntriesD(nx,":",X);
               SRFMatrix->GetEntriesD(nx,":",SRFlux);

               fprintf(f,"%s %e ",Tag,real(Omega));
               if (kBloch) 
                fprintf(f,"%e %e ",kBloch[0],kBloch[1]);
               fprintf(f,"%e %e %e %i ",X[0],X[1],X[2],nss);
               for(int nfc=0; nfc<NUMSRFLUX; nfc++)
                { int Index=GetSRQIndex(SNEQD, nt, nss, nx, nfc); 
                  Flux[Index]=SRFMatrix->GetEntryD(nx,nfc);
                  fprintf(f,"%e ",Flux[Index]);
                };
               fprintf(f,"\n");
             };
            fclose(f);
          };

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
