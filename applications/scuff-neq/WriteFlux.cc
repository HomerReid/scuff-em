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
/* Compute the dressed Rytov matrix for sources contained in   */
/* SourceSurface. The matrix is stored in the DRMatrix         */
/* field of the SNEQD structure.                               */
/***************************************************************/
#define RYTOVPF (-4.0/M_PI)
void ComputeDRMatrix(SNEQData *SNEQD, int SourceSurface)
{
  Log("...computing DR matrix");

  RWGGeometry *G  = SNEQD->G;
  HMatrix *M      = SNEQD->M;
  HMatrix *DR     = SNEQD->DRMatrix;

  int NBFS        = G->Surfaces[SourceSurface]->NumBFs;
  int OffsetS     = G->BFIndexOffset[SourceSurface];
  HMatrix *TInt   = SNEQD->TInt[SourceSurface];

  /***************************************************************/
  /* stamp Sym(T_s) = (T_s + T_s^\dagger) / 2                    */
  /* into the sth diagonal block of DRMatrix,                    */
  /* undoing the SCUFF matrix transformation along the way.      */
  /***************************************************************/
  DR->Zero();
  for(int a=0; a<(NBFS/2); a++)
   for(int b=a; b<(NBFS/2); b++)
    {
      cdouble TEEab = ZVAC * TInt->GetEntry(2*a+0, 2*b+0);
      cdouble TEMab = -1.0 * TInt->GetEntry(2*a+0, 2*b+1);
      cdouble TMEab =        TInt->GetEntry(2*a+1, 2*b+0);
      cdouble TMMab = -1.0 * TInt->GetEntry(2*a+1, 2*b+1) / ZVAC;

      cdouble TEEba = ZVAC * TInt->GetEntry(2*b+0, 2*a+0);
      cdouble TEMba = -1.0 * TInt->GetEntry(2*b+0, 2*a+1);
      cdouble TMEba =        TInt->GetEntry(2*b+1, 2*a+0);
      cdouble TMMba = -1.0 * TInt->GetEntry(2*b+1, 2*a+1) / ZVAC;

      cdouble SymTEE = 0.5*RYTOVPF*( TEEab + conj(TEEba) );
      cdouble SymTEM = 0.5*RYTOVPF*( TEMab + conj(TMEba) );
      cdouble SymTME = 0.5*RYTOVPF*( TMEab + conj(TEMba) );
      cdouble SymTMM = 0.5*RYTOVPF*( TMMab + conj(TMMba) );

      DR->SetEntry(OffsetS + 2*a+0, OffsetS + 2*b+0, SymTEE );
      DR->SetEntry(OffsetS + 2*b+0, OffsetS + 2*a+0, conj(SymTEE) );

      DR->SetEntry(OffsetS + 2*a+0, OffsetS + 2*b+1, SymTEM );
      DR->SetEntry(OffsetS + 2*b+1, OffsetS + 2*a+0, conj(SymTEM) );

      DR->SetEntry(OffsetS + 2*a+1, OffsetS + 2*b+0, SymTME );
      DR->SetEntry(OffsetS + 2*b+0, OffsetS + 2*a+1, conj(SymTME) );

      DR->SetEntry(OffsetS + 2*a+1, OffsetS + 2*b+1, SymTMM );
      DR->SetEntry(OffsetS + 2*b+1, OffsetS + 2*a+1, conj(SymTMM) );
    };

  /***************************************************************/
  /* set DR = W * DR * W' ****************************************/
  /* by computing DR = M \ (M \ DR)'       ***********************/
  /***************************************************************/
  M->LUSolve(DR);
  DR->Adjoint();
  M->LUSolve(DR);

  Log("...done with DR matrix");
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
/* evaluate trace formulas for the contribution of sources     */
/* inside SourceSurface to the fluxes of all spatially-        */
/* integrated quantities on all surfaces.                      */
/* return 0 if the calculation was skipped for one reason or   */
/* another, or 1 if the calculation was done.                  */
/***************************************************************/
int GetSIFlux(SNEQData *SNEQD, int SourceSurface, cdouble Omega,
               int PFTMethod, HMatrix *PFTMatrix)
{
  RWGGeometry *G      = SNEQD->G;
  int NumSurfaces     = G->NumSurfaces;
  bool OmitSelfTerms  = SNEQD->OmitSelfTerms;
  HMatrix *DRMatrix   = SNEQD->DRMatrix;
  int DestOnly        = SNEQD->DestOnly;

  if (     PFTMethod==SCUFF_PFT_DSI
       && !DoDSIAtThisFrequency(SNEQD, Omega)
     ) return 0;

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
        if ( (DestOnly!=-1 && DestSurface!=DestOnly) )
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

  return 1;

} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WriteFlux(SNEQData *SNEQD, cdouble Omega, double *kBloch)
{
  /***************************************************************/
  /* extract fields from SNEQData structure **********************/
  /***************************************************************/
  RWGGeometry *G      = SNEQD->G;
  HMatrix *M          = SNEQD->M;
  HMatrix **TExt      = SNEQD->TExt;
  HMatrix **TInt      = SNEQD->TInt; 
  HMatrix **U         = SNEQD->U;
  int NS              = SNEQD->G->NumSurfaces;
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

     if ( !(G->Surfaces[ns]->IsPEC) )
      { G->RegionMPs[ G->Surfaces[ns]->RegionIndices[1] ]->UnZero();
        G->AssembleBEMMatrixBlock(ns, ns, Omega, kBloch, TInt[ns]);
        G->RegionMPs[ G->Surfaces[ns]->RegionIndices[1] ]->Zero();
      };

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
     /*- transform the geometry -------------------------------------*/
     /*--------------------------------------------------------------*/
     char *Tag=SNEQD->GTCs[nt]->Tag;
     G->Transform(SNEQD->GTCs[nt]);
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
     /*- stamp all blocks into the BEM matrix and LU-factorize      -*/
     /*--------------------------------------------------------------*/
     for(int nb=0, ns=0; ns<NS; ns++)
      { 
        int RowOffset=G->BFIndexOffset[ns];
        M->InsertBlock(TExt[ns], RowOffset, RowOffset);
        if( !(G->Surfaces[ns]->IsPEC) )
         M->AddBlock(TInt[ns], RowOffset, RowOffset);

        for(int nsp=ns+1; nsp<NS; nsp++, nb++)
         { int ColOffset=G->BFIndexOffset[nsp];
           M->InsertBlock(U[nb], RowOffset, ColOffset);
           M->InsertBlockTranspose(U[nb], ColOffset, RowOffset);
         };
      };
     UndoSCUFFMatrixTransformation(M);
     Log("LU factorizing...");
     M->LUFactorize();

     /*--------------------------------------------------------------*/
     /*- compute the requested quantities for all objects           -*/
     /*- note: nss = 'num surface, source'                          -*/
     /*-       nsd = 'num surface, destination'                     -*/
     /*--------------------------------------------------------------*/
     int SourceOnly            = SNEQD->SourceOnly;
     int NumPFTMethods         = SNEQD->NumPFTMethods;
     int *PFTMethods           = SNEQD->PFTMethods;
     HMatrix *PFTMatrix        = SNEQD->PFTMatrix;
     HMatrix *PFTByRegion      = SNEQD->PFTByRegion;
     HMatrix **RegionRegionPFT = SNEQD->RegionRegionPFT;
     int NR = G->NumRegions;
     if (RegionRegionPFT)
      for(int npm=0; npm<NumPFTMethods; npm++)
       RegionRegionPFT[npm]->Zero();
     for(int nss=0; nss<NS; nss++)
      {
        if ( SourceOnly!=-1 && nss!=SourceOnly )
         continue;
        // PEC bodies do not act as thermal sources
        if (G->Surfaces[nss]->IsPEC) 
         continue;

        // compute the "dressed Rytov" matrix for this source
        ComputeDRMatrix(SNEQD, nss);

        // compute spatially-integrated flux quantities for
        // all destination objects using all requested 
        // calculation methods
        for(int npm=0; npm<NumPFTMethods; npm++)
         { 
           int Status=GetSIFlux(SNEQD, nss, Omega,
                                PFTMethods[npm], PFTMatrix);
           if (Status==0)
            continue;

           FILE *f=vfopen(SNEQD->SIFluxFileNames[npm],"a");
           for(int nsd=0; nsd<NS; nsd++)
            { 
              fprintf(f,"%s %e ",Tag,real(Omega));
              if (kBloch) fprintVec(f,kBloch,G->LDim);
              fprintf(f,"%i%i ",nss+1,nsd+1);
              for(int nq=0; nq<NUMPFT; nq++)
               fprintf(f,"%+.8e ",PFTMatrix->GetEntryD(nsd,nq));
              fprintf(f,"\n");
            };
           fclose(f);

           if (PFTByRegion)
            { GetPFTByRegion(G, PFTMatrix, PFTByRegion);
              int nsr1 = G->Surfaces[nss]->RegionIndices[0]; // source region 1
              int nsr2 = G->Surfaces[nss]->RegionIndices[1]; // source region 2
              for(int ndr=0; ndr<NR; ndr++) // ndr = destination region
               for(int nq=0; nq<NUMPFT; nq++)
                { double PFT = PFTByRegion->GetEntryD(ndr, nq);
                  if (nsr1!=0)
                   RegionRegionPFT[npm]->AddEntry( (nsr1+1)*(NR+1) + ndr+1, nq, PFT);
                  if (nsr2!=0)
                   RegionRegionPFT[npm]->AddEntry( (nsr2+1)*(NR+1) + ndr+1, nq, PFT);
                };
            };
         };

        // compute spatially-resolved flux quantities for
        // all evaluation points
        if (SNEQD->SRXMatrix)
          { 
            HMatrix *SRXMatrix = SNEQD->SRXMatrix;
            HMatrix *SRFMatrix = SNEQD->SRFMatrix;
            HMatrix *DRMatrix  = SNEQD->DRMatrix;
            GetSRFluxTrace(G, SRXMatrix, Omega, DRMatrix, SRFMatrix);

            FILE *f=vfopen("%s.SRFlux","a",FileBase);
            for(int nx=0; nx<SRXMatrix->NR; nx++)
             {
               double X[3], SRFlux[NUMSRFLUX];
               SRXMatrix->GetEntriesD(nx,":",X);
               SRFMatrix->GetEntriesD(nx,":",SRFlux);

               fprintf(f,"%s %e ",Tag,real(Omega));
               if (kBloch) 
                fprintVec(f, kBloch, G->LDim);
               fprintf(f,"%e %e %e %i ",X[0],X[1],X[2],nss);
               for(int nfc=0; nfc<NUMSRFLUX; nfc++)
                fprintf(f,"%e ",SRFMatrix->GetEntryD(nx,nfc));
               fprintf(f,"\n");
             };
            fclose(f);
          };

      }; // for(int nss=0; nss<NS; nss++)

     if (PFTByRegion)
      for(int npm=0; npm<NumPFTMethods; npm++)
       { 
         for(int nsr=0; nsr<NR; nsr++)
          for(int ndr=0; ndr<NR; ndr++)
           for(int nq=0; nq<NUMPFT; nq++)
            { double PFT=RegionRegionPFT[npm]->GetEntryD( (nsr+1)*(NR+1) + ndr+1, nq);
              RegionRegionPFT[npm]->AddEntry( 0*(NR+1) + ndr+1, nq, PFT );
              RegionRegionPFT[npm]->AddEntry( 0*(NR+1) +     0, nq, PFT );
            };

         FILE *f=vfopen("%s.byRegion","a",SNEQD->SIFluxFileNames[npm]);
         for(int nsr=0; nsr<=NR; nsr++) // ndr = number of destination region
          for(int ndr=0; ndr<=NR; ndr++) // ndr = number of destination region
           { fprintf(f,"%s %e ",Tag,real(Omega));
             if (kBloch) fprintVec(f,kBloch,G->LDim);
             fprintf(f,"%i%i ",nsr,ndr);
             for(int nq=0; nq<NUMPFT; nq++)
              fprintf(f,"%+.8e ",RegionRegionPFT[npm]->GetEntryD( nsr*(NR+1) + ndr, nq));
             fprintf(f,"\n");
           };
          fclose(f);
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
