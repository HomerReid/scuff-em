/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This fMAXQUANTITIESile is part of SCUFF-EM.
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
void LUInvert2(HMatrix *W, HMatrix *Scratch)
{ 
  Scratch->Zero();
  for(int n=0; n<W->NR; n++)
   Scratch->SetEntry(n,n,1.0);
  W->LUSolve(Scratch);
  W->Copy(Scratch);
}

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

/***************************************************************/
/***************************************************************/
/***************************************************************/
double GetTrace2(SNEQData *SNEQD, int QIndex, int SourceSurface, int DestSurface)
{
  RWGGeometry *G      = SNEQD->G;
  HMatrix *W          = SNEQD->W;
  HMatrix **SymG      = SNEQD->SymG;
  SMatrix ***SArray   = SNEQD->SArray;

  int OffsetS         = G->BFIndexOffset[DestSurface];
  SMatrix *OMatrix    = SArray[DestSurface][ 1 + QIndex ];

  int OffsetD         = G->BFIndexOffset[SourceSurface];
   
  int NS = G->Surfaces[SourceSurface]->NumBFs;
  int ND = G->Surfaces[DestSurface]->NumBFs;

  // set S1 = Sym(G_S)
  HMatrix *S1=new HMatrix(NS, NS, LHM_COMPLEX, LHM_NORMAL, SNEQD->Buffer[0]);
  S1->Zero();
  HMatrix *Gs = SymG[SourceSurface];
  for(int nr=0; nr<NS; nr++)
   for(int nc=0; nc<NS; nc++)
    S1->SetEntry(nr, nc, 0.5*(Gs->GetEntry(nr,nc) + conj(Gs->GetEntry(nc,nr)) ) ); 

  // set S2 = Sym(G_S) * W_{SD}
  HMatrix *WSD = new HMatrix(NS, ND, LHM_COMPLEX, LHM_NORMAL, SNEQD->Buffer[1]);
  W->ExtractBlock(OffsetS, OffsetD, WSD);
  HMatrix *S2 = new HMatrix(NS, NS, LHM_COMPLEX, LHM_NORMAL, SNEQD->Buffer[2]);
  S1->Multiply(WSD, S2);
  delete S1;

  // set S1 = W_{SD}^\dagger Sym(G_S) * W_{SD}
  S1=new HMatrix(ND, ND, LHM_COMPLEX, LHM_NORMAL, SNEQD->Buffer[0]);
  WSD->Adjoint();
  WSD->Multiply(S2, S1);

  // compute tr(O*S1)
  int *qValues;
  cdouble *O1Entries;
  cdouble FMPTrace=0.0; //'four-matrix-product trace'
  for(int p=0; p<ND; p++)
   { 
     int nnzq=OMatrix->GetRow(p, &qValues, (void **)&O1Entries);
     for(int nq=0, q=qValues[0]; nq<nnzq; q=qValues[++nq] )
      FMPTrace +=  O1Entries[nq] * S1->GetEntry(q, p);
   }; // for (p=0... 
  FMPTrace *= (-1.0/16.0);

 delete S1;
 delete S2;
 delete WSD;
 return real(FMPTrace);

} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetSelfForces(SNEQData *SNEQD, int ns, double *SelfForce)
{
  RWGGeometry *G      = SNEQD->G;
  HMatrix **SymG      = SNEQD->SymG;
  SMatrix ***SArray   = SNEQD->SArray;

  int N = G->Surfaces[ns]->NumBFs;

  /*--------------------------------------------------------------*/
  /*- this code snippet sets S1 = W_0^\dagger Sym(G_S) W_0        */
  /*--------------------------------------------------------------*/
  // set S1 = Sym(G_S)
  HMatrix *S1=new HMatrix(N, N, LHM_COMPLEX, LHM_NORMAL, SNEQD->Buffer[0]);
  S1->Zero();
  HMatrix *Gs = SymG[ns];
  for(int nr=0; nr<N; nr++)
   for(int nc=0; nc<N; nc++)
    S1->SetEntry(nr, nc, 0.5*(Gs->GetEntry(nr,nc) + conj(Gs->GetEntry(nc,nr)) ) ); 

  // set S2 = Sym(G_S) * W_0
  HMatrix *W0 = new HMatrix(N, N, LHM_COMPLEX, LHM_NORMAL, SNEQD->Buffer[1]);
  W0->Copy(SNEQD->T[ns]);
  UndoSCUFFMatrixTransformation(W0);
  W0->LUFactorize();
  W0->LUInvert();
  HMatrix *S2 = new HMatrix(N, N, LHM_COMPLEX, LHM_NORMAL, SNEQD->Buffer[2]);
  S1->Multiply(W0, S2);

  // set S1 = W_0^\dagger Sym(G_S) * W_0
  W0->Adjoint();
  W0->Multiply(S2, S1);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  for(int nf=0; nf<3; nf++)
   { if (SNEQD->NeedMatrix[2+nf] == false) 
      SelfForce[nf]=0.0;
     else 
      { SMatrix *OMatrix    = SArray[ns][ 2 + nf ];

        // compute tr(O*S1)
        int *qValues;
        cdouble *O1Entries;
        cdouble FMPTrace=0.0; //'four-matrix-product trace'
        for(int p=0; p<N; p++)
         { 
          int nnzq=OMatrix->GetRow(p, &qValues, (void **)&O1Entries);
          for(int nq=0, q=qValues[0]; nq<nnzq; q=qValues[++nq] )
           FMPTrace +=  O1Entries[nq] * S1->GetEntry(q, p);
         }; // for (p=0... 
        SelfForce[nf] = real(FMPTrace) * (-1.0/16.0);
      };
   };

 delete S1;
 delete S2; 
 delete W0; 

} 

/***************************************************************/
/* the computed quantities are ordered in the output vector    */
/* like this:                                                  */
/*                                                             */
/*  FI[ nt*NS2NQ + ns*NSNQ + nsp*NQ + nq ]                     */
/*   = contribution of sources inside surface #nsp to flux of  */
/*     quantity #nq into surface #ns, all at transformation #nt*/
/*                                                             */
/*  where    NQ = number of quantities (1--4)                  */
/*  where  NSNQ = number of surface * NQ                       */
/*  where NS2NQ = (number of surfaces)^2* NQ                   */
/***************************************************************/
void GetFlux(SNEQData *SNEQD, cdouble Omega, double *FI)
{
  Log("Computing neq quantities at omega=%s...",z2s(Omega));

  /***************************************************************/
  /* extract fields from SNEQData structure ************************/
  /***************************************************************/
  RWGGeometry *G      = SNEQD->G;
  HMatrix *W          = SNEQD->W;
  HMatrix **T         = SNEQD->T;
  HMatrix **SymG      = SNEQD->SymG;
  HMatrix **U         = SNEQD->U;
  int QuantityFlags   = SNEQD->QuantityFlags;
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
  double SelfForce[50][3];
  for(int ns=0; ns<NS; ns++)
   if ( G->Mate[ns]!=-1 )
    memcpy(SelfForce[ns], SelfForce[G->Mate[ns]], 3*sizeof(double));
   else
    { GetSelfForces(SNEQD, ns, SelfForce[ns]);
      Log("SelfForce(%e,%i)=(%e,%e,%e)",real(Omega),ns,
           SelfForce[ns][0],SelfForce[ns][1],SelfForce[ns][2]);
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
     for(int nss=0; nss<NS; nss++)
      for(int nsd=0; nsd<NS; nsd++)
       {
         fprintf(f,"%e %s %i%i ",real(Omega),Tag,nss+1,nsd+1);

         int nq=0;
         if ( QuantityFlags & QFLAG_POWER )
          { int i = GetIndex(SNEQD, nt, nss, nsd, nq++);
            FI[i] = GetTrace2(SNEQD, QINDEX_POWER, nss, nsd);
            fprintf(f,"%.8e ",FI[i]);
          };
         if ( QuantityFlags & QFLAG_XFORCE )
          { int i = GetIndex(SNEQD, nt, nss, nsd, nq++);
            FI[i] = GetTrace2(SNEQD, QINDEX_XFORCE, nss, nsd);
            if (nss==nsd) FI[i] -= SelfForce[nss][0];
            fprintf(f,"%.8e ",FI[i]);
          }
         if ( QuantityFlags & QFLAG_YFORCE )
          { int i = GetIndex(SNEQD, nt, nss, nsd, nq++);
            FI[i] = GetTrace2(SNEQD, QINDEX_YFORCE, nss, nsd);
            if (nss==nsd) FI[i] -= SelfForce[nss][1];
            fprintf(f,"%.8e ",FI[i]);
          }
         if ( QuantityFlags & QFLAG_ZFORCE )
          { int i = GetIndex(SNEQD, nt, nss, nsd, nq++);
            FI[i] = GetTrace2(SNEQD, QINDEX_ZFORCE, nss, nsd);
            if (nss==nsd) FI[i] -= SelfForce[nss][2];
            fprintf(f,"%.8e ",FI[i]);
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
