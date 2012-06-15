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
 * FrequencyIntegrand.cc -- evaluate spectral density of power/momentum 
 *                       -- transfer at a single frequency
 *
 * homer reid            -- 2/2012
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
/***************************************************************/
/***************************************************************/
#define OVERLAP_OVERLAP     0
#define OVERLAP_CROSS       1
#define OVERLAP_XBULLET     2
#define OVERLAP_XNABLANABLA 3
#define OVERLAP_XTIMESNABLA 4
#define OVERLAP_YBULLET     5
#define OVERLAP_YNABLANABLA 6
#define OVERLAP_YTIMESNABLA 7
#define OVERLAP_ZBULLET     8
#define OVERLAP_ZNABLANABLA 9
#define OVERLAP_ZTIMESNABLA 10
void AssembleOverlapMatrices(SNEQData *SNEQD, cdouble Omega)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGGeometry *G=SNEQD->G;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  SMatrix **OMatrices = SNEQD->OMatrices;
  SMatrix *PFMatrix;
  SMatrix *xMFMatrix;
  SMatrix *yMFMatrix;
  SMatrix *zMFMatrix;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble Eps, Mu;
  G->ExteriorMP->GetEpsMu(Omega,&Eps,&Mu);
  cdouble K2 = Eps*Mu*Omega*Omega;
  cdouble Z = ZVAC * sqrt(Mu/Eps);
 
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int no, neAlpha, neBeta;
  RWGObject *O;
  double Overlaps[11];
  cdouble PFEntry, MFEntry1, MFEntry2, MFEntry3;
  for(no=0; no<G->NumObjects; no++)
   { 
     O=G->Objects[no];

     PFMatrix  = OMatrices[ no*MAXQUANTITIES + QINDEX_POWER  ];
     xMFMatrix = OMatrices[ no*MAXQUANTITIES + QINDEX_XFORCE ];
     yMFMatrix = OMatrices[ no*MAXQUANTITIES + QINDEX_YFORCE ];
     zMFMatrix = OMatrices[ no*MAXQUANTITIES + QINDEX_ZFORCE ];
 
     for(neAlpha=0; neAlpha<O->NumEdges; neAlpha++)
      for(neBeta=0; neBeta<O->NumEdges; neBeta++)
       { 
         O->GetOverlaps(neAlpha, neBeta, Overlaps);
         if (Overlaps[0]==0.0) continue; 

         PFEntry=Overlaps[OVERLAP_CROSS];
         PFMatrix->SetEntry(2*neAlpha+0, 2*neBeta+1, PFEntry);
         PFMatrix->SetEntry(2*neAlpha+1, 2*neBeta+0, PFEntry);

         MFEntry1 = Overlaps[OVERLAP_XBULLET] - Overlaps[OVERLAP_XNABLANABLA]/K2;
         MFEntry2 = 2.0*Overlaps[OVERLAP_XTIMESNABLA] / (II*Omega);
         xMFMatrix->SetEntry(2*neAlpha+0, 2*neBeta+0, Z*MFEntry1);
         xMFMatrix->SetEntry(2*neAlpha+0, 2*neBeta+1, MFEntry2);
         xMFMatrix->SetEntry(2*neAlpha+1, 2*neBeta+0, -MFEntry2);
         xMFMatrix->SetEntry(2*neAlpha+1, 2*neBeta+1, MFEntry1/Z);

         MFEntry1 = Overlaps[OVERLAP_YBULLET] - Overlaps[OVERLAP_YNABLANABLA]/K2;
         MFEntry2 = 2.0*Overlaps[OVERLAP_YTIMESNABLA] / (II*Omega);
         yMFMatrix->SetEntry(2*neAlpha+0, 2*neBeta+0, Z*MFEntry1);
         yMFMatrix->SetEntry(2*neAlpha+0, 2*neBeta+1, MFEntry2);
         yMFMatrix->SetEntry(2*neAlpha+1, 2*neBeta+0, -MFEntry2);
         yMFMatrix->SetEntry(2*neAlpha+1, 2*neBeta+1, MFEntry1/Z);

         MFEntry1 = Overlaps[OVERLAP_ZBULLET] - Overlaps[OVERLAP_ZNABLANABLA]/K2;
         MFEntry2 = 2.0*Overlaps[OVERLAP_ZTIMESNABLA] / (II*Omega);
         zMFMatrix->SetEntry(2*neAlpha+0, 2*neBeta+0, Z*MFEntry1);
         zMFMatrix->SetEntry(2*neAlpha+0, 2*neBeta+1, MFEntry2);
         zMFMatrix->SetEntry(2*neAlpha+1, 2*neBeta+0, -MFEntry2);
         zMFMatrix->SetEntry(2*neAlpha+1, 2*neBeta+1, MFEntry1/Z);
       };
         
   };

}

/***************************************************************/
/* evaluate the four-matrix trace formula for the contribution */
/* of fluctuations within SourceObject to the flux of quantity */
/* QIndex into DestObject.                                     */
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
/*                                                             */
/***************************************************************/
double GetTrace(SNEQData *SNEQD, int QIndex, 
                int SourceObject, int DestObject, 
                FILE *ByOmegaFile)
{
  RWGGeometry *G      = SNEQD->G;
  HMatrix *W          = SNEQD->W;
  SMatrix **OMatrices = SNEQD->OMatrices;

  int Offset1       = G->BFIndexOffset[DestObject];
  SMatrix *OMatrix1 = OMatrices[ DestObject*MAXQUANTITIES + QIndex ];

  int Offset2       = G->BFIndexOffset[SourceObject];
  SMatrix *OMatrix2 = OMatrices[ SourceObject*MAXQUANTITIES + QINDEX_POWER ];

  int p, q, r, s; 
  int nnzq, nq, nnzs, ns;
  int qValues[11], sValues[11];
  double O1Entries[10], O2Entries[10];

  cdouble FMPTrace=0.0; //'four-matrix-product trace'

  for(p=0; p<OMatrix1->NR; p++)
   { 
     nnzq=OMatrix1->GetRow(p, qValues, O1Entries);

     for(r=0; r<OMatrix2->NR; r++)
      { 
        nnzs=OMatrix2->GetRow(r, sValues, O2Entries);

        for(nq=0, q=qValues[0]; nq<nnzq; q=qValues[++nq] )
         for(ns=0, s=sValues[0]; ns<nnzs; s=sValues[++ns] )
          FMPTrace +=  O1Entries[nq]
                      *W->GetEntry( Offset1+q, Offset2+r )
                      *O2Entries[ns]
                      *conj( W->GetEntry( Offset1+p, Offset2+s ));
      };
   };
  FMPTrace *= (-1.0/16.0);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (ByOmegaFile)
   fprintf(ByOmegaFile,"%e ",real(FMPTrace));

 return real(FMPTrace);

} 

/***************************************************************/
/* the computed quantities are ordered in the output vector    */
/* like this:                                                  */
/*                                                             */
/*  FI[ nt*NO2NQ + no*NONQ + nop*NQ + nq ]                     */
/*   = contribution of sources in object #nop to flux of       */
/*     quantity #nq into object #no                            */
/*                                                             */
/*  where    NQ = number of quantities (1--4)                  */
/*  where  NONQ = number of objects * NQ                       */
/*  where NO2NQ = (number of objects)^2* NQ                    */
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
  int QuantityFlags   = SNEQD->QuantityFlags;

  /***************************************************************/
  /* preinitialize an argument structure for the BEM matrix      */
  /* block assembly routine                                      */
  /***************************************************************/
  ABMBArgStruct MyABMBArgStruct, *Args=&MyABMBArgStruct;
  InitABMBArgs(Args);
  Args->G         = G;
  Args->Omega     = Omega;

  /***************************************************************/
  /* before entering the loop over transformations, we first     */
  /* assemble the (transformation-independent) T matrix blocks.  */
  /***************************************************************/
  int no, nop, nb, NO=G->NumObjects;
  for(no=0; no<NO; no++)
   { 
     Log(" Assembling self contributions to T(%i)...",no);

     Args->Oa = Args->Ob = G->Objects[no];
     Args->B = T[no];
     Args->Symmetric=1;
     AssembleBEMMatrixBlock(Args);
   };

  /***************************************************************/
  /* also before entering the loop over transformations, we      */
  /* pause to assemble the overlap matrices.                     */
  /***************************************************************/
  AssembleOverlapMatrices(SNEQD, Omega);

  /***************************************************************/
  /* now loop over transformations.                              */
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
     UndoSCUFFMatrixTransformation(W);
     W->LUFactorize();
     W->LUInvert();

     /*--------------------------------------------------------------*/
     /*- compute the requested quantities for all objects -----------*/
     /*--------------------------------------------------------------*/
     FILE *f=0;
     int nfi=0;
     for(no=0; no<NO; no++)
      for(nop=0; nop<NO; nop++)
       {
         if (SNEQD->ByOmegaFileNames)
          { f=fopen(SNEQD->ByOmegaFileNames[no*NO+nop],"a");
            fprintf(f,"%e %s ",real(Omega),Tag);
          };

         if ( QuantityFlags & QFLAG_POWER )
          FI[nfi++] = GetTrace(SNEQD, QINDEX_POWER,  no, nop, f);
         if ( QuantityFlags & QFLAG_XFORCE )
          FI[nfi++] = GetTrace(SNEQD, QINDEX_XFORCE, no, nop, f);
         if ( QuantityFlags & QFLAG_YFORCE )
          FI[nfi++] = GetTrace(SNEQD, QINDEX_YFORCE, no, nop, f);
         if ( QuantityFlags & QFLAG_ZFORCE )
          FI[nfi++] = GetTrace(SNEQD, QINDEX_ZFORCE, no, nop, f);

         if (f)
          { fprintf(f,"\n");
            fclose(f);
          };
      };

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
