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
void SpyPlot(HMatrix *M, char *Title, char *FileName)
{
  if (!FileName)
   FileName="/tmp/.spyplot.dat";
   
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
printf("george jetson\n");
 // pclose(f);
printf("doomatage\n");
  
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
  SMatrix ***SArray   = SNEQD->SArray;

  int Offset1       = G->BFIndexOffset[DestObject];
  SMatrix *OMatrix1 = SArray[SourceObject][ 1 + QIndex ];

  int Offset2       = G->BFIndexOffset[SourceObject];
  SMatrix *OMatrix2 = SArray[DestObject][ 1 + QINDEX_POWER];

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
  int *NeedMatrix     = SNEQD->NeedMatrix;
  SMatrix ***SArray   = SNEQD->SArray;

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
  for(no=0; no<NO; no++)
   G->Objects[no]->GetOverlapMatrices(NeedMatrix, SArray[no], Omega, G->ExteriorMP);

#if 0
//void *pCC=HMatrix::OpenMATLABContext("SMatrices.out");
char buffer1[100];
char buffer2[100];
for(no=0; no<NO; no++)
 for(int nom=0; nom<SCUFF_NUM_OMATRICES; nom++)
  { if (NeedMatrix[nom])
     { HMatrix *M=new HMatrix(SArray[no][nom]);
  //     M->ExportToMATLAB(pCC,"M_%i_%i",no,nom);
sprintf(buffer1,"M_%i_%i",no,nom);
sprintf(buffer2,"M_%i_%i.dat",no,nom);
SpyPlot(M,buffer1,buffer2);
printf("Suppatage foryaf %i %i \n",no,nom);
       delete M;
     };
  };
//HMatrix::CloseMATLABContext(pCC);
#endif
       

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
cdouble Eps1, Eps2, Mu;
G->Objects[0]->MP->GetEpsMu(Omega,&Eps1,&Mu);
G->Objects[1]->MP->GetEpsMu(Omega,&Eps2,&Mu);
fprintf(f,"%e %e %e %e ",real(Eps1),imag(Eps1),real(Eps2),imag(Eps2));
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
