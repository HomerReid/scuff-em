/* 
 * XQIntegrand.cc -- calculate casimir energy and/or force integrand at a
 *                -- single (xi,q) point 
 *
 * homer reid     -- 10/2008 -- 10/2010
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdarg.h>
#include <unistd.h>

#include <libhrutil.h>
#include "libTDRT.h"
#include "scuff-cas2D.h"

extern "C" {
/* Subroutine */ int zgetrf_(int *m, int *n, cdouble *a, 
	int *lda, int *ipiv, int *info);
}

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void *mypCC;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/


/***************************************************************/
/* compute \log \det \{ M^{-1} MInfinity \}                    */
/*                                                             */
/* ResultPtr is a pointer into an output buffer. We write      */
/* one or two values into this buffer and update ResultPtr     */ 
/* accordingly.                                                */ 
/*                                                             */
/* Thus, suppose on entry *ResultPtr points to a buffer Buf.   */
/* Then, on return,                                            */
/*   if                                                        */
/*      W->TETM==0                                             */
/*   then                                                      */
/*      Buf[0] = full (TE+TM) result                           */
/*   and                                                       */
/*      ResultPtr has been incremented by 1                    */
/*                                                             */
/*   whereas                                                   */
/*   if                                                        */
/*      W->TETM==1                                             */
/*   then                                                      */
/*      Buf[0] = TE contribution to result                     */
/*      Buf[1] = TM contribution to result                     */
/*   and                                                       */
/*      ResultPtr has been incremented by 2.                   */
/*                                                             */
/* For the TE/TM decomposition, we note that parallel-type     */
/* electric currents and z-type magnetic currents, i.e. the    */
/* 1st and 4th elements in each 4-element chunk of diagonal    */
/* entries, correspond to TE modes, while z-type electric      */
/* currents and parallel-type magnetic currents (2nd and 3rd   */
/* entries in each 4-entry block) correspond to TM modes.      */
/***************************************************************/
void GetLNDetMInvMInf(C2DWorkspace *W, double **ResultPtr)
{ 
  HMatrix *M; 
  double *DRMInf;
  int n, N;
  double LNDetTE, LNDetTM;
  double *p;

  M=W->M;
  DRMInf=W->DRMInf;
  N=W->N;

  LNDetTE=LNDetTM=0.0;
  if (W->G->AllPEC)
   { 
     for(n=0; n<N/2; n++)
      { LNDetTE +=  log( DRMInf[2*n+0] / abs(M->GetEntry(2*n+0,2*n+0)) );
        LNDetTM +=  log( DRMInf[2*n+1] / abs(M->GetEntry(2*n+1,2*n+1)) );
      };
   }
  else
   { 
     for(n=0; n<N/4; n++) 
      { LNDetTE +=  log( DRMInf[4*n+0] / abs(M->GetEntry(4*n+0,4*n+0)) )
                   +log( DRMInf[4*n+3] / abs(M->GetEntry(4*n+3,4*n+3)) );
        LNDetTM +=  log( DRMInf[4*n+1] / abs(M->GetEntry(4*n+1,4*n+1)) )
                   +log( DRMInf[4*n+2] / abs(M->GetEntry(4*n+2,4*n+2)) );
      };
   };

  LNDetTE *= -1.0/(2.0*M_PI*M_PI);
  LNDetTM *= -1.0/(2.0*M_PI*M_PI);

  p=*ResultPtr;
  if (W->TETM)
   { (*p++) = LNDetTE;
     (*p++) = LNDetTM;
   }
  else
    (*p++) = (LNDetTE + LNDetTM);
  *ResultPtr=p;

} 

/***************************************************************/
/* compute force integrand                                     */
/* see comments to previous routine for details of calling     */
/* convention                                                  */
/***************************************************************/
void GetTraceMInvdM(C2DWorkspace *W, char XY, double **ResultPtr)
{
  TDRTGeometry *G;
  HMatrix *M, *dM, **dU0b;
  int n, N1, niv, nop, RowOffset;
  double TraceTE,TraceTM;
  double *p;
   
  if (TDRTGeometry::LogLevel>=2)
   Log(" Computing %c force ...",XY);
  Tic();

  /***************************************************************/
  /* extract fields from workspace structure *********************/
  /***************************************************************/
  G=W->G;
  M=W->M;
  dM=W->dM;
  dU0b=( XY=='X' ? W->dU0bdX : W->dU0bdY ); 

  /***************************************************************/
  /* stamp dU blocks into dM matrix ******************************/
  /***************************************************************/
  RowOffset=G->Objects[0]->NumBFs;
  dM->Zero();
  for(nop=1; nop<G->NumObjects; nop++)
   { dM->InsertBlockTranspose(dU0b[nop],RowOffset,0);
     RowOffset += G->Objects[nop]->NumBFs;
   };

  /***************************************************************/
  /* stamp in dT0dY block if a ground plane is present and we    */
  /* are doing a y-force calculation                             */
  /***************************************************************/
  if ( W->GroundPlane && XY=='Y' )
   dM->InsertBlock(W->dT0dY, 0, 0);

  /***************************************************************/
  /* compute the first N1 (=NBF) columns of M\dMdZ               */
  /***************************************************************/
  M->LUSolve(dM);
  if (mypCC)
   dM->ExportToMATLAB(mypCC,"dMd%c%s",XY,W->CurrentTag);

  /***************************************************************/
  /* compute the trace of the upper-leftmost block of dMdZ.      */
  /***************************************************************/
  TraceTE=TraceTM=0.0;
  if (W->NumContributingIVs>0)
   { 
     if (W->G->AllPEC)
      { 
        for(niv=0; niv<W->NumContributingIVs; niv++)
         { n=W->ContributingIVIndices[niv];
           TraceTE +=  dM->GetEntryD(2*n, 2*n);
           TraceTM +=  dM->GetEntryD(2*n+1, 2*n+1);
         };
      }
     else
      { 
        for(niv=0; niv<W->NumContributingIVs; niv++)
         { n=W->ContributingIVIndices[niv];
           TraceTE +=  dM->GetEntryD(4*n, 4*n) + dM->GetEntryD(4*n+3, 4*n+3);
           TraceTM +=  dM->GetEntryD(4*n+1, 4*n+1) + dM->GetEntryD(4*n+2,4*n+2);
         };
      };
   }
  else
   { 
     N1=G->Objects[0]->NumBFs;
     if (W->G->AllPEC)
      { 
        for(n=0; n<N1/2; n++)
         { TraceTE +=  dM->GetEntryD(2*n, 2*n);
           TraceTM +=  dM->GetEntryD(2*n+1, 2*n+1);
         };
      }
     else
      { 
        for(n=0; n<N1/4; n++) 
         { TraceTE +=  dM->GetEntryD(4*n, 4*n) + dM->GetEntryD(4*n+3, 4*n+3);
           TraceTM +=  dM->GetEntryD(4*n+1, 4*n+1) + dM->GetEntryD(4*n+2,4*n+2);
         };
      };
   };
  TraceTE *= 2.0;  /* because the contribution of the first N1 entries */
  TraceTM *= 2.0;  /* gives us exactly 1/2 the full trace              */

  TraceTE *= -1.0/(2.0*M_PI*M_PI);
  TraceTM *= -1.0/(2.0*M_PI*M_PI);

  p=*ResultPtr;
  if (W->TETM)
   { (*p++) = TraceTE;
     (*p++) = TraceTM;
   }
  else
   (*p++)=TraceTE + TraceTM;
  *ResultPtr=p;
}

/***************************************************************/
/* assemble full BEM matrix from its constituent T and U       */
/* blocks, then attempt to compute its LU-factorization.       */
/* returns 0 on failure or nonzero on success.                 */
/***************************************************************/
int Factorize(C2DWorkspace *W)
{ 
  TDRTGeometry *G;
  HMatrix *M;
  int no, nop, NO, Offset, RowOffset, ColOffset;
  int info;

  if (TDRTGeometry::LogLevel>=2)
   Log(" Factorizing the M matrix...");
  Tic();

  /***************************************************************/ 
  /* unpack fields from workspace structure **********************/
  /***************************************************************/
  G=W->G;
  NO=G->NumObjects;
  M=W->M;

  /***************************************************************/
  /* stamp T and U blocks into R matrix  *************************/
  /***************************************************************/
  M->Zero();

  for(no=0, Offset=0; no<NO; Offset=G->BFIndexOffset[++no])
   M->InsertBlock(W->T[no], Offset, Offset);

  for(no=0, RowOffset=0; no<NO; RowOffset=G->BFIndexOffset[++no])
   for(nop=no+1, ColOffset=G->BFIndexOffset[nop]; nop<NO; ColOffset=G->BFIndexOffset[++nop])
    { M->InsertBlock(W->Uab[no][nop], RowOffset, ColOffset);
      M->InsertBlockTranspose(W->Uab[no][nop], ColOffset, RowOffset);
    };

  /***************************************************************/
  /* add T0 image blocks if groundplane is present.              */
  /* unfortunately i don't have an AddBlock() method in the      */
  /* HMatrix class, so this has to be written out explicitly.    */
  /***************************************************************/
  if ( W->GroundPlane )
   { int NBF, nbfa, nbfb;
     for (no=0; no<NO; no++)
      { Offset=G->BFIndexOffset[no];
        NBF=G->Objects[no]->NumBFs;
        for(nbfa=0; nbfa<NBF; nbfa++)
         for(nbfb=0; nbfb<NBF; nbfb++)
          M->AddEntry(Offset+nbfa, Offset+nbfb, W->TI[no]->GetEntry(nbfa,nbfb));
      };
   };

  /***************************************************************/
  /* attempt to compute LU-factorization  ************************/
  /***************************************************************/
  info=M->LUFactorize();

  if (info)
   { if (TDRTGeometry::LogLevel>=2) 
      LogC("FAILED with error code %i",info);
     return 0;
   };

  if (TDRTGeometry::LogLevel>=2)
   LogC("%.1f ms",1.0e3*Toc());
  return 1;

} 

/***************************************************************/
/* CacheRead: attempt to bypass an entire calculation of the   */
/* Xi-q integrand (or of the full q integral) by reading       */
/* results from a cache file.                                  */ 
/*                                                             */
/* If q==-1.0, we look in the .byXi file for a set of data     */
/* lines whose Xi values match the specified value of Xi.      */
/*                                                             */
/* Otherwise, we look in the .byXQ file for a set of data      */
/* lines whose (Xi,q) values match the specified values of Xi  */
/* and q.                                                      */
/*                                                             */
/* Returns 1 if successful (which means the values of the      */
/* energy/force integrand for ALL transformations at this      */
/* value of Xi,q were successfully read from the file) or 0    */
/* on failure).                                                */
/***************************************************************/
int CacheRead(C2DWorkspace *W, double Xi, double q, double *EF)
{ 
  FILE *f;
  double Q[3];
  char Line[1000], fTag[1000];
  double fXi, fq;
  int nt, ntnq, nRead, NumQuantitiesRead, FoundFirst;

  /*----------------------------------------------------------*/
  /* 0. try to open the cache file. --------------------------*/
  /*----------------------------------------------------------*/
  if ( q==-1.0 )
   f=fopen(W->ByXiFileName,"r");
  else 
   f=fopen(W->ByXQFileName,"r");

  if (!f)
   return 0;

  /*----------------------------------------------------------*/
  /* 1. skip down through the cache file until we find a line */
  /*    whose Xi, q, Tag values match.                        */
  /*----------------------------------------------------------*/
  FoundFirst=0;
  while( !FoundFirst && fgets(Line,1000,f) )
   { 
     if ( Line[0]=='#' || Line[0]=='\n' )
      continue;

     sscanf(Line,"%s %le %le",fTag,&fXi,&fq);

     if (    ( fabs(fXi-Xi) <= 1.0e-6*Xi)
          && ( q==-1.0 || fabs(fq-q) <= 1.0e-6*q)
          && ( !strcmp(fTag,W->Tags[0])) 
        )
      FoundFirst=1;
   };
  if ( !FoundFirst ) 
   { 
     if (TDRTGeometry::LogLevel>=2)
      { if (q==-1.0)
         Log(" did not find (Tag,Xi)=(%s,%e) in .byXi file.",W->Tags[0],Xi);
        else 
         Log(" did not find (Tag,Xi,q)=(%s,%e,%e) in .byXQ file.",W->Tags[0],Xi,q);
      };
     fclose(f); 
     return 0;
   };

  if (TDRTGeometry::LogLevel>=2)
   { if (q==-1.0)
      Log(" found (Tag,Xi)=(%s,%e) in .byXi file...",fTag,Xi);
     else
      Log(" found (Tag,Xi,q)=(%s,%e,%e) in .byXQ file...",fTag,Xi,q);
   };

  /*----------------------------------------------------------*/
  /* 2. verify that the line we just read from the cache file */
  /*    contains data for all the quantities we need          */
  /*----------------------------------------------------------*/
  if (q==-1.0)
   { nRead=sscanf(Line,"%s %le %le %le %le",fTag,&fXi,Q,Q+1,Q+2);
     NumQuantitiesRead=nRead-2;
   }
  else
   { nRead=sscanf(Line,"%s %le %le %le %le %le",fTag,&fXi,&fq,Q,Q+1,Q+2);
     NumQuantitiesRead=nRead-3;
   };
  
  if ( NumQuantitiesRead < W->NumQuantities )
   { if (TDRTGeometry::LogLevel>=2)
      Log(" ...but number of quantities is wrong (skipping)");
     fclose(f);
     return 0; 
   };
  memcpy(EF,Q,W->NumQuantities*sizeof(double));
  ntnq=W->NumQuantities;

  /*----------------------------------------------------------*/
  /* 3. ok, since that worked, now keep going ----------------*/
  /*----------------------------------------------------------*/
  for(nt=1; nt<W->NumTransforms; nt++)
   { 
     /* check for premature end of file */
     if ( !fgets(Line,1000,f) )
      { if (TDRTGeometry::LogLevel>=2)
         Log(" ...but data for some transforms were missing (skipping)");
        fclose(f);
        return 0; 
      };

     if (q==-1.0)
      { nRead=sscanf(Line,"%s %le %le %le %le",fTag,&fXi,Q,Q+1,Q+2);
        NumQuantitiesRead=nRead-2;
      }
     else
      { nRead=sscanf(Line,"%s %le %le %le %le %le",fTag,&fXi,&fq,Q,Q+1,Q+2);
        NumQuantitiesRead=nRead-3;
      };
  
     /* check for incorrect number of quantities */
     if ( NumQuantitiesRead < W->NumQuantities )
      { if (TDRTGeometry::LogLevel>=2)
         Log(" ...but number of quantities is wrong (skipping)");
        fclose(f);
        return 0; 
      };

     /* check for tag and/or Xi,q mismatch */
     if (     ( fabs(fXi-Xi)>1.0e-6*Xi)
           || ( q!=-1.0 && fabs(fq-q)>1.0e-6*q )
           || (strcmp(fTag,W->Tags[nt])) 
        )
      { if (TDRTGeometry::LogLevel>=2)
         Log(" ...but tag #%i did not match (%s != %s) (skipping)",nt,fTag,W->Tags[nt]);
        fclose(f);
        return 0; 
      };

     memcpy(EF+ntnq,Q,W->NumQuantities*sizeof(double));
     ntnq+=W->NumQuantities;

   };

  if (TDRTGeometry::LogLevel>=1)
   { if ( q==-1.0 )
      Log(" (Xi)=(%g): successfully read all data from .byXi file",Xi);
     else 
      Log(" (Xi,q)=(%g,%g): successfully read all data from .byXQ file",Xi,q);
   };

  fclose(f);
  return 1;

} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
void XQIntegrand(C2DWorkspace *W, double Xi, double q, double *EF)
{ 
  TDRTGeometry *G;
  int no, nop, NO;
  int nbf, nbfp, NBF;
  int nt, nq, ntnq;
  int N;
  int AllConverged;
  int info;
  int Offset, MateOffset;
  double EpsOut, MuOut;
  double *EFPtr;
  void *pCC;

  char *Tag, MovedString[100];

  FILE *f;

  G=W->G;
  N=G->TotalBFs;
  NO=G->NumObjects;

  if (W->WriteHDF5)
   pCC=HMatrix::OpenMATLABContext("%s_%g_%g",GetFileBase(G->GeoFileName),Xi,q);
  else
   pCC=0;


  /***************************************************************/
  /* attempt to bypass calculation by looking up results in cache*/
  /***************************************************************/
  if ( W->WriteHDF5==0 && CacheRead(W, Xi, q, EF) )
   return;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (TDRTGeometry::LogLevel>=1)
   Log("Computing Casimir quantities at (Xi,q)=(%e,%e)",Xi,q);

  /***************************************************************/
  /* get the permittivity and permeability of the external medium*/
  /* at this frequency                                           */
  /***************************************************************/
  cdouble zEps, zMu;
  G->MP->GetEpsMu( cdouble(0,Xi), &zEps, &zMu);
  EpsOut = real(zEps);
  MuOut = real(zMu);

  /***************************************************************/
  /* assemble T matrices                                         */
  /***************************************************************/
  for(no=0; no<G->NumObjects; no++)
   { 
     /* skip if this object is identical to a previous object */
     if ( (nop=G->Mate[no]) !=-1 )
      { if (TDRTGeometry::LogLevel>=2)
         Log("Object %i is identical to object %i (skipping)",no,nop);
        continue;
      };

     if (TDRTGeometry::LogLevel>=2)
      { Log("Assembling T%i at (Xi,Q)=(%g,%g)...",no+1,Xi,q);
        Tic();
      };

     AssembleT(W->G->Objects[no],Xi,q,EpsOut,MuOut,
               W->NumThreads, W->TSSSIDataTables[no], W->T[no]);

     if (TDRTGeometry::LogLevel>=2)
      LogC("%.1f ms",1.0e3*Toc());

     if (pCC)
      W->T[no]->ExportToMATLAB(pCC,"T%i",no);

   }; // for(no=0; no<G->NumObjects; no++)

  /***************************************************************/
  /* if an energy calculation was requested, precompute the      */
  /* diagonal of the cholesky factor of the MInfinity matrix.    */
  /***************************************************************/
  if (W->WhichQuantities & QUANTITY_ENERGY)
   {
     for(Offset=no=0; no<G->NumObjects; no++)
      { 
        NBF=G->Objects[no]->NumBFs;

        if ( (nop=G->Mate[no]) != -1 )
         { MateOffset=G->BFIndexOffset[nop];
           memcpy(W->DRMInf+Offset,W->DRMInf+MateOffset,NBF*sizeof(double));
           Offset+=NBF;
         }
        else
         { 
           /* KIND OF HACKY: we use the data buffer inside M as temporary  */
           /* storage for the content of T; this means that we have to     */
           /* call the lapack routines directly instead of using the nice  */
           /* wrappers provided by libhmat                                 */
           if (TDRTGeometry::LogLevel>=2)
            Log("LU-factorizing T%i at (Xi,Q)=(%g,%g)...",no+1,Xi,q);
           for(nbf=0; nbf<NBF; nbf++)
            for(nbfp=nbf; nbfp<NBF; nbfp++)
             { W->M->ZM[nbf + nbfp*NBF]=W->T[no]->GetEntry(nbf,nbfp);
               if (nbf!=nbfp) 
                W->M->ZM[nbfp + nbf*NBF]=conj(W->T[no]->GetEntry(nbf,nbfp));
             };
           Tic();
           zgetrf_(&NBF,&NBF,W->M->ZM,&NBF,W->ipiv,&info);
           if (TDRTGeometry::LogLevel>=2)
            LogC("%.1f ms",1.0e3*Toc());

           for(nbf=0; nbf<NBF; nbf++)
            W->DRMInf[Offset++] = abs(W->M->ZM[nbf+nbf*NBF]);
         };
      };

   };
      
  /***************************************************************/
  /* for each line in the TransFile, apply the specified         */
  /* transformation, then calculate all quantities requested.    */
  /***************************************************************/
  f=fopen(W->ByXQFileName,"a");
  EFPtr=EF;
  for(ntnq=nt=0; nt<W->NumTransforms; nt++)
   { 
     Tag=W->CurrentTag=W->Tags[nt];

     /******************************************************************/
     /* skip if all quantities are already converged at this transform */
     /******************************************************************/
     AllConverged=1;
     for(nq=0; AllConverged==1 && nq<W->NumQuantities; nq++)
      if ( !W->Converged[ ntnq + nq ] )
       AllConverged=0;

     if (AllConverged)
      { if (TDRTGeometry::LogLevel>=2)
         Log("All quantities already converged at Tag %s",Tag);
        fprintf(f,"%s %.15e %.15e ",Tag,Xi,q);
        for(nq=W->NumQuantities; nq>0; nq--)
         fprintf(f,"%.15e ",0.0);
        fprintf(f,"\n");
        memset(EFPtr,0,W->NumQuantities*sizeof(double));
        EFPtr+=W->NumQuantities;
        continue;
      };

     /***************************************************************/
     /* apply transform to geometry *********************************/
     /***************************************************************/
     G->Transform( W->TransLines[nt] );
     for(no=0; no<G->NumObjects; no++)
      MovedString[no] = G->ObjectMoved[no] ? '+' : '-';
     MovedString[no]=0;
     if (TDRTGeometry::LogLevel>=2)
      Log("Tag %s: %s",Tag,MovedString);

     /* apply transform to image objects */
     if (W->GroundPlane)
      { for(no=0; no<G->NumObjects; no++)
         W->ImageObjects[no]->Displace(  G->Objects[no]->Displacement[0],
                                        -G->Objects[no]->Displacement[1]);
      };

     /***************************************************************/
     /* assemble TI blocks at this transform if a ground plane is   */
     /* present                                                     */
     /***************************************************************/
     if (W->GroundPlane)  
      { for(no=0; no<G->NumObjects; no++)
         {  
           if ( nt>0 && G->ObjectMoved[no]==0 )
            continue;

           if (TDRTGeometry::LogLevel>=2)
            { Log(" Assembling TI(%i)...",no);
              Tic();
            };
   
           AssembleTI(G->Objects[no], W->ImageObjects[no], Xi, q,
                      EpsOut, MuOut, W->NumThreads, 0,
                      W->TI[no], no==0 ? W->dT0dY : 0);
   
           if (pCC)
            W->TI[no]->ExportToMATLAB(pCC,"TI%i_%s",no,Tag);

           if (TDRTGeometry::LogLevel>=2)
            LogC("%.1f ms",1.0e3*Toc());
         };
      };
   
     /***************************************************************/
     /* assemble U_{ab} blocks at this transformation               */
     /***************************************************************/
     for(no=0; no<G->NumObjects; no++)
      for(nop=no+1; nop<G->NumObjects; nop++)
       { 
         if ( nt>0 && G->ObjectMoved[no]==0 && G->ObjectMoved[nop]==0)
          continue;

         if (TDRTGeometry::LogLevel>=2)
          Log(" Assembling U(%i,%i) at Tag %s...",no,nop,Tag);
   
         Tic();

         if (no==0)
          AssembleU(G->Objects[no], G->Objects[nop], Xi, q, 
                    EpsOut, MuOut, W->NumThreads,
                    W->USSSIDataTables[nt][no+nop*NO],
                    W->Uab[no][nop], W->dU0bdX[nop], W->dU0bdY[nop]);
         else
          AssembleU(G->Objects[no], G->Objects[nop], Xi, q, 
                    EpsOut, MuOut, W->NumThreads,
                    W->USSSIDataTables[nt][no+nop*NO],
                    W->Uab[no][nop], 0, 0);

         if (TDRTGeometry::LogLevel>=2)
          LogC("%.1f ms",1.0e3*Toc());
   
         if (W->GroundPlane)
          { 
            if (TDRTGeometry::LogLevel>=2)
             Log("   Adding contributions of image objects...");
            Tic();
            if (no==0)
             AddImageContributionToU(G->Objects[no], W->ImageObjects[nop], Xi, q, 
                                     EpsOut, MuOut, W->NumThreads, 0,
                                     W->Uab[no][nop], W->dU0bdX[nop], W->dU0bdY[nop]);
            else
             AddImageContributionToU(G->Objects[no], W->ImageObjects[nop], Xi, q, 
                                     EpsOut, MuOut, W->NumThreads, 0,
                                     W->Uab[no][nop], 0, 0);
            if (TDRTGeometry::LogLevel>=2)
             LogC("%.1f ms",1.0e3*Toc());
          };

         if (pCC)
          { W->Uab[no][nop]->ExportToMATLAB(pCC,"U%i%i_%s",no,nop,Tag);
            if (no==0 && W->dU0bdX[nop])
             W->dU0bdX[nop]->ExportToMATLAB(pCC,"dU0%idX_%s",nop,Tag);
            if (no==0 && W->dU0bdY[nop])
             W->dU0bdY[nop]->ExportToMATLAB(pCC,"dU0%idY_%s",nop,Tag);
          };

       }; // for(no=0; ... for(nop=no+1...)

     /***************************************************************/
     /* attempt to cholesky-factorize M matrix and compute all      */
     /* requested quantities if successful                          */
     /***************************************************************/
     if ( Factorize(W) )
      { 
        if ( W->WhichQuantities & QUANTITY_ENERGY )
         GetLNDetMInvMInf(W,&EFPtr);
        if ( W->WhichQuantities & QUANTITY_XFORCE )
         GetTraceMInvdM(W,'X',&EFPtr);
        if ( W->WhichQuantities & QUANTITY_YFORCE )
         GetTraceMInvdM(W,'Y',&EFPtr);
      }
     else
      { memset(EFPtr,0,W->NumQuantities*sizeof(double));
        EFPtr+=W->NumQuantities;
      };

     /***************************************************************/
     /* write results to .byxq file                                 */
     /***************************************************************/
     fprintf(f,"%s %.15e %.15e ",Tag,Xi,q);
     for(nq=W->NumQuantities; nq>0; nq--)
      fprintf(f,"%.15e ",EF[ntnq++]);
     fprintf(f,"\n");
     fflush(f);
  
     /* undo transform */
     G->UnTransform();

     /* undo transform on image objects if present */
     if (W->GroundPlane)
      for(no=0; no<G->NumObjects; no++)
       W->ImageObjects[no]->UnDisplace();

   }; // for(ntnq=nt=0; nt<W->NumTransforms; nt++)
  fclose(f);

  if (pCC)
   HMatrix::CloseMATLABContext(pCC);

}
