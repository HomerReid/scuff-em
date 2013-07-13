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
 * CreateSNEQData.cc -- a utility function to initialize a 
 *                   -- scuff-neq data structure for a given
 *                   -- run of the code
 *
 * homer reid      -- 2/2012
 *
 */

#include "scuff-neq.h"

#define MAXSTR 1000

/***************************************************************/
/***************************************************************/
/***************************************************************/
SNEQData *CreateSNEQData(char *GeoFile, char *TransFile, int QuantityFlags, 
                         int PlotFlux, char *pFileBase, 
                         bool SymGSource, bool SymGDest)
{

  SNEQData *SNEQD=(SNEQData *)mallocEC(sizeof(*SNEQD));

  SNEQD->WriteCache=0;
  SNEQD->PlotFlux=PlotFlux;
  SNEQD->SymGSource=SymGSource;
  SNEQD->SymGDest=SymGDest;

  /*--------------------------------------------------------------*/
  /*-- try to create the RWGGeometry -----------------------------*/
  /*--------------------------------------------------------------*/
  RWGGeometry *G=new RWGGeometry(GeoFile);
  G->SetLogLevel(SCUFF_VERBOSELOGGING);
  SNEQD->G=G;
  
  if (pFileBase)
   SNEQD->FileBase = strdup(pFileBase);
  else
   SNEQD->FileBase = strdup(GetFileBase(G->GeoFileName));

  /*--------------------------------------------------------------*/
  /*- this code does not make sense if any of the objects are PEC */
  /*--------------------------------------------------------------*/
  int ns, nsp;
  for(ns=0; ns<G->NumSurfaces; ns++)
   if ( G->Surfaces[ns]->IsPEC ) 
    ErrExit("%s: object %s: PEC objects are not allowed in scuff-neq", G->GeoFileName,G->Surfaces[ns]->Label);

  /*--------------------------------------------------------------*/
  /*- read the transformation file if one was specified and check */
  /*- that it plays well with the specified geometry file.        */
  /*- note if TransFile==0 then this code snippet still works; in */
  /*- this case the list of GTComplices is initialized to contain */
  /*- a single empty GTComplex and the check automatically passes.*/
  /*--------------------------------------------------------------*/
  SNEQD->GTCList=ReadTransFile(TransFile, &(SNEQD->NumTransformations));
  char *ErrMsg=G->CheckGTCList(SNEQD->GTCList, SNEQD->NumTransformations);
  if (ErrMsg)
   ErrExit("file %s: %s",TransFile,ErrMsg);

  /*--------------------------------------------------------------*/
  /*- figure out which quantities were specified                 -*/
  /*--------------------------------------------------------------*/
  SNEQD->QuantityFlags=QuantityFlags;
  SNEQD->NQ=0;
  if ( QuantityFlags & QFLAG_POWER  ) SNEQD->NQ++;
  if ( QuantityFlags & QFLAG_XFORCE ) SNEQD->NQ++;
  if ( QuantityFlags & QFLAG_YFORCE ) SNEQD->NQ++;
  if ( QuantityFlags & QFLAG_ZFORCE ) SNEQD->NQ++;
  if ( QuantityFlags & QFLAG_XTORQUE ) SNEQD->NQ++;
  if ( QuantityFlags & QFLAG_YTORQUE ) SNEQD->NQ++;
  if ( QuantityFlags & QFLAG_ZTORQUE ) SNEQD->NQ++;

  SNEQD->NSNQ = G->NumSurfaces * SNEQD->NQ; 
  SNEQD->NTNSNQ = SNEQD->NumTransformations * SNEQD->NSNQ;

  /*--------------------------------------------------------------*/
  /*- allocate arrays of matrix subblocks that allow us to reuse -*/
  /*- chunks of the BEM matrices for multiple geometrical        -*/
  /*- transformations.                                           -*/
  /*--------------------------------------------------------------*/
  int nb, NS=G->NumSurfaces, NBF, NBFp;
  SNEQD->T = (HMatrix **)mallocEC(NS*sizeof(HMatrix *));
  SNEQD->TSelf = (HMatrix **)mallocEC(NS*sizeof(HMatrix *));
  SNEQD->U = (HMatrix **)mallocEC( ((NS*(NS-1))/2)*sizeof(HMatrix *));
  Log("Before T, U blocks: mem=%3.1f GB",GetMemoryUsage()/1.0e9);
  for(nb=ns=0; ns<G->NumSurfaces; ns++)
   { 
     NBF=G->Surfaces[ns]->NumBFs;

     if (G->Mate[ns]==-1)
      SNEQD->T[ns] = new HMatrix(NBF, NBF, LHM_COMPLEX, LHM_SYMMETRIC);
     else
      SNEQD->T[ns] = SNEQD->T[ G->Mate[ns] ];

     for(nsp=ns+1; nsp<G->NumSurfaces; nsp++, nb++)
      { NBFp=G->Surfaces[nsp]->NumBFs;
        SNEQD->U[nb] = new HMatrix(NBF, NBFp, LHM_COMPLEX);
      };

     if ( SNEQD->SymGSource || SNEQD->SymGDest )
      { if (G->Mate[ns]==-1)
         SNEQD->TSelf[ns] = new HMatrix(NBF, NBF, LHM_COMPLEX);//, LHM_SYMMETRIC);
        else
         SNEQD->TSelf[ns] = SNEQD->TSelf[ G->Mate[ns] ];
      };
   };
  Log("After T, U blocks: mem=%3.1f GB",GetMemoryUsage()/1.0e9);

  /*--------------------------------------------------------------*/
  /*- allocate BEM matrix ----------------------------------------*/
  /*--------------------------------------------------------------*/
  SNEQD->W = new HMatrix(G->TotalBFs, G->TotalBFs, LHM_COMPLEX );
  Log("After W: mem=%3.1f GB",GetMemoryUsage()/1.0e9);

  /*--------------------------------------------------------------*/
  /*- Buffer[0..2] are data storage buffers with enough room to  -*/
  /*- hold MaxBFs^2 cdoubles, where MaxBFs is the maximum number -*/
  /*- of basis functions on any object, i.e. the max dimension   -*/
  /*- of any BEM matrix subblock.                                -*/
  /*--------------------------------------------------------------*/
  int MaxBFs=G->Surfaces[0]->NumBFs;
  for(ns=1; ns<G->NumSurfaces; ns++)
   if (G->Surfaces[ns]->NumBFs > MaxBFs) 
    MaxBFs = G->Surfaces[ns]->NumBFs;
  
  SNEQD->Buffer[0] = malloc(MaxBFs*MaxBFs*sizeof(cdouble));
  SNEQD->Buffer[1] = malloc(MaxBFs*MaxBFs*sizeof(cdouble));
  SNEQD->Buffer[2] = malloc(MaxBFs*MaxBFs*sizeof(cdouble));

  /*--------------------------------------------------------------*/
  /*- allocate sparse matrices to store the various overlap      -*/
  /*- matrices. note that all overlap matrices have 10 nonzero   -*/
  /*- entries per row.                                           -*/
  /*-                                                            -*/
  /*- SArray[ns] is an array of SCUFF_NUM_OMATRICES pointers to  -*/
  /*- SMatrix structures for object #ns. SArray[ns][nom] is only -*/
  /*- non-NULL if we need the nomth type of overlap matrix. (Here-*/ 
  /*- nom=1..7 for power, xyz-force, xyz-torque, as defined in    */ 
  /*- libscuff.h).                                                */ 
  /*-                                                            -*/
  /*--------------------------------------------------------------*/
  bool *NeedMatrix=SNEQD->NeedMatrix;
  memset(NeedMatrix, 0, SCUFF_NUM_OMATRICES*sizeof(bool));
  NeedMatrix[SCUFF_OMATRIX_OVERLAP]  = 0;
  NeedMatrix[SCUFF_OMATRIX_POWER  ]  = true;
  NeedMatrix[SCUFF_OMATRIX_XFORCE ]  = QuantityFlags & QFLAG_XFORCE;
  NeedMatrix[SCUFF_OMATRIX_YFORCE ]  = QuantityFlags & QFLAG_YFORCE;
  NeedMatrix[SCUFF_OMATRIX_ZFORCE ]  = QuantityFlags & QFLAG_ZFORCE;
  NeedMatrix[SCUFF_OMATRIX_XTORQUE ] = QuantityFlags & QFLAG_XTORQUE;
  NeedMatrix[SCUFF_OMATRIX_YTORQUE ] = QuantityFlags & QFLAG_YTORQUE;
  NeedMatrix[SCUFF_OMATRIX_ZTORQUE ] = QuantityFlags & QFLAG_ZTORQUE;

  SNEQD->SArray=(SMatrix ***)mallocEC(NS*sizeof(SMatrix **));
  for(ns=0; ns<NS; ns++)
   { 
     SNEQD->SArray[ns]=(SMatrix **)mallocEC(SCUFF_NUM_OMATRICES*sizeof(SMatrix *));

     for(int nom=0; nom<SCUFF_NUM_OMATRICES; nom++)
      if (NeedMatrix[nom]) 
       SNEQD->SArray[ns][nom] = new SMatrix(G->Surfaces[ns]->NumBFs,G->Surfaces[ns]->NumBFs,LHM_COMPLEX);
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  time_t MyTime;
  struct tm *MyTm;
  char TimeString[30];
  MyTime=time(0);
  MyTm=localtime(&MyTime);
  strftime(TimeString,30,"%D::%T",MyTm);
  FILE *f=vfopen("%s.flux","a",SNEQD->FileBase);
  fprintf(f,"\n");
  fprintf(f,"# scuff-neq run on %s (%s)\n",GetHostName(),TimeString);
  fprintf(f,"# data file columns: \n");
  fprintf(f,"# 1 omega \n");
  fprintf(f,"# 2 transform tag\n");
  fprintf(f,"# 3 (sourceObject,destObject) \n");
  int nq=4;
  if (SNEQD->QuantityFlags & QFLAG_POWER) 
   fprintf(f,"# %i power flux spectral density\n",nq++);
  if (SNEQD->QuantityFlags & QFLAG_XFORCE) 
   fprintf(f,"# %i x-force flux spectral density\n",nq++);
  if (SNEQD->QuantityFlags & QFLAG_YFORCE) 
   fprintf(f,"# %i y-force flux spectral density\n",nq++);
  if (SNEQD->QuantityFlags & QFLAG_ZFORCE) 
   fprintf(f,"# %i z-force flux spectral density\n",nq++);
  if (SNEQD->QuantityFlags & QFLAG_XTORQUE) 
   fprintf(f,"# %i x-torque flux spectral density\n",nq++);
  if (SNEQD->QuantityFlags & QFLAG_YTORQUE) 
   fprintf(f,"# %i y-torque flux spectral density\n",nq++);
  if (SNEQD->QuantityFlags & QFLAG_ZTORQUE) 
   fprintf(f,"# %i z-torque flux spectral density\n",nq++);
  fclose(f);
  

  Log("After CreateSNEQData: mem=%3.1f GB",GetMemoryUsage()/1.0e9);
  return SNEQD;

}
