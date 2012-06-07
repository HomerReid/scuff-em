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
SNEQData *CreateSNEQData(char *GeoFile, char *TransFile,
                         int QuantityFlags, int PlotFlux)
{

  SNEQData *SNEQD=(SNEQData *)mallocEC(sizeof(*SNEQD));

  SNEQD->WriteCache=0;
  SNEQD->PlotFlux=PlotFlux;

  /*--------------------------------------------------------------*/
  /*-- try to create the RWGGeometry -----------------------------*/
  /*--------------------------------------------------------------*/
  RWGGeometry *G=new RWGGeometry(GeoFile);
  G->SetLogLevel(SCUFF_VERBOSELOGGING);
  SNEQD->G=G;

  /*--------------------------------------------------------------*/
  /*- this code does not make sense if any of the objects are PEC */
  /*--------------------------------------------------------------*/
  int no, nop;
  for(no=0; no<G->NumObjects; no++)
   if ( G->Objects[no]->MP->IsPEC() ) 
    ErrExit("%s: object %s: PEC objects are not allowed in scuff-neq", G->GeoFileName,G->Objects[no]->Label);

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

  SNEQD->NONQ = G->NumObjects * SNEQD->NQ; 
  SNEQD->NTNONQ = SNEQD->NumTransformations * SNEQD->NONQ;

  /*--------------------------------------------------------------*/
  /*- set the name of the .byOmega output file -------------------*/
  /*--------------------------------------------------------------*/
#if 0
  if (ByOmegaFile)
   SNEQD->ByOmegaFile = ByOmegaFile;
  else if (PlotFlux)
   SNEQD->ByOmegaFile = 0;
  else
   { SNEQD->ByOmegaFile = vstrdup("%s.byOmega",GetFileBase(GeoFile));
     char MyFileName[MAXSTR];
     FILE *f=CreateUniqueFile(SNEQD->ByOmegaFile, 1, MyFileName);
     fclose(f);
     SNEQD->ByOmegaFile=strdup(MyFileName);
   };
#endif

  /*--------------------------------------------------------------*/
  /*- allocate arrays of matrix subblocks that allow us to reuse -*/
  /*- chunks of the BEM matrices for multiple geometrical        -*/
  /*- transformations.                                           -*/
  /*--------------------------------------------------------------*/
  int nb, nq, NO=G->NumObjects, NBF, NBFp;
  SNEQD->T = (HMatrix **)mallocEC(NO*sizeof(HMatrix *));
  SNEQD->U = (HMatrix **)mallocEC( ((NO*(NO-1))/2)*sizeof(HMatrix *));
  for(nb=no=0; no<G->NumObjects; no++)
   { NBF=G->Objects[no]->NumBFs;
     SNEQD->T[no] = new HMatrix(NBF, NBF, LHM_COMPLEX, LHM_SYMMETRIC);
     for(nop=no+1; nop<G->NumObjects; nop++, nb++)
      { NBFp=G->Objects[nop]->NumBFs;
        SNEQD->U[nb] = new HMatrix(NBF, NBFp, LHM_COMPLEX);
      };
   };

  /*--------------------------------------------------------------*/
  /*- allocate BEM matrix ----------------------------------------*/
  /*--------------------------------------------------------------*/
  SNEQD->W = new HMatrix(G->TotalBFs, G->TotalBFs, LHM_COMPLEX );

  /*--------------------------------------------------------------*/
  /*- allocate sparse matrices to store the various overlap      -*/
  /*- matrices. note that all overlap matrices have 10 nonzero   -*/
  /*- entries per row.                                           -*/
  /*-                                                            -*/
  /*- Also note: we allocate space for all overlap matrices even -*/
  /*- though they may not all be required depending on which     -*/
  /*- quantities the user requested. (For example, if only       -*/
  /*- --power and --zforce were specified, then the x- and y-    -*/
  /*- momentum flux matrices are not needed.) We could save some -*/
  /*- small amount of memory by allocating space for only the    -*/
  /*- matrices we will actually need.                            -*/
  /*--------------------------------------------------------------*/
  SNEQD->OMatrices=(SMatrix **)mallocEC(MAXQUANTITIES*NO*sizeof(SMatrix *));
  for(no=0; no<NO; no++)
   for(nq=0; nq<MAXQUANTITIES; nq++)
    SNEQD->OMatrices[ no*MAXQUANTITIES + nq ]=new SMatrix(G->Objects[no]->NumBFs,10);

  /*--------------------------------------------------------------*/
  /*- create frequency-resolved output files for each object in  -*/
  /*- the geometry and write a file header to each file.         -*/
  /*--------------------------------------------------------------*/
  SNEQD->ByOmegaFileNames=0;
  int WriteByOmegaFiles=1;
  if (WriteByOmegaFiles)
   { 
     SNEQD->ByOmegaFileNames=(char **)mallocEC(NO*NO*sizeof(char *));
     FILE *f;
     for(no=0; no<NO; no++)
      for(nop=0; nop<NO; nop++)
       { 
         SNEQD->ByOmegaFileNames[no*NO+nop] = vstrdup("From%sTo%s.byOmega",
                                                       G->Objects[no]->Label,
                                                       G->Objects[nop]->Label);
   
         f=fopen(SNEQD->ByOmegaFileNames[no*NO+nop],"a");
         if (!f)
          ErrExit("could not create file %s",SNEQD->ByOmegaFileNames[no]);
         fprintf(f,"# data file columns: \n");
         fprintf(f,"# 1: angular frequency in units of 3e14 rad/sec \n");
         fprintf(f,"# 2: transformation tag \n");
        
         nq=3;
         if (QuantityFlags && QFLAG_POWER)
          fprintf(f,"# %i: spectral density of power flux \n",nq++);
         if (QuantityFlags && QFLAG_XFORCE)
          fprintf(f,"# %i: spectral density of x-momentum flux\n",nq++);
         if (QuantityFlags && QFLAG_YFORCE)
          fprintf(f,"# %i: spectral density of y-momentum flux\n",nq++);
         if (QuantityFlags && QFLAG_ZFORCE)
          fprintf(f,"# %i: spectral density of z-momentum flux\n",nq++);
   
         fclose(f);
       }; // for ( no = ...) for (nop = ...)

   }; //if (WriteByOmegaFiles)

  return SNEQD;

}
