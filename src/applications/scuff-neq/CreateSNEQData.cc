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
SNEQData *CreateSNEQData(char *GeoFile, char *TransFile, char *ByOmegaFile, 
                         int QuantityFlags, int PlotFlux)
{

  SNEQData *SNEQD=(SNEQData *)mallocEC(sizeof(*SNEQD));

  SNEQD->PlotFlux=PlotFlux;
  SNEQD->WriteCache=0;

  /*--------------------------------------------------------------*/
  /*-- try to create the RWGGeometry -----------------------------*/
  /*--------------------------------------------------------------*/
  RWGGeometry *G=new RWGGeometry(GeoFile);
  G->SetLogLevel(SCUFF_VERBOSELOGGING);
  SNEQD->G=G;

  /*--------------------------------------------------------------*/
  /*- this code does not make sense if any of the objects are PEC */
  /*--------------------------------------------------------------*/
  int no;
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
  SNEQD->NumQuantities=0;
  if ( QuantityFlags & QFLAG_POWER  ) SNEQD->NumQuantities++;
  if ( QuantityFlags & QFLAG_XFORCE ) SNEQD->NumQuantities++;
  if ( QuantityFlags & QFLAG_YFORCE ) SNEQD->NumQuantities++;
  if ( QuantityFlags & QFLAG_ZFORCE ) SNEQD->NumQuantities++;
  SNEQD->NTNQ = (SNEQD->NumQuantities * SNEQD->NumTransformations);

  /*--------------------------------------------------------------*/
  /*- set the name of the .byOmega output file -------------------*/
  /*--------------------------------------------------------------*/
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

  /*--------------------------------------------------------------*/
  /*- allocate arrays of matrix subblocks that allow us to reuse -*/
  /*- chunks of the BEM matrices for multiple geometrical        -*/
  /*- transformations.                                           -*/
  /*--------------------------------------------------------------*/
  int no, nop, nb, NO=G->NumObjects, NBF, NBFp;
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

}
