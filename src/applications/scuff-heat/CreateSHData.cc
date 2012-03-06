/*
 * CreateSHData.cc -- a utility function to initialize a 
 *                 -- scuff-heat data structure for a given
 *                 -- run of the code
 *
 * homer reid      -- 2/2012
 *
 */

#include "scuff-heat.h"

/***************************************************************/
/***************************************************************/
/***************************************************************/
SHData *CreateSHData(char *GeoFile, char *TransFile, int PlotFlux, 
                     char *ByOmegaFile, int nThread)
{
  SHData *SHD=(SHData *)malloc(sizeof(*SHD));

  /*--------------------------------------------------------------*/
  /*-- try to create the RWGGeometry ------------------------------*/
  /*--------------------------------------------------------------*/
  RWGGeometry *G=new RWGGeometry(GeoFile);
  G->SetLogLevel(SCUFF_VERBOSELOGGING);
  SHD->G=G;

  /*--------------------------------------------------------------*/
  /*- read the transformation file if one was specified and check */
  /*- that it plays well with the specified geometry file.        */
  /*- note if TransFile==0 then this code snippet still works (in */
  /*- this case the list of GTComplices is initialized to contain */
  /*- a single empty GTComplex and the check automatically passes)*/
  /*--------------------------------------------------------------*/
  SHD->GTCList=ReadTransFile(TransFile, &(SHD->NumGTCs));
  char *ErrMsg=G->CheckGTCList(SHD->GTCList, SHD->NumGTCs);
  if (ErrMsg)
   ErrExit("file %s: %s",TransFile,ErrMsg);

  /*--------------------------------------------------------------*/
  /*- the DV field is only needed for generating flux plots      -*/
  /*--------------------------------------------------------------*/
  SHD->PlotFlux=PlotFlux;
  if (PlotFlux)
   SHD->DV=new HVector(G->TotalBFs); // diagonal vector (note real-valued)

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  SHD->nThread=nThread;
  if (SHD->nThread==0)
   SHD->nThread=GetNumThreads();

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (ByOmegaFile)
   SHD->ByOmegaFile = ByOmegaFile;
  else 
   { SHD->ByOmegaFile = vstrdup("%s.byOmega",GetFileBase(GeoFile));
     char MyFileName[MAXSTR];
     FILE *f=CreateUniqueFile(SHD->ByOmegaFile, 1, MyFileName);
     fclose(f);
     SHD->ByOmegaFile=strdup(MyFileName);
   };

  /*--------------------------------------------------------------*/
  /*- allocate arrays of matrix subblocks that allow us to reuse -*/
  /*- chunks of the BEM matrices for multiple geometrical        -*/
  /*- transformations                                            -*/
  /*--------------------------------------------------------------*/
  int no, nop, nb, NO=G->NumObjects, NBF, NBFp;

  // SHD->TBlocks[no] = (no,no) block of matrix 
  SHD->TBlocks = (HMatrix **)malloc(NO*sizeof(HMatrix *));
  for(no=0; no<G->NumObjects; no++)
   { NBF=G->Objects[no]->NumBFs;
     SHD->TBlocks[no] = new HMatrix(NBF, NBF, LHM_COMPLEX, LHM_SYMMETRIC);
   };

  // SHD->UBlocks[0]    = 0,1    block
  // SHD->UBlocks[1]    = 0,2    block
  //             ...    = ...
  // SHD->UBlocks[NO-1] = 0,NO-1 block
  // SHD->UBlocks[NO]   = 1,2    block
  // etc.                                         
  SHD->UBlocks = (HMatrix **)malloc( ( NO*(NO-1)/2)*sizeof(HMatrix *));
  for(nb=0, no=0; no<G->NumObjects; no++)
   for(nop=no+1; nop<G->NumObjects; nop++, nb++)
    { NBF=G->Objects[no]->NumBFs;
      NBFp=G->Objects[nop]->NumBFs;
      SHD->UBlocks[nb] = new HMatrix(NBF, NBFp, LHM_COMPLEX);
    };

  SHD->M0 = SHD->G->AllocateBEMMatrix();
  SHD->M1 = SHD->G->AllocateBEMMatrix();
  SHD->M2 = SHD->G->AllocateBEMMatrix();
}
