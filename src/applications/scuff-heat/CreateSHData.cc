/*
 * CreateSHData.cc -- a utility function to initialize a 
 *                 -- scuff-heat data structure for a given
 *                 -- run of the code
 *
 * homer reid      -- 2/2012
 *
 */

#include "scuff-heat.h"

#define MAXSTR 1000

/***************************************************************/
/***************************************************************/
/***************************************************************/
SHData *CreateSHData(char *GeoFile, char *TransFile, int PlotFlux, 
                     char *ByOmegaFile, int nThread)
{
  SHData *SHD=(SHData *)mallocEC(sizeof(*SHD));

  /*--------------------------------------------------------------*/
  /*-- try to create the RWGGeometry -----------------------------*/
  /*--------------------------------------------------------------*/
  RWGGeometry *G=new RWGGeometry(GeoFile);
  G->SetLogLevel(SCUFF_VERBOSELOGGING);
  SHD->G=G;

  /*--------------------------------------------------------------*/
  /*- read the transformation file if one was specified and check */
  /*- that it plays well with the specified geometry file.        */
  /*- note if TransFile==0 then this code snippet still works; in */
  /*- this case the list of GTComplices is initialized to contain */
  /*- a single empty GTComplex and the check automatically passes.*/
  /*--------------------------------------------------------------*/
  SHD->GTCList=ReadTransFile(TransFile, &(SHD->NumTransformations));
  char *ErrMsg=G->CheckGTCList(SHD->GTCList, SHD->NumTransformations);
  if (ErrMsg)
   ErrExit("file %s: %s",TransFile,ErrMsg);

  /*--------------------------------------------------------------*/
  /*- the DV field is only needed for generating flux plots.     -*/
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

  SHD->WriteCache=0;

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

  // SHD->TSelf[no]   = contribution of object #no to (no,no) block of matrix 
  // SHD->TMedium[no] = contribution of external medium to (no,no) block of matrix 
  SHD->TSelf= (HMatrix **)mallocEC(NO*sizeof(HMatrix *));
  SHD->TMedium= (HMatrix **)mallocEC(NO*sizeof(HMatrix *));
  for(no=0; no<G->NumObjects; no++)
   { NBF=G->Objects[no]->NumBFs;
     SHD->TSelf[no] = new HMatrix(NBF, NBF, LHM_COMPLEX, LHM_SYMMETRIC);
     SHD->TMedium[no] = new HMatrix(NBF, NBF, LHM_COMPLEX, LHM_SYMMETRIC);
   };

  // SHD->UMedium[0]    = 0,1    block
  // SHD->UMedium[1]    = 0,2    block
  //             ...    = ...
  // SHD->UMedium[NO-1] = 0,NO-1 block
  // SHD->UMedium[NO]   = 1,2    block
  // etc.                                         
  SHD->UMedium = (HMatrix **)mallocEC( ( NO*(NO-1)/2)*sizeof(HMatrix *));
  for(nb=0, no=0; no<NO; no++)
   for(nop=no+1; nop<NO; nop++, nb++)
    { NBF=G->Objects[no]->NumBFs;
      NBFp=G->Objects[nop]->NumBFs;
      SHD->UMedium[nb] = new HMatrix(NBF, NBFp, LHM_COMPLEX);
    };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int N = SHD->G->TotalBFs;
  int N1 = SHD->N1 = SHD->G->Objects[0]->NumBFs;
  int N2 = SHD->N2 = N - N1;
  SHD->SymG1      = new HMatrix(N1, N1, LHM_COMPLEX );
  SHD->SymG2      = new HMatrix(N2, N2, LHM_COMPLEX );
  SHD->W          = new HMatrix(N,  N,  LHM_COMPLEX );
  SHD->W21        = new HMatrix(N2, N1, LHM_COMPLEX );
  SHD->W21SymG1   = new HMatrix(N2, N1, LHM_COMPLEX );
  SHD->W21DSymG2  = new HMatrix(N1, N2, LHM_COMPLEX );
  SHD->Scratch    = new HMatrix(N,  N1, LHM_COMPLEX );

  return SHD;

}
