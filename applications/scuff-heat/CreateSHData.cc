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

  SHD->PlotFlux=PlotFlux;

  SHD->nThread=nThread;
  if (SHD->nThread==0)
   SHD->nThread=GetNumThreads();

  SHD->WriteCache=0;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (ByOmegaFile)
   SHD->ByOmegaFile = ByOmegaFile;
  else if (PlotFlux)
   SHD->ByOmegaFile = 0;
  else
   { SHD->ByOmegaFile = vstrdup("%s.byOmega",GetFileBase(GeoFile));
     char MyFileName[MAXSTR];
     FILE *f=CreateUniqueFile(SHD->ByOmegaFile, 1, MyFileName);
     fclose(f);
     SHD->ByOmegaFile=strdupEC(MyFileName);
   };

  /*--------------------------------------------------------------*/
  /*- allocate arrays of matrix subblocks that allow us to reuse -*/
  /*- chunks of the BEM matrices for multiple geometrical        -*/
  /*- transformations                                            -*/
  /*--------------------------------------------------------------*/
  int ns, nsp, nb, NS=G->NumSurfaces, NBF, NBFp;

  // SHD->TSelf[ns]   = contribution of medium inside surface #ns to (ns,ns) block of matrix 
  // SHD->TMedium[ns] = contribution of external medium to (ns,ns) block of matrix 
  // note that the T blocks of the BEM matrix as computed by libscuff are symmetric,
  // but here (in contrast to the case in scuff-cas3D), we CANNOT use packed-storage
  // symmetric matrices to store them, because after assembling them via libscuff
  // we flip the signs of the magnetic columns, which kills the symmetry.
  SHD->TSelf= (HMatrix **)mallocEC(NS*sizeof(HMatrix *));
  SHD->TMedium= (HMatrix **)mallocEC(NS*sizeof(HMatrix *));
  for(ns=0; ns<G->NumSurfaces; ns++)
   { NBF=G->Surfaces[ns]->NumBFs;
     SHD->TSelf[ns] = new HMatrix(NBF, NBF, LHM_COMPLEX);
     SHD->TMedium[ns] = new HMatrix(NBF, NBF, LHM_COMPLEX);
   };

  // SHD->UMedium[0]    = 0,1    block
  // SHD->UMedium[1]    = 0,2    block
  //             ...    = ...
  // SHD->UMedium[NS-1] = 0,NS-1 block
  // SHD->UMedium[NS]   = 1,2    block
  // etc.                                         
  SHD->UMedium = (HMatrix **)mallocEC( ( NS*(NS-1)/2)*sizeof(HMatrix *));
  for(nb=0, ns=0; ns<NS; ns++)
   for(nsp=ns+1; nsp<NS; nsp++, nb++)
    { NBF=G->Surfaces[ns]->NumBFs;
      NBFp=G->Surfaces[nsp]->NumBFs;
      SHD->UMedium[nb] = new HMatrix(NBF, NBFp, LHM_COMPLEX);
    };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int N, N1, N2;
  N = SHD->G->TotalBFs;
  N1 = SHD->N1 = SHD->G->Surfaces[0]->NumBFs;
  if (NS==1)
   N2 = N1;
  else
   N2 = SHD->N2 = N - N1;

  SHD->SymG1      = new HMatrix(N1, N1, LHM_COMPLEX );
  SHD->SymG2      = new HMatrix(N2, N2, LHM_COMPLEX );
  SHD->W          = new HMatrix(N,  N,  LHM_COMPLEX );
  SHD->W21        = new HMatrix(N2, N1, LHM_COMPLEX );
  SHD->W21SymG1   = new HMatrix(N2, N1, LHM_COMPLEX );
  SHD->W21DSymG2  = new HMatrix(N1, N2, LHM_COMPLEX );
  SHD->Scratch    = new HMatrix(N,  N1, LHM_COMPLEX );

  SHD->DV         = new HVector(N2, LHM_REAL);

  return SHD;

}
