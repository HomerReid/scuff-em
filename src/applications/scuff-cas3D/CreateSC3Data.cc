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
 * CreateSC3Data.cc -- a utility function to initialize a 
 *                 -- scuff-cas3D data structure for a given
 *                 -- run of the code
 *
 * homer reid      -- 2/2012
 *
 */

#include "scuff-cas3D.h"
#include <libhrutil.h>

#define MAXSTR 1000

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
SC3Data *CreateSC3Data(RWGGeometry *G, char *TransFile,
                       int WhichQuantities, int NumQuantities,
                       int NumTorqueAxes, double TorqueAxes[9])
{
  SC3Data *SC3D=(SC3Data *)mallocEC(sizeof(*SC3D));
  SC3D->G = G;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  SC3D->WhichQuantities=WhichQuantities;
  SC3D->NumQuantities=NumQuantities;
  SC3D->NumTorqueAxes=NumTorqueAxes;
  if (SC3D->NumTorqueAxes)
   memcpy(SC3D->TorqueAxes, TorqueAxes, 3*NumTorqueAxes*sizeof(cdouble));

  /*--------------------------------------------------------------*/
  /*- read the transformation file if one was specified and check */
  /*- that it plays well with the specified geometry file.        */
  /*- note if TransFile==0 then this code snippet still works; in */
  /*- this case the list of GTComplices is initialized to contain */
  /*- a single empty GTComplex and the check automatically passes.*/
  /*--------------------------------------------------------------*/
  SC3D->GTCList=ReadTransFile(TransFile, &(SC3D->NumTransformations));
  char *ErrMsg=G->CheckGTCList(SC3D->GTCList, SC3D->NumTransformations);
  if (ErrMsg)
   ErrExit("file %s: %s",TransFile,ErrMsg);

  SC3D->NTNQ = SC3D->NumTransformations * SC3D->NumQuantities;

  SC3D->Converged = (int *)mallocEC( (SC3D->NTNQ) * sizeof(int) );

  /*--------------------------------------------------------------*/
  /*- allocate arrays of matrix subblocks that allow us to reuse  */
  /*- chunks of the BEM matrices for multiple geometrical         */
  /*- transformations.                                            */
  /*-                                                             */
  /*- TBlocks[ns]       = (ns,ns) (diagonal) block                */
  /*- UBlocks[nb]       = nbth above-diagonal block               */
  /*- dUBlocks[3*nb+Mu] = Mu derivative of nbth above-diagonal blk*/
  /*--------------------------------------------------------------*/
  int ns, nsp, nb, NBF, NBFp;
  int NS=G->NumSurfaces; 
  SC3D->TBlocks  = (HMatrix **)mallocEC(NS*sizeof(HMatrix *));
  for(ns=0; ns<G->NumSurfaces; ns++)
   { NBF=G->Surfaces[ns]->NumBFs;
     SC3D->TBlocks[ns] = new HMatrix(NBF, NBF);
   };

  // SC3D->UBlocks[0]    = 0,1    block
  // SC3D->UBlocks[1]    = 0,2    block
  //             ...    = ...
  // SC3D->UBlocks[NS-1] = 0,NS   block
  // SC3D->UBlocks[NS]   = 1,2    block
  // etc.                                         
  int NumBlocks   = NS*(NS-1)/2; // number of above-diagonal blocks 
  SC3D->UBlocks   = (HMatrix **)mallocEC(   NumBlocks * sizeof(HMatrix *));
  SC3D->dUBlocks  = (HMatrix **)mallocEC( 3*NumBlocks * sizeof(HMatrix *));

  for(nb=0, ns=0; ns<NS; ns++)
   for(nsp=ns+1; nsp<NS; nsp++, nb++)
    { NBF=G->Surfaces[ns]->NumBFs;
      NBFp=G->Surfaces[nsp]->NumBFs;
      SC3D->UBlocks[nb] = new HMatrix(NBF, NBFp);
      if (WhichQuantities & QUANTITY_XFORCE)
       SC3D->dUBlocks[3*nb+0] = new HMatrix(NBF, NBFp);
      if (WhichQuantities & QUANTITY_YFORCE)
       SC3D->dUBlocks[3*nb+1] = new HMatrix(NBF, NBFp);
      if (WhichQuantities & QUANTITY_YFORCE)
       SC3D->dUBlocks[3*nb+2] = new HMatrix(NBF, NBFp);
    };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int N  = SC3D->N  = SC3D->G->TotalBFs;
  int N1 = SC3D->N1 = SC3D->G->Surfaces[0]->NumBFs;
  SC3D->M          = new HMatrix(N,  N);
  SC3D->dM         = new HMatrix(N,  N1);

  if (WhichQuantities & QUANTITY_ENERGY)
   { SC3D->MInfLUDiagonal = new HVector(G->TotalBFs);
     SC3D->ipiv = (int *)mallocEC(N*sizeof(int));
   }
  else
   { SC3D->MInfLUDiagonal=0;
     SC3D->ipiv=0;
   };

  return SC3D;

}
