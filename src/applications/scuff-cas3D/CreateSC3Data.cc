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

#define MAXSTR 1000

/***************************************************************/
/***************************************************************/
/***************************************************************/
SC3Data *CreateSC3Data(RWGGeometry *G, char *TransFile,
                       char *ByXiFile, char *ByXikBlochFile)
{
  SC3Data *SC3D=(SC3Data *)mallocEC(sizeof(*SC3D));
  SC3D->G = G;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  SC3D->WhichQuantities=WhichQuantities;
  SC3D->NumQuantities=NumQuantities;
  SC3D->NumTorqueAxes=0;
  if (NumQuantities & QUANTITY_TORQUE1)
   SC3D->NumTorqueAxes++;
  if (NumQuantities & QUANTITY_TORQUE2)
   SC3D->NumTorqueAxes++;
  if (NumQuantities & QUANTITY_TORQUE3)
   SC3D->NumTorqueAxes++;
  if (SC3D->NumTorqueAxes)
   memcpy(SC3D->TorqueAxes, TorqueArgs, 3*NumTorqueAxes*sizeof(cdouble));

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

  /*--------------------------------------------------------------*/
  /*- FIXME ------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if ( TransFile && G->NumLatticeBasisVectors>0 )
   ErrExit("--TransFiles are not yet supported for PBC geometries");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (ByXiFile)
   SC3D->ByXiFile = ByOmegaFile;
  else 
   { SC3D->ByXiFile = vstrdup("%s.byXi",GetFileBase(GeoFile));
     char MyFileName[MAXSTR];
     FILE *f=CreateUniqueFile(SC3D->ByXiFile, 1, MyFileName);
     fclose(f);
     SC3D->ByXiFile=strdupEC(MyFileName);
   };

  if (ByXikBlochFile)
   SC3D->ByXiFile = ByXikBlochFile;
  else 
   ByXikBlochFile=0;

  /*--------------------------------------------------------------*/
  /*- allocate arrays of matrix subblocks that allow us to reuse -*/
  /*- chunks of the BEM matrices for multiple geometrical        -*/
  /*- transformations                                            -*/
  /*--------------------------------------------------------------*/
  int no, nop, nb, NO=G->NumObjects, NBF, NBFp;

  // SC3D->TSelf[no]   = contribution of object #no to (no,no) block of matrix 
  // SC3D->TMedium[no] = contribution of external medium to (no,no) block of matrix 
  SC3D->TSelf= (HMatrix **)mallocEC(NO*sizeof(HMatrix *));
  SC3D->TMedium= (HMatrix **)mallocEC(NO*sizeof(HMatrix *));
  for(no=0; no<G->NumObjects; no++)
   { NBF=G->Objects[no]->NumBFs;
     SC3D->TSelf[no] = new HMatrix(NBF, NBF, LHM_COMPLEX, LHM_SYMMETRIC);
     SC3D->TMedium[no] = new HMatrix(NBF, NBF, LHM_COMPLEX, LHM_SYMMETRIC);
   };

  // SC3D->UMedium[0]    = 0,1    block
  // SC3D->UMedium[1]    = 0,2    block
  //             ...    = ...
  // SC3D->UMedium[NO-1] = 0,NO-1 block
  // SC3D->UMedium[NO]   = 1,2    block
  // etc.                                         
  SC3D->UMedium = (HMatrix **)mallocEC( ( NO*(NO-1)/2)*sizeof(HMatrix *));
  for(nb=0, no=0; no<NO; no++)
   for(nop=no+1; nop<NO; nop++, nb++)
    { NBF=G->Objects[no]->NumBFs;
      NBFp=G->Objects[nop]->NumBFs;
      SC3D->UMedium[nb] = new HMatrix(NBF, NBFp, LHM_COMPLEX);
    };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int N = SC3D->G->TotalBFs;
  int N1 = SC3D->N1 = SC3D->G->Objects[0]->NumBFs;
  int N2 = SC3D->N2 = N - N1;
  SC3D->SymG1      = new HMatrix(N1, N1, LHM_COMPLEX );
  SC3D->SymG2      = new HMatrix(N2, N2, LHM_COMPLEX );
  SC3D->W          = new HMatrix(N,  N,  LHM_COMPLEX );
  SC3D->W21        = new HMatrix(N2, N1, LHM_COMPLEX );
  SC3D->W21SymG1   = new HMatrix(N2, N1, LHM_COMPLEX );
  SC3D->W21DSymG2  = new HMatrix(N1, N2, LHM_COMPLEX );
  SC3D->Scratch    = new HMatrix(N,  N1, LHM_COMPLEX );

  return SC3D;

}
