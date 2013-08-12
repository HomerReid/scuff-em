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
#include <time.h>
#include <sys/time.h>
#include <libhrutil.h>

#define MAXSTR 1000

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
SC3Data *CreateSC3Data(RWGGeometry *G, char *TransFile,
                       int WhichQuantities, int NumQuantities,
                       int NumTorqueAxes, double TorqueAxes[9],
                       bool NewEnergyMethod)
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
   { memcpy(SC3D->TorqueAxes, TorqueAxes, 3*NumTorqueAxes*sizeof(cdouble));
     for(int nta=0; nta<NumTorqueAxes; nta++)
      CreateGammaMatrix(TorqueAxes+3*nta, SC3D->GammaMatrix+9*nta);
   };

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
  /*- compute a basis for the reciprocal lattice -----------------*/
  /*--------------------------------------------------------------*/
  if (G->NumLatticeBasisVectors==1)
   ErrExit("1D lattice periodicity not yet supported");
  if (G->NumLatticeBasisVectors==2)
   { 
     double *L1 = G->LatticeBasisVectors[0];
     double *L2 = G->LatticeBasisVectors[1];

     double PreFac = 2.0*M_PI / fabs( L1[0]*L2[1] - L1[1]*L2[0]);

     SC3D->RLBasisVectors[0][0] = PreFac*L2[1];
     SC3D->RLBasisVectors[0][1] = -PreFac*L2[0];

     SC3D->RLBasisVectors[1][0] = -PreFac*L1[1];
     SC3D->RLBasisVectors[1][1] = PreFac*L1[0];
   };

  /*--------------------------------------------------------------*/
  /*- allocate arrays of matrix subblocks that allow us to reuse  */
  /*- chunks of the BEM matrices for multiple geometrical         */
  /*- transformations.                                            */
  /*-                                                             */
  /*- TBlocks[ns]       = (ns,ns) (diagonal) block                */
  /*- UBlocks[nb]       = nbth above-diagonal block               */
  /*- dUBlocks[6*nb+Mu] = Mu derivative of nbth above-diagonal blk*/
  /*-                     where Mu=0,1,2, for x,y,z-displacement  */
  /*-                           Mu=3,4,5, for axis 1,2,3 rotation */
  /*--------------------------------------------------------------*/
  int ns, nsp, nb, NBF, NBFp;
  int NS=G->NumSurfaces; 
  SC3D->TBlocks  = (HMatrix **)mallocEC(NS*sizeof(HMatrix *));
  for(ns=0; ns<G->NumSurfaces; ns++)
   { if ( (nsp=G->Mate[ns])!=-1 )
      SC3D->TBlocks[ns] = SC3D->TBlocks[nsp];
     else
      { NBF=G->Surfaces[ns]->NumBFs;
        SC3D->TBlocks[ns] = new HMatrix(NBF, NBF, LHM_REAL, LHM_SYMMETRIC);
      };
   };

  // SC3D->UBlocks[0]    = 0,1    block
  // SC3D->UBlocks[1]    = 0,2    block
  //             ...    = ...
  // SC3D->UBlocks[NS-1] = 0,NS   block
  // SC3D->UBlocks[NS]   = 1,2    block
  // etc.                                         
  int NumBlocks   = NS*(NS-1)/2; // number of above-diagonal blocks 
  SC3D->UBlocks   = (HMatrix **)mallocEC(   NumBlocks * sizeof(HMatrix *));
  SC3D->dUBlocks  = (HMatrix **)mallocEC( 6*NumBlocks * sizeof(HMatrix *));

  for(nb=0, ns=0; ns<NS; ns++)
   for(nsp=ns+1; nsp<NS; nsp++, nb++)
    { NBF=G->Surfaces[ns]->NumBFs;
      NBFp=G->Surfaces[nsp]->NumBFs;
      SC3D->UBlocks[nb] = new HMatrix(NBF, NBFp);
      if (WhichQuantities & QUANTITY_XFORCE)
       SC3D->dUBlocks[6*nb+0] = new HMatrix(NBF, NBFp);
      if (WhichQuantities & QUANTITY_YFORCE)
       SC3D->dUBlocks[6*nb+1] = new HMatrix(NBF, NBFp);
      if (WhichQuantities & QUANTITY_ZFORCE)
       SC3D->dUBlocks[6*nb+2] = new HMatrix(NBF, NBFp);
      if (WhichQuantities & QUANTITY_TORQUE1)
       SC3D->dUBlocks[6*nb+3] = new HMatrix(NBF, NBFp);
      if (WhichQuantities & QUANTITY_TORQUE2)
       SC3D->dUBlocks[6*nb+4] = new HMatrix(NBF, NBFp);
      if (WhichQuantities & QUANTITY_TORQUE3)
       SC3D->dUBlocks[6*nb+5] = new HMatrix(NBF, NBFp);
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

  /*--------------------------------------------------------------*/
  /*- 20130427 for the 'new energy method' we need to allocate   -*/
  /*- storage for an extra full-size BEM matrix.                 -*/
  /*--------------------------------------------------------------*/
  SC3D->NewEnergyMethod = NewEnergyMethod;
  if (NewEnergyMethod && (WhichQuantities & QUANTITY_ENERGY) )
   SC3D->MM1MInf = new HMatrix(N, N);
  else
   SC3D->MM1MInf = 0;

  return SC3D;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WriteFilePreamble(FILE *f, SC3Data *SC3D, int PreambleType)
{
  if (f==0) return;

  char DateStr[40];
  time_t MyTime = time(0);
  struct tm *MyTm=localtime(&MyTime);
  strftime(DateStr,30,"%D::%T",MyTm);
  fprintf(f,"# scuff-cas3D run on %s at %s",GetHostName(),DateStr);
  fprintf(f,"# data file columns: \n");

  int nc=1;

  fprintf(f,"#%i: transform tag\n",nc++);

  if ( PreambleType==PREAMBLE_BYXI || PreambleType==PREAMBLE_BYXIK) 
   fprintf(f,"#%i: imaginary angular frequency\n",nc++);

  if ( PreambleType==PREAMBLE_BYXIK )
   { fprintf(f,"#%i,%i: bloch wavevector k_x,k_y\n",nc,nc+1);
     nc+=2;
   };
  
  if (PreambleType == PREAMBLE_OUT) 
   {
     if ( SC3D->WhichQuantities & QUANTITY_ENERGY )
      { fprintf(f,"#%i: energy \n",nc++);
        fprintf(f,"#%i: energy error \n",nc++);
      };
     if ( SC3D->WhichQuantities & QUANTITY_XFORCE )
      { fprintf(f,"#%i: x-force \n",nc++);
        fprintf(f,"#%i: x-force error \n",nc++);
      };
     if ( SC3D->WhichQuantities & QUANTITY_YFORCE )
      { fprintf(f,"#%i: y-force \n",nc++);
        fprintf(f,"#%i: y-force error \n",nc++);
      };
     if ( SC3D->WhichQuantities & QUANTITY_ZFORCE )
      { fprintf(f,"#%i: z-force \n",nc++);
        fprintf(f,"#%i: z-force error \n",nc++);
      };
     if ( SC3D->WhichQuantities & QUANTITY_TORQUE1 )
      { fprintf(f,"#%i: 1-torque \n",nc++);
        fprintf(f,"#%i: 1-torque error \n",nc++);
      };
     if ( SC3D->WhichQuantities & QUANTITY_TORQUE2 )
      { fprintf(f,"#%i: 2-torque \n",nc++);
        fprintf(f,"#%i: 2-torque error \n",nc++);
      };
     if ( SC3D->WhichQuantities & QUANTITY_TORQUE3 )
      { fprintf(f,"#%i: 3-torque \n",nc++);
        fprintf(f,"#%i: 3-torque error \n",nc++);
      };
   }
  else
   {
     if ( SC3D->WhichQuantities & QUANTITY_ENERGY )
      fprintf(f,"#%i: energy integrand \n",nc++);
     if ( SC3D->WhichQuantities & QUANTITY_XFORCE )
      fprintf(f,"#%i: x-force integrand \n",nc++);
     if ( SC3D->WhichQuantities & QUANTITY_YFORCE )
      fprintf(f,"#%i: y-force integrand \n",nc++);
     if ( SC3D->WhichQuantities & QUANTITY_ZFORCE )
      fprintf(f,"#%i: z-force integrand \n",nc++);
     if ( SC3D->WhichQuantities & QUANTITY_TORQUE1 )
      fprintf(f,"#%i: 1-torque integrand \n",nc++);
     if ( SC3D->WhichQuantities & QUANTITY_TORQUE2 )
      fprintf(f,"#%i: 2-torque integrand \n",nc++);
     if ( SC3D->WhichQuantities & QUANTITY_TORQUE3 )
      fprintf(f,"#%i: 3-torque integrand \n",nc++);
   };

}
