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
 *                  -- scuff-cas3D data structure for a given
 *                  -- run of the code
 *
 * homer reid       -- 2/2012
 *
 */

#include "scuff-cas3D.h"
#include <time.h>
#include <sys/time.h>
#include <libhrutil.h>
#include <libTriInt.h>

#define MAXSTR 1000

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
SC3Data *CreateSC3Data(RWGGeometry *G, char *TransFile,
                       int WhichQuantities, int NumQuantities,
                       int NumTorqueAxes, double TorqueAxes[9],
                       bool NewEnergyMethod, char *FileBase)
{
  SC3Data *SC3D=(SC3Data *)mallocEC(sizeof(*SC3D));
  SC3D->G = G;
  int LDim = G->LBasis ? G->LBasis->NC : 0;
  bool PBC = (LDim>0);

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
  SC3D->GTCs=ReadTransFile(TransFile);
  SC3D->NumTransformations = SC3D->GTCs.size();
  char *ErrMsg=G->CheckGTCList(SC3D->GTCs);
  if (ErrMsg)
   ErrExit("file %s: %s",TransFile,ErrMsg);

  SC3D->NTNQ = SC3D->NumTransformations * SC3D->NumQuantities;

  SC3D->XiConverged = (bool *)mallocEC( (SC3D->NTNQ) * sizeof(bool) );
  if (LDim==0)
   SC3D->BZConverged = 0;
  else
   SC3D->BZConverged = (bool *)mallocEC( (SC3D->NTNQ) * sizeof(bool) );

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
  int NS=G->NumSurfaces;
  SC3D->TBlocks  = (HMatrix **)mallocEC(NS*sizeof(HMatrix *));
  SC3D->TAccelerators = PBC ? (void **)mallocEC(NS*sizeof(void *)) : 0;
  bool NeedDMDZ = (WhichQuantities & QUANTITY_ZFORCE);
  for(int ns=0; ns<G->NumSurfaces; ns++)
   { int nsp = G->Mate[ns];
     if ( nsp!=-1 )
      SC3D->TBlocks[ns] = SC3D->TBlocks[nsp];
     else
      { int NBF=G->Surfaces[ns]->NumBFs;
        if (PBC)
         { SC3D->TBlocks[ns] = new HMatrix(NBF, NBF, LHM_COMPLEX);
           SC3D->TAccelerators[ns] = G->CreateABMBAccelerator(ns, ns, true, NeedDMDZ);
         }
        else
         SC3D->TBlocks[ns] = new HMatrix(NBF, NBF, LHM_REAL, LHM_SYMMETRIC);
      };
   };

  // SC3D->UBlocks[0]    = 0,1    block
  // SC3D->UBlocks[1]    = 0,2    block
  //             ...    = ...
  // SC3D->UBlocks[NS-1] = 0,NS   block
  // SC3D->UBlocks[NS]   = 1,2    block
  // etc.                                         
  int NumBlocks       = NS*(NS-1)/2; // number of above-diagonal blocks
  SC3D->UBlocks       = (HMatrix **)mallocEC(   NumBlocks * sizeof(HMatrix *));
  SC3D->dUBlocks      = (HMatrix **)mallocEC( 6*NumBlocks * sizeof(HMatrix *));
  int RealComplex     = PBC ? LHM_COMPLEX : LHM_REAL;
  for(int ns=0, nb=0; ns<NS; ns++)
   for(int nsp=ns+1; nsp<NS; nsp++, nb++)
    { int NBF=G->Surfaces[ns]->NumBFs;
      int NBFp=G->Surfaces[nsp]->NumBFs;
      SC3D->UBlocks[nb] = new HMatrix(NBF, NBFp, RealComplex);
      if (WhichQuantities & QUANTITY_XFORCE)
       SC3D->dUBlocks[6*nb+0] = new HMatrix(NBF, NBFp, RealComplex);
      if (WhichQuantities & QUANTITY_YFORCE)
       SC3D->dUBlocks[6*nb+1] = new HMatrix(NBF, NBFp, RealComplex);
      if (WhichQuantities & QUANTITY_ZFORCE)
       SC3D->dUBlocks[6*nb+2] = new HMatrix(NBF, NBFp, RealComplex);
      if (WhichQuantities & QUANTITY_TORQUE1)
       SC3D->dUBlocks[6*nb+3] = new HMatrix(NBF, NBFp, RealComplex);
      if (WhichQuantities & QUANTITY_TORQUE2)
       SC3D->dUBlocks[6*nb+4] = new HMatrix(NBF, NBFp, RealComplex);
      if (WhichQuantities & QUANTITY_TORQUE3)
       SC3D->dUBlocks[6*nb+5] = new HMatrix(NBF, NBFp, RealComplex);
    };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int N  = SC3D->N  = SC3D->G->TotalBFs;
  int N1 = SC3D->N1 = SC3D->G->Surfaces[0]->NumBFs;
  SC3D->M           = new HMatrix(N,  N,  RealComplex);
  SC3D->dM          = new HMatrix(N,  N1, RealComplex);
  SC3D->NewEnergyMethod  = NewEnergyMethod;

  if (WhichQuantities & QUANTITY_ENERGY)
   { SC3D->MInfLUDiagonal = new HVector(G->TotalBFs);
     SC3D->ipiv = (int *)mallocEC(N*sizeof(int));
     if (NewEnergyMethod)
      SC3D->MM1MInf = new HMatrix(N, N, RealComplex);
   }
  else
   { SC3D->MInfLUDiagonal=0;
     SC3D->ipiv=0;
     SC3D->MM1MInf = 0;
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (!FileBase)
   FileBase = strdup(GetFileBase(G->GeoFileName));

  SC3D->FileBase = strdup(FileBase);

  SC3D->OutFileName=vstrdup("%s.out",FileBase);

  SC3D->ByXiFileName=vstrdup("%s.byXi",FileBase);
  WriteFilePreamble(SC3D, PREAMBLE_BYXI);

  if (LDim>0)
   { SC3D->ByXiKFileName=vstrdup("%s.byXikBloch",FileBase);
     WriteFilePreamble(SC3D, PREAMBLE_BYXIK);
   }

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  SC3D->UAccelerators = 0;
  if (PBC)
   { int NT = SC3D->NumTransformations;
     SC3D->UAccelerators = (void ***)mallocEC(NT*sizeof(void *));
     for(int nt=0; nt<NT; nt++)
      { SC3D->UAccelerators[nt] = (void **)mallocEC(NumBlocks*sizeof(void *));
        for(int ns=0, nb=0; ns<NS; ns++)
         for(int nsp=ns+1; nsp<NS; nsp++, nb++)
          SC3D->UAccelerators[nt][nb] 
           = G->CreateABMBAccelerator(ns, nsp, true, NeedDMDZ && ns==0);
      };
   };

  return SC3D;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WriteFilePreamble(SC3Data *SC3D, int PreambleType)
{
  char *FileName=0;
  switch(PreambleType)
   { case PREAMBLE_OUT:   FileName=SC3D->OutFileName;   break;
     case PREAMBLE_BYXI:  FileName=SC3D->ByXiFileName;  break;
     case PREAMBLE_BYXIK: FileName=SC3D->ByXiKFileName; break;
     default: ErrExit("%s:%i: internal error",__FILE__,__LINE__);
   };

  FILE *f=fopen(FileName,"a");
  if (f==0) 
   ErrExit("could not open file %s",FileName);

  /*--------------------------------------------------------------*/
  /*- write some common basic header information -----------------*/
  /*--------------------------------------------------------------*/
  char DateStr[40];
  time_t MyTime = time(0);
  struct tm *MyTm=localtime(&MyTime);
  strftime(DateStr,30,"%D::%T",MyTm);
  fprintf(f,"# scuff-cas3D " VERSIONSTRING " run on %s at %s\n",GetHostName(),DateStr);
  fprintf(f,"# data file columns: \n");
  fprintf(f,"#1: transform tag\n");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  const char *ErrorString;
  const char *IntegrandString;
  int nc=2;
  int LDim=SC3D->G->LBasis ? SC3D->G->LBasis->NC : 0;

  if (PreambleType == PREAMBLE_OUT)
   { 
     IntegrandString="";
     ErrorString="error due to numerical Xi integration";
   }
  else if (PreambleType == PREAMBLE_BYXI )
   { 
     fprintf(f,"#%i: imaginary angular frequency\n",nc++);

     IntegrandString="Xi integrand";
     ErrorString = (LDim==0) ? 0 :
                   "error due to numerical Brillouin-zone integration";
   }
  else  // PREAMBLE_BYXIK
   { 
     fprintf(f,"#%i: imaginary angular frequency\n",nc++);
     fprintf(f,"#%i: bloch wavevector kx \n",nc++);
     if (LDim>=2)
      fprintf(f,"#%i: bloch wavevector ky\n",nc++);

     IntegrandString="Brillouin-zone integrand";
     ErrorString=0;
   };
  
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if ( SC3D->WhichQuantities & QUANTITY_ENERGY )
   { fprintf(f,"#%i: energy %s\n",nc++,IntegrandString);
     if (ErrorString) fprintf(f,"#%i: energy %s \n",nc++,ErrorString);
   };

  if ( SC3D->WhichQuantities & QUANTITY_XFORCE )
   { fprintf(f,"#%i: x-force %s\n",nc++,IntegrandString);
     if (ErrorString) fprintf(f,"#%i: x-force %s \n",nc++,ErrorString);
   };

  if ( SC3D->WhichQuantities & QUANTITY_YFORCE )
   { fprintf(f,"#%i: y-force %s\n",nc++,IntegrandString);
     if (ErrorString) fprintf(f,"#%i: y-force %s \n",nc++,ErrorString);
   };

  if ( SC3D->WhichQuantities & QUANTITY_ZFORCE )
   { fprintf(f,"#%i: z-force %s\n",nc++,IntegrandString);
     if (ErrorString) fprintf(f,"#%i: z-force %s \n",nc++,ErrorString);
   };

  if ( SC3D->WhichQuantities & QUANTITY_TORQUE1 )
   { fprintf(f,"#%i: x-torque %s\n",nc++,IntegrandString);
     if (ErrorString) fprintf(f,"#%i: x-torque %s \n",nc++,ErrorString);
   };

  if ( SC3D->WhichQuantities & QUANTITY_TORQUE2 )
   { fprintf(f,"#%i: y-torque %s\n",nc++,IntegrandString);
     if (ErrorString) fprintf(f,"#%i: y-torque %s \n",nc++,ErrorString);
   };

  if ( SC3D->WhichQuantities & QUANTITY_TORQUE3 )
   { fprintf(f,"#%i: z-torque %s\n",nc++,IntegrandString);
     if (ErrorString) fprintf(f,"#%i: z-torque %s \n",nc++,ErrorString);
   };

  fclose(f);

  SC3D->UTIntegralBuffer[0]=SC3D->UTIntegralBuffer[1]=0;

  SC3D->BZIArgs=0;

}
