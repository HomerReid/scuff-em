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
const char *QuantityNames[NUMPFT]=
 { "Power",
   "XForce",  "YForce",  "ZForce",
   "XTorque", "YTorque", "ZTorque"
 };

/***************************************************************/
/***************************************************************/
/***************************************************************/
SNEQData *CreateSNEQData(char *GeoFile, char *TransFile,
                         int QuantityFlags, char *EPFile,
                         bool PlotFlux, char *pFileBase)
{

  SNEQData *SNEQD=(SNEQData *)mallocEC(sizeof(*SNEQD));
  SNEQD->WriteCache=0;

  /*--------------------------------------------------------------*/
  /*-- try to create the RWGGeometry -----------------------------*/
  /*--------------------------------------------------------------*/
  RWGGeometry *G=new RWGGeometry(GeoFile);
  SNEQD->G=G;
  
  if (pFileBase)
   SNEQD->FileBase = strdup(pFileBase);
  else
   SNEQD->FileBase = strdup(GetFileBase(G->GeoFileName));

  /*--------------------------------------------------------------*/
  /*- this code does not make sense if any of the objects are PEC */
  /*--------------------------------------------------------------*/
  for(int ns=0; ns<G->NumSurfaces; ns++)
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
  bool *NeedQuantity=SNEQD->NeedQuantity;
  memset(NeedQuantity, 0, NUMPFT*sizeof(bool));
  int NQ=0;
  if ( QuantityFlags & QFLAG_POWER  ) 
   { NeedQuantity[QINDEX_POWER]   = true; NQ++; };
  if ( QuantityFlags & QFLAG_XFORCE ) 
   { NeedQuantity[QINDEX_XFORCE]  = true; NQ++; };
  if ( QuantityFlags & QFLAG_YFORCE ) 
   { NeedQuantity[QINDEX_YFORCE]  = true; NQ++; };
  if ( QuantityFlags & QFLAG_ZFORCE ) 
   { NeedQuantity[QINDEX_ZFORCE]  = true; NQ++; };
  if ( QuantityFlags & QFLAG_XTORQUE ) 
   { NeedQuantity[QINDEX_XTORQUE] = true; NQ++; };
  if ( QuantityFlags & QFLAG_YTORQUE ) 
   { NeedQuantity[QINDEX_YTORQUE] = true; NQ++; };
  if ( QuantityFlags & QFLAG_ZTORQUE ) 
   { NeedQuantity[QINDEX_ZTORQUE] = true; NQ++; };
  
  SNEQD->NQ = NQ;

  int NT = SNEQD->NumTransformations;
  int NS = SNEQD->G->NumSurfaces;
  SNEQD->NumSIQs = NT*NS*NS*NQ;

  /*--------------------------------------------------------------*/
  /*- read the list of evaluation points for spatially-resolved  -*/
  /*- quantities                                                 -*/
  /*--------------------------------------------------------------*/
  SNEQD->XPoints = 0;
  SNEQD->NX      = 0; 
  if (EPFile)
   { SNEQD->XPoints = new HMatrix(EPFile);
     if (SNEQD->XPoints->ErrMsg)
      ErrExit(SNEQD->XPoints->ErrMsg);
     SNEQD->NX = SNEQD->XPoints->NR;
     ErrExit("--EPFile option is not yet supported");
   };
  SNEQD->NumSRQs = NT*NS*(SNEQD->NX)*NQ;

  /*--------------------------------------------------------------*/
  /*- allocate arrays of matrix subblocks that allow us to reuse -*/
  /*- chunks of the BEM matrices for multiple geometrical        -*/
  /*- transformations.                                           -*/
  /*--------------------------------------------------------------*/
  SNEQD->T = (HMatrix **)mallocEC(NS*sizeof(HMatrix *));
  SNEQD->TSelf = (HMatrix **)mallocEC(NS*sizeof(HMatrix *));
  SNEQD->U = (HMatrix **)mallocEC( ((NS*(NS-1))/2)*sizeof(HMatrix *));
  Log("Before T, U blocks: mem=%3.1f GB",GetMemoryUsage()/1.0e9);
  for(int nb=0, ns=0; ns<NS; ns++)
   { 
     int NBF=G->Surfaces[ns]->NumBFs;

     if (G->Mate[ns]==-1)
      { SNEQD->T[ns]     = new HMatrix(NBF, NBF, LHM_COMPLEX);
        SNEQD->TSelf[ns] = new HMatrix(NBF, NBF, LHM_COMPLEX);
      }
     else
      { SNEQD->T[ns]     = SNEQD->T[ G->Mate[ns] ];
        SNEQD->TSelf[ns] = SNEQD->TSelf[ G->Mate[ns] ];
      };

     for(int nsp=ns+1; nsp<G->NumSurfaces; nsp++, nb++)
      { int NBFp=G->Surfaces[nsp]->NumBFs;
        SNEQD->U[nb] = new HMatrix(NBF, NBFp, LHM_COMPLEX);
      };
   };
  Log("After T, U blocks: mem=%3.1f GB",GetMemoryUsage()/1.0e9);

  /*--------------------------------------------------------------*/
  /*- allocate BEM matrix and Sigma matrix -----------------------*/
  /*--------------------------------------------------------------*/
  SNEQD->W     = new HMatrix(G->TotalBFs, G->TotalBFs, LHM_COMPLEX );
  SNEQD->Sigma = new HMatrix(G->TotalBFs, G->TotalBFs, LHM_COMPLEX );
  Log("After W, Sigma: mem=%3.1f GB",GetMemoryUsage()/1.0e9);

  /*--------------------------------------------------------------*/
  /*- Buffer[0..nBuffer-1] are data storage buffers with enough  -*/
  /*- room to hold MaxBFs^2 cdoubles, where MaxBFs is the maximum-*/
  /*- number of basis functions on any object, i.e. the max      -*/
  /*- dimension of any BEM matrix subblock.                      -*/
  /*--------------------------------------------------------------*/
  int MaxBFs=G->Surfaces[0]->NumBFs;
  for(int ns=1; ns<G->NumSurfaces; ns++)
   if (G->Surfaces[ns]->NumBFs > MaxBFs) 
    MaxBFs = G->Surfaces[ns]->NumBFs;
  
  int nBuffer = 3;
  int BufSize = MaxBFs * MaxBFs * sizeof(cdouble);
  SNEQD->Buffer[0] = mallocEC(nBuffer*BufSize);
  for(int nb=1; nb<nBuffer; nb++)
   SNEQD->Buffer[nb] = (void *)( (char *)SNEQD->Buffer[nb-1] + BufSize);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int fdim = SNEQD->NumSIQs +  SNEQD->NumSRQs;
  SNEQD->OmegaConverged = (bool *)mallocEC(fdim*sizeof(bool));

  /*--------------------------------------------------------------*/
  /*- ByEdge[nd][nq][ne] = contribution of edge #ne on surface   -*/
  /*- #nd to the flux of quantity #nq. This is used to produce   -*/
  /*- plots of spatially-resolved flux density.                  -*/
  /*--------------------------------------------------------------*/
  SNEQD->ByEdge=0;
  if (PlotFlux)
   { double ***ByEdge = (double ***)mallocEC(NS*sizeof(double **));
     for(int nds=0; nds<NS; nds++)
      { ByEdge[nds] = (double **)mallocEC(NUMPFT*sizeof(double **));
        int NE=G->Surfaces[nds]->NumEdges;
        for(int nq=0; nq<NUMPFT; nq++)
         { if (NeedQuantity[nq])
            ByEdge[nds][nq]=(double *)mallocEC(NE*sizeof(double));
         };
      };
     SNEQD->ByEdge=ByEdge;
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  FILE *f=vfopen("%s.SIFlux","a",SNEQD->FileBase);
  fprintf(f,"\n");
  fprintf(f,"# scuff-neq run on %s (%s)\n",GetHostName(),GetTimeString());
  fprintf(f,"# data file columns: \n");
  fprintf(f,"# 1 omega \n");
  fprintf(f,"# 2 transform tag\n");
  int nq=3;
  if (G->LDim==1)
   fprintf(f,"# %i kBloch_x \n",nq++);
  else if (G->LDim==2)
   { fprintf(f,"# %i kBloch_x \n",nq++);
     fprintf(f,"# %i kBloch_x \n",nq++);
   };
  fprintf(f,"# %i (sourceObject,destObject) \n",nq++);
  for(int nPFT=0; nPFT<NUMPFT; nPFT++)
   if (NeedQuantity[nPFT])
    fprintf(f,"# %i %s flux spectral density\n",nq++,QuantityNames[nPFT]);
  fclose(f);

  Log("After CreateSNEQData: mem=%3.1f GB",GetMemoryUsage()/1.0e9);
  return SNEQD;

}
