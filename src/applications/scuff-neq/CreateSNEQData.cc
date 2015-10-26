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
 { "PAbs", "PRad",
   "XForce",  "YForce",  "ZForce",
   "XTorque", "YTorque", "ZTorque"
 };

/***************************************************************/
/***************************************************************/
/***************************************************************/
SNEQData *CreateSNEQData(char *GeoFile, char *TransFile,
                         char **TempStrings, int nTempStrings,
                         int QuantityFlags, char *EPFile,
                         char *pFileBase, bool JDEPFT)
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
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  InitPFTOptions( &(SNEQD->PFTOpts) );
  SNEQD->NeedQuantity = SNEQD->PFTOpts.NeedQuantity;

  /*******************************************************************/
  /* process --temperature options ***********************************/
  /*******************************************************************/
  SNEQD->TEnvironment=0.0;
  SNEQD->TSurfaces=(double *)mallocEC(G->NumSurfaces*sizeof(double));
  for(int nts=0; nts<nTempStrings; nts++)
   { 
     double TTemp;
     int WhichSurface;
     if ( 1!=sscanf(TempStrings[2*nts+1],"%le",&TTemp) )
      ErrExit("invalid temperature (%s) passed for --temperature option",TempStrings[2*nts+1]);

     if (    !strcasecmp(TempStrings[2*nts],"MEDIUM")
          || !strcasecmp(TempStrings[2*nts],"EXTERIOR")
          || !strcasecmp(TempStrings[2*nts],"ENVIRONMENT")
        )
      { SNEQD->TEnvironment=TTemp;
        Log("Setting environment temperature to %g kelvin.",TTemp);
        printf("Setting environment temperature to %g kelvin.\n",TTemp);
      }
     else if ( G->GetSurfaceByLabel(TempStrings[2*nts],&WhichSurface) )
      { SNEQD->TSurfaces[WhichSurface]=TTemp;
        Log("Setting temperature of object %s to %g kelvin.",TempStrings[2*nts],TTemp);
        printf("Setting temperature of object %s to %g kelvin.\n",TempStrings[2*nts],TTemp);
      }
     else 
      ErrExit("unknown surface/region %s in --temperature specification",TempStrings[2*nts]);
   };

  /*--------------------------------------------------------------*/
  /*- figure out which quantities were specified                 -*/
  /*--------------------------------------------------------------*/
  bool *NeedQuantity=SNEQD->NeedQuantity;
  memset(NeedQuantity, 0, NUMPFT*sizeof(bool));
  int NQ=0;
  if ( QuantityFlags & QFLAG_PABS   ) 
   { NeedQuantity[QINDEX_PABS]    = true; NQ++; };
  if ( QuantityFlags & QFLAG_PRAD   ) 
   { NeedQuantity[QINDEX_PRAD]    = true; NQ++; };
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
  
  SNEQD->NPFT = NQ;
  int NT = SNEQD->NumTransformations;
  int NS = SNEQD->G->NumSurfaces;
  SNEQD->NumSIQs = NT*NS*NS*NQ;

  /*--------------------------------------------------------------*/
  /*- read the list of evaluation points for spatially-resolved  -*/
  /*- quantities                                                 -*/
  /*--------------------------------------------------------------*/
  SNEQD->SRXMatrix = 0;
  SNEQD->SRFMatrix = 0;
  SNEQD->NX        = 0; 
  SNEQD->NumSRQs   = 0;
  if (EPFile)
   { SNEQD->SRXMatrix = new HMatrix(EPFile);
     if (SNEQD->SRXMatrix->ErrMsg)
      ErrExit(SNEQD->SRXMatrix->ErrMsg);
     int NX = SNEQD->NX = SNEQD->SRXMatrix->NR;
     SNEQD->NumSRQs = NT*NS*NX*NUMSRFLUX;
     SNEQD->SRFMatrix = new HMatrix(NX, NUMSRFLUX);
   };

  /*--------------------------------------------------------------*/
  /*- allocate arrays of matrix subblocks that allow us to reuse -*/
  /*- chunks of the BEM matrices for multiple geometrical        -*/
  /*- transformations.                                           -*/
  /*--------------------------------------------------------------*/
  SNEQD->TExt = (HMatrix **)mallocEC(NS*sizeof(HMatrix *));
  SNEQD->TInt = (HMatrix **)mallocEC(NS*sizeof(HMatrix *));
  SNEQD->U = (HMatrix **)mallocEC( ((NS*(NS-1))/2)*sizeof(HMatrix *));
  Log("Before T, U blocks: mem=%3.1f GB",GetMemoryUsage()/1.0e9);
  for(int nb=0, ns=0; ns<NS; ns++)
   { 
     int NBF=G->Surfaces[ns]->NumBFs;

     if (G->Mate[ns]==-1)
      { SNEQD->TExt[ns]  = new HMatrix(NBF, NBF, LHM_COMPLEX);
        SNEQD->TInt[ns]  = new HMatrix(NBF, NBF, LHM_COMPLEX);
      }
     else
      { SNEQD->TExt[ns] = SNEQD->TExt[ G->Mate[ns] ];
        SNEQD->TInt[ns] = SNEQD->TInt[ G->Mate[ns] ];
      };

     for(int nsp=ns+1; nsp<G->NumSurfaces; nsp++, nb++)
      { int NBFp=G->Surfaces[nsp]->NumBFs;
        SNEQD->U[nb] = new HMatrix(NBF, NBFp, LHM_COMPLEX);
      };
   };
  Log("After T, U blocks: mem=%3.1f GB",GetMemoryUsage()/1.0e9);

  /*--------------------------------------------------------------*/
  /*- allocate BEM matrix and Rytov matrix -----------------------*/
  /*--------------------------------------------------------------*/
  SNEQD->W     = new HMatrix(G->TotalBFs, G->TotalBFs, LHM_COMPLEX );
  SNEQD->Rytov = new HMatrix(G->TotalBFs, G->TotalBFs, LHM_COMPLEX );
  Log("After W, Rytov: mem=%3.1f GB",GetMemoryUsage()/1.0e9);

  /*--------------------------------------------------------------*/
  /*- Buffer[0..nBuffer-1] are data storage buffers with enough  -*/
  /*- room to hold MaxBFs^2 cdoubles, where MaxBFs is the maximum-*/
  /*- number of basis functions on any object, i.e. the max      -*/
  /*- dimension of any BEM matrix subblock.                      -*/
  /*--------------------------------------------------------------*/
  size_t MaxBFs=G->Surfaces[0]->NumBFs;
  for(int ns=1; ns<G->NumSurfaces; ns++)
   if ( (size_t)G->Surfaces[ns]->NumBFs > MaxBFs) 
    MaxBFs = G->Surfaces[ns]->NumBFs;
  
  int nBuffer = 3;
  size_t BufSize = MaxBFs * MaxBFs * sizeof(cdouble);
  SNEQD->Buffer[0] = mallocEC(nBuffer*BufSize);
  for(int nb=1; nb<nBuffer; nb++)
   SNEQD->Buffer[nb] = (void *)( (char *)SNEQD->Buffer[nb-1] + BufSize);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int fdim = SNEQD->NumSIQs + SNEQD->NumSRQs;
  SNEQD->OmegaConverged = (bool *)mallocEC(fdim*sizeof(bool));

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  SNEQD->JDEPFT = JDEPFT;
  SNEQD->SIFluxFileName
   = vstrdup("%s.%s",SNEQD->FileBase, JDEPFT ? "JDEPFT" : "SIFlux");
  if (SNEQD->NumSIQs>0)
   { FILE *f=fopen(SNEQD->SIFluxFileName,"a");
     fprintf(f,"\n");
     fprintf(f,"# scuff-neq run on %s (%s)\n",GetHostName(),GetTimeString());
     fprintf(f,"# data file columns: \n");
     fprintf(f,"# 1 transform tag\n");
     fprintf(f,"# 2 omega \n");
     int nq=3;
     if (G->LDim>=1)
      fprintf(f,"# %i kBloch_x \n",nq++);
     if (G->LDim==2)
      fprintf(f,"# %i kBloch_y \n",nq++);
     fprintf(f,"# %i (sourceObject,destObject) \n",nq++);
     for(int nPFT=0; nPFT<NUMPFT; nPFT++)
      if (NeedQuantity[nPFT])
       fprintf(f,"# %i %s flux spectral density\n",nq++,QuantityNames[nPFT]);
     fclose(f);
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (SNEQD->NumSRQs>0)
   { FILE *f=vfopen("%s.SRFlux","a",SNEQD->FileBase);
     fprintf(f,"\n");
     fprintf(f,"# scuff-neq run on %s (%s)\n",GetHostName(),GetTimeString());
     fprintf(f,"# data file columns: \n");
     fprintf(f,"# 1 transform tag\n");
     fprintf(f,"# 2 omega \n");
     int nq=3;
     if (G->LDim>=1)
      fprintf(f,"# %i kBloch_x \n",nq++);
     if (G->LDim==2)
      fprintf(f,"# %i kBloch_y \n",nq++);
     fprintf(f,"# %i, %i, %i x,y,z (coordinates of eval point)\n",nq,nq+1,nq+2);
     nq+=3;
     fprintf(f,"# %i sourceObject \n",nq++);
     fprintf(f,"# %2i %2i %2i Px,    Py,    Pz  \n",nq,nq+1,nq+2); nq+=3;
     fprintf(f,"# %2i %2i %2i Txx,   Txy,   Txz \n",nq,nq+1,nq+2); nq+=3;
     fprintf(f,"# %2i %2i %2i Tyx,   Tyy,   Tyz \n",nq,nq+1,nq+2); nq+=3;
     fprintf(f,"# %2i %2i %2i Tzx,   Tzy,   Tzz \n",nq,nq+1,nq+2); nq+=3;
     fclose(f);
   };

  Log("After CreateSNEQData: mem=%3.1f GB",GetMemoryUsage()/1.0e9);
  return SNEQD;

}
