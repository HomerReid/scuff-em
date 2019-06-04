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
void WriteSIFluxFilePreamble(SNEQData *SNEQD, char *FileName, bool ByRegion=false)
{
  FILE *f = ByRegion ? vfopen("%s.byRegion","a",FileName) : fopen(FileName,"a");
  fprintf(f,"\n");
  fprintf(f,"# scuff-neq " VERSIONSTRING " run on ");
  fprintf(f,"%s (%s)\n",GetHostName(),GetTimeString());
  fprintf(f,"# data file columns: \n");
  fprintf(f,"# 1 transform tag\n");
  fprintf(f,"# 2 omega \n");
  int nq=3;
  int LDim = SNEQD->G->LBasis ? SNEQD->G->LBasis->NC : 0;
  if (LDim>=1)
   fprintf(f,"# %i kBloch_x \n",nq++);
  if (LDim>=2)
   fprintf(f,"# %i kBloch_y \n",nq++);
  if (ByRegion)
   fprintf(f,"# %i (sourceRegion,destRegion) \n",nq++);
  else
   fprintf(f,"# %i (sourceSurface,destSurface) \n",nq++);
  for(int nPFT=0; nPFT<NUMPFT; nPFT++)
   fprintf(f,"# %i %s flux spectral density\n",
  nq++,QuantityNames[nPFT]);

  fclose(f);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
SNEQData *CreateSNEQData(char *GeoFile, char *TransFile,
                         int *PFTMethods, int NumPFTMethods,
                         char *EPFile, char *pFileBase)
{

  SNEQData *SNEQD=(SNEQData *)mallocEC(sizeof(*SNEQD));
  SNEQD->WriteCache=0;
  SNEQD->SourceOnly=-1;
  SNEQD->DestOnly=-1;

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
  /*- read the transformation file if one was specified and check */
  /*- that it plays well with the specified geometry file.        */
  /*- note if TransFile==0 then this code snippet still works; in */
  /*- this case the list of GTComplices is initialized to contain */
  /*- a single empty GTComplex and the check automatically passes.*/
  /*--------------------------------------------------------------*/
  SNEQD->GTCs=ReadTransFile(TransFile);
  SNEQD->NumTransformations = SNEQD->GTCs.size();
  char *ErrMsg=G->CheckGTCList(SNEQD->GTCs);
  if (ErrMsg)
   ErrExit("file %s: %s",TransFile,ErrMsg);
  if (SNEQD->NumTransformations==1 && !strcmp(SNEQD->GTCs[0]->Tag,"DEFAULT"))
   sprintf(SNEQD->GTCs[0]->Tag,"0.0");

  /*--------------------------------------------------------------*/
  /*- figure out which PFT methods were requested and write       */
  /*- SIFlux file preambles                                       */
  /*--------------------------------------------------------------*/
  SNEQD->PFTMatrix = new HMatrix(G->NumSurfaces, NUMPFT);
  InitPFTOptions( &(SNEQD->PFTOpts) );
  SNEQD->NumPFTMethods = NumPFTMethods;
  SNEQD->DSIOmegaPoints=0;
  for(int npm=0; npm<NumPFTMethods; npm++)
   { 
     SNEQD->PFTMethods[npm] = PFTMethods[npm];

     char PFTName[20];
     if (PFTMethods[npm] == -SCUFF_PFT_DSI)
      { SNEQD->PFTMethods[npm] = SCUFF_PFT_DSI;
        snprintf(PFTName,20,"DSIMesh");
      }
     else if (PFTMethods[npm]>20)
      { SNEQD->PFTMethods[npm]=SCUFF_PFT_DSI;
        snprintf(PFTName,20,"DSI%i",PFTMethods[npm]);
      }
     else if (PFTMethods[npm]==SCUFF_PFT_OVERLAP)
      sprintf(PFTName,"OPFT");
     else if (PFTMethods[npm]==SCUFF_PFT_EMT_EXTERIOR)
      sprintf(PFTName,"EMTPFT");

     SNEQD->SIFluxFileNames[npm]
      = vstrdup("%s.SIFlux.%s",SNEQD->FileBase,PFTName);
     WriteSIFluxFilePreamble(SNEQD, SNEQD->SIFluxFileNames[npm]);
   };
  
  int NT = SNEQD->NumTransformations;
  int NS = SNEQD->G->NumSurfaces;
  int NR = SNEQD->G->NumRegions;

  /*--------------------------------------------------------------*/
  /*- read the list of evaluation points for spatially-resolved  -*/
  /*- quantities                                                 -*/
  /*--------------------------------------------------------------*/
  SNEQD->SRXMatrix = 0;
  SNEQD->SRFMatrix = 0;
  SNEQD->NX        = 0; 
  SNEQD->NumSRQs   = 0;
  if (EPFile)
   { 
     // if the EPFile has extension .msh, attempt to
     // read it in as a GMSH msh and take the vertices as
     // the eval-point coordinates
     if(!strcasecmp(GetFileExtension(EPFile),"msh"))
      { RWGSurface *S=new RWGSurface(EPFile);
        if(S->ErrMsg) ErrExit(S->ErrMsg);
        SNEQD->SRXMatrix = new HMatrix(S->NumVertices, 3);
        for(int nv=0; nv<S->NumVertices; nv++)
         SNEQD->SRXMatrix->SetEntriesD(nv,"0:2",S->Vertices + 3*nv);
        Log("Read %i SR eval points from mesh %s",S->NumVertices,EPFile);
        delete S;
      }
     else
      {
        SNEQD->SRXMatrix = new HMatrix(EPFile);
        if (SNEQD->SRXMatrix->ErrMsg)
         ErrExit(SNEQD->SRXMatrix->ErrMsg); 
      }
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
  /*- allocate BEM matrix and dressed Rytov matrix ---------------*/
  /*--------------------------------------------------------------*/
  SNEQD->M        = new HMatrix(G->TotalBFs, G->TotalBFs, LHM_COMPLEX );
  SNEQD->DRMatrix = new HMatrix(G->TotalBFs, G->TotalBFs, LHM_COMPLEX );
  Log("After W, Rytov: mem=%3.1f GB",GetMemoryUsage()/1.0e9);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (SNEQD->NumSRQs>0)
   { FILE *f=vfopen("%s.SRFlux","a",SNEQD->FileBase);
     fprintf(f,"\n");
     fprintf(f,"# scuff-neq " VERSIONSTRING " run on %s (%s)\n",GetHostName(),GetTimeString());
     fprintf(f,"# data file columns: \n");
     fprintf(f,"# 1 transform tag\n");
     fprintf(f,"# 2 omega \n");
     int nq=3;
     int LDim = SNEQD->G->LBasis ? SNEQD->G->LBasis->NC : 0;
     if (LDim>=1)
      fprintf(f,"# %i kBloch_x \n",nq++);
     if (LDim==2)
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

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  SNEQD->PFTByRegion=0;
  SNEQD->RegionRegionPFT=0;
  if( NS > (NR-1) ) 
   { SNEQD->PFTByRegion = new HMatrix(NR, NUMPFT);
     SNEQD->RegionRegionPFT = (HMatrix **)mallocEC( NumPFTMethods*sizeof(HMatrix *));
     for(int npm=0; npm<NumPFTMethods; npm++)
      { SNEQD->RegionRegionPFT[npm]=new HMatrix( (NR+1)*(NR+1), NUMPFT);
        WriteSIFluxFilePreamble(SNEQD, SNEQD->SIFluxFileNames[npm], true);
      };
   };

  Log("After CreateSNEQData: mem=%3.1f GB",GetMemoryUsage()/1.0e9);
  return SNEQD;

}
