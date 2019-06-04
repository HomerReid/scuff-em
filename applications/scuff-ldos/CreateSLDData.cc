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
 * CreateSLDData -- utility routine to initialize an Data
 *               -- structure containing all information passed
 *               -- around between the various scuff-ldos routines
 *
 * homer reid      -- 3/2015
 *
 */

#include "scuff-ldos.h"

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WriteFilePreamble(char *FileName, int FileType, int LDim,
                       bool HaveGTCList, bool TwoPointDGF)
{
  FILE *f=fopen(FileName,"a");

  fprintf(f,"# scuff-ldos " VERSIONSTRING " run on %s ",GetHostName());
  fprintf(f,"%s\n",GetTimeString());
  fprintf(f,"# columns: \n");
  
  int nc;
  if (TwoPointDGF)
   { fprintf(f,"# 1 2 3: (x, y, z) (destination point) \n");
     fprintf(f,"# 4 5 6: (x, y, z) (source point) \n");
     nc=7;
   }
  else
   { fprintf(f,"# 1 2 3:\n");
     nc=4;
   };

  fprintf(f,"# %i %i : real(Omega) imag(Omega)\n",nc,nc+1);
  nc+=2;

  if (HaveGTCList) 
   fprintf(f,"# %i: transform label ",nc++);

  if (FileType==FILETYPE_BYK && LDim==1)
   fprintf(f,"# %i: kx\n",nc++);
  else if (FileType==FILETYPE_BYK && LDim==2)
   { fprintf(f,"# %i,%i: kx ky\n",nc,nc+1);
     nc+=2;
   };

  fprintf(f,"# %2i: electric LDOS\n",nc++);
  fprintf(f,"# %2i: magnetic LDOS\n",nc++);
  fprintf(f,"# %2i, %2i: re, im GE_{xx} \n",nc,nc+1); nc+=2;
  fprintf(f,"# %2i, %2i: re, im GE_{xy} \n",nc,nc+1); nc+=2;
  fprintf(f,"# %2i, %2i: re, im GE_{xz} \n",nc,nc+1); nc+=2;
  fprintf(f,"# %2i, %2i: re, im GE_{yx} \n",nc,nc+1); nc+=2;
  fprintf(f,"# %2i, %2i: re, im GE_{yy} \n",nc,nc+1); nc+=2;
  fprintf(f,"# %2i, %2i: re, im GE_{yz} \n",nc,nc+1); nc+=2;
  fprintf(f,"# %2i, %2i: re, im GE_{zx} \n",nc,nc+1); nc+=2;
  fprintf(f,"# %2i, %2i: re, im GE_{zy} \n",nc,nc+1); nc+=2;
  fprintf(f,"# %2i, %2i: re, im GE_{zz} \n",nc,nc+1); nc+=2;
  fprintf(f,"# %2i, %2i: re, im GM_{xx} \n",nc,nc+1); nc+=2;
  fprintf(f,"# %2i, %2i: re, im GM_{xy} \n",nc,nc+1); nc+=2;
  fprintf(f,"# %2i, %2i: re, im GM_{xz} \n",nc,nc+1); nc+=2;
  fprintf(f,"# %2i, %2i: re, im GM_{yx} \n",nc,nc+1); nc+=2;
  fprintf(f,"# %2i, %2i: re, im GM_{yy} \n",nc,nc+1); nc+=2;
  fprintf(f,"# %2i, %2i: re, im GM_{yz} \n",nc,nc+1); nc+=2;
  fprintf(f,"# %2i, %2i: re, im GM_{zx} \n",nc,nc+1); nc+=2;
  fprintf(f,"# %2i, %2i: re, im GM_{zy} \n",nc,nc+1); nc+=2;
  fprintf(f,"# %2i, %2i: re, im GM_{zz} \n",nc,nc+1); nc+=2;

  fclose(f);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
SLDData *CreateSLDData(char *GeoFile, char *TransFile,
                       char **EPFiles, int nEPFiles)
{
  SetDefaultCD2SFormat("%.8e %.8e");

  SLDData *Data=(SLDData *)mallocEC(sizeof(*Data));
  Data->RelTol      = 1.0e-2;
  Data->AbsTol      = 1.0e-8;
  Data->MaxEvals    = 100000;
  Data->HalfSpaceMP = 0;
  Data->GroundPlane = false;

  /***************************************************************/
  /* read in geometry and allocate BEM matrix and RHS vector     */
  /***************************************************************/
  RWGGeometry *G = Data->G = new RWGGeometry(GeoFile);
  Data->M = G->AllocateBEMMatrix();

  /***************************************************************/
  /* read in geometrical transformation file if any **************/
  /***************************************************************/
  Data->GTCs=ReadTransFile(TransFile);
  Data->NumTransforms = Data->GTCs.size();
  char *ErrMsg=G->CheckGTCList(Data->GTCs);
  if (ErrMsg)
   ErrExit("file %s: %s",TransFile,ErrMsg);
  bool HaveGTCList = (TransFile!=0);

  Data->TBlocks=Data->UBlocks=0;
  if (HaveGTCList)
   { int NS=G->NumSurfaces;
     int NADB = NS*(NS-1)/2; // number of above-diagonal blocks
     HMatrix **TBlocks = Data->TBlocks = (HMatrix **)mallocEC(NS*sizeof(HMatrix *));
     HMatrix **UBlocks = Data->UBlocks = (HMatrix **)mallocEC(NADB*sizeof(HMatrix *));
     for(int ns=0, nb=0; ns<NS; ns++)
      { 
        int nsMate = G->Mate[ns];
        if ( nsMate!=-1 )
         TBlocks[ns] = TBlocks[nsMate];
        else
         { int NBF=G->Surfaces[ns]->NumBFs;
           TBlocks[ns] = new HMatrix(NBF, NBF, Data->M->RealComplex);
         };

        for(int nsp=ns+1; nsp<NS; nsp++, nb++)
         { int NBF=G->Surfaces[ns]->NumBFs;
           int NBFp=G->Surfaces[nsp]->NumBFs;
           UBlocks[nb] = new HMatrix(NBF, NBFp, Data->M->RealComplex);
         };
      };
   };

  /***************************************************************/
  /* read in lists of evaluation points **************************/
  /***************************************************************/
  HMatrix **XMatrices=0;
  char **EPFileBases=0;
  int NumXMatrices=0, TotalEvalPoints=0;
  for(int n=0; n<nEPFiles; n++)
   { 
     HMatrix *XMatrix = new HMatrix(EPFiles[n]);
     if (XMatrix->ErrMsg)
      { Warn("%s (skipping file)",XMatrix->ErrMsg);
        continue;
      };
     int NX = XMatrix->NR;

     XMatrices = (HMatrix **)reallocEC(XMatrices, (NumXMatrices+1)*sizeof(HMatrix *));
     EPFileBases = (char **)reallocEC(EPFileBases, (NumXMatrices+1)*sizeof(char *));

     XMatrices[NumXMatrices]   = XMatrix;
     EPFileBases[NumXMatrices] = strdup(GetFileBase(EPFiles[n]));
   
     NumXMatrices++;
     TotalEvalPoints+=NX;
     
     Log("Read %i evaluation points from file %s.",NX,EPFiles[n]);
   };
  Log("%i total evaluation points from %i files.",TotalEvalPoints,NumXMatrices);
  Data->XMatrices       = XMatrices;
  Data->EPFileBases     = EPFileBases;
  Data->NumXMatrices    = NumXMatrices;
  Data->TotalEvalPoints = TotalEvalPoints;

  Data->WrotePreamble[0] = (bool *)mallocEC(NumXMatrices*sizeof(bool));
  Data->WrotePreamble[1] = (bool *)mallocEC(NumXMatrices*sizeof(bool));
  for(int nm=0; nm<NumXMatrices; nm++)
   Data->WrotePreamble[0][nm]=Data->WrotePreamble[1][nm]=false;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NumTransforms = Data->NumTransforms;
  int NumGMatrices = NumXMatrices * NumTransforms;
  HMatrix **GMatrices = 
   (HMatrix **)mallocEC(NumGMatrices*sizeof(HMatrix *));
  for(int nt=0; nt<NumTransforms; nt++)
   for(int nm=0; nm<NumXMatrices; nm++)
    GMatrices[nt*NumXMatrices + nm]
     = new HMatrix(XMatrices[nm]->NR, 18, LHM_COMPLEX);
  Data->GMatrices = GMatrices;

  /***************************************************************/
  /* For PBC geometries we need to do some preliminary setup     */
  /***************************************************************/
  Data->ABMBCache=0;
  if (G->LDim>0)
   { 
     int NS = G->NumSurfaces;
     int NB = NS*(NS+1)/2;
     Data->ABMBCache = (void **)malloc(NB*sizeof(void *));
     for(int nsa=0, nb=0; nsa<NS; nsa++)
      for(int nsb=nsa; nsb<NS; nsb++, nb++)
       Data->ABMBCache[nb]=G->CreateABMBAccelerator(nsa, nsb, false, false);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  return Data;
}

