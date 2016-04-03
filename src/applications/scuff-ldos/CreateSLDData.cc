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
void WriteFilePreamble(char *FileName, int FileType, int LDim)
{
  FILE *f=fopen(FileName,"a");

  fprintf(f,"# scuff-ldos run on %s ",GetHostName());
  fprintf(f,"%s\n",GetTimeString());
  fprintf(f,"# columns: \n");
  fprintf(f,"# 1 2 3 4 5 : x y z re(Omega) im(Omega)\n");
  int nc=6;

  if (FileType==FILETYPE_BYK && LDim==1)
   { fprintf(f,"# %i: kx\n",nc);
     nc+=1;
   }
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
SLDData *CreateSLDData(char *GeoFile, char **EPFiles, int nEPFiles)
{
  SetDefaultCD2SFormat("%.8e %.8e");

  SLDData *Data=(SLDData *)mallocEC(sizeof(*Data));
  Data->RelTol      = 1.0e-3;
  Data->MaxEvals    = 10000;
  Data->HalfSpaceMP = 0;

  /***************************************************************/
  /* read in geometry and allocate BEM matrix and RHS vector     */
  /***************************************************************/
  RWGGeometry *G = Data->G = new RWGGeometry(GeoFile);
  Data->M = G->AllocateBEMMatrix();

  /***************************************************************/
  /* read in lists of evaluation points **************************/
  /***************************************************************/
  HMatrix **XMatrices=0, **GMatrices=0;
  char **EPFileBases=0;
  int NumXGMatrices=0, TotalEvalPoints=0;
  for(int n=0; n<nEPFiles; n++)
   { 
     HMatrix *XMatrix = new HMatrix(EPFiles[n]);
     if (XMatrix->ErrMsg)
      { Warn("%s (skipping file)",XMatrix->ErrMsg);
        continue;
      };
     int NX = XMatrix->NR;

     XMatrices = (HMatrix **)reallocEC(XMatrices, (NumXGMatrices+1)*sizeof(HMatrix *));
     GMatrices = (HMatrix **)reallocEC(GMatrices, (NumXGMatrices+1)*sizeof(HMatrix *));
     EPFileBases = (char **)reallocEC(EPFileBases, (NumXGMatrices+1)*sizeof(char *));

     XMatrices[NumXGMatrices]   = XMatrix;
     GMatrices[NumXGMatrices]   = new HMatrix(NX, 18, LHM_COMPLEX);
     EPFileBases[NumXGMatrices] = strdup(GetFileBase(EPFiles[n]));
   
     NumXGMatrices++;
     TotalEvalPoints+=NX;
     
     Log("Read %i evaluation points from file %s.",NX,EPFiles[n]);
   };
  Log("%i total evaluation points from %i files.",TotalEvalPoints,NumXGMatrices);
  Data->XMatrices       = XMatrices;
  Data->GMatrices       = GMatrices;
  Data->EPFileBases     = EPFileBases;
  Data->NumXGMatrices   = NumXGMatrices;
  Data->TotalEvalPoints = TotalEvalPoints;

  Data->WrotePreamble[0] = (bool *)mallocEC(NumXGMatrices*sizeof(bool));
  Data->WrotePreamble[1] = (bool *)mallocEC(NumXGMatrices*sizeof(bool));
  for(int nm=0; nm<NumXGMatrices; nm++)
   Data->WrotePreamble[0][nm]=Data->WrotePreamble[1][nm]=false;

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

