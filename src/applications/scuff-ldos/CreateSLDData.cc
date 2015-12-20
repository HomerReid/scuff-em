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
  fprintf(f,"# 1 2 3 4: x y z Omega\n");
  int nc=5;

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
  fprintf(f,"# %2i, %2i: re, im GM_{zz} \n",nc,nc+1); nc+=2;
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
SLDData *CreateSLDData(char *GeoFile, char *EPFile)
{
  SetDefaultCD2SFormat("%.8e %.8e");

  SLDData *Data=(SLDData *)mallocEC(sizeof(*Data));
  Data->RelTol      = 1.0e-2;
  Data->MaxEvals    = 1000;
  Data->HalfSpaceMP = 0;

  /***************************************************************/
  /* read in geometry and allocate BEM matrix and RHS vector     */
  /***************************************************************/
  RWGGeometry *G = Data->G = new RWGGeometry(GeoFile);
  Data->M = G->AllocateBEMMatrix();

  /***************************************************************/
  /* read in list of evaluation points ***************************/
  /***************************************************************/
  Data->XMatrix = new HMatrix(EPFile);
  if (Data->XMatrix->ErrMsg)
   ErrExit(Data->XMatrix->ErrMsg); 

  Data->GMatrix = new HMatrix(Data->XMatrix->NR, 18, LHM_COMPLEX);

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

