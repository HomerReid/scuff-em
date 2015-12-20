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
 * CreateSCPData.cc  -- create the SCPD data structure that is 
 *                   -- passed around among scuff-caspol routines
 *
 * homer reid        -- 10/2006 -- 8/2013
 *
 */
#include "scuff-caspol.h"
#include <time.h>

/***************************************************************/
/***************************************************************/
/***************************************************************/
SCPData *CreateSCPData(char *GeoFile,
                       char **Atoms, int NumBIAtoms,
                       char **Particles, int NumParticles,
                       char *EPFile, char *FileBase)
{
  SCPData *SCPD=(SCPData *)mallocEC(sizeof(SCPData));

  /***************************************************************/
  /* create the RWGGeometry, allocate BEM matrix and RHS vector, and */
  /* set up the data structure passed to the computational routines  */
  /***************************************************************/
  if (GeoFile)
   { SCPD->G  = new RWGGeometry(GeoFile, SCUFF_TERSELOGGING);
     SCPD->M  = SCPD->G->AllocateBEMMatrix(SCUFF_PUREIMAGFREQ);
     SCPD->KN = SCPD->G->AllocateRHSVector(SCUFF_PUREIMAGFREQ);
   }
  else
   { SCPD->G  = 0; // in this case we take the
     SCPD->M  = 0; // geometry to be a PEC plate in the xy plane
     SCPD->KN = 0;
     GeoFile = strdup("PECPlate");
   };

  /*******************************************************************/
  /* create polarizability models for all specified atoms            */
  /*******************************************************************/
  int NumAtoms = SCPD->NumAtoms = NumBIAtoms + NumParticles;
  SCPD->PolModels    = (PolModel **)malloc(NumAtoms * sizeof(PolModel *));
  SCPD->Alphas       = (HMatrix **) malloc(NumAtoms * sizeof(HMatrix *) );

  int na=0;
  for(int nbi=0; nbi<NumBIAtoms; nbi++)
   { SCPD->PolModels[na] = new PolModel(Atoms[nbi], PM_BUILTIN);
     if (SCPD->PolModels[na]->ErrMsg)
      ErrExit(SCPD->PolModels[na]->ErrMsg);
     SCPD->Alphas[na] = new HMatrix(3,3);
     na++;
   };
  for(int np=0; np<NumParticles; np++)
   { SCPD->PolModels[na] = new PolModel(Particles[np], PM_FILE);
     if (SCPD->PolModels[na]->ErrMsg)
      ErrExit(SCPD->PolModels[na]->ErrMsg);
     SCPD->Alphas[na] = new HMatrix(3,3);
     na++;
   };

  /*******************************************************************/
  /* process list of evaluation points *******************************/
  /*******************************************************************/
  SCPD->EPMatrix = new HMatrix(EPFile);

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  if (FileBase==0)
   FileBase=strdup(GetFileBase(GeoFile));
  SCPD->ByXiFileName=vstrdup("%s.byXi",FileBase);

  /***************************************************************/
  /* For PBC geometries we need to do some preliminary setup     */
  /***************************************************************/
  SCPD->ABMBCache=0;
  SCPD->ByXikFileName=0;
  SCPD->GBZIArgs=0;
  RWGGeometry *G=SCPD->G;
  if (G && G->LDim>0)
   { 
     /*--------------------------------------------------------------*/
     /*- create output file for k-resolved data ---------------------*/
     /*--------------------------------------------------------------*/
     SCPD->ByXikFileName=vstrdup("%s.byXikFile",FileBase);

     /*--------------------------------------------------------------*/
     /* allocate accelerator for matrix block assembly -------------*/
     /*--------------------------------------------------------------*/
     int NS = G->NumSurfaces;
     int NB = NS*(NS+1)/2;
     SCPD->ABMBCache = (void **)malloc(NB*sizeof(void *));
     for(int nsa=0, nb=0; nsa<NS; nsa++)
      for(int nsb=nsa; nsb<NS; nsb++, nb++)
       SCPD->ABMBCache[nb]=G->CreateABMBAccelerator(nsa, nsb, false, false);
  
     /*--------------------------------------------------------------*/
     /*- initialize options structure for Brillouin-zone integration-*/
     /*--------------------------------------------------------------*/
     SCPD->GBZIArgs              = CreateGetBZIArgs(G->LBasis);
     SCPD->GBZIArgs->BZIFunc     = GetCPIntegrand;
     SCPD->GBZIArgs->UserData    = (void *)SCPD;
     SCPD->GBZIArgs->FDim        = NumAtoms * (SCPD->EPMatrix)->NR;

   };
     /*--------------------------------------------------------------*/

  return SCPD;
}
