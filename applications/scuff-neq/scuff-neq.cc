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
 * scuff-neq   -- a standalone code within the scuff-em suite
 *             -- for implementing the fluctuating-surface-current
 *             -- approach to nonequilibrium phenomena (more 
 *             -- specifically, for computing heat radiation, 
 *             -- heat transfer, and nonequilibrium casimir forces) 
 *
 * homer reid  -- 5/2012
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <BZIntegration.h>
#include "scuff-neq.h"
#include <libhrutil.h>

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define MAXFREQ  10    // max number of frequencies 
#define MAXCACHE 10    // max number of cache files for preload

#define MAXSTR   1000

#define II cdouble(0.0,1.0)

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  InstallHRSignalHandler();
  InitializeLog(argv[0]);

  /***************************************************************/
  /* pre-process command-line arguments to extract arguments     */
  /* any relevant for Brillouin-zone integration                 */
  /***************************************************************/
  //GetBZIArgStruct *BZIArgs=InitBZIArgs(argc, argv);

  /***************************************************************/
  /* process options *********************************************/
  /***************************************************************/
  char *GeoFile=0;
  char *TransFile=0;
  char *SourceObject=0;
  char *DestObject=0;

  /*--------------------------------------------------------------*/
  char *EPFile=0;

  /*--------------------------------------------------------------*/
  cdouble OmegaVals[MAXFREQ];        int nOmegaVals;
  char *OmegaFile;                   int nOmegaFiles;
  char *OmegaKBFile=0;

  /*--------------------------------------------------------------*/
  bool EMTPFT      = false;
  bool OPFT        = false;
  bool DSIPFT      = false;
  int DSIPoints    = 0;
  int DSIPoints2   = 0;
  char *DSIMesh    = 0;
  double DSIRadius = 0.0;
  bool DSIFarField = false;
  char *DSIOmegaFile = 0;

  /*--------------------------------------------------------------*/
  char *FileBase=0;

  /*--------------------------------------------------------------*/
  bool PlotFlux=false;
  bool OmitSelfTerms=false;

  /*--------------------------------------------------------------*/
  char *Cache=0;
  char *ReadCache[MAXCACHE];         int nReadCache;
  char *WriteCache=0;

  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { 
     {"Geometry",       PA_STRING,  1, 1,       (void *)&GeoFile,    0,             "geometry file"},
     {"TransFile",      PA_STRING,  1, 1,       (void *)&TransFile,  0,             "list of geometrical transformation"},
/**/
     {"SourceObject",   PA_STRING,  1, 1,       (void *)&SourceObject, 0,           "label of source object"},
     {"DestObject",     PA_STRING,  1, 1,       (void *)&DestObject, 0,             "label of destination object"},
/**/     
     {"EPFile",         PA_STRING,  1, 1,       (void *)&EPFile, 0,             "list of evaluation points for spatially-resolved flux"},
/**/     
     {"Omega",          PA_CDOUBLE, 1, MAXFREQ, (void *)OmegaVals,   &nOmegaVals,   "(angular) frequency"},
     {"OmegaFile",      PA_STRING,  1, 1,       (void *)&OmegaFile,  &nOmegaFiles,  "list of (angular) frequencies"},
     {"OmegaKBFile",    PA_STRING,  1, 1,       (void *)&OmegaKBFile, 0,             "list of (Omega, kx, ky) points"},
/**/
     {"FileBase",       PA_STRING,  1, 1,       (void *)&FileBase,   0,             "prefix for names of .log, .flux, .out files"},
/**/     
     {"PlotFlux",       PA_BOOL,    0, 1,       (void *)&PlotFlux,   0,             "write spatially-resolved flux data"},
     {"OmitSelfTerms",  PA_BOOL,    0, 1,       (void *)&OmitSelfTerms,   0,             "write spatially-resolved flux data"},
/**/
     {"EMTPFT",         PA_BOOL,    0, 1,       (void *)&EMTPFT,     0,             "compute SIFlux using EMT method"},
     {"OPFT",           PA_BOOL,    0, 1,       (void *)&OPFT,       0,             "compute SIFlux using overlap method"},
     {"DSIPFT",         PA_BOOL,    0, 1,       (void *)&DSIPFT,     0,             "compute SIFlux using displaced surface integral method"},
     {"DSIRadius",      PA_DOUBLE,  1, 1,       (void *)&DSIRadius,  0,             "bounding-sphere radius for DSIPFT"},
     {"DSIPoints",      PA_INT,     1, 1,       (void *)&DSIPoints,  0,             "number of cubature points for DSIPFT over sphere (6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810)"},
     {"DSIPoints2",     PA_INT,     1, 1,       (void *)&DSIPoints2, 0,             "number of cubature points for DSIPFT second opinion"},
     {"DSIMesh",        PA_STRING,  1, 1,       (void *)&DSIMesh,    0,             "bounding surface .msh file for DSIPFT"},
     {"DSIOmegaFile",   PA_STRING,  1, 1,       (void *)&DSIOmegaFile,  0,          "list of frequencies at which to perform DSI calculations"},
/**/
     {"Cache",          PA_STRING,  1, 1,       (void *)&Cache,      0,             "read/write cache"},
     {"ReadCache",      PA_STRING,  1, MAXCACHE,(void *)ReadCache,   &nReadCache,   "read cache"},
     {"WriteCache",     PA_STRING,  1, 1,       (void *)&WriteCache, 0,             "write cache"},
/**/     
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  if (GeoFile==0)
   OSUsage(argv[0], VERSION, OSArray, "--geometry option is mandatory");
  if (!FileBase)
   FileBase=vstrdup(GetFileBase(GeoFile));

  if ( Cache!=0 && WriteCache!=0 )
   ErrExit("--cache and --writecache options are mutually exclusive");

  if (DSIPFT && DSIPoints==0)
   DSIPoints=230;

  /*******************************************************************/
  /* determine which PFT methods were requested       ****************/
  /*******************************************************************/
  int NumPFTMethods=0, PFTMethods[MAXPFTMETHODS];
  if(EMTPFT)
   PFTMethods[NumPFTMethods++] = SCUFF_PFT_EMT;
  if(OPFT)
   PFTMethods[NumPFTMethods++] = SCUFF_PFT_OVERLAP;
  if(DSIMesh)
   PFTMethods[NumPFTMethods++] = -SCUFF_PFT_DSI;
  if(DSIPoints)
   PFTMethods[NumPFTMethods++] = DSIPoints;
  if(DSIPoints2)
   PFTMethods[NumPFTMethods++] = DSIPoints2;

  if (NumPFTMethods==0 && EPFile==0)
   PFTMethods[NumPFTMethods++] = SCUFF_PFT_EMT;

  if (DSIMesh) 
   PlotFlux=true;

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  HMatrix *OmegaKBPoints=0;
  if (OmegaKBFile)
   { OmegaKBPoints = new HMatrix(OmegaKBFile,LHM_TEXT,"--ncol 3");
     if (OmegaKBPoints->ErrMsg)
      ErrExit(OmegaKBPoints->ErrMsg);
   };

  /*******************************************************************/
  /* process frequency-related options to construct a list of        */
  /* frequencies at which to run simulations                         */
  /*******************************************************************/
  HVector *OmegaPoints=0, *OmegaPoints0;
  int NumFreqs=0;
  if (nOmegaFiles==1) // first process --OmegaFile option if present
   { 
     OmegaPoints=new HVector(OmegaFile,LHM_TEXT);
     if (OmegaPoints->ErrMsg)
      ErrExit(OmegaPoints->ErrMsg);
     NumFreqs=OmegaPoints->N;
     Log("Read %i frequencies from file %s.",NumFreqs,OmegaFile);
   };

  // now add any individually specified --Omega options
  if (nOmegaVals>0)
   { 
     NumFreqs += nOmegaVals;
     OmegaPoints0=OmegaPoints;
     OmegaPoints=new HVector(NumFreqs, LHM_COMPLEX);
     int nFreq=0;
     if (OmegaPoints0)
      { for(nFreq=0; nFreq<OmegaPoints0->N; nFreq++)
         OmegaPoints->SetEntry(nFreq, OmegaPoints0->GetEntry(nFreq));
        delete OmegaPoints0;
      };
     for(int nOV=0; nOV<nOmegaVals; nOV++)
      OmegaPoints->SetEntry(nFreq+nOV, OmegaVals[nOV]);
     Log("Read %i frequencies from command line.",nOmegaVals);
   };

  /*******************************************************************/
  /* create the SNEQData structure that contains all the info needed */
  /* to evaluate the neq transfer at a single frequency              */
  /*******************************************************************/
  SNEQData *SNEQD=CreateSNEQData(GeoFile, TransFile,
                                 PFTMethods, NumPFTMethods,
                                 EPFile, FileBase);
  RWGGeometry *G=SNEQD->G;
  SNEQD->PlotFlux                = PlotFlux;
  SNEQD->OmitSelfTerms           = OmitSelfTerms;
  SNEQD->PFTOpts.DSIMesh         = DSIMesh;
  SNEQD->PFTOpts.DSIRadius       = DSIRadius;
  SNEQD->PFTOpts.DSIFarField     = DSIFarField;
  SNEQD->DSIOmegaPoints          = DSIOmegaFile ? new HVector(DSIOmegaFile) : 0;

  if (OmegaKBPoints && !G->LBasis)
   ErrExit("--OmegaKBPoints may only be used with extended geometries");
  else if (G->LBasis && !OmegaKBPoints==0)
   ErrExit("--OmegaKBPoints is required for extended geometries");

  SNEQD->SourceOnly=-1;
  if (SourceObject)
   { G->GetSurfaceByLabel(SourceObject, &(SNEQD->SourceOnly) );
     if (SNEQD->SourceOnly == -1)
      ErrExit("geometry contains no object with label %s",SourceObject);
     Log("Computing only quantities sourced by object %s.",SourceObject);
   };
  SNEQD->DestOnly=-1;
  if (DestObject)
   { G->GetSurfaceByLabel(DestObject, &(SNEQD->DestOnly) );
     if (SNEQD->DestOnly == -1)
      ErrExit("geometry contains no object with label %s",DestObject);
     Log("Computing only PFT quantities for object %s.",DestObject);
   };
         
  /*******************************************************************/
  /* preload the scuff cache with any cache preload files the user   */
  /* may have specified                                              */
  /*******************************************************************/
  for (int nrc=0; nrc<nReadCache; nrc++)
   PreloadCache( ReadCache[nrc] );
  if (Cache)
   PreloadCache( Cache );

  if (Cache) WriteCache=Cache;
  SNEQD->WriteCache = WriteCache;

  /*******************************************************************/
  /* now switch off based on the requested frequency behavior to     */
  /* perform the actual calculations.                                */
  /*******************************************************************/
  if (OmegaKBPoints)
   { for (int nok=0; nok<OmegaKBPoints->NR; nok++)
      {  
        cdouble Omega; 
        double kBloch[2];
        Omega     = OmegaKBPoints->GetEntryD(nok, 0);
        kBloch[0] = OmegaKBPoints->GetEntryD(nok, 1);
        kBloch[1] = OmegaKBPoints->GetEntryD(nok, 2);
        WriteFlux(SNEQD, Omega, kBloch);
      };
   }
  else
   for (int nFreq=0; nFreq<NumFreqs; nFreq++)
    WriteFlux(SNEQD, OmegaPoints->GetEntry(nFreq));

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  printf("Thank you for your support.\n");

}
