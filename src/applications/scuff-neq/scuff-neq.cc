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
#include <BZIntegration.h>
#include "scuff-neq.h"
#include <libhrutil.h>

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define MAXFREQ  10    // max number of frequencies 
#define MAXCACHE 10    // max number of cache files for preload

#define MAXDSI   10

#define MAXSTR   1000

#define II cdouble(0.0,1.0)

/***************************************************************/
/***************************************************************/
/***************************************************************/
DSIPFTData *CreateDSIPFTData(RWGGeometry *G, char *DSIMesh,
                             char *SurfaceLabel=0,
                             char *DSIRadius=0, char *DSIPoints=0)
{ 
  DSIPFTData *Data=(DSIPFTData *)mallocEC(sizeof(DSIPFTData));
  Data->DSIMesh=DSIMesh;
  Data->DSIRadius=0.0;
  Data->DSIPoints=0;

  if (SurfaceLabel)
   { Data->TrackingSurface = G->GetSurfaceByLabel(SurfaceLabel);
     if ( !(Data->TrackingSurface))
      ErrExit("geometry has no tracking surface %s",SurfaceLabel);
   };
  if (DSIRadius) sscanf(DSIRadius,"%le",&(Data->DSIRadius));
  if (DSIPoints) sscanf(DSIPoints,"%i",&(Data->DSIPoints));

  if (!SurfaceLabel)
   { Data->Label=vstrdup("%s",GetFileBase(DSIMesh));
     Log("Computing DSIPFT with fixed mesh %s.",DSIMesh);
   }
  else if (DSIMesh)
   { Data->Label=vstrdup("%s.%s",SurfaceLabel,GetFileBase(DSIMesh));
     Log("Computing DSIPFT for %s (mesh %s)",SurfaceLabel,DSIMesh);
   }
  else 
   { Data->Label=vstrdup("%s.R%g.N%i",SurfaceLabel,Data->DSIRadius,Data->DSIPoints);
     Log("Computing DSIPFT for %s (radius %g, %i points)",SurfaceLabel,Data->DSIRadius,Data->DSIPoints);
   }

  return Data;
}

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
  char *SourceRegion=0;

  /*--------------------------------------------------------------*/
  char *EPFile=0;

  /*--------------------------------------------------------------*/
  cdouble OmegaVals[MAXFREQ];        int nOmegaVals;
  char *OmegaFile;                   int nOmegaFiles;
  char *OmegaKBFile=0;

  /*--------------------------------------------------------------*/
  bool EMTPFT      = false;

  /*--------------------------------------------------------------*/
  char *FDSIMeshes[MAXDSI];    int nFDSIMeshes=0;
  char *TDSIMeshes[2*MAXDSI];  int nTDSIMeshes=0;
  char *TDSISpheres[3*MAXDSI]; int nTDSISpheres=0;
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
     {"Geometry",       PA_STRING,  1, 1,       (void *)&GeoFile,      0,             "geometry file"},
     {"TransFile",      PA_STRING,  1, 1,       (void *)&TransFile,    0,             "list of geometrical transformation"},
/**/
     {"SourceRegion",   PA_STRING,  1, 1,       (void *)&SourceRegion, 0,             "compute only quantities due to sources in region xx"},
/**/     
     {"EPFile",         PA_STRING,  1, 1,       (void *)&EPFile,       0,             "list of evaluation points for spatially-resolved flux"},
/**/     
     {"Omega",          PA_CDOUBLE, 1, MAXFREQ, (void *)OmegaVals,     &nOmegaVals,   "(angular) frequency"},
     {"OmegaFile",      PA_STRING,  1, 1,       (void *)&OmegaFile,    &nOmegaFiles,  "list of (angular) frequencies"},
     {"OmegaKBFile",    PA_STRING,  1, 1,       (void *)&OmegaKBFile,  0,             "list of (Omega, kx, ky) points"},
/**/
     {"FileBase",       PA_STRING,  1, 1,       (void *)&FileBase,     0,             "prefix for names of .log, .flux, .out files"},
/**/     
     {"PlotFlux",       PA_BOOL,    0, 1,       (void *)&PlotFlux,     0,             "write spatially-resolved flux data"},
/**/
     {"EMTPFT",         PA_BOOL,    0, 1,       (void *)&EMTPFT,       0,             "compute SIFlux using EMT method"},
/**/
     {"FixedDSIMesh",      PA_STRING, 1, MAXDSI,  (void *)FDSIMeshes,  &nFDSIMeshes,  "fixed DSIPFT mesh"},
     {"TrackingDSIMesh",   PA_STRING, 2, MAXDSI,  (void *)TDSIMeshes,  &nTDSIMeshes,  "tracking DSIPFT object, mesh (e.g. --TrackingDSIMesh MyObject MyMesh.msh)"},
     {"TrackingDSISphere", PA_STRING, 3, MAXDSI,  (void *)TDSISpheres, &nTDSISpheres, "tracking DSIPFT object, radius, #points (e.g. --TrackingDSISphere MyObject 3.0 302) (#points=6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810)"},
     {"DSIOmegaFile",   PA_STRING,  1, 1,       (void *)&DSIOmegaFile,  0,             "list of frequencies at which to perform DSI calculations"},
/**/
     {"Cache",          PA_STRING,  1, 1,       (void *)&Cache,      0,                "read/write cache"},
     {"ReadCache",      PA_STRING,  1, MAXCACHE,(void *)ReadCache,   &nReadCache,      "read cache"},
     {"WriteCache",     PA_STRING,  1, 1,       (void *)&WriteCache, 0,                "write cache"},
/**/     
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  if (GeoFile==0)
   OSUsage(argv[0], OSArray, "--geometry option is mandatory");
  if (!FileBase)
   FileBase=vstrdup(GetFileBase(GeoFile));

  if ( Cache!=0 && WriteCache!=0 )
   ErrExit("--cache and --writecache options are mutually exclusive");

  RWGGeometry *G=new RWGGeometry(GeoFile);

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
  /* process DSIPFT-related options **********************************/
  /*******************************************************************/
  if (EMTPFT)
   Log("Computing EMTPFT");
  DSIPFTDataList DSIPFTs;
  for(int n=0; n<nFDSIMeshes; n++)
   DSIPFTs.push_back(CreateDSIPFTData(G,FDSIMeshes[n]));
  for(int n=0; n<nTDSIMeshes; n++)
   DSIPFTs.push_back(CreateDSIPFTData(G,TDSIMeshes[2*n+1],TDSIMeshes[2*n+0]));
  for(int n=0; n<nTDSISpheres; n++)
   DSIPFTs.push_back(CreateDSIPFTData(G,0,TDSISpheres[3*n+0],
                                          TDSISpheres[3*n+1],
                                          TDSISpheres[3*n+2]));

  /*******************************************************************/
  /* create the SNEQData structure that contains all the info needed */
  /* to evaluate the neq transfer at a single frequency              */
  /*******************************************************************/
  SNEQData *SNEQD=CreateSNEQData(G, TransFile, EMTPFT, DSIPFTs, EPFile, FileBase);
  SNEQD->PlotFlux                = PlotFlux;
  SNEQD->OmitSelfTerms           = OmitSelfTerms;
  SNEQD->DSIOmegaPoints          = DSIOmegaFile ? new HVector(DSIOmegaFile) : 0;

  if (OmegaKBPoints && !G->LBasis)
   ErrExit("--OmegaKBPoints may only be used with extended geometries");
  else if (G->LBasis && !OmegaKBPoints==0)
   ErrExit("--OmegaKBPoints is required for extended geometries");

  SNEQD->SourceRegion=-1;
  if (SourceRegion)
   { SNEQD->SourceRegion=G->GetRegionByLabel(SourceRegion);
     if (SNEQD->SourceRegion== -1)
      ErrExit("geometry contains no region with label %s",SourceRegion);
     Log("Computing only quantities sourced by region %s.",SourceRegion);
   }
         
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
