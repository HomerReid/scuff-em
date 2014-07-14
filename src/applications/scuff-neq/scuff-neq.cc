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

#include "scuff-neq.h"
#include <libhrutil.h>

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define MAXFREQ  10    // max number of frequencies 
#define MAXCACHE 10    // max number of cache files for preload

#define MAXTEMPS 10    // max number of objects for which temperatures may be set

#define MAXSTR   1000

#define II cdouble(0.0,1.0)

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  InstallHRSignalHandler();

  /***************************************************************/
  /* process options *********************************************/
  /***************************************************************/
  char *GeoFile=0;
  char *TransFile=0;

  /*--------------------------------------------------------------*/
  int Power=0;
  int XForce=0;
  int YForce=0;
  int ZForce=0;
  int XTorque=0;
  int YTorque=0;
  int ZTorque=0;

  /*--------------------------------------------------------------*/
  cdouble OmegaVals[MAXFREQ];        int nOmegaVals;
  char *OmegaFile;                   int nOmegaFiles;
  char *OmegaKFile=0;
  double OmegaMin=0.00;              int nOmegaMin;
  double OmegaMax=-1.0;              int nOmegaMax;

  /*--------------------------------------------------------------*/
  char *TempStrings[2*MAXTEMPS];     int nTempStrings;

  /*--------------------------------------------------------------*/
  double AbsTol=0.0;
  double RelTol=5.0e-2;
  int Intervals=25;

  /*--------------------------------------------------------------*/
  char *FileBase=0;

  /*--------------------------------------------------------------*/
  int PlotFlux=0;

  /*--------------------------------------------------------------*/
  char *Cache=0;
  char *ReadCache[MAXCACHE];         int nReadCache;
  char *WriteCache=0;

  /*--------------------------------------------------------------*/
  bool SymGPower=false;

  /*--------------------------------------------------------------*/
  bool UseExistingData=false;
  bool SubtractSelfTerms=false;
  bool Visualize=false;
   
  /*--------------------------------------------------------------*/
  double SIRadius    = 100.0;
  double SINumPoints = 31;

  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { 
     {"Geometry",       PA_STRING,  1, 1,       (void *)&GeoFile,    0,             "geometry file"},
     {"TransFile",      PA_STRING,  1, 1,       (void *)&TransFile,  0,             "list of geometrical transformation"},
/**/     
     {"Power",          PA_BOOL,    0, 1,       (void *)&Power,      0,             "compute power transfer"},
     {"XForce",         PA_BOOL,    0, 1,       (void *)&XForce,     0,             "compute X-force"},
     {"YForce",         PA_BOOL,    0, 1,       (void *)&YForce,     0,             "compute Y-force"},
     {"ZForce",         PA_BOOL,    0, 1,       (void *)&ZForce,     0,             "compute Z-force"},
     {"XTorque",        PA_BOOL,    0, 1,       (void *)&XTorque,    0,             "compute X-torque"},
     {"YTorque",        PA_BOOL,    0, 1,       (void *)&YTorque,    0,             "compute Y-torque"},
     {"ZTorque",        PA_BOOL,    0, 1,       (void *)&ZTorque,    0,             "compute Z-torque"},
/**/     
     {"Omega",          PA_CDOUBLE, 1, MAXFREQ, (void *)OmegaVals,   &nOmegaVals,   "(angular) frequency"},
     {"OmegaFile",      PA_STRING,  1, 1,       (void *)&OmegaFile,  &nOmegaFiles,  "list of (angular) frequencies"},
     {"OmegaMin",       PA_DOUBLE,  1, 1,       (void *)&OmegaMin,   &nOmegaMin,    "lower integration limit"},
     {"OmegaMax",       PA_DOUBLE,  1, 1,       (void *)&OmegaMax,   &nOmegaMax,    "upper integration limit"},
/**/     
     {"OmegaKFile",     PA_STRING,  1, 1,       (void *)&OmegaKFile, 0,             "list of (Omega, kx, ky) points"},
/**/     
     {"Temperature",    PA_STRING,  2, MAXTEMPS, (void *)TempStrings, &nTempStrings,  "set object xx temperature to xx"},
/**/     
     {"FileBase",       PA_STRING,  1, 1,       (void *)&FileBase,   0,             "prefix for names of .log, .flux, .out files"},
/**/     
     {"PlotFlux",       PA_BOOL,    0, 1,       (void *)&PlotFlux,   0,             "write spatially-resolved flux data"},
/**/     
     {"AbsTol",         PA_DOUBLE,  1, 1,       (void *)&AbsTol,     0,             "absolute tolerance for frequency quadrature"},
     {"RelTol",         PA_DOUBLE,  1, 1,       (void *)&RelTol,     0,             "relative tolerance for frequency quadrature"},
     {"Intervals",      PA_INT,     1, 1,       (void *)&Intervals,  0,             "number of intervals for frequency quadrature"},
/**/
     {"SymGPower",      PA_BOOL,    0, 1,       (void *)&SymGPower,  0,             "use Sym(G) instead of overlap matrix for power computation"},
/**/
     {"SIRadius",       PA_DOUBLE,  1, 1,       (void *)&SIRadius,   0,             "bounding-sphere radius for SIPFT"},
     {"SINumPoints",    PA_INT,     1, 1,       (void *)&SINumPoints,0,             "number of quadrature points for SIPFT"},
/**/
     {"UseExistingData", PA_BOOL,   0, 1,       (void *)&UseExistingData, 0,        "read existing data from .flux files"},
     {"Visualize",      PA_BOOL,    0, 1,       (void *)&Visualize,  0,             "visualize flux profiles"},
/**/
     {"Cache",          PA_STRING,  1, 1,       (void *)&Cache,      0,             "read/write cache"},
     {"ReadCache",      PA_STRING,  1, MAXCACHE,(void *)ReadCache,   &nReadCache,   "read cache"},
     {"WriteCache",     PA_STRING,  1, 1,       (void *)&WriteCache, 0,             "write cache"},
/**/     
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  if (GeoFile==0)
   OSUsage(argv[0], OSArray, "--geometry option is mandatory");

  if ( Cache!=0 && WriteCache!=0 )
   ErrExit("--cache and --writecache options are mutually exclusive");

  if (FileBase)
   SetLogFileName("%s.log",FileBase);
  else
   SetLogFileName("scuff-neq.log");

  Log("scuff-neq running on %s",GetHostName());

  /*******************************************************************/
  /* determine which output quantities were requested ****************/
  /*******************************************************************/
  int QuantityFlags=0;
  if (Power)  QuantityFlags|=QFLAG_POWER;
  if (XForce) QuantityFlags|=QFLAG_XFORCE;
  if (YForce) QuantityFlags|=QFLAG_YFORCE;
  if (ZForce) QuantityFlags|=QFLAG_ZFORCE;
  if (XTorque) QuantityFlags|=QFLAG_XTORQUE;
  if (YTorque) QuantityFlags|=QFLAG_YTORQUE;
  if (ZTorque) QuantityFlags|=QFLAG_ZTORQUE;
  if (QuantityFlags==0)
   ErrExit("you must specify at least one quantity to compute");

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  HMatrix *OmegaKPoints=0; 
  if (OmegaKFile)
   { OmegaKPoints = new HMatrix(OmegaKFile,LHM_TEXT,"--ncol 3");
     if (OmegaKPoints->ErrMsg)
      ErrExit(OmegaKPoints->ErrMsg);
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
  /* check that the user didn't simultaneously ask for a discrete set*/
  /* of frequencies and a frequency range over which to integrate;   */
  /* if a range was specified check that it makes sense              */
  /*******************************************************************/
  if ( OmegaKPoints || OmegaPoints ) 
   { if ( nOmegaMin>0 || nOmegaMax>0 )
      ErrExit("--OmegaMin/--OmegaMax options may not be used with --Omega/--OmegaFile");
     if ( nTempStrings>0 )
      ErrExit("--Temperature option may not be used with --Omega/--OmegaFile");
     Log("Computing spectral density at %i frequencies.",NumFreqs);
   }
  else
   { 
     if ( nOmegaMin==1 && OmegaMin<0.0 )
      ErrExit("invalid value specified for --OmegaMin");
     if ( nOmegaMax==1 && OmegaMax<OmegaMin )
      ErrExit("invalid value specified for --OmegaMax");

     if ( OmegaMax==-1.0 )
      Log("Integrating over range Omega=(%g,infinity).",OmegaMin);
     else
      Log("Integrating over range Omega=(%g,%g).",OmegaMin,OmegaMax);
   };

  /*******************************************************************/
  /* create the SNEQData structure that contains all the info needed*/
  /* to evaluate the neq transfer at a single frequency              */
  /*******************************************************************/
  SNEQData *SNEQD=CreateSNEQData(GeoFile, TransFile, QuantityFlags, 
                                 PlotFlux, FileBase, SymGPower);
  RWGGeometry *G=SNEQD->G;
  SNEQD->UseExistingData   = UseExistingData;
  SNEQD->SubtractSelfTerms = SubtractSelfTerms;
  SNEQD->Visualize         = Visualize;
  SNEQD->SIRadius          = SIRadius;  
  SNEQD->SINumPoints       = SINumPoints;

  if (OmegaKPoints && G->NumLatticeBasisVectors==0)
   ErrExit("--OmegaKPoints may only be used with extended geometries");
  else if (G->NumLatticeBasisVectors!=0 && OmegaKPoints==0)
   ErrExit("--OmegaKPoints is required for extended geometries");

  /*******************************************************************/
  /* process any temperature specifications **************************/
  /*******************************************************************/
  double TEnvironment=0.0;
  double *TSurfaces=(double *)malloc(G->NumSurfaces*sizeof(double));
  memset(TSurfaces, 0, G->NumSurfaces*sizeof(double));
  for(int nts=0; nts<nTempStrings; nts++)
   { 
     double TTemp;
     int WhichSurface;
     if ( 1!=sscanf(TempStrings[2*nts+1],"%le",&TTemp) )
      ErrExit("invalid temperature (%s) passed for --temperature option",TempStrings[2*nts+1]);

     if (    !strcasecmp(TempStrings[2*nts],"MEDIUM")
          || !strcasecmp(TempStrings[2*nts],"ENVIRONMENT")
        )
      { TEnvironment=TTemp;
        Log("Setting environment temperature to %g kelvin.\n",TTemp);
        printf("Setting environment temperature to %g kelvin.\n",TTemp);
      }
     else if ( G->GetSurfaceByLabel(TempStrings[2*nts],&WhichSurface) )
      { TSurfaces[WhichSurface]=TTemp;
        Log("Setting temperature of object %s to %g kelvin.\n",TempStrings[2*nts],TTemp);
        printf("Setting temperature of object %s to %g kelvin.\n",TempStrings[2*nts],TTemp);
      }
     else 
      ErrExit("unknown surface/region %s in --temperature specification",TempStrings[2*nts]);

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
  /* perform the actual calculations                                 */
  /*******************************************************************/
  int OutputVectorLength = SNEQD->NumTransformations * G->NumSurfaces * G->NumSurfaces * SNEQD->NQ;
  double *I = new double[ OutputVectorLength ];
  if (OmegaKPoints)
   { for (int nok=0; nok<OmegaKPoints->NR; nok++)
      {  
        cdouble Omega; 
        double kBloch[2];
        Omega     = OmegaKPoints->GetEntryD(nok, 0);
        kBloch[0] = OmegaKPoints->GetEntryD(nok, 1);
        kBloch[1] = OmegaKPoints->GetEntryD(nok, 2);
        GetFlux(SNEQD, Omega, kBloch, I);
      };
   }
  else if (NumFreqs>0)
   { for (int nFreq=0; nFreq<NumFreqs; nFreq++)
      GetFlux(SNEQD, OmegaPoints->GetEntry(nFreq), I);
   }
  else
   { 
      double *E = new double[ OutputVectorLength ];
      EvaluateFrequencyIntegral2(SNEQD, OmegaMin, OmegaMax, 
                                 TSurfaces, TEnvironment, 
                                 Intervals, I, E);
      delete[] E;
   };
  delete[] I;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  printf("Thank you for your support.\n");

}
