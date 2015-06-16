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
  bool PAbs=false;
  bool PRad=false;
  bool XForce=false;
  bool YForce=false;
  bool ZForce=false;
  bool XTorque=false;
  bool YTorque=false;
  bool ZTorque=false;

  /*--------------------------------------------------------------*/
  char *EPFile=0;

  /*--------------------------------------------------------------*/
  cdouble OmegaVals[MAXFREQ];        int nOmegaVals;
  char *OmegaFile;                   int nOmegaFiles;
  char *OmegaKFile=0;
  double OmegaMin=0.00;              int nOmegaMin;
  double OmegaMax=-1.0;              int nOmegaMax;
  char *OQString=0;   

  /*--------------------------------------------------------------*/
  char *TempStrings[2*MAXTEMPS];     int nTempStrings;

  /*--------------------------------------------------------------*/
  double AbsTol=1.0e-8;
  double RelTol=5.0e-2;
  int Intervals=25.0;                int nIntervals=0;
   
  /*--------------------------------------------------------------*/
  char *DSIMesh    = 0;
  double DSIRadius = 10.0;
  int DSIPoints    = 302;
  bool DSICCQ      = false;
  bool DSIFarField = false;

  /*--------------------------------------------------------------*/
  bool OmitSelfTerms = false;
  bool OmitZeroTemperatureFlux = false;
  bool ForceDSI      = false;

  /*--------------------------------------------------------------*/
  char *FileBase=0;

  /*--------------------------------------------------------------*/
  bool PlotFlux=false;
  int PlotRytovVectors=0;

  /*--------------------------------------------------------------*/
  char *Cache=0;
  char *ReadCache[MAXCACHE];         int nReadCache;
  char *WriteCache=0;

  /*--------------------------------------------------------------*/
  bool UseExistingData=false;

  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { 
     {"Geometry",       PA_STRING,  1, 1,       (void *)&GeoFile,    0,             "geometry file"},
     {"TransFile",      PA_STRING,  1, 1,       (void *)&TransFile,  0,             "list of geometrical transformation"},
/**/     
     {"Power",          PA_BOOL,    0, 1,       (void *)&PAbs,       0,             "compute power transfer"},
     {"PAbs",           PA_BOOL,    0, 1,       (void *)&PAbs,       0,             "(synonym for --power)"},
     {"PRad",           PA_BOOL,    0, 1,       (void *)&PRad,       0,             "compute radiated power"},
     {"XForce",         PA_BOOL,    0, 1,       (void *)&XForce,     0,             "compute X-force"},
     {"YForce",         PA_BOOL,    0, 1,       (void *)&YForce,     0,             "compute Y-force"},
     {"ZForce",         PA_BOOL,    0, 1,       (void *)&ZForce,     0,             "compute Z-force"},
     {"XTorque",        PA_BOOL,    0, 1,       (void *)&XTorque,    0,             "compute X-torque"},
     {"YTorque",        PA_BOOL,    0, 1,       (void *)&YTorque,    0,             "compute Y-torque"},
     {"ZTorque",        PA_BOOL,    0, 1,       (void *)&ZTorque,    0,             "compute Z-torque"},
/**/     
     {"EPFile",         PA_STRING,  1, 1,       (void *)&EPFile, 0,             "list of evaluation points for spatially-resolved flux"},
/**/     
     {"Omega",          PA_CDOUBLE, 1, MAXFREQ, (void *)OmegaVals,   &nOmegaVals,   "(angular) frequency"},
     {"OmegaFile",      PA_STRING,  1, 1,       (void *)&OmegaFile,  &nOmegaFiles,  "list of (angular) frequencies"},
     {"OmegaKFile",     PA_STRING,  1, 1,       (void *)&OmegaKFile, 0,             "list of (Omega, kx, ky) points"},
     {"OmegaMin",       PA_DOUBLE,  1, 1,       (void *)&OmegaMin,   &nOmegaMin,    "lower integration limit"},
     {"OmegaMax",       PA_DOUBLE,  1, 1,       (void *)&OmegaMax,   &nOmegaMax,    "upper integration limit"},
     {"OmegaQuadrature",PA_STRING,  1, 1,       (void *)&OQString,   0,             "quadrature method for omega integral"},
     {"AbsTol",         PA_DOUBLE,  1, 1,       (void *)&AbsTol,     0,             "absolute tolerance for adaptive/cliff frequency quadrature"},
     {"RelTol",         PA_DOUBLE,  1, 1,       (void *)&RelTol,     0,             "relative tolerance for adaptive/cliff frequency quadrature"},
     {"Intervals",      PA_INT,     1, 1,       (void *)&Intervals,  &nIntervals,   "number of intervals for TrapSimp frequency quadrature"},
/**/
     {"Temperature",    PA_STRING,  2, MAXTEMPS, (void *)TempStrings, &nTempStrings,  "set object xx temperature to xx"},
/**/     
     {"FileBase",       PA_STRING,  1, 1,       (void *)&FileBase,   0,             "prefix for names of .log, .flux, .out files"},
/**/     
     {"PlotFlux",       PA_BOOL,    0, 1,       (void *)&PlotFlux,   0,             "write spatially-resolved flux data"},
/**/
     {"DSIMesh",        PA_STRING,  1, 1,       (void *)&DSIMesh,    0,             "bounding surface .msh file for DSIPFT"},
     {"DSIRadius",      PA_DOUBLE,  1, 1,       (void *)&DSIRadius,  0,             "bounding-sphere radius for DSIPFT"},
     {"DSIPoints",      PA_INT,     1, 1,       (void *)&DSIPoints,  0,             "number of quadrature points for DSIPFT"},
     {"DSICCQ",         PA_BOOL,    0, 1,       (void *)&DSICCQ,    0,              "use Clenshaw-Curtis cubature for DSIPFT"},
     {"DSIFarField",    PA_BOOL,    0, 1,       (void *)&DSIFarField,   0,          "retain only far-field contributions to DSIPFT"},
/**/
     {"OmitSelfTerms",  PA_BOOL,    0, 1,       (void *)&OmitSelfTerms, 0,          "omit the calculation of self terms"},
     {"ForceDSI",       PA_BOOL,    0, 1,       (void *)&ForceDSI,   0,             "use DSIPFT instead of OPFT/EPPFT"},
/**/
     {"UseExistingData", PA_BOOL,   0, 1,       (void *)&UseExistingData, 0,        "read existing data from .flux files"},
     {"PlotRytovVectors",PA_INT,    1, 1,       (void *)&PlotRytovVectors, 0, "plot the first N Rytov surface-current vectors"},
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
   OSUsage(argv[0], OSArray, "--geometry option is mandatory");
  if (!FileBase)
   FileBase=vstrdup(GetFileBase(GeoFile));
  SetLogFileName("%s.log",FileBase);
  Log("scuff-neq running on %s",GetHostName());

  if ( Cache!=0 && WriteCache!=0 )
   ErrExit("--cache and --writecache options are mutually exclusive");

  /*******************************************************************/
  /* determine which output quantities were requested ****************/
  /*******************************************************************/
  int QuantityFlags=0;
  if (PAbs)   QuantityFlags|=QFLAG_PABS;
  if (PRad)   QuantityFlags|=QFLAG_PRAD;
  if (XForce) QuantityFlags|=QFLAG_XFORCE;
  if (YForce) QuantityFlags|=QFLAG_YFORCE;
  if (ZForce) QuantityFlags|=QFLAG_ZFORCE;
  if (XTorque) QuantityFlags|=QFLAG_XTORQUE;
  if (YTorque) QuantityFlags|=QFLAG_YTORQUE;
  if (ZTorque) QuantityFlags|=QFLAG_ZTORQUE;
  if (QuantityFlags==0 && EPFile==0)
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

  int OQMethod = QMETHOD_CLIFF;
  if ( OQString )
   { if (OmegaPoints)
      ErrExit("--OmegaQuadrature may not be used with --Omega or --OmegaFile");
     if ( !strcasecmp(OQString,"adaptive") )
      OQMethod = QMETHOD_ADAPTIVE;
     else if ( !strcasecmp(OQString,"cliff") )
      OQMethod = QMETHOD_CLIFF;
     else if ( !strcasecmp(OQString,"trapsimp") )
      OQMethod = QMETHOD_TRAPSIMP;
     else
      ErrExit("unknown --OmegaQuadrature method %s",OQString);    
   };
  if ( nIntervals!=0 && OQMethod!=QMETHOD_TRAPSIMP )
   ErrExit("--Intervals may only be used with --OmegaQuadrature TrapSimp");

  /*******************************************************************/
  /* create the SNEQData structure that contains all the info needed */
  /* to evaluate the neq transfer at a single frequency              */
  /*******************************************************************/
  SNEQData *SNEQD=CreateSNEQData(GeoFile, TransFile,
                                 TempStrings, nTempStrings,
                                 QuantityFlags, EPFile, FileBase);
  RWGGeometry *G=SNEQD->G;
  SNEQD->UseExistingData         = UseExistingData;
  SNEQD->PlotFlux                = PlotFlux;
  SNEQD->OmitSelfTerms           = OmitSelfTerms;
  SNEQD->OmitZeroTemperatureFlux = OmitZeroTemperatureFlux;
  SNEQD->ForceDSI                = ForceDSI;
  SNEQD->AbsTol                  = AbsTol;
  SNEQD->RelTol                  = RelTol;
  SNEQD->PFTOpts.DSIMesh         = DSIMesh;
  SNEQD->PFTOpts.DSIRadius       = DSIRadius;
  SNEQD->PFTOpts.DSIPoints       = DSIPoints;
  SNEQD->PFTOpts.DSIFarField     = DSIFarField;
  SNEQD->PlotRytovVectors        = PlotRytovVectors;

  if (SNEQD->NumSIQs==0 && SNEQD->NumSRQs==0)
   ErrExit("you must specify at least one quantity to compute");

  if (OmegaKPoints && G->LDim==0)
   ErrExit("--OmegaKPoints may only be used with extended geometries");
  else if (G->LDim!=0 && OmegaKPoints==0)
   ErrExit("--OmegaKPoints is required for extended geometries");
         
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
  int OutputVectorLength = SNEQD->NumSIQs + SNEQD->NumSRQs;
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
   { 
     for (int nFreq=0; nFreq<NumFreqs; nFreq++)
      GetFlux(SNEQD, OmegaPoints->GetEntry(nFreq), I);
   }
  else
   { 
      SNEQD->OmitZeroTemperatureFlux=true;

      double *E = new double[ OutputVectorLength ];

      switch(OQMethod) 
       { case QMETHOD_ADAPTIVE:
          GetOmegaIntegral_Adaptive(SNEQD, OmegaMin, OmegaMax, I, E);

         case QMETHOD_TRAPSIMP:
          GetOmegaIntegral_TrapSimp(SNEQD, OmegaMin, OmegaMax, Intervals, I, E);

         case QMETHOD_CLIFF:
          GetOmegaIntegral_Cliff(SNEQD, OmegaMin, OmegaMax, I, E);
       };

      if (SNEQD->NumSIQs > 0) 
       WriteSIOutputFile(SNEQD, I, E);

      if (SNEQD->NumSRQs > 0) 
       WriteSROutputFile(SNEQD, I, E);

      delete[] E;
   };
  delete[] I;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  printf("Thank you for your support.\n");

}
