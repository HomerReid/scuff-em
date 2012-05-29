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

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define MAXFREQ  10    // max number of frequencies 
#define MAXCACHE 10    // max number of cache files for preload

#define MAXTEMPS 10    // max number of objects for which temperatures may be set

#define MAXSTR   1000

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

  int Power=0;
  int xForce=0;
  int yForce=0;
  int zForce=0;

  cdouble OmegaVals[MAXFREQ];        int nOmegaVals;
  char *OmegaFile;                   int nOmegaFiles;
  cdouble OmegaMin;                  int nOmegaMin;
  cdouble OmegaMax;                  int nOmegaMax;

  char *TempStrings[2*MAXTEMPS];     int nTempStrings;

  char *OutputFile=0;
  char *ByOmegaFile=0;
  char *LogFile=0;

  int PlotFlux=0;

  char *Cache=0;
  char *ReadCache[MAXCACHE];         int nReadCache;
  char *WriteCache=0;

  int nThread=0;

  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { 
     {"Geometry",       PA_STRING,  1, 1,       (void *)&GeoFile,    0,             "geometry file"},
     {"TransFile",      PA_STRING,  1, 1,       (void *)&TransFile,  0,             "list of geometrical transformation"},
/**/     
     {"Power",          PA_BOOL,    0, 1,       (void *)&Power,      0,             "compute power transfer"},
     {"xForce",         PA_BOOL,    0, 1,       (void *)&xForce,     0,             "compute X-momentum transfer"},
     {"yForce",         PA_BOOL,    0, 1,       (void *)&yForce,     0,             "compute Y-momentum transfer"},
     {"zForce",         PA_BOOL,    0, 1,       (void *)&zForce,     0,             "compute Z-momentum transfer"},
/**/     
     {"Omega",          PA_CDOUBLE, 1, MAXFREQ, (void *)OmegaVals,   &nOmegaVals,   "(angular) frequency"},
     {"OmegaFile",      PA_STRING,  1, 1,       (void *)&OmegaFile,  &nOmegaFiles,  "list of (angular) frequencies"},
     {"OmegaMin",       PA_CDOUBLE, 1, 1,       (void *)&OmegaMin,   &nOmegaMin,    "lower integration limit"},
     {"OmegaMax",       PA_CDOUBLE, 1, 1,       (void *)&OmegaMax,   &nOmegaMax,    "upper integration limit"},
/**/     
     {"Temperature",    PA_STRING,  2, MAXTEMPS, (void *)TempStrings, &nTempStrings,  "set object xx temperature to xx"},
/**/     
     {"OutputFile",     PA_STRING,  1, 1,       (void *)&OutputFile, 0,             "name of frequency-integrated output file"},
     {"ByOmegaFile",    PA_STRING,  1, 1,       (void *)&ByOmegaFile,0,             "name of frequency-resolved output file"},
     {"LogFile",        PA_STRING,  1, 1,       (void *)&LogFile,    0,             "name of log file"},
/**/     
     {"PlotFlux",       PA_BOOL,    0, 1,       (void *)&PlotFlux,   0,             "write spatially-resolved flux data"},
/**/     
     {"Cache",          PA_STRING,  1, 1,       (void *)&Cache,      0,             "read/write cache"},
     {"ReadCache",      PA_STRING,  1, MAXCACHE,(void *)ReadCache,   &nReadCache,   "read cache"},
     {"WriteCache",     PA_STRING,  1, 1,       (void *)&WriteCache, 0,             "write cache"},
/**/     
     {"nThread",        PA_INT,     1, 1,       (void *)&nThread,    0,             "number of CPU threads to use"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  if (GeoFile==0)
   OSUsage(argv[0], OSArray, "--geometry option is mandatory");

  if ( Cache!=0 && WriteCache!=0 )
   ErrExit("--cache and --writecache options are mutually exclusive");

  if (PlotFlux && ByOmegaFile!=0 )
   ErrExit("--PlotFlux and --ByOmegaFile options are mutually exclusive");

  if (LogFile)
   SetLogFileName(LogFile);
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
  if (QuantityFlags==0)
   ErrExit("you must specify at least one quantity to compute");

  /*******************************************************************/
  /* process frequency-related options to construct a list of        */
  /* frequencies at which to run simulations                         */
  /*******************************************************************/
  HVector *OmegaList=0, *OmegaList0;
  int nFreq, nOV, NumFreqs=0;
  if (nOmegaFiles==1) // first process --OmegaFile option if present
   { 
     OmegaList=new HVector(OmegaFile,LHM_TEXT);
     if (OmegaList->ErrMsg)
      ErrExit(OmegaList->ErrMsg);
     NumFreqs=OmegaList->N;
     Log("Read %i frequencies from file %s.",NumFreqs,OmegaFile);
   };

  // now add any individually specified --Omega options
  if (nOmegaVals>0)
   { 
     NumFreqs += nOmegaVals;
     OmegaList0=OmegaList;
     OmegaList=new HVector(NumFreqs, LHM_COMPLEX);
     nFreq=0;
     if (OmegaList0)
      { for(nFreq=0; nFreq<OmegaList0->N; nFreq++)
         OmegaList->SetEntry(nFreq, OmegaList0->GetEntry(nFreq));
        delete OmegaList0;
      };
     for(nOV=0; nOV<nOmegaVals; nOV++)
      OmegaList->SetEntry(nFreq+nOV, OmegaVals[nOV]);
     Log("Read %i frequencies from command line.",nOmegaVals);
   };

  // check that the user didn't simultaneously ask for a discrete
  // list of frequencies and a frequency range over which to integrate;
  // if a range was specified check that it makes sense 
  if ( NumFreqs>0 ) 
   { if ( nOmegaMin>0 || nOmegaMax>0 )
      ErrExit("--OmegaMin/--OmegaMax options may not be used with --Omega/--OmegaFile");
     if ( nTempStrings>0 )
      ErrExit("--Temperature option may not be used with --Omega/--OmegaFile");
     Log("Computing spectral density at %i frequencies.",NumFreqs);
   }
  else if (NumFreqs==0)
   { 
     // if --OmegaMin and/or --OmegaMax values were specified,
     // check that they make sense
     if ( nOmegaMin==1 && (real(OmegaMin)<0.0 || imag(OmegaMin)!=0.0) )
      ErrExit("invalid value specified for --OmegaMin");
     if ( nOmegaMax==1 && (real(OmegaMax)<real(OmegaMin) || imag(OmegaMax)!=0.0 ) )
      ErrExit("invalid value specified for --OmegaMax");

     if ( OmegaMax==-1.0 )
      Log("Integrating over range Omega=(%g,infinity).",real(OmegaMin));
     else
      Log("Integrating over range Omega=(%g,%g).",real(OmegaMin),real(OmegaMax));
   };

  /*******************************************************************/
  /* create the SHNEQData structure that contains all the info needed*/
  /* to evaluate the neq transfer at a single frequency              */
  /*******************************************************************/
  SHNEQData *SHNEQD=CreateSHNEQData(GeoFile, TransFile, ByOmegaFile, QuantityFlags, PlotFlux);
  RWGGeometry *G=SHNEQD->G;

  /*******************************************************************/
  /* process any temperature specifications **************************/
  /*******************************************************************/
  double TEnvironment=0.0;
  double *TObjects=(double *)malloc(G->NumObjects*sizeof(double));
  memset(TObjects, 0, G->NumObjects*sizeof(double));
  if (nTempStrings)
   { 
     RWGObject *O;
     int WhichObject;
     double TTemp;

     for(int nts=0; nts<nTempStrings; nts++)
      { O=G->GetObjectByLabel(TempStrings[2*nts],&WhichObject);

        if(WhichObject==-2)
         ErrExit("unknown object (%s) passed for --temperature option",TempStrings[2*nts]);
        if ( 1!=sscanf(TempStrings[2*nts+1],"%le",&TTemp) )
         ErrExit("invalid temperature (%s) passed for --temperature option",TempStrings[2*nts+1]);

        if(WhichObject==-1)
         { TEnvironment=TTemp;
           Log("Setting environment temperature to %g kelvin.\n",TTemp);
           printf("Setting environment temperature to %g kelvin.\n",TTemp);
         }
        else
         { TObjects[WhichObject]=TTemp;
           Log("Setting temperature of object %s to %g kelvin.\n",TempStrings[2*nts],TTemp);
           printf("Setting temperature of object %s to %g kelvin.\n",TempStrings[2*nts],TTemp);
         };
      };
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
  SHNEQD->WriteCache = WriteCache;

  /*******************************************************************/
  /* now switch off based on the requested frequency behavior to     */
  /* perform the actual calculations                                 */
  /*******************************************************************/
  double *I = new double[SHNEQD->NTNQ];
  if (NumFreqs>0)
   { for (nFreq=0; nFreq<NumFreqs; nFreq++)
      GetFrequencyIntegrand(SHNEQD, OmegaList->GetEntry(nFreq), I);
   }
  else
   { // frequency integration not yet implemented 
     ErrExit("frequency integration is not yet implemented");
   };
  delete[] I;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  printf("Thank you for your support.\n");

}
