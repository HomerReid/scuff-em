/*
 * scuff-cas3D.cc  -- a standalone code within the scuff-EM suite for 
 *                 -- computing the Casimir energy, force, and/or
 *                 -- torque for a collection of interacting objects
 *
 * homer reid     -- 3/2007 -- 3/2012
 *
 * --------------------------------------------------------------
 *
 * this program has a large number of command-line options, which
 * subdivide into a number of categories as described below.
 *
 * in addition to the command line, options may also be specified
 * in an input file, piped into standard input with one option-value
 * pair per line; thus, if --option1 and --option2 are command-line 
 * options as described below, then you may create an input file 
 * (call it 'myOptions') with the content 
 *                  
 *   ...
 *   option1 value1
 *   option2 value2
 *   ...
 * 
 * and then running 
 * 
 *  scuff-cas3D < myOptions 
 * 
 * is equivalent to 
 * 
 *  scuff-cas3D --option1 value1 --option value2.
 * 
 * (if any options are specified both on standard input and 
 *  on the command line, the values given on the command line take
 *  precedence.)
 * 
 *       -------------------------------------------------
 * 
 * a. options describing the geometry and the (optional)
 *    list of geometrical transformations applied to it
 * 
 *     --geometry  MyGeometry.scuffgeo
 *     --transfile MyGeometry.trans
 * 
 *       -------------------------------------------------
 * 
 * b. options describing the quantities to compute
 * 
 *     --energy 
 *     --xforce 
 *     --yforce 
 *     --zforce 
 *     --torque ABOUT nx ny nz 
 * 
 *       -------------------------------------------------
 * 
 * c. options specifying frequency behavior
 * 
 *     --Xi xx
 *
 *         Specify a single imaginary frequency at which to
 *         evaluate the spectral density of contributions
 *         to the requested Casimir quantities.
 *
 *         Note that --Xi may be specified more than once.
 * 
 *     --XiFile MyXiFile
 * 
 *         Specify a file containing a list of --Xi values.
 * 
 *     --Temperature T
 * 
 *         Specify a temperature (in Kelvin) at which to  
 *         evaluate the Matsubara sum for the requested Casimir
 *         quantities.
 * 
 *     Note: if no frequency options are specified, the 
 *           default behavior is to integrate over the  
 *           entire positive imaginary frequency axis to compute
 *           the full 
 * 
 *       -------------------------------------------------
 * 
 * c. options specifying output file names 
 * 
 *     --ByXiFile MyFileName.byXi
 * 
 *         Set the name of the frequency-resolved output  
 *         file. If this option is not specified, the 
 *         frequency-resolved output file will be called
 *         Geometry.byXi, where Geometry.scuffgeo is the
 *         geometry file specified using the --geometry option.
 * 
 *     --OutputFile MyFile.Out
 * 
 *         Set the name of the output file (which is only generated 
 *         if your command-line options call for scuff-cas3D to 
 *         evaluate a frequency integral or Matsubara sum.)
 *         If this option is not specified, the output file will 
 *         be called Geometry.out, where Geometry.scuffgeo is the
 *         geometry file specified using the --geometry option.
 * 
 *     --LogFile MyFile.log
 * 
 *         Set the name of the log file. If this option is not
 *         specified, the log file will be named scuff-cas3D.log.
 * 
 *       -------------------------------------------------
 * 
 * d. options specifying scuff caches 
 * 
 *     -- ReadCache MyReadCache.scuffcache
 * 
 *         Specify a scuff cache file to be preloaded.
 *         This option may be specified more than once. 
 * 
 *     -- WriteCache MyWriteCache.scuffcache
 * 
 *         Specify the name of a scuff cache file to which
 *         scuff-cas3D will dump the contents of its cache
 *         after completing all computations.
 * 
 *     -- Cache MyCache.scuffcache
 * 
 *         Specify a cache file that scuff-cas3D will both
 *         (a) preload before starting its computations, and
 *         (b) overwrite after completing its computations.
 *         Specifying this option is equivalent to setting
 *         --ReadCache and --WriteCache both to MyCache.scuffcache.
 * 
 * e. other options 
 * 
 *     --nThread xx   (use xx computational threads)
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

#include "scuff-cas3D.h"

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define MAXFREQ  10    // max number of frequencies 
#define MAXCACHE 10    // max number of cache files for preload

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
  int Energy=0;
  int XForce=0;
  int YForce=0;
  int ZForce=0;
  char *TorqueArgs[TORQUEARGS*MAXTORQUE];   int nTorque;
  double XiVals[MAXFREQ];	   int nXiVals;
  char *XiFile;			   int nXiFiles;
  double Temperature=0.0;          int nTemperature;
  char *ByXiFile=0;
  char *OutputFile=0;
  char *LogFile=0;
  char *Cache=0;
  char *ReadCache[MAXCACHE];         int nReadCache;
  char *WriteCache=0;
  int nThread=0;
  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { {"Geometry",       PA_STRING,  1, 1,       (void *)&GeoFile,    0,             "geometry file"},
     {"TransFile",      PA_STRING,  1, 1,       (void *)&TransFile,  0,             "list of geometrical transformation"},
     {"Energy",         PA_BOOL,    0, 1,       (void *)&Energy,     0,             "compute Casimir energy"},
     {"XForce",         PA_BOOL,    0, 1,       (void *)&XForce,     0,             "compute Casimir energy"},
     {"YForce",         PA_BOOL,    0, 1,       (void *)&YForce,     0,             "compute Casimir energy"},
     {"ZForce",         PA_BOOL,    0, 1,       (void *)&ZForce,     0,             "compute Casimir energy"},
     {"Torque",         PA_STRING,  4, 3,       (void *)TorqueArgs,  &nTorque,      "compute Casimir torque"},
     {"Xi",             PA_DOUBLE,  1, MAXFREQ, (void *)XiVals,      &nXiVals,      "imaginary frequency"},
     {"XiFile",         PA_STRING,  1, 1,       (void *)&XiFile,     &nXiFiles,     "list of --Xi values"},
     {"Temperature",    PA_DOUBLE,  1, 1,       (void *)&Temperature, &nTemperature, "temperature in Kelvin"},
     {"OutputFile",     PA_STRING,  1, 1,       (void *)&OutputFile, 0,             "name of frequency-integrated output file"},
     {"ByXiFile",       PA_STRING,  1, 1,       (void *)&ByXiFile,   0,             "name of frequency-resolved output file"},
     {"LogFile",        PA_STRING,  1, 1,       (void *)&LogFile,    0,             "name of log file"},
     {"Cache",          PA_STRING,  1, 1,       (void *)&Cache,      0,             "read/write cache"},
     {"ReadCache",      PA_STRING,  1, MAXCACHE,(void *)ReadCache,   &nReadCache,   "read cache"},
     {"WriteCache",     PA_STRING,  1, 1,       (void *)&WriteCache, 0,             "write cache"},
     {"nThread",        PA_INT,     1, 1,       (void *)&nThread,    0,             "number of CPU threads to use"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  if (GeoFile==0)
   OSUsage(argv[0], OSArray, "--geometry option is mandatory");

  if (LogFile)
   SetLogFileName(LogFile);
  else
   SetLogFileName("scuff-cas3D.log");

  Log("scuff-cas3D running on %s",GetHostName());

  /*******************************************************************/
  /* process any frequency-related options to construct a list of    */
  /* frequencies at which to compute Casimir quantities              */
  /*******************************************************************/
  HVector *Xiist=0, *XiList0;
  int nFreq, nOV, NumFreqs=0;
  if (nXiFiles==1) // first process --XiFile option if present
   { 
     XiList=new HVector(XiFile,LHM_TEXT);
     if (XiList->ErrMsg)
      ErrExit(XiList->ErrMsg);
     NumFreqs=XiList->N;
     Log("Read %i frequencies from file %s.",NumFreqs,XiFile);
   };

  // now add any individually specified --Xi options
  if (nXiVals>0)
   { 
     NumFreqs += nXiVals;
     XiList0=XiList;
     XiList=new HVector(NumFreqs);
     nFreq=0;
     if (XiList0)
      { for(nFreq=0; nFreq<XiList0->N; nFreq++)
         XiList->SetEntry(nFreq, XiList0->GetEntry(nFreq));
        delete XiList0;
      };
     for(nOV=0; nOV<nXiVals; nOV++)
      XiList->SetEntry(nFreq+nOV, XiVals[nOV]);
     Log("Read %i frequencies from command line.",nXiVals);
   };

  // if the user specified a temperature, check that the
  // temperature makes sense and that the user didn't also 
  // specify a list of frequencies
  if ( nTemperature!=0 )
   { if (Temperature<0.0)
      ErrExit("negative value specified for --Temperature");
     if ( NumFreqs>0 )
      ErrExit("--Temperature may not be used with --Xi/--XiList");
   };

  /*******************************************************************/
  /* create the SC3Data structure that contains all the info needed  */
  /* to evaluate the contributions of a single frequency to the      */
  /* Casimir quantities                                              */
  /*******************************************************************/

  SC3Data *SC3D=(SC3Data *)malloc(sizeof(SC3Data));

  RWGGeometry *G = new RWGGeometry(GeoFile);
  G->SetLogLevel(SCUFF_VERBOSELOGGING);
  SC3D->G=G;

  if (ByXiFile)
   SC3Data->ByXiFile = ByXiFile;a
  else
   { SHD->ByXiFile = vstrdup("%s.byXi",GetFileBase(GeoFile));
     char MyFileName[MAXSTR];
     FILE *f=CreateUniqueFile(SHD->ByOmegaFile, 1, MyFileName);
     fclose(f);
     SHD->ByOmegaFile=strdup(MyFileName);
   };

  /*******************************************************************/
  /* preload the scuff cache with any cache preload files the user   */
  /* may have specified                                              */
  /*******************************************************************/
  if ( Cache!=0 && WriteCache!=0 )
   ErrExit("--cache and --writecache options are mutually exclusive");
  for (int nrc=0; nrc<nReadCache; nrc++)
   PreloadGlobalFIPPICache( ReadCache[nrc] );
  if (Cache)
   PreloadGlobalFIPPICache( Cache );

  if (Cache) WriteCache=Cache;
  SC3D->WriteCache = WriteCache;

  /*******************************************************************/
  /* now switch off based on the requested frequency behavior to     */
  /* perform the actual calculations                                 */
  /*******************************************************************/
  double *I = new double[SC3D->NumTransformations];
  if (NumFreqs>0)
   { for (nFreq=0; nFreq<NumFreqs; nFreq++)
      GetFrequencyIntegrand(SC3D, XiList->GetEntry(nFreq), I);
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
