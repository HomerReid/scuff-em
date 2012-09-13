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
 *     --XikBlochFile MyXikBlochFile
 * 
 *         Specify a file containing a list of --Xi and --kBloch
 *         values (see below)
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
 * d. options specifying bloch wavevector (for periodically 
 *    extended geometries only)                                 
 * 
 *  --kBloch xx yy 
 *
 *         Specify a single bloch wavevector at which to 
 *         evaluate the spectral density of contributions
 *         to the requested Casimir quantities.
 *
 *  --kBlochFile MykBlochFile
 * 
 *         Specify a file containing a list of --kBloch
 *         values.
 * 
 *     Note: if none of the options --kBloch, --kBlochFile, or
 *           --XikBlochFile are specified, the default is
 *           to integrate over the Brillouin zone. 
 * 
 *       -------------------------------------------------
 * 
 * e. options specifying output file names 
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
 * f. options specifying scuff caches 
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
 * g. other options 
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

/*******************************************************************/
/* subroutine to process command-line options related to           */
/* frequencies and bloch wavevectors to figure out the list of     */
/* (Xi, kx, ky) points at which we will compute Casimir quantities.*/
/*                                                                 */
/* On return, one of the following is true.                        */
/*                                                                 */
/*  a) *pXikBlochPoints is non-NULL and points to an Nx3 HMatrix   */
/*     whose rows are the (Xi, kx, ky) points at which we compute  */
/*     Casimir quantities.                                         */
/*                                                                 */
/*  b) *pXiPoints is non-NULL and points to an N-dimensional       */
/*      HVector whose entries are the Xi points at which we        */
/*      compute Casimir quantities.                                */
/*      In this case, if a lattice is present, *pBZIMethod tells   */
/*      us which method to use for integrating over the Brillouin  */
/*      zone at each frequency.                                    */
/*                                                                 */
/*  c) *pXikBlochPoints and *pXiPoints are both NULL. This means   */
/*      the user wants us to integrate over the entire positive    */
/*      Xi axis (if no temperature was specified) or else to       */
/*      evaluate the full Matsubara sum at the given temperature.  */
/*      In this case, if a lattice is present, *pBZIMethod tells   */
/*      us which method to use for integrating over the Brillouin  */
/*      zone at each frequency.                                    */
/*                                                                 */
/*  Note that if *pBZIMethod = BZIMETHOD_SINGLEPOINT, then kBloch  */
/*  stores the components of a single Bloch vector at which we do  */
/*  calculations at each frequency.                                */
/*******************************************************************/
void ProcessFrequencyOptions(bool HaveLattice,
                             double *XiVals, int nXiVals,
                             double *kBloch, int nkBloch,
                             char *XiFile, char *XikBlochFile, char *BZIString,
                             double Temperature, int nTemperature,
                             HMatrix **pXikBlochPoints, HVector **pXiPoints,
                             int *pBZIMethod)
{
  /*--------------------------------------------------------------*/
  /*- a deluge of sanity checks ----------------------------------*/
  /*--------------------------------------------------------------*/
  if ( HaveLattice==0 )
   { if ( (nkBLoch!=0) )
       ErrExit("--kBloch may only be specified for .scuffgeo files with LATTICE statements");
     if ( (XikBlochFile!=0) )
       ErrExit("--XikBlochFile may only be specified for .scuffgeo files with LATTICE statements");
     if ( (BZIString!=0) )
       ErrExit("--BZIMethod may only be specified for .scuffgeo files with LATTICE statements");
   };

  if ( BZIString!=0 && nkBloch!=0 )
   ErrExit("--BZIMethod and --kBloch options are mutually exclusive");
  if ( (nkBLoch!=0) && XikBlochFile!=0 )
   ErrExit("--kBloch and --XikBlochFile options are mutually exclusive");
  if ( (BZIString!=0) && XikBlochFile!=0 )
   ErrExit("--BZIMethod and --XikBlochFile options are mutually exclusive");
  if ( (nXiVals!=0) && XikBlochFile!=0 )
   ErrExit("--Xi and --XikBlochFile options are mutually exclusive");

  if ( nTemperature > 0 )
   { if ( nXiVals!=0 )
      ErrExit("--Xi and --Temperature options are mutually exclusive");
     if ( XiFile!=0 )
      ErrExit("--XiFile and --Temperature options are mutually exclusive");
     if ( XikBlochFile!=0 )
      ErrExit("--XikBlochFile and --Temperature options are mutually exclusive");
     if ( Temperature < 0.0 )
      ErrExit("negative value specified for the --Temperature option")
   };

  /*--------------------------------------------------------------*/
  /*- process BZIString argument ---------------------------------*/
  /*--------------------------------------------------------------*/
  *pBZIMethod=BZIMETHOD_MP7;              // default if nothing specified
  if ( BZIString )
   { if ( !strcasecmp(BZIString,"MP7") )
      { Log("Using 7-point monkhorst-pack scheme for Brillouin zone integration.");
        *pBZIMethod=BZIMETHOD_MP7;
      }
     else if ( !strcasecmp(BZIString,"MP15") )
      { Log("Using 15-point monkhorst-pack scheme for Brillouin zone integration.");
        *pBZIMethod=BZIMETHOD_MP15;
      }
     else
      { fprintf(f," error: unknown Brillouin-zone integration scheme %s (aborting)\n",BZIString);
        fprintf(f," \n");
        fprintf(f," Available schemes for Brillouin-zone integration:\n");
        fprintf(f,"  --bziMethod MP7   (7-point monkhorst-pack scheme)\n");
        fprintf(f,"  --bziMethod MP15  (15-point monkhorst-pack scheme)\n");
        exit(1);
      };
   }
  else if ( nkBloch!=0 )
   *pBZIMethod = BZIMETHOD_SINGLEPOINT;
     
  /*--------------------------------------------------------------*/
  /*- look for a XikBloch file                                    */
  /*--------------------------------------------------------------*/
  if (XikBlochFile)
   { *pXikBlochPoints=new HMatrix(XikBlochFile);
     if (*pXikBlochPoints->ErrMsg) 
      ErrExit(*pXikBlochPoints->ErrMsg);
     Log("Read %i (Xi, kBloch) points from file %s.",*pXikBlockPoints->NR, XikBlochFile);
     *pXiPoints=0;
     return;
   };

  /*--------------------------------------------------------------*/
  /*- otherwise, the user may have given just a list of Xi points,*/
  /*- on the command line and/or in a --XiFile                    */
  /*--------------------------------------------------------------*/
  HVector *XiList1 = nXiVals==0 ? 0 : new HVector(XiVals, nXiVals);
  HVector *XiList2 =  XiFile==0 ? 0 : new HVector(XiFile);
  if (XiList2 && XiList2->ErrMsg)
   ErrExit(XiList2->ErrMsg);
  *pXiPoints = Concat(XiList1, XiList2); // note the return value may be NULL


}
  

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
  int AllTorque=0;
  char *TorqueArgs[TORQUEARGS*MAXTORQUE];   int nTorque;
  double XiVals[MAXFREQ];	            int nXiVals;
  double kBloch[2];                         int nkBloch;
  char *XiFile=0;
  char *XikBlochFile=0;
  double Temperature=0.0;                   int nTemperature;
  char *BZIMethod;      
  char *OutputFile=0;
  char *ByXiFile=0;
  char *ByXikBlochFile=0;
  char *LogFile=0;
  char *Cache=0;
  char *ReadCache[MAXCACHE];                int nReadCache;
  char *WriteCache=0;
  int nThread=0;
  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { {"Geometry",       PA_STRING,  1, 1,       (void *)&GeoFile,       0,             "geometry file"},
     {"TransFile",      PA_STRING,  1, 1,       (void *)&TransFile,     0,             "list of geometrical transformation"},
     {"Energy",         PA_BOOL,    0, 1,       (void *)&Energy,        0,             "compute Casimir energy"},
     {"XForce",         PA_BOOL,    0, 1,       (void *)&XForce,        0,             "compute x-directed Casimir force"},
     {"YForce",         PA_BOOL,    0, 1,       (void *)&YForce,        0,             "compute y-directed Casimir force"},
     {"ZForce",         PA_BOOL,    0, 1,       (void *)&ZForce,        0,             "compute z-directed Casimir force"},
     {"Torque",         PA_STRING,  3, 3,       (void *)TorqueArgs,     &nTorque,      "compute Casimir torque about a given axis"},
     {"AllTorque",      PA_BOOL,    0, 1,       (void *)AllTorque,      0,             "compute all three Casimir torque components"},
     {"Xi",             PA_DOUBLE,  1, MAXFREQ, (void *)XiVals,         &nXiVals,      "imaginary frequency"},
     {"kBloch",         PA_STRING,  2, 1,       (void *)kBloch,         &nkBloch,      "bloch wavevector"},
     {"XiFile",         PA_STRING,  1, 1,       (void *)&XiFile,        &nXiFiles,     "list of --Xi values"},
     {"XikBlochFile",   PA_STRING,  1, 1,       (void *)&XikBlochFile,  0,             "list of (--Xi, --kBloch) values"},
     {"Temperature",    PA_DOUBLE,  1, 1,       (void *)&Temperature,   &nTemperature, "temperature in Kelvin"},
     {"BZIMethod",      PA_STRING,  1, 1,       (void *)&BZIString,     0,             "Brillouin-zone integration method"},
     {"OutputFile",     PA_STRING,  1, 1,       (void *)&OutputFile,    0,             "name of frequency-integrated output file"},
     {"ByXiFile",       PA_STRING,  1, 1,       (void *)&ByXiFile,      0,             "name of frequency-resolved output file"},
     {"ByXikBlochFile", PA_STRING,  1, 1,       (void *)&ByXikBlochFile, 0,            "name of frequency- and kBloch-resolved output file"},
     {"LogFile",        PA_STRING,  1, 1,       (void *)&LogFile,       0,             "name of log file"},
     {"Cache",          PA_STRING,  1, 1,       (void *)&Cache,         0,             "read/write cache"},
     {"ReadCache",      PA_STRING,  1, MAXCACHE,(void *)ReadCache,      &nReadCache,   "read cache"},
     {"WriteCache",     PA_STRING,  1, 1,       (void *)&WriteCache,    0,             "write cache"},
     {"nThread",        PA_INT,     1, 1,       (void *)&nThread,       0,             "number of CPU threads to use"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  if (LogFile)
   SetLogFileName(LogFile);
  else
   SetLogFileName("scuff-cas3D.log");
  Log("scuff-cas3D running on %s",GetHostName());

  /***************************************************************/
  /* try to create the geometry  *********************************/
  /***************************************************************/
  if (GeoFile==0)
   OSUsage(argv[0], OSArray, "--geometry option is mandatory");
  RWGGeometry *G = new RWGGeometry(GeoFile);
  G->SetLogLevel(SCUFF_VERBOSELOGGING);

  /***************************************************************/
  /* process frequency- and kBloch-related options               */
  /***************************************************************/
  HMatrix *XikBlochFile;
  HVector *XiPoints;
  ProcessFrequencyOptions(G->NumLatticeBasisVectors > 0 ,
                          XiVals, nXiVals, kBloch, nkBloch, 
                          XiFile, XikBlochFile, BZIString,
                          Temperature, nTemperature,
                          &XikBlochPoints, &XiPoints, &BZIMethod);

  /*******************************************************************/
  /* figure out which quantities to compute **************************/
  /*******************************************************************/
  int WhichQuantities=0, NumQuantities=0;
  if (Energy)    { NumQuantities++; WhichQuantities |= QUANTITY_ENERGY; };
  if (XForce)    { NumQuantities++; WhichQuantities |= QUANTITY_XFORCE; };
  if (YForce)    { NumQuantities++; WhichQuantities |= QUANTITY_YFORCE; };
  if (ZForce)    { NumQuantities++; WhichQuantities |= QUANTITY_ZFORCE; };
  if (nTorque>0) { NumQuantities++; WhichQuantities |= QUANTITY_TORQUE1; };
  if (nTorque>1) { NumQuantities++; WhichQuantities |= QUANTITY_TORQUE2; };
  if (nTorque>2) { NumQuantities++; WhichQuantities |= QUANTITY_TORQUE3; };

  /*******************************************************************/
  /* create the SC3Data structure that contains all the info needed  */
  /* to evaluate the contributions of a single frequency and kBloch  */
  /* point to the Casimir quantities                                 */
  /*******************************************************************/
  SC3Data *SC3D=CreateSC3Data(G, TransFile, ByXiFile, ByXikBlochFile, 
                              WhichQuantities, NumQuantities, Torque);

  /*******************************************************************/
  /* preload the scuff cache with any cache preload files the user   */
  /* may have specified                                              */
  /*******************************************************************/
  if ( Cache!=0 && WriteCache!=0 )
   ErrExit("--cache and --writecache options are mutually exclusive");
  for (int nrc=0; nrc<nReadCache; nrc++)
   PreloadCache( ReadCache[nrc] );
  if (Cache)
   PreloadCache( Cache );

  if (Cache) WriteCache=Cache;
  SC3D->WriteCache = WriteCache;

  /*******************************************************************/
  /* now switch off based on the requested frequency behavior to     */
  /* perform the actual calculations                                 */
  /*******************************************************************/
  double *I = new double[SC3D->NumTransformations];
  if ( XikBlochPoints )
   { 
     for(int nr=0; nr<XikBlochPoints->NR; nr++)
      { 
        Xi        = XikBlochPoints->GetEntry(nr, 0);
        kBloch[0] = XikBlochPoints->GetEntry(nr, 1);
        kBloch[1] = XikBlochPoints->GetEntry(nr, 2);
        GetCasimirIntegrand(SC3D, Xi, kBloch, I); 
      };
   }
  else if ( XiPoints > 0)
   { 
     for (int nr=0; nr<XiPoints->NR; nr++)
      GetFrequencyIntegrand(SC3D, XiList->GetEntry(nr), I);
   }
  else if ( Temperature > 0.0) 
   { 
     GetMatsubaraSum(SC3D, Temperature, I);
   }
  else
   { 
     GetXiIntegral(SC3D, I);
   }
  delete[] I;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  printf("Thank you for your support.\n");

et 
