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
 *           the full zero-temperature Casimir quantities.
 * 
 *       -------------------------------------------------
 * 
 * d. options specifying output file names 
 * 
 *     --ByXiFile MyFileName.byXi
 * 
 *         Set the name of the frequency-resolved output  
 *         file. If this option is not specified, the 
 *         frequency-resolved output file will be called
 *         Geometry.byXi, where Geometry.scuffgeo is the
 *         geometry file specified using the --geometry option.
 * 
 *     --ByXikBlochFile MyFileName.byXikBloch
 * 
 *         Set the name of the frequency- and Bloch-vector-resolved 
 *         output file. (For periodic geometries only.) If this option 
 *         is not specified, this output file will not be created.
 * 
 *     --OutputFile MyFile.Out
 * 
 *         Set the name of the output file. (This file is only generated 
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
 *     --BZIMethod xx 
 * 
 *       Use method xx for integrating over the Brillouin zone.
 * 
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
//
  bool Energy=false;
  bool XForce=false;
  bool YForce=false;
  bool ZForce=false;
  double TorqueAxes[9];                     int nTorque;
  bool AllTorque=false;
//
  double Temperature=0.0;                   int nTemperature;
  double XiVals[MAXFREQ];	            int nXiVals;
  char *XiFile=0;
  char *XikBlochFile=0;
  char *BZIString=0;
//
  int MaxXiPoints=10000;
  double AbsTol=0.0;
  double RelTol=1.0e-2;
  int Intervals=50;
//
  char *OutputFile=0;
  char *ByXiFile=0;
  char *ByXikBlochFile=0;
  char *LogFile=0;
//
  char *Cache=0;
  char *ReadCache[MAXCACHE];                int nReadCache;
  char *WriteCache=0;
//
  bool NewEnergyMethod = false;
//
  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { {"Geometry",       PA_STRING,  1, 1,       (void *)&GeoFile,       0,             "geometry file"},
     {"TransFile",      PA_STRING,  1, 1,       (void *)&TransFile,     0,             "list of geometrical transformation"},
//
     {"Energy",         PA_BOOL,    0, 1,       (void *)&Energy,        0,             "compute Casimir energy"},
     {"XForce",         PA_BOOL,    0, 1,       (void *)&XForce,        0,             "compute x-directed Casimir force"},
     {"YForce",         PA_BOOL,    0, 1,       (void *)&YForce,        0,             "compute y-directed Casimir force"},
     {"ZForce",         PA_BOOL,    0, 1,       (void *)&ZForce,        0,             "compute z-directed Casimir force"},
     {"Torque",         PA_DOUBLE,  3, 3,       (void *)TorqueAxes,     &nTorque,      "compute Casimir torque about a given axis"},
     {"AllTorque",      PA_BOOL,    0, 1,       (void *)&AllTorque,     0,             "compute all three Casimir torque components"},
//
     {"Temperature",    PA_DOUBLE,  1, 1,       (void *)&Temperature,   &nTemperature, "temperature in Kelvin"},
     {"Xi",             PA_DOUBLE,  1, MAXFREQ, (void *)XiVals,         &nXiVals,      "imaginary frequency"},
     {"XiFile",         PA_STRING,  1, 1,       (void *)&XiFile,        0,             "list of --Xi values"},
     {"XikBlochFile",   PA_STRING,  1, 1,       (void *)&XikBlochFile,  0,             "list of (--Xi, --kBloch) values"},
     {"BZIMethod",      PA_STRING,  1, 1,       (void *)&BZIString,     0,             "Brillouin-zone integration method"},
//
     {"OutputFile",     PA_STRING,  1, 1,       (void *)&OutputFile,    0,             "name of frequency-integrated output file"},
     {"ByXiFile",       PA_STRING,  1, 1,       (void *)&ByXiFile,      0,             "name of frequency-resolved output file"},
     {"ByXikBlochFile", PA_STRING,  1, 1,       (void *)&ByXikBlochFile, 0,            "name of frequency- and kBloch-resolved output file"},
     {"LogFile",        PA_STRING,  1, 1,       (void *)&LogFile,       0,             "name of log file"},
//
     {"MaxXiPoints",    PA_INT,     1, 1,       (void *)&MaxXiPoints,   0,             "maximum number of Xi integrand evaluations "},
     {"AbsTol",         PA_DOUBLE,  1, 1,       (void *)&AbsTol,        0,             "absolute tolerance for sums and integrations"},
     {"RelTol",         PA_DOUBLE,  1, 1,       (void *)&RelTol,        0,             "relative tolerance for sums and integrations"},
     {"Intervals",      PA_INT,     1, 1,       (void *)&Intervals,     0,             "number of subintervals for frequency quadrature"},
//
     {"Cache",          PA_STRING,  1, 1,       (void *)&Cache,         0,             "read/write cache"},
     {"ReadCache",      PA_STRING,  1, MAXCACHE,(void *)ReadCache,      &nReadCache,   "read cache"},
     {"WriteCache",     PA_STRING,  1, 1,       (void *)&WriteCache,    0,             "write cache"},
//
     {"NewEnergyMethod", PA_BOOL,   0, 1,       (void *)&NewEnergyMethod, 0,           "use alternative method for energy calculation"},
//
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
  //G->SetLogLevel(SCUFF_TERSELOGGING);
  G->SetLogLevel(SCUFF_VERBOSELOGGING);

  /***************************************************************/
  /* process frequency- and kBloch-related options               */
  /***************************************************************/
  HMatrix *XikBlochPoints=0;
  HMatrix *XiPoints=0;

  if ( XikBlochFile )
   { if (XiFile)       ErrExit("--XikBlochFile and --XiFile options are mutually exclusive");
     if (nXiVals>0)    ErrExit("--XikBlochFile and --Xi options are mutually exclusive");
     if (nTemperature) ErrExit("--XikBlochFile and --Temperature options are mutually exclusive");
     if (G->NumLatticeBasisVectors==0) ErrExit("--XikBlochFile may only be used for periodic geometries");
     XikBlochPoints=new HMatrix(XikBlochFile,LHM_TEXT,"--nc 3 --strict");
     if (XikBlochPoints->ErrMsg)
      ErrExit(XikBlochPoints->ErrMsg);
     Log("Read %i (Xi, kBloch) points from file %s.",XikBlochPoints->NR, XikBlochFile);
   }
  else if ( XiFile )
   { if (nXiVals>0)    ErrExit("--XiFile and --Xi options are mutually exclusive");
     if (nTemperature) ErrExit("--XiFile and --Temperature options are mutually exclusive");
     XiPoints = new HMatrix(XiFile,LHM_TEXT,"--nc 1 --strict");
     if (XiPoints->ErrMsg)
      ErrExit(XiPoints->ErrMsg);
     Log("Read %i Xi points from file %s.",XiPoints->NR, XiFile);
   }
  else if ( nXiVals>0 )
   { if (nTemperature) ErrExit("--Xi and --Temperature options are mutually exclusive");
     XiPoints = new HMatrix(nXiVals, 1);
     for(int nxv=0; nxv<nXiVals; nxv++) 
      XiPoints->SetEntry(nxv, 0, XiVals[nxv]);
     Log("Performing Casimir calculations at %i command-line Xi points.",nXiVals);
   }
  else if ( nTemperature==1 )
   Log("Computing full Matsubara-summed Casimir quantities at T=%e Kelvin.",Temperature);
  else
   Log("Computing full zero-temperature Casimir quantities.");

  int BZIMethod=BZIMETHOD_DEFAULT;
  if ( BZIString )
   { if (G->NumLatticeBasisVectors==0) ErrExit("--BZIMethod option may only be used for periodic geometries");
     if ( !strcasecmp(BZIString,"MP3") )
      { Log("Using 3-point monkhorst-pack scheme for Brillouin zone integration.");
        BZIMethod=BZIMETHOD_MP3;
      }
     else if ( !strcasecmp(BZIString,"MP6") )
      { Log("Using 6-point monkhorst-pack scheme for Brillouin zone integration.");
        BZIMethod=BZIMETHOD_MP6;
      }
     else if ( !strcasecmp(BZIString,"MP10") )
      { Log("Using 10-point monkhorst-pack scheme for Brillouin zone integration.");
        BZIMethod=BZIMETHOD_MP10;
      }
     else if ( !strcasecmp(BZIString,"MP15") )
      { Log("Using 15-point monkhorst-pack scheme for Brillouin zone integration.");
        BZIMethod=BZIMETHOD_MP15;
      }
     else if ( !strcasecmp(BZIString,"adaptive") )
      { Log("Using adaptive cubature scheme for Brillouin zone integration.");
        BZIMethod=BZIMETHOD_ADAPTIVE;
      };
   };

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
  if (AllTorque)
   { if (nTorque>0) ErrExit("--AllTorque and --Torque options are mutually exclusive");
     NumQuantities+=3; 
     WhichQuantities |= (QUANTITY_TORQUE1 + QUANTITY_TORQUE2 + QUANTITY_TORQUE3);
     nTorque=3;
     memset(TorqueAxes, 0, 9*sizeof(double));
     TorqueAxes[0]=TorqueAxes[4]=TorqueAxes[8]=1.0;
   };

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

  /*******************************************************************/
  /* create the SC3Data structure that contains all the info needed  */
  /* to evaluate the contributions of a single frequency and kBloch  */
  /* point to the Casimir quantities                                 */
  /*******************************************************************/
  SC3Data *SC3D=CreateSC3Data(G, TransFile, WhichQuantities, NumQuantities, 
                              nTorque, TorqueAxes, NewEnergyMethod);

  SC3D->BZIMethod          = BZIMethod;
  SC3D->ByXiFileName       = ByXiFile; 
  SC3D->ByXikBlochFileName = ByXikBlochFile;
  SC3D->WriteCache         = WriteCache;
  SC3D->AbsTol             = AbsTol;
  SC3D->RelTol             = RelTol;
  SC3D->MaxXiPoints        = MaxXiPoints;

  if (SC3D->ByXiFileName==0)
   SC3D->ByXiFileName=vstrdup("%s.byXi",GetFileBase(G->GeoFileName));

  if (G->NumLatticeBasisVectors==0 && SC3D->ByXikBlochFileName!=0)
   { Warn("--byXikBlochFile option only makes sense for periodic geometries (disabling)"); 
     SC3D->ByXikBlochFileName=0;
   };

  /*******************************************************************/
  /* now switch off based on the requested frequency behavior to     */
  /* perform the actual calculations                                 */
  /*******************************************************************/
  double *EFT = new double[SC3D->NTNQ];
  double *Error=0; 
  if ( XikBlochPoints )
   { 
     double Xi, kBloch[2];
     for(int nr=0; nr<XikBlochPoints->NR; nr++)
      { 
        Xi        = XikBlochPoints->GetEntryD(nr, 0);
        kBloch[0] = XikBlochPoints->GetEntryD(nr, 1);
        kBloch[1] = XikBlochPoints->GetEntryD(nr, 2);
        GetCasimirIntegrand(SC3D, Xi, kBloch, EFT);
      };
   }
  else if ( XiPoints )
   { 
     for (int nr=0; nr<XiPoints->NR; nr++)
      GetXiIntegrand(SC3D, XiPoints->GetEntryD(nr,0), EFT);
   }
  else if ( Temperature > 0.0)
   { 
     Error = new double[SC3D->NTNQ];
     GetMatsubaraSum(SC3D, Temperature, EFT, Error);
   }
  else
   { 
     Error = new double[SC3D->NTNQ];
     GetXiIntegral2(SC3D, Intervals, EFT, Error);
   };

  /***************************************************************/
  /* write output file if we computed summed or integrated quantities */
  /***************************************************************/
  if (Error)
   { if (!OutputFile)
      OutputFile=vstrdup("%s.out",GetFileBase(G->GeoFileName));
     FILE *OutFile=CreateUniqueFile(OutputFile,1);
     int ntnq=0;
     for(int nt=0; nt<SC3D->NumTransformations; nt++)
      { fprintf(OutFile,"%s ",SC3D->GTCList[nt]->Tag);
        for(int nq=0; nq<SC3D->NumQuantities; nq++, ntnq++)
         fprintf(OutFile,"%e %e ",EFT[ntnq],Error[ntnq]);
        fprintf(OutFile,"\n");
      };
     fclose(OutFile); 
     delete[] Error;
   };

  delete[] EFT;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  printf("Thank you for your support.\n");

}
