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
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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
  InitializeLog(argv[0]);

  /***************************************************************/
  /* pre-process command-line arguments to extract arguments     */
  /* any relevant for Brillouin-zone integration                 */
  /***************************************************************/
  GetBZIArgStruct *BZIArgs=InitBZIArgs(argc, argv);

  /***************************************************************/
  /* process options *********************************************/
  /***************************************************************/
  char *GeoFile=0;
  char *TransFile=0;

  //
  // options specifying quantities computed
  //
  bool Energy=false;
  bool XForce=false;
  bool YForce=false;
  bool ZForce=false;
  double TorqueAxes[9];                     int nTorque;
  bool AllTorque=false;

  //
  // options affecting frequency integration or sampling
  //
  double Temperature=0.0;                   int nTemperature;
  double XiVals[MAXFREQ];	            int nXiVals;
  char *XiFile=0;
  char *XikBlochFile=0;
  char *XiQuadrature=0;
  double XiMin=0.001;
  int MaxXiPoints=10000;
  int Intervals=50;
  double AbsTol=0.0;
  double RelTol=1.0e-2;

  //
  // option allowing user to override default output file names
  //
  char *FileBase=0;

  //
  // cache options
  //
  char *Cache=0;
  char *ReadCache[MAXCACHE];                int nReadCache;
  char *WriteCache=0;

  //
  // other miscellaneous flags
  //
  bool UseExistingData = false;
  bool NewEnergyMethod = false;
  bool WriteHDF5Files  = false;

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
     {"XiFile",         PA_STRING,  1, 1,       (void *)&XiFile,        0,             "file containing Xi values"},
     {"XikBlochFile",   PA_STRING,  1, 1,       (void *)&XikBlochFile,       0,             "file containing (Xi, kx, ky) values"},
     {"XiQuadrature",   PA_STRING,  1, 1,       (void *)&XiQuadrature,  0,             "quadrature method for Xi integration"},
     {"XiMin",          PA_DOUBLE,  1, 1,       (void *)&XiMin,         0,             "assume Xi integrand constant below this value"},
     {"MaxXiPoints",    PA_INT,     1, 1,       (void *)&MaxXiPoints,   0,             "maximum number of Xi integrand evaluations"},
     {"Intervals",      PA_INT,     1, 1,       (void *)&Intervals,     0,             "number of subintervals for frequency quadrature"},
     {"AbsTol",         PA_DOUBLE,  1, 1,       (void *)&AbsTol,        0,             "absolute tolerance for sums and integrations"},
     {"RelTol",         PA_DOUBLE,  1, 1,       (void *)&RelTol,        0,             "relative tolerance for sums and integrations"},
//
     {"FileBase",       PA_STRING,  1, 1,       (void *)&FileBase,      0,             "base filename for output files"},
//
     {"Cache",          PA_STRING,  1, 1,       (void *)&Cache,         0,             "read/write cache"},
     {"ReadCache",      PA_STRING,  1, MAXCACHE,(void *)ReadCache,      &nReadCache,   "read cache"},
     {"WriteCache",     PA_STRING,  1, 1,       (void *)&WriteCache,    0,             "write cache"},
//
     {"UseExistingData", PA_BOOL,   0, 1,       (void *)&UseExistingData, 0,           "reuse data from existing .byXi files"},
//
     {"NewEnergyMethod", PA_BOOL,   0, 1,       (void *)&NewEnergyMethod, 0,           "use alternative method for energy calculation"},
//
     {"WriteHDF5Files", PA_BOOL,    1, 1,       (void *)&WriteHDF5Files,0,             "write BEM matrices to .hdf5 files"},
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

  /***************************************************************/
  /* try to create the geometry  *********************************/
  /***************************************************************/
  RWGGeometry *G = new RWGGeometry(GeoFile);
  //G->SetLogLevel(SCUFF_TERSELOGGING);
  G->SetLogLevel(SCUFF_VERBOSELOGGING);

  /***************************************************************/
  /* process frequency- and kBloch-related options               */
  /***************************************************************/
  HMatrix *XiPoints=0;
  HMatrix *XiKPoints=0;

  if ( XikBlochFile )
   { if (XiFile)       ErrExit("--XikBlochFile and --XiFile options are mutually exclusive");
     if (nXiVals>0)    ErrExit("--XikBlochFile and --Xi options are mutually exclusive");
     if (nTemperature) ErrExit("--XikBlochFile and --Temperature options are mutually exclusive");
     if (G->LDim==0)   ErrExit("--XikBlochFile may only be used for periodic geometries");
     XiKPoints=new HMatrix(XikBlochFile,LHM_TEXT,"--nc 3 --strict");
     if (XiKPoints->ErrMsg)
      ErrExit(XiKPoints->ErrMsg);
     Log("Read %i (Xi, kBloch) points from file %s.",XiKPoints->NR, XikBlochFile);
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

  /*******************************************************************/
  /* process integration options *************************************/
  /*******************************************************************/
  if ( XiQuadrature && XiPoints )
   ErrExit("--XiQuadrature is incompatible with --Xi or --XiFile");
  if ( XiQuadrature && XiKPoints )
   ErrExit("--XiQuadrature is incompatible with --XikBlochFile ");
  if ( XiQuadrature && XiKPoints )
   ErrExit("--XiQuadrature is incompatible with --Temperature ");

  int XiQMethod = QMETHOD_CLIFF;
  if (XiQuadrature)
   { if ( !strcasecmp(XiQuadrature,"CLIFF") )
      { XiQMethod = QMETHOD_CLIFF;
        Log("Using cliff integration method for Xi quadrature.");
      }
     if ( !strcasecmp(XiQuadrature,"TRAPSIMP") )
      { XiQMethod = QMETHOD_TRAPSIMP;
        Log("Using trap/simp integration method (%i intervals) for Xi quadrature.",Intervals);
      }
     else if ( !strcasecmp(XiQuadrature,"ADAPTIVE") )
      { XiQMethod = QMETHOD_ADAPTIVE;
        Log("Using adaptive integration method for Xi quadrature.");
      }
     else
      ErrExit("unknown value %s specified for --XiQuadrature",XiQuadrature);
   };

  /*******************************************************************/
  /* figure out which quantities to compute **************************/
  /*******************************************************************/
  if (G->LDim>0)
   { if ( XForce || YForce || nTorque>0 || AllTorque)
      ErrExit("only --energy and --zforce are available for extended geometries");
   };

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
                              nTorque, TorqueAxes, NewEnergyMethod, FileBase);

  SC3D->WriteCache         = WriteCache;
  SC3D->WriteHDF5Files     = WriteHDF5Files;
  SC3D->AbsTol             = AbsTol;
  SC3D->RelTol             = RelTol;
  SC3D->UseExistingData    = UseExistingData;
  SC3D->MaxXiPoints        = MaxXiPoints;
  SC3D->XiMin              = XiMin;

  if (G->LDim>=1)
   { UpdateBZIArgs(BZIArgs, G->RLBasis, G->RLVolume);
     BZIArgs->BZIFunc  = GetCasimirIntegrand;
     BZIArgs->UserData = (void *)SC3D;
     BZIArgs->FDim     = SC3D->NTNQ;
     SC3D->BZIArgs     = BZIArgs;
   };

  /*******************************************************************/
  /* now switch off based on the requested frequency behavior to     */
  /* perform the actual calculations                                 */
  /*******************************************************************/
  double *EFT = new double[SC3D->NTNQ];
  double *Error=0; 
  if ( XiKPoints )
   { 
     for(int nr=0; nr<XiKPoints->NR; nr++)
      { 
        double Xi, kBloch[2]={0.0, 0.0};
        Xi        = XiKPoints->GetEntryD(nr, 0);
        kBloch[0] = XiKPoints->GetEntryD(nr, 1);
        if (G->LDim>=2)
         kBloch[1] = XiKPoints->GetEntryD(nr, 2);
        GetCasimirIntegrand(SC3D, cdouble(0.0,Xi), kBloch, EFT);
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
  else // full Xi integration for zero-temperature casimir quantities
   { 
     Error = new double[SC3D->NTNQ];
     if (XiQMethod == QMETHOD_ADAPTIVE)
      GetXiIntegral_Adaptive(SC3D, EFT, Error);
     else if (XiQMethod == QMETHOD_TRAPSIMP)
      GetXiIntegral_TrapSimp(SC3D, Intervals, EFT, Error);
     else
      GetXiIntegral_Cliff(SC3D, EFT, Error);
   };

  /***************************************************************/
  /* write output file if we computed summed or integrated quantities */
  /***************************************************************/
  if (Error)
   { 
     WriteFilePreamble(SC3D, PREAMBLE_OUT);
     FILE *OutFile=fopen(SC3D->OutFileName,"a");
     int ntnq=0;
     for(int nt=0; nt<SC3D->NumTransformations; nt++)
      { fprintf(OutFile,"%s ",SC3D->GTCs[nt]->Tag);
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
