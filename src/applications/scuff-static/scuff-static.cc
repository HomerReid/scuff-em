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
 * scuff-static.cc   -- a standalone code within the scuff-EM suite for 
 *                   -- solving electrostatics problems 
 *
 * homer reid        -- 10/2006 -- 5/2013
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

#include <libscuff.h>
#include <SSSolver.h>
#include <cmatheval.h>

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define MAXEPF   10    // max number of evaluation-point files
#define MAXCACHE 10    // max number of cache files for preload

#define MAXSTR   1000

/***************************************************************/
/* routines in OutputModules.cc ********************************/
/***************************************************************/
void GetPolarizabilities(SSSolver *SSS, HMatrix *M,
                         HVector *Sigma, char *FileName);

void GetCapacitanceMatrix(SSSolver *SSS, HMatrix *M,
                          HVector *Sigma, char *FileName);

void GetCMatrix(SSSolver *SSS, HMatrix *M,
                HVector *Sigma, char *FileName);

void DoFieldCalculation(SSSolver *SSS, HMatrix *M, HVector *Sigma,
                        char *PotFile, char *PhiExt, int ConstFieldDirection,
                        char *PlotFile, char **EPFiles, int nEPFiles);

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  /***************************************************************/
  /* process options *********************************************/
  /***************************************************************/
  InstallHRSignalHandler();
  char *GeoFile     = 0;
  char *PolFile     = 0;
  char *CapFile     = 0;
  char *CMatrixFile = 0;
  char *PotFile     = 0;
  char *PhiExt      = 0;
  char *EPFiles[MAXEPF];             int nEPFiles;
  char *PlotFile    = 0;
  char *Cache       = 0;
  char *ReadCache[MAXCACHE];         int nReadCache;
  char *WriteCache  = 0;
  char *ConstField  = 0;
  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { 
     {"Geometry",       PA_STRING,  1, 1,       (void *)&GeoFile,    0,             "geometry file"},
/**/
     {"PolFile",        PA_STRING,  1, 1,       (void *)&PolFile,    0,             "polarizability output file"},
/**/
     {"CapFile",        PA_STRING,  1, 1,       (void *)&CapFile,    0,             "capacitance output file"},
/**/
     {"CMatrixFile",    PA_STRING,  1, 1,       (void *)&CMatrixFile, 0,            "c-matrix file"},
/**/
     {"PotentialFile",  PA_STRING,  1, 1,       (void *)&PotFile,    0,             "list of conductor potentials"},
/**/
     {"PhiExternal",    PA_STRING,  1, 1,       (void *)&PhiExt,     0,             "functional form of external potential"},
     {"ConstField",     PA_STRING,  1, 1,       (void *)&ConstField, 0,             "direction of constant unit-strength E field (x,y,z) "},
/**/
     {"EPFile",         PA_STRING,  1, MAXEPF,  (void *)EPFiles,     &nEPFiles,     "list of evaluation points"},
     {"PlotFile",       PA_STRING,  1, 1,       (void *)&PlotFile,   0,             "surface charge visualization output file"},
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

  /*******************************************************************/
  /* sanity check on input arguments *********************************/
  /*******************************************************************/
  if ( (PolFile || CapFile || CMatrixFile) && (nEPFiles>0 || PotFile!=0 || PhiExt!=0 ) )
   ErrExit("(--EPFile,--PotFile,--PhiExternal) may not be used with (--polfile, --capfile)");
  if ( nEPFiles==0 && !PolFile && !CapFile && !CMatrixFile && !PlotFile )
   OSUsage(argv[0], OSArray, "you have not selected any type of calculation");

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  int ConstFieldDirection=-1;
  if (ConstField)
   { switch(tolower(ConstField[0]))
      { case 'x': ConstFieldDirection=0; break;
        case 'y': ConstFieldDirection=1; break;
        case 'z': ConstFieldDirection=2; break;
        default:  ErrExit("invalid --ConstField specification");
      };
   };

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  SetLogFileName("scuff-static.log");
  Log("scuff-static running on %s",GetHostName());

  /*******************************************************************/
  /* create the ScuffStaticGeometry **********************************/
  /*******************************************************************/
  SSSolver *SSS   = new SSSolver(GeoFile);
  RWGGeometry *G  = SSS->G;
  HMatrix *M      = SSS->AllocateBEMMatrix();
  HVector *Sigma  = SSS->AllocateRHSVector();

  /*******************************************************************/
  /* preload the scuff cache with any cache preload files the user   */
  /* may have specified                                              */
  /*******************************************************************/
  if ( Cache!=0 && WriteCache!=0 )
   ErrExit("--cache and --writecache options are mutually exclusive");
  if (Cache) 
   WriteCache=Cache;
  for (int nrc=0; nrc<nReadCache; nrc++)
   PreloadCache( ReadCache[nrc] );
  if (Cache)
   PreloadCache( Cache );

  /*******************************************************************/
  /* assemble and factorize the BEM matrix                           */
  /*******************************************************************/
  SSS->AssembleBEMMatrix(M);
  M->LUFactorize();

  /*******************************************************************/
  /* now switch off depending on the type of calculation the user    */
  /* requested                                                       */
  /*******************************************************************/
  if (PolFile)
   GetPolarizabilities(SSS, M, Sigma, PolFile);
  if (CapFile)
   GetCapacitanceMatrix(SSS, M, Sigma, CapFile);
  if (CMatrixFile)
   GetCMatrix(SSS, M, Sigma, CMatrixFile);
  if (nEPFiles>0 || PlotFile )
   DoFieldCalculation(SSS, M, Sigma, 
                      PotFile, PhiExt, ConstFieldDirection, 
                      PlotFile, EPFiles, nEPFiles);

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  printf("Thank you for your support.\n");

}
