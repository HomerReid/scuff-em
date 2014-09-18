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
 * scuff-scatter.cc  -- a standalone code within the scuff-EM suite for 
 *                   -- solving problems involving the scattering of 
 *                   -- electromagnetic radiation from an arbitrary 
 *                   -- compact object
 *
 * homer reid        -- 10/2006 -- 1/2012
 *
 * documentation at: http://homerreid.com/scuff-em/
 *
 */
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>

#include "scuff-scatter.h"

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define MAXPW    1     // max number of plane waves
#define MAXGB    1     // max number of gaussian beams
#define MAXPS    10    // max number of point sources
#define MAXFREQ  10    // max number of frequencies 
#define MAXEPF   10    // max number of evaluation-point files
#define MAXFM    10    // max number of flux meshes
#define MAXCACHE 10    // max number of cache files for preload

#define MAXSTR   1000

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{

  /***************************************************************/
  /* process options *********************************************/
  /***************************************************************/
  InstallHRSignalHandler();
  char *GeoFile=0;
  double pwDir[3*MAXPW];             int npwDir;
  cdouble pwPol[3*MAXPW];            int npwPol;
  double gbDir[3*MAXGB];             int ngbDir;
  cdouble gbPol[3*MAXGB];            int ngbPol;
  double gbCenter[3*MAXGB];          int ngbCenter;
  double gbWaist[MAXGB];             int ngbWaist;
  double psLoc[3*MAXPS];             int npsLoc;
  cdouble psStrength[3*MAXPS];       int npsStrength;
  cdouble OmegaVals[MAXFREQ];        int nOmegaVals;
  char *OmegaFile;                   int nOmegaFiles;
  char *EPFiles[MAXEPF];             int nEPFiles;
  char *OPFTFile=0;
  char *EPPFTFile=0;
  char *SIPFTFile=0;
  double SIRadius = 100.0;
  int SIPoints = 11;
  char *SIMeshFile=0;
  char *PSDFile=0;
  char *MomentFile=0;
  char *FluxMeshes[MAXFM];           int nFluxMeshes;
  int PlotSurfaceCurrents=0;
  int nThread=0;
  char *HDF5File=0;
  char *Cache=0;
  char *ReadCache[MAXCACHE];         int nReadCache;
  char *WriteCache=0;
  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { 
     {"geometry",       PA_STRING,  1, 1,       (void *)&GeoFile,    0,             "geometry file"},
/**/
     {"Omega",          PA_CDOUBLE, 1, MAXFREQ, (void *)OmegaVals,   &nOmegaVals,   "(angular) frequency"},
     {"OmegaFile",      PA_STRING,  1, 1,       (void *)&OmegaFile,  &nOmegaFiles,  "list of (angular) frequencies"},
/**/
     {"pwDirection",    PA_DOUBLE,  3, MAXPW,   (void *)pwDir,       &npwDir,       "plane wave direction"},
     {"pwPolarization", PA_CDOUBLE, 3, MAXPW,   (void *)pwPol,       &npwPol,       "plane wave polarization"},
/**/
     {"gbDirection",    PA_DOUBLE,  3, MAXGB,   (void *)gbDir,       &ngbDir,       "gaussian beam direction"},
     {"gbPolarization", PA_CDOUBLE, 3, MAXGB,   (void *)gbPol,       &ngbPol,       "gaussian beam polarization"},
     {"gbCenter",       PA_DOUBLE,  3, MAXGB,   (void *)gbCenter,    &ngbCenter,    "gaussian beam center"},
     {"gbWaist",        PA_DOUBLE,  1, MAXGB,   (void *)gbWaist,     &ngbWaist,     "gaussian beam waist"},
/**/
     {"psLocation",     PA_DOUBLE,  3, MAXPS,   (void *)psLoc,       &npsLoc,       "point source location"},
     {"psStrength",     PA_CDOUBLE, 3, MAXPS,   (void *)psStrength,  &npsStrength,  "point source strength"},
/**/
     {"EPFile",         PA_STRING,  1, MAXEPF,  (void *)EPFiles,     &nEPFiles,     "list of evaluation points"},
     {"FluxMesh",       PA_STRING,  1, MAXFM,   (void *)FluxMeshes,  &nFluxMeshes,  "flux mesh"},
/**/
     {"OPFTFile",       PA_STRING,  1, 1,       (void *)&OPFTFile,   0,             "name of overlap PFT output file"},
/**/
     {"EPPFTFile",      PA_STRING,  1, 1,       (void *)&EPPFTFile,  0,             "name of equivalence-principle PFT file"},
     {"PFTFile",        PA_STRING,  1, 1,       (void *)&OPFTFile,   0,             "name of equivalence-principle PFT file"},
/**/
     {"SIPFTFile",      PA_STRING,  1, 1,       (void *)&SIPFTFile,  0,             "name of surface-integral PFT output file"},
     {"SIRadius",       PA_DOUBLE,  1, 1,       (void *)&SIRadius,   0,             "radius of bounding sphere for surface-integral PFT"},
     {"SIPoints",       PA_INT,     1, 1,       (void *)&SIPoints,   0,             "number of quadrature points for surface-integral PFT"},
     {"SIMesh",         PA_STRING,  1, 1,       (void *)&SIMeshFile, 0,             "mesh file for surface-integral PFT"},
/**/
     {"MomentFile",     PA_STRING,  1, 1,       (void *)&MomentFile, 0,             "name of dipole moment output file"},
     {"PSDFile",        PA_STRING,  1, 1,       (void *)&PSDFile,    0,             "name of panel source density file"},
     {"PlotSurfaceCurrents", PA_BOOL, 0, 1,     (void *)&PlotSurfaceCurrents,  0,   "generate surface current visualization files"},
/**/
     {"Cache",          PA_STRING,  1, 1,       (void *)&Cache,      0,             "read/write cache"},
     {"ReadCache",      PA_STRING,  1, MAXCACHE,(void *)ReadCache,   &nReadCache,   "read cache"},
     {"WriteCache",     PA_STRING,  1, 1,       (void *)&WriteCache, 0,             "write cache"},
/**/
     {"nThread",        PA_INT,     1, 1,       (void *)&nThread,    0,             "number of CPU threads to use"},
     {"HDF5File",       PA_STRING,  0, 1,       (void *)&HDF5File,   0,             "name of HDF5 file for BEM matrix/vector export"},
/**/
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  if (GeoFile==0)
   OSUsage(argv[0], OSArray, "--geometry option is mandatory");

  if (nThread!=0)
   SetNumThreads(nThread);

  /*******************************************************************/
  /* process frequency-related options to construct a list of        */
  /* frequencies at which to run calculations                        */
  /*******************************************************************/
  HVector *OmegaList=0, *OmegaList0;
  int nFreq, nOV, NumFreqs=0;
  if (nOmegaFiles==1) // first process --OmegaFile option if present
   { 
     OmegaList=new HVector(OmegaFile,LHM_TEXT);
     if (OmegaList->ErrMsg)
      ErrExit(OmegaList->ErrMsg);
     NumFreqs=OmegaList->N;
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
   };

  if ( !OmegaList || OmegaList->N==0)
   OSUsage(argv[0], OSArray, "you must specify at least one frequency");

  /*******************************************************************/
  /* process incident-field-related options to construct the data    */
  /* used to generate the incident field in our scattering problem   */
  /*******************************************************************/
  if ( npwPol != npwDir )
   ErrExit("numbers of --pwPolarization and --pwDirection options must agree");
  if ( ngbPol != ngbDir || ngbDir!=ngbCenter || ngbCenter!=ngbWaist )
   ErrExit("numbers of --gbPolarization, --gbDirection, --gbCenter, and --gbWaist options must agree ");
  if ( npsLoc!=npsStrength )
   ErrExit("numbers of --psLocation and --psStrength options must agree");

  IncField *IFDList=0, *IFD;
  int npw, ngb, nps;
  for(npw=0; npw<npwPol; npw++)
   { IFD=new PlaneWave(pwPol + 3*npw, pwDir + 3*npw);
     IFD->Next=IFDList;
     IFDList=IFD;
   };
  for(ngb=0; ngb<ngbCenter; ngb++)
   { IFD=new GaussianBeam(gbCenter + 3*ngb, gbDir + 3*ngb, gbPol + 3*ngb, gbWaist[ngb]);
     IFD->Next=IFDList;
     IFDList=IFD;
   };
  for(nps=0; nps<npsLoc; nps++)
   { IFD=new PointSource(psLoc + 3*nps, psStrength + 3*nps);
     IFD->Next=IFDList;
     IFDList=IFD;
   };

  /*******************************************************************/
  /* sanity check to make sure the user specified an incident field  */
  /* if one is required for the outputs the user requested           */
  /*******************************************************************/
  bool NeedIncidentField = (    MomentFile!=0
                             || OPFTFile!=0
                             || EPPFTFile!=0
                             || SIPFTFile!=0
                             || nEPFiles>0
                             || nFluxMeshes>0
                             || PlotSurfaceCurrents
                           );
  if ( NeedIncidentField && IFDList==0 )
   ErrExit("you must specify at least one incident field source");

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  SetLogFileName("scuff-scatter.log");
  Log("scuff-scatter running on %s",GetHostName());

  /*******************************************************************/
  /* create the SSData structure containing everything we need to    */
  /* execute scattering calculations                                 */
  /*******************************************************************/
  SSData MySSData, *SSD=&MySSData;

  RWGGeometry *G = SSD->G = new RWGGeometry(GeoFile);
  G->SetLogLevel(SCUFF_VERBOSELOGGING);
  HMatrix *M = SSD->M =G->AllocateBEMMatrix();
  SSD->RHS = G->AllocateRHSVector();
  HVector *KN = SSD->KN =G->AllocateRHSVector();
  SSD->IF=IFDList;

  char GeoFileBase[MAXSTR];
  strncpy(GeoFileBase, GetFileBase(GeoFile), MAXSTR);

  /*******************************************************************/
  /* sanity check: for now (20120924), calculations involving        */
  /* extended geometries must have only a single incident field      */
  /* source, which must be a plane wave, and the bloch wavevector    */
  /* is extracted from the plane wave direction                      */
  /*******************************************************************/
  double kBlochBuffer[3];
  if (G->LDim>0)
   { if ( npwPol!=1 || ngbCenter!=0 || npsLoc!=0 )
      ErrExit("for extended geometries, the incident field must be a single plane wave");
     SSD->kBloch = kBlochBuffer;
   }
  else
   SSD->kBloch=0;

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
  /*******************************************************************/
  /*******************************************************************/
  void *HDF5Context=0;
  if (HDF5File)
   HDF5Context=HMatrix::OpenHDF5Context(HDF5File);

  /*******************************************************************/
  /* loop over frequencies *******************************************/
  /*******************************************************************/
  char OmegaStr[MAXSTR];
  cdouble Omega;
  cdouble Eps, Mu;
  for(nFreq=0; nFreq<NumFreqs; nFreq++)
   { 
     Omega = OmegaList->GetEntry(nFreq);
     z2s(Omega, OmegaStr);
     Log("Working at frequency %s...",OmegaStr);

     /*******************************************************************/
     /* assemble the BEM matrix at this frequency                       */
     /*******************************************************************/
     if ( G->LDim==0 )
      G->AssembleBEMMatrix(Omega, M);
     else
      { cdouble EpsExterior, MuExterior;
        G->RegionMPs[0]->GetEpsMu(Omega, &EpsExterior, &MuExterior);
        double kExterior = real( csqrt2(EpsExterior*MuExterior) * Omega );
        SSD->kBloch[0] = kExterior*pwDir[0];
        SSD->kBloch[1] = kExterior*pwDir[1];
        SSD->kBloch[2] = 0.0;
        G->AssembleBEMMatrix(Omega, SSD->kBloch, M);
      };

     /*******************************************************************/
     /* dump the scuff cache to a cache storage file if requested. note */
     /* we do this only once per execution of the program, after the    */
     /* assembly of the BEM matrix at the first frequency, since at that*/
     /* point all cache elements that are to be computed will have been */
     /* computed and the cache will not grow any further for the rest   */
     /* of the program run.                                             */
     /*******************************************************************/
     if (WriteCache)
      { StoreCache( WriteCache );
        WriteCache=0;       
      };

     /*******************************************************************/
     /* export BEM matrix to a binary .hdf5 file if that was requested  */
     /*******************************************************************/
     if (HDF5Context)
      M->ExportToHDF5(HDF5Context,"M_%s",OmegaStr);

     /*******************************************************************/
     /* if the user requested no output options (for example, if she   **/
     /* just wanted to export the matrix to a binary file), don't      **/
     /* bother LU-factorizing the matrix or assembling the RHS vector. **/
     /*******************************************************************/
     if ( !NeedIncidentField )
      continue;

     /*******************************************************************/
     /* LU-factorize the BEM matrix to prepare for solving scattering   */
     /* problems                                                        */
     /*******************************************************************/
     Log("  LU-factorizing BEM matrix...");
     M->LUFactorize();

     /***************************************************************/
     /* set up the incident field profile and assemble the RHS vector */
     /***************************************************************/
     Log("  Assembling the RHS vector..."); 
     G->AssembleRHSVector(Omega, IFDList, KN);
     SSD->RHS->Copy(SSD->KN); // copy RHS vector for later 

     /***************************************************************/
     /* solve the BEM system*****************************************/
     /***************************************************************/
     Log("  Solving the BEM system...");
     M->LUSolve(KN);

     if (HDF5Context)
      { SSD->RHS->ExportToHDF5(HDF5Context,"RHS_%s",OmegaStr);
        SSD->KN->ExportToHDF5(HDF5Context,"KN_%s",OmegaStr);
      };

     /***************************************************************/
     /* now process all requested outputs                           */
     /***************************************************************/
     SSD->Omega=Omega;

     /*--------------------------------------------------------------*/
     /*- overlap PFT ------------------------------------------------*/
     /*--------------------------------------------------------------*/
     if (OPFTFile)
      WriteOPFTFile(SSD, OPFTFile);

     /*--------------------------------------------------------------*/
     /*- equivalence-principle PFT ----------------------------------*/
     /*--------------------------------------------------------------*/
     if (EPPFTFile)
      WriteEPPFTFile(SSD, EPPFTFile);

     /*--------------------------------------------------------------*/
     /*- surface-integral PFT           -----------------------------*/
     /*--------------------------------------------------------------*/
     if (SIPFTFile)
      WriteSIPFTFile(SSD, SIPFTFile, SIRadius, SIPoints, SIMeshFile);

     /*--------------------------------------------------------------*/
     /*- panel source densities -------------------------------------*/
     /*--------------------------------------------------------------*/
     if (PSDFile)
      WritePSDFile(SSD, PSDFile);
 
     /*--------------------------------------------------------------*/
     /*- scattered fields at user-specified points ------------------*/
     /*--------------------------------------------------------------*/
     int nepf;
     for(nepf=0; nepf<nEPFiles; nepf++)
      ProcessEPFile(SSD, EPFiles[nepf]);

     /*--------------------------------------------------------------*/
     /*- induced dipole moments       -------------------------------*/
     /*--------------------------------------------------------------*/
     if (MomentFile)
      GetMoments(SSD, MomentFile);

     /*--------------------------------------------------------------*/
     /*- surface current visualization-------------------------------*/
     /*--------------------------------------------------------------*/
     if (PlotSurfaceCurrents)
      G->PlotSurfaceCurrents(KN, Omega, "%s.%s.pp",GetFileBase(GeoFile),z2s(Omega));

     /*--------------------------------------------------------------*/
     /*- flux meshes ------------------------------------------------*/
     /*--------------------------------------------------------------*/
     int nfm;
     for(nfm=0; nfm<nFluxMeshes; nfm++)
      CreateFluxPlot(SSD, FluxMeshes[nfm]);

   }; //  for(nFreq=0; nFreq<NumFreqs; nFreqs++)

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (HDF5Context)
   HMatrix::CloseHDF5Context(HDF5Context);
  printf("Thank you for your support.\n");

}
