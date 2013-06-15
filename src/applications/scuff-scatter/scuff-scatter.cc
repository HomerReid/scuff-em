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
 *  scuff-scatter < myOptions 
 * 
 * is equivalent to 
 * 
 *  scuff-scatter --option1 value1 --option value2.
 * 
 * (if any options are specified both on standard input and 
 *  on the command line, the values given on the command line take
 *  precedence.)
 * 
 *       -------------------------------------------------
 * 
 * a. options describing the scatterer
 * 
 *     --geometry MyGeometry.scuffgeo
 * 
 *       -------------------------------------------------
 * 
 * b. options specifying the incident field: 
 * 
 *     --pwPolarization Ex Ey Ez 
 *     --pwDirection    Nx Ny Nz 
 * 
 *          incident field is a plane wave with E-field 
 *          polarization (Ex,Ey,Ez) and propagating in the
 *          (Nx,Ny,Nz) direction.
 * 
 *     --gbCenter xx yy zz
 *     --gbDirection Nx Ny Nz
 *     --gbWaist W
 * 
 *          incident field is a gaussian beam, traveling in 
 *          the (Nx, Ny, Nz) direction, with beam center point
 *          (xx,yy,zz) and beam waist W. 
 * 
 *     --psLocation    xx yy zz
 *     --psStrength    Px Py Pz
 * 
 *          incident field is the field of a unit-strength point 
 *          electric dipole source at coordinates (xx,yy,zz) and
 *          strength (dipole moment) (Px,Py,Pz)
 *          
 *          (note that --psLocation and --psStrength may be 
 *           specified more than once (up to 10 times) to 
 *           represent multiple point sources.) 
 * 
 *       -------------------------------------------------
 * 
 * c. options specifying the output: 
 * 
 *     --EPFile MyEPFile
 * 
 *         a file containing a list of points at which the scattered 
 *         field is to be evaluated.
 * 
 *         each line of EPFile should contain 3 numbers, the cartesian
 *         components of the scattered field. (blank lines and comments,
 *         i.e. lines beginning with a '#', are skipped)
 * 
 *         (Note that --epfile may be specified more than once.)
 * 
 *     --PFTFile MyPFTFile
 * 
 *        Requests that the power, force, and torque delivered to 
 *        the scattering objects from the incident field be written
 *        to the file MyPFTFile.
 * 
 *     --PSDFile MyPSDFile
 * 
 *        Requests that panel source densities be written to the file 
 *        MyPSDFile.
 * 
 *     --FluxMesh  MyFluxMesh.msh
 *        
 *        This specification will cause field visualization data
 *        to be written to the file MyFluxMesh.pp.
 * 
 *     --MomentFile MyMoments.dat
 * 
 *        if the --MomentFile option is present, then the electric
 *        and magnetic dipole moments induced on each object in the 
 *        geometry will be computed and reported in the file MyMoments.dat. 
 * 
 *     --PlotSurfaceCurrents
 *
 *        If this flag is specified, scuff-scatter will generate 
 *        GMSH graphical data files that you can open in GMSH 
 *        to visualize the electric and magnetic surface currents 
 *        induced by the incident fields on the scattering objects.
 *
 *       -------------------------------------------------
 * 
 * d. options describing the frequency
 * 
 *     --Omega xx    
 *
 *         Specify the frequency at which to conduct the 
 *         simulation. 
 *
 *         The specified frequency may be complex, i.e.
 *         valid specifications include 
 *
 *           --Omega 3.4i
 *           --Omega 0.2e3-4.5e2I
 *
 *         Note that --Omega may be specified more than 
 *         once.
 * 
 *     --OmegaFile MyOmegaFile
 * 
 *         Specify a file containing a list of frequencies
 *         at which to conduct simulations.               
 * 
 *       -------------------------------------------------
 * 
 * e. other options 
 * 
 *     --nThread xx   (use xx computational threads)
 *     --ExportMatrix (export the BEM matrix to an .hdf5 data file)
 * 
 * --------------------------------------------------------------
 *
 * if this program terminates successfully, the following output 
 * will have been generated: 
 * 
 *  (a) if the --EPFile option was specified, you will have
 *      files named 'MyGeometry.scattered' and 'MyGeometry.total'
 *      (where MyGeometry.rwggeo was the option you passed to 
 *       --geometry) tabulating values of the scattered and total 
 *       E and H fields at each point specified in the MyEPFile 
 *      (where MyEPFile was the option you passed to --EPFile).
 *
 *      each of these files will contain one line for each line in 
 *      MyEPFile. each line will contain 15 numbers: the 3 coordinates
 *      of the evaluation point, the 3 cartesian components of the 
 *      E field at that point (real and imaginary components), and 
 *      the 3 cartesian components of the scattered H field at that 
 *      point (real and imaginary components.)
 * 
 *  (b) if one or more --FluxMesh options was specified, you will 
 *      have files named 'MyFluxMesh.pp,' 'MySecondFluxMesh.pp,'
 *      etc. that may be opened in GMSH.
 *      
 *      each of these files will contain several data sets:
 *      
 *       -- poynting flux (of both the scattered and total fields)
 *          plotted as a scalar field over the meshed surface
 *       -- scattered and total E and H fields (real and imag parts)
 *          plotted as normalized arrows at the centroid of each
 *          triangle in the meshed surface
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
  char *PFTFile=0;
  char *PSDFile=0;
  char *MomentFile=0;
  char *FluxMeshes[MAXFM];           int nFluxMeshes;
  int PlotSurfaceCurrents=0;
  int nThread=0;
  int ExportMatrix=0;
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
     {"PFTFile",        PA_STRING,  1, 1,       (void *)&PFTFile,    0,             "name of power/force/torque output file"},
     {"MomentFile",     PA_STRING,  1, 1,       (void *)&MomentFile, 0,             "name of dipole moment output file"},
     {"PSDFile",        PA_STRING,  1, 1,       (void *)&PSDFile,    0,             "name of panel source density file"},
     {"PlotSurfaceCurrents", PA_BOOL, 0, 1,     (void *)&PlotSurfaceCurrents,  0,   "generate surface current visualization files"},
/**/
     {"Cache",          PA_STRING,  1, 1,       (void *)&Cache,      0,             "read/write cache"},
     {"ReadCache",      PA_STRING,  1, MAXCACHE,(void *)ReadCache,   &nReadCache,   "read cache"},
     {"WriteCache",     PA_STRING,  1, 1,       (void *)&WriteCache, 0,             "write cache"},
/**/
     {"nThread",        PA_INT,     1, 1,       (void *)&nThread,    0,             "number of CPU threads to use"},
     {"ExportMatrix",   PA_BOOL,    0, 1,       (void *)&ExportMatrix, 0,           "export BEM matrix to file"},
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
  if ( (MomentFile!=0 || PFTFile!=0 || nEPFiles>0 || nFluxMeshes>0 || PlotSurfaceCurrents) && IFDList==0 )
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
  if (G->NumLatticeBasisVectors>0)
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
     if ( G->NumLatticeBasisVectors==0 )
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
     /* export BEM matrix to a binary file if that was requested        */
     /*******************************************************************/
     if (ExportMatrix)
      { void *pCC=HMatrix::OpenMATLABContext("%s_%s",GeoFileBase,OmegaStr);
        M->ExportToMATLAB(pCC,"M");
        HMatrix::CloseMATLABContext(pCC);
      };

     /*******************************************************************/
     /* if the user requested no output options (for example, if she   **/
     /* just wanted to export the matrix to a binary file), don't      **/
     /* bother LU-factorizing the matrix or assembling the RHS vector. **/
     /*******************************************************************/
     if ( PFTFile==0 && PSDFile==0 && MomentFile==0 && nEPFiles==0 && nFluxMeshes==0 && PlotSurfaceCurrents==0 )
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

     /***************************************************************/
     /* now process all requested outputs                           */
     /***************************************************************/
     SSD->Omega=Omega;

     /*--------------------------------------------------------------*/
     /*- scattered and absorbed power -------------------------------*/
     /*--------------------------------------------------------------*/
     if (PFTFile)
      WritePFTFile(SSD, PFTFile);

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
  printf("Thank you for your support.\n");

}
