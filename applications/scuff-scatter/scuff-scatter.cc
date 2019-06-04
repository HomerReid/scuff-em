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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "scuff-scatter.h"

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define MAXPW    10    // max number of plane waves
#define MAXGB    1     // max number of gaussian beams
#define MAXPS    10    // max number of point sources
#define MAXFREQ  10    // max number of frequencies
#define MAXEPF   10    // max number of evaluation-point files
#define MAXFVM   10    // max number of field visualization meshes
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
  InitializeLog(argv[0]);

//
  char *GeoFile=0;
//
  cdouble OmegaVals[MAXFREQ];        int nOmegaVals;
  char *OmegaFile=0;
  cdouble LambdaVals[MAXFREQ];       int nLambdaVals;
  char *LambdaFile=0;
//
  double pwDir[3*MAXPW];             int npwDir;
  cdouble pwPol[3*MAXPW];            int npwPol;
//
  double gbDir[3*MAXGB];             int ngbDir;
  cdouble gbPol[3*MAXGB];            int ngbPol;
  double gbCenter[3*MAXGB];          int ngbCenter;
  double gbWaist[MAXGB];             int ngbWaist;
//
  double psLoc[3*MAXPS];             int npsLoc;
  cdouble psStrength[3*MAXPS];       int npsStrength;
//
  char *IFFile=0;
//
  char *TransFile=0;                 
//
  char *EPFiles[MAXEPF];             int nEPFiles;
  char *FileBase=0;
//
  char *FVMeshes[MAXFVM];            int nFVMeshes;
  char *FVMeshTransFiles[MAXFVM];    int nFVMeshTransFiles;
  char *FVFuncs[MAXFVM];             int nFVFuncs;
  memset(FVMeshTransFiles, 0, MAXFVM*sizeof(char *));
  memset(FVFuncs,          0, MAXFVM*sizeof(char *));
//
  char *PFTFile=0;
  char *OPFTFile=0;
  char *MomentPFTFile=0;
  char *EMTPFTFile=0;
  char *DSIPFTFile = 0;
  double DSIRadius = 0.0;
  int DSIPoints    = 302;
  int DSIPoints2   = 0;
  char *DSIMesh    = 0;
  bool DSIFarField = false;
//
  bool GetRegionPFTs = false;
//
  bool PlotPFTFlux=false;
//
  char *MomentFile=0;
  char *PSDFile=0;
  bool PlotSurfaceCurrents=false;
//
  char *HDF5File=0;
  char *Cache=0;
  char *ReadCache[MAXCACHE];         int nReadCache;
  char *WriteCache=0;
  char *LogLevel=0;
  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { 
     {"geometry",       PA_STRING,  1, 1,       (void *)&GeoFile,    0,             "geometry file\n"},
/**/
     {"Omega",          PA_CDOUBLE, 1, MAXFREQ, (void *)OmegaVals,   &nOmegaVals,   "(angular) frequency"},
     {"OmegaFile",      PA_STRING,  1, 1,       (void *)&OmegaFile,  0,             "file listing angular frequencies"},
     {"Lambda",         PA_CDOUBLE, 1, MAXFREQ, (void *)LambdaVals,  &nLambdaVals,  "wavelength"},
     {"LambdaFile",     PA_STRING,  1, 1,       (void *)&LambdaFile, 0,             "file listing wavelengths\n"},
/**/
     {"pwDirection",    PA_DOUBLE,  3, MAXPW,   (void *)pwDir,       &npwDir,       "plane wave direction"},
     {"pwPolarization", PA_CDOUBLE, 3, MAXPW,   (void *)pwPol,       &npwPol,       "plane wave polarization\n"},
/**/
     {"gbDirection",    PA_DOUBLE,  3, MAXGB,   (void *)gbDir,       &ngbDir,       "gaussian beam direction"},
     {"gbPolarization", PA_CDOUBLE, 3, MAXGB,   (void *)gbPol,       &ngbPol,       "gaussian beam polarization"},
     {"gbCenter",       PA_DOUBLE,  3, MAXGB,   (void *)gbCenter,    &ngbCenter,    "gaussian beam center"},
     {"gbWaist",        PA_DOUBLE,  1, MAXGB,   (void *)gbWaist,     &ngbWaist,     "gaussian beam waist\n"},
/**/
     {"psLocation",     PA_DOUBLE,  3, MAXPS,   (void *)psLoc,       &npsLoc,       "point source location"},
     {"psStrength",     PA_CDOUBLE, 3, MAXPS,   (void *)psStrength,  &npsStrength,  "point source strength\n"},
/**/
     {"IFFile",         PA_STRING,  1, 1,       (void *)&IFFile,     0,             "list of incident fields\n"},
/**/
     {"TransFile",      PA_STRING,  1, 1,       (void *)&TransFile,  0,             "list of geometrical transformations\n"},
/**/
     {"EPFile",         PA_STRING,  1, MAXEPF,  (void *)EPFiles,     &nEPFiles,     "list of evaluation points"},
     {"FileBase",       PA_STRING,  1, 1,       (void *)&FileBase,   0,             "base filename for EPFile output\n\n"},
/**/
     {"FVMesh",         PA_STRING,  1, MAXFVM,  (void *)FVMeshes,    &nFVMeshes,    "field visualization mesh"},
     {"FVMeshTransFile", PA_STRING,  1, MAXFVM,  (void *)FVMeshTransFiles,    &nFVMeshTransFiles,    "list of geometrical transformations applied to FVMesh"},
     {"FVFunctions",    PA_STRING,  1, MAXFVM,  (void *)FVFuncs,     &nFVFuncs,     "comma-separated list of field functions to visualize\n"},
/**/
     {"PFTFile",        PA_STRING,  1, 1,       (void *)&PFTFile,    0,             "name of power, force, and torque output file\n"},
/**/
     {"OPFTFile",       PA_STRING,  1, 1,       (void *)&OPFTFile,   0,             "name of overlap PFT output file"},
     {"MomentPFTFile",  PA_STRING,  1, 1,       (void *)&MomentPFTFile,   0,        "name of multipole-moment PFT output file"},
     {"EMTPFTFile",     PA_STRING,  1, 1,       (void *)&EMTPFTFile, 0,             "name of energy/momentum-transfer PFT output file"},
     {"DSIPFTFile",     PA_STRING,  1, 1,       (void *)&DSIPFTFile, 0,             "name of displaced surface-integral PFT output file"},
     {"DSIMesh",        PA_STRING,  1, 1,       (void *)&DSIMesh,    0,             "mesh file for surface-integral PFT"},
     {"DSIRadius",      PA_DOUBLE,  1, 1,       (void *)&DSIRadius,  0,             "radius of bounding sphere for surface-integral PFT"},
     {"DSIPoints",      PA_INT,     1, 1,       (void *)&DSIPoints,  0,             "number of quadrature points for surface-integral PFT (6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810)"},
     {"DSIPoints2",     PA_INT,     1, 1,       (void *)&DSIPoints2, 0,             "number of quadrature points for DSIPFT second opinion"},
     {"DSIFarField",    PA_BOOL,    0, 1,       (void *)&DSIFarField, 0,            "retain only far-field contributions to DSIPFT\n"},
/**/
     {"GetRegionPFTs",  PA_BOOL,    0, 1,       (void *)&GetRegionPFTs,   0,          "report PFTs for each region"},
     {"PlotPFTFlux",    PA_BOOL,    0, 1,       (void *)&PlotPFTFlux,   0,          "generate plots of spatially-resolved PFT flux\n"},
/**/
     {"MomentFile",     PA_STRING,  1, 1,       (void *)&MomentFile, 0,             "name of dipole moment output file"},
     {"PSDFile",        PA_STRING,  1, 1,       (void *)&PSDFile,    0,             "name of panel source density file"},
     {"PlotSurfaceCurrents", PA_BOOL, 0, 1,     (void *)&PlotSurfaceCurrents,  0,   "generate surface current visualization files\n"},
/**/
     {"HDF5File",       PA_STRING,  1, 1,       (void *)&HDF5File,   0,             "name of HDF5 file for BEM matrix/vector export\n"},
/**/
     {"LogLevel",       PA_STRING,  1, 1,       (void *)&LogLevel,   0,             "none | terse | verbose | verbose2\n"},
/**/
     {"Cache",          PA_STRING,  1, 1,       (void *)&Cache,      0,             "read/write cache"},
     {"ReadCache",      PA_STRING,  1, MAXCACHE,(void *)ReadCache,   &nReadCache,   "read cache"},
     {"WriteCache",     PA_STRING,  1, 1,       (void *)&WriteCache, 0,             "write cache"},
/**/
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  if (GeoFile==0)
   OSUsage(argv[0], VERSION, OSArray, "--geometry option is mandatory");

  if (FileBase==0) 
   FileBase=vstrdup(GetFileBase(GeoFile));

  /*******************************************************************/
  /* process frequency-related options                               */
  /*******************************************************************/
  HVector *OmegaList=GetOmegaList(OmegaFile, OmegaVals, nOmegaVals,
                                  LambdaFile, LambdaVals, nLambdaVals);
  if (OmegaList==0)
   OSUsage(argv[0], VERSION, OSArray, "you must specify at least one frequency");

  /*******************************************************************/
  /* process incident-field-related options to construct the data    */
  /* used to generate the incident field in our scattering problem.  */
  /*******************************************************************/
  if ( npwPol != npwDir )
   ErrExit("numbers of --pwPolarization and --pwDirection options must agree");
  if ( ngbPol != ngbDir || ngbDir!=ngbCenter || ngbCenter!=ngbWaist )
   ErrExit("numbers of --gbPolarization, --gbDirection, --gbCenter, and --gbWaist options must agree ");
  if ( npsLoc!=npsStrength )
   ErrExit("numbers of --psLocation and --psStrength options must agree");

  IncField *IF=0;
  for(int npw=0; npw<npwPol; npw++)
   { PlaneWave *PW=new PlaneWave(pwPol + 3*npw, pwDir + 3*npw);
     PW->Next=IF;
     IF=PW;
   };
  for(int ngb=0; ngb<ngbCenter; ngb++)
   { GaussianBeam *GB=new GaussianBeam(gbCenter + 3*ngb, gbDir + 3*ngb, gbPol + 3*ngb, gbWaist[ngb]);
     GB->Next=IF;
     IF=GB;
   };
  for(int nps=0; nps<npsLoc; nps++)
   { PointSource *PS=new PointSource(psLoc + 3*nps, psStrength + 3*nps);
     PS->Next=IF;
     IF=PS;
   };

  IncFieldList *IFList=0;
  if (IFList!=0 && IF!=0)
   ErrExit("--IFFile is incompatible with other incident-field specifications");
  else if (IFFile!=0 && IF==0)
   IFList = ReadIncFieldList(IFFile);
  else if (IFFile==0 && IF!=0)
   IFList = AddIncFieldToList(IF,const_cast<char *>("Default"));

  /*******************************************************************/
  /* make --PFTFile synonymous with --EMTPFTFile *********************/
  /*******************************************************************/
  if (PFTFile)
   { if (EMTPFTFile)
      Warn("specifying both --PFTFile and --EMTPFTFile is redundant (ignoring --PFTFile)");
     else 
      EMTPFTFile=PFTFile;
   };

  /*******************************************************************/
  /* sanity check to make sure the user specified an incident field  */
  /* if one is required for the outputs the user requested           */
  /*******************************************************************/
  bool NeedIncidentField = (    MomentFile!=0
                             || OPFTFile!=0
                             || MomentPFTFile!=0
                             || EMTPFTFile!=0
                             || DSIPFTFile!=0
                             || nEPFiles>0
                             || nFVMeshes>0
                             || PlotSurfaceCurrents
                           );
  if ( NeedIncidentField && IFList==0 )
   ErrExit("you must specify at least one incident field source");

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  if (nFVMeshTransFiles!=0 && nFVMeshes!=nFVMeshTransFiles)
   ErrExit("--FVMeshTransFile must be specified for every --FVMesh");
  if (nFVFuncs!=0 && nFVMeshes!=nFVFuncs)
   ErrExit("--FVFuncs must be specified for every --FVMesh");

  /*******************************************************************/
  /* PFT options *****************************************************/
  /*******************************************************************/
  PFTOptions MyPFTOpts, *PFTOpts=&MyPFTOpts;
  InitPFTOptions(PFTOpts);
  PFTOpts->DSIMesh       = DSIMesh;
  PFTOpts->DSIRadius     = DSIRadius;
  PFTOpts->DSIPoints     = DSIPoints;
  PFTOpts->DSIFarField   = DSIFarField;
  PFTOpts->GetRegionPFTs = GetRegionPFTs;

  char *DSIPFTFile2 = 0;
  if (DSIPFTFile && DSIPoints2)
   DSIPFTFile2=vstrdup("%s.DSI%i",GetFileBase(DSIPFTFile),DSIPoints2);

  /*******************************************************************/
  /* create the SSData structure containing everything we need to    */
  /* execute scattering calculations                                 */
  /*******************************************************************/
  SSData MySSData, *SSD=&MySSData;

  RWGGeometry *G      = SSD->G   = new RWGGeometry(GeoFile);
  HMatrix *M          = SSD->M   = G->AllocateBEMMatrix();
  HVector *RHS        = SSD->RHS = G->AllocateRHSVector();
  HVector *KN         = SSD->KN  = G->AllocateRHSVector();
  double *kBloch      = SSD->kBloch = 0;
  SSD->IF             = 0;
  SSD->TransformLabel = 0;
  SSD->IFLabel        = 0;
  SSD->FileBase       = FileBase;

  if (LogLevel) G->SetLogLevel(LogLevel);

  /*--------------------------------------------------------------*/
  /*- read the transformation file if one was specified and check */
  /*- that it plays well with the specified geometry file.        */
  /*- note if TransFile==0 then this code snippet still works; in */
  /*- this case the list of GTComplices is initialized to contain */
  /*- a single empty GTComplex and the check automatically passes.*/
  /*--------------------------------------------------------------*/
  GTCList GTCs=ReadTransFile(TransFile);
  int NumTransformations=GTCs.size();
  char *ErrMsg=G->CheckGTCList(GTCs);
  if (ErrMsg)
   ErrExit("file %s: %s",TransFile,ErrMsg);

  /*******************************************************************/
  /* for periodic geometries, all incident field sources that are    */
  /* active at a given time must involve  single incident field      */
  /* source, which must be a plane wave, and the bloch wavevector    */
  /* is extracted from the plane wave direction                      */
  /*******************************************************************/
  double kBlochBuffer[3];
  if (G->LDim>0)
   { if ( npwPol!=1 || ngbCenter!=0 || npsLoc!=0 )
      ErrExit("for extended geometries, the incident field must be a single plane wave");
     kBloch = SSD->kBloch = kBlochBuffer;
   };

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
  /* if we have more than one geometrical transformation,            */
  /* allocate storage for BEM matrix blocks                          */
  /*******************************************************************/
  HMatrix **TBlocks=0, **UBlocks=0;
  int NS=G->NumSurfaces;
  if (NumTransformations>1)
   { int NADB = NS*(NS-1)/2; // number of above-diagonal blocks
     TBlocks  = (HMatrix **)mallocEC(NS*sizeof(HMatrix *));
     UBlocks  = (HMatrix **)mallocEC(NADB*sizeof(HMatrix *));
     for(int ns=0, nb=0; ns<NS; ns++)
      { 
        int nsMate = G->Mate[ns];
        if ( nsMate!=-1 )
         TBlocks[ns] = TBlocks[nsMate];
        else
         { int NBF=G->Surfaces[ns]->NumBFs;
           TBlocks[ns] = new HMatrix(NBF, NBF, M->RealComplex);
         };

        for(int nsp=ns+1; nsp<NS; nsp++, nb++)
         { int NBF=G->Surfaces[ns]->NumBFs;
           int NBFp=G->Surfaces[nsp]->NumBFs;
           UBlocks[nb] = new HMatrix(NBF, NBFp, M->RealComplex);
         };
      };
   };

  /*******************************************************************/
  /* loop over frequencies *******************************************/
  /*******************************************************************/
  for(int nFreq=0; nFreq<OmegaList->N; nFreq++)
   { 
     cdouble Omega = OmegaList->GetEntry(nFreq);
     SSD->Omega    = Omega;

     char OmegaStr[MAXSTR];
     z2s(Omega, OmegaStr);
     Log("Working at frequency %s...",OmegaStr);

     /*******************************************************************/
     /* for periodic geometries, extract bloch wavevector from incident field */
     /*******************************************************************/
     if (G->LDim>0)
      { cdouble EpsExterior, MuExterior;
        G->RegionMPs[0]->GetEpsMu(Omega, &EpsExterior, &MuExterior);
        double kExterior = real( csqrt2(EpsExterior*MuExterior) * Omega );
        kBloch[0] = kExterior*pwDir[0];
        kBloch[1] = kExterior*pwDir[1];
        kBloch[2] = 0.0;
      };

     /*******************************************************************/
     /* if we have more than one transformation, pre-assemble diagonal  */
     /* matrix blocks at this frequency; otherwise just assemble the    */
     /* whole matrix                                                    */
     /*******************************************************************/
     if (NumTransformations==1)
      G->AssembleBEMMatrix(Omega, kBloch, M);
     else
      for(int ns=0; ns<G->NumSurfaces; ns++)
       if (G->Mate[ns]==-1)
        G->AssembleBEMMatrixBlock(ns, ns, Omega, kBloch, TBlocks[ns]);

     /*******************************************************************/
     /* dump the scuff cache to a cache storage file if requested. note */
     /* we do this only once per execution of the program, after the    */
     /* assembly of the diagonal BEM matrix blocks at first frequency,  */
     /* since at that point all cache elements that are to be computed  */
     /* will have been computed and the cache will not grow any further */
     /* for the rest of the program run.                                */
     /*******************************************************************/
     if (WriteCache)
      { StoreCache( WriteCache );
        WriteCache=0;       
      };

     /*******************************************************************/
     /*******************************************************************/
     /*******************************************************************/
     for(int nt=0; nt<NumTransformations; nt++)
      {
        char TransformStr[100]="";
        if (TransFile)
         { G->Transform(GTCs[nt]);
           SSD->TransformLabel=GTCs[nt]->Tag;
           Log("Working at transformation %s...",SSD->TransformLabel);
           snprintf(TransformStr,100,"_%s",SSD->TransformLabel);
         };

        /*******************************************************************/
        /* assemble and insert off-diagonal blocks as necessary ************/
        /*******************************************************************/
        if (NumTransformations>1)
         { for(int ns=0, nb=0; ns<G->NumSurfaces; ns++)
            for(int nsp=ns+1; nsp<G->NumSurfaces; nsp++, nb++)
             G->AssembleBEMMatrixBlock(ns, nsp, Omega, kBloch, UBlocks[nb]);

           for(int ns=0, nb=0; ns<G->NumSurfaces; ns++)
            { int RowOffset=G->BFIndexOffset[ns];
              M->InsertBlock(TBlocks[ns], RowOffset, RowOffset);
              for(int nsp=ns+1; nsp<G->NumSurfaces; nsp++, nb++)
               { int ColOffset=G->BFIndexOffset[nsp];
                 M->InsertBlock(UBlocks[nb], RowOffset, ColOffset);
                 M->InsertBlockAdjoint(UBlocks[nb], ColOffset, RowOffset);
               };
            };
         };

        /*******************************************************************/
        /* export BEM matrix to a binary .hdf5 file if that was requested  */
        /*******************************************************************/
        if (HDF5Context)
         M->ExportToHDF5(HDF5Context,"M_%s%s",OmegaStr,TransformStr);

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
        /* loop over incident fields                                   */
        /***************************************************************/
        for(int nIF=0; nIF<IFList->NumIFs; nIF++)
         { 
           IF = SSD->IF = IFList->IFs[nIF];
           SSD->IFLabel = IFFile ? IFList->Labels[nIF] : 0;
           if (SSD->IFLabel)
            Log("  Processing incident field %s...",SSD->IFLabel);

           char IFStr[100]="";
           if (SSD->IFLabel)
            snprintf(IFStr,100,"_%s",SSD->IFLabel);
   
           /***************************************************************/
           /* assemble RHS vector and solve BEM system*********************/
           /***************************************************************/
           Log("  Assembling RHS vector...");
           G->AssembleRHSVector(Omega, kBloch, IF, KN);
           RHS->Copy(KN); // copy RHS vector for later 
           Log("  Solving the BEM system...");
           M->LUSolve(KN);
   
           if (HDF5Context)
            { RHS->ExportToHDF5(HDF5Context,"RHS_%s%s%s",OmegaStr,TransformStr,IFStr);
              KN->ExportToHDF5(HDF5Context,"KN_%s%s%s",OmegaStr,TransformStr,IFStr);
            };
   
           /***************************************************************/
           /* now process all requested outputs                           */
           /***************************************************************/
   
           /*--------------------------------------------------------------*/
           /*- power, force, torque by various methods --------------------*/
           /*--------------------------------------------------------------*/
           if (OPFTFile)
            WritePFTFile(SSD, PFTOpts, SCUFF_PFT_OVERLAP, PlotPFTFlux, OPFTFile);
      
           if (MomentPFTFile)
            WritePFTFile(SSD, PFTOpts, SCUFF_PFT_MOMENTS, PlotPFTFlux, MomentPFTFile);
      
           if (DSIPFTFile)
            { PFTOpts->DSIPoints=DSIPoints;
              WritePFTFile(SSD, PFTOpts, SCUFF_PFT_DSI, PlotPFTFlux, DSIPFTFile);
            };
      
           if (DSIPFTFile2)
            { PFTOpts->DSIPoints=DSIPoints2;
              WritePFTFile(SSD, PFTOpts, SCUFF_PFT_DSI, PlotPFTFlux, DSIPFTFile2);
            };
      
           if (EMTPFTFile)
            WritePFTFile(SSD, PFTOpts, SCUFF_PFT_EMT, PlotPFTFlux, EMTPFTFile);
      
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
            G->PlotSurfaceCurrents(KN, Omega, SSD->kBloch, "%s.pp", FileBase);
      
           /*--------------------------------------------------------------*/
           /*- field visualization meshes ---------------------------------*/
           /*--------------------------------------------------------------*/
           int nfm;
           for(nfm=0; nfm<nFVMeshes; nfm++)
            VisualizeFields(SSD, FVMeshes[nfm], FVMeshTransFiles[nfm], FVFuncs[nfm]);

         }; // for(int nIF=0; nIF<IFList->NumIFs; nIF++
      
        /*******************************************************************/
        /*******************************************************************/
        /*******************************************************************/
        G->UnTransform();

      }; // for(int nt=0; nt<NumTransformations; nt++)

   }; //  for(nFreq=0; nFreq<NumFreqs; nFreqs++)
  
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (HDF5Context)
   HMatrix::CloseHDF5Context(HDF5Context);
  printf("Thank you for your support.\n");
   
}
