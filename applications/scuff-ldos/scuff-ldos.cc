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
#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <libhrutil.h>
#include <libscuff.h>
#include <libscuffInternals.h>
#include <libSGJC.h>
#include <BZIntegration.h>
#include "scuff-ldos.h"

#define II cdouble(0.0,1.0)
#define MAXEPFILES 100
#define MAXFREQ    10

using namespace scuff;

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
  /* process remaining command-line arguments ********************/
  /***************************************************************/
  char *GeoFile;
/**/
  cdouble OmegaVals[MAXFREQ];	int nOmegaVals;   int nLambdaVals;
  char *OmegaFile=0;
  char *LambdaFile=0;
  char *OkBFile=0;
/**/
  char *EPFiles[MAXEPFILES];   int nEPFiles=0;
/**/
  char *TransFile=0;
/**/
  bool GroundPlane=false;
  char *HalfSpace=0;
  bool SkipBZIntegration=false;
/**/
  double RelTol=1.0e-2;
  double AbsTol=1.0e-8;
  int MaxEvals=10000;
/**/
  char *FileBase=0;
  bool LDOSOnly=false;
  bool FullTPDGF=false;
/**/
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"geometry",    PA_STRING,  1, 1, (void *)&GeoFile,      0,  ".scuffgeo file"},
//
     {"HalfSpace",   PA_STRING,  1, 1, (void *)&HalfSpace,    0,  "simulate an infinite half-space for z<0 with the given material"},
     {"GroundPlane", PA_BOOL,    1, 1, (void *)&GroundPlane,  0,  "simulate an infinite PEC ground plane at z=0"},
     {"SkipBZIntegration", PA_BOOL,  1, 1, (void *)&SkipBZIntegration, 0,  "bypass BZ integration for analytical half-space/ground-plane calculation"},
//
     {"EPFile",      PA_STRING,  1, MAXEPFILES, (void *)EPFiles, &nEPFiles,  "list of evaluation points"},
//
     {"TransFile",   PA_STRING,  1, 1, (void *)&TransFile,    0,  "list of geometrical transformations"},
//
     {"Omega",       PA_CDOUBLE, 1, MAXFREQ, (void *)OmegaVals,  &nOmegaVals,  "angular frequency"},
     {"Lambda",      PA_CDOUBLE, 1, MAXFREQ, (void *)OmegaVals,  &nLambdaVals, "(free-space) wavelength"},
     {"OmegaFile",   PA_STRING,  1, 1, (void *)&OmegaFile,    0,  "list of omega points "},
     {"LambdaFile",  PA_STRING,  1, 1,       (void *)&LambdaFile,   0,       "list of (free-space) wavelengths"},
     {"OmegakBlochFile",   PA_STRING,  1, 1, (void *)&OkBFile,    0,  "list of (omega,kBloch) points "},
//
     {"RelTol",      PA_DOUBLE,  1, 1, (void *)&RelTol,        0,  "relative tolerance"},
     {"AbsTol",      PA_DOUBLE,  1, 1, (void *)&AbsTol,        0,  "absolute tolerance"},
     {"MaxEvals",    PA_INT,     1, 1, (void *)&MaxEvals,      0,  "maximum number of integrand/summand evaluation"},
//
     {"FileBase",    PA_STRING,  1, 1, (void *)&FileBase,      0,  "base name for output files"},
     {"LDOSOnly",    PA_BOOL,    0, 1, (void *)&LDOSOnly,      0,  "omit DGF components from Brillouin-zone integration"},
     {"FullTPDGF",   PA_BOOL,    0, 1, (void *)&FullTPDGF,     0,  "compute full (bare+scattered) two-point DGF (default is scattering part only)"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFile==0 && (SkipBZIntegration==false) )
   OSUsage(argv[0], VERSION, OSArray,"--geometry option is mandatory");
  if (nEPFiles==0)
   OSUsage(argv[0], VERSION, OSArray,"you must specify at least one --EPFile");
  if (!FileBase)
   FileBase = vstrdup(GetFileBase(GeoFile));
  if (HalfSpace)
   { char *NewFileBase=vstrdup("%s.%s",FileBase,HalfSpace);
     FileBase=NewFileBase;
   };
  if (GroundPlane)
   { char *NewFileBase=vstrdup("%s.GroundPlane",FileBase);
     FileBase=NewFileBase;
   };

  /***************************************************************/
  /* process --Omega, --OmegaFile, --OmegakBlochFile arguments   */
  /***************************************************************/
  HVector *OmegaPoints=0;
  HMatrix *OkBPoints=0;
  if (OkBFile)
   { OkBPoints = new HMatrix(OkBFile);
     if (OkBPoints->ErrMsg)
      ErrExit(OkBPoints->ErrMsg);
   }
  else
   { OmegaPoints=GetOmegaList(OmegaFile, OmegaVals, nOmegaVals, LambdaFile, OmegaVals, nLambdaVals);
     if (!OmegaPoints)
      ErrExit("you must specify at least one frequency");
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (SkipBZIntegration)
   { if (OkBFile)
      ErrExit("--SkipBZIntegration is incompatible with --OmegakBlochFile");
     if ( !GroundPlane && HalfSpace==0 )
      ErrExit("--SkipBZIntegration requires either --HalfSpace or --GroundPlane");
     ProcessHalfSpaceDGFs(OmegaPoints, EPFiles, nEPFiles,
                          HalfSpace, RelTol, AbsTol, MaxEvals);
     exit(1);
   };  

  /***************************************************************/
  /* create the SLDData structure that will be passed to all     */
  /* computational routines                                      */
  /***************************************************************/
  SLDData *Data     = CreateSLDData(GeoFile, TransFile, 
                                    EPFiles, nEPFiles);
  Data->RelTol      = RelTol;
  Data->AbsTol      = AbsTol;
  Data->MaxEvals    = MaxEvals;
  Data->FileBase    = FileBase;
  Data->LDOSOnly    = LDOSOnly;
  Data->ScatteringOnly = !FullTPDGF;
  Data->GroundPlane = GroundPlane;
  Data->HalfSpaceMP = HalfSpace ? new MatProp(HalfSpace) : 0;

  // set LDOSOnly = false if any EPFiles have 6 coordinates
  // (for two-point DGF calculations)
  for(int n=0; n<nEPFiles; n++)
   if (Data->XMatrices[n]->NC==6)
    Data->LDOSOnly=false;

  int LDim = Data->G->LDim;
  if (HalfSpace && LDim!=2)
   OSUsage(argv[0], VERSION, OSArray,
       "--HalfSpace requires a 2D-periodic geometry unless you also say --SkipBZIntegration");

  int NX         = Data->TotalEvalPoints;
  int NFun       = Data->LDOSOnly ? 2 : 38; // # outputs per eval pt
  int FDim       = NX*NFun;
  double *Result = (double *)mallocEC(FDim*sizeof(double));

  /***************************************************************/
  /* now switch off to figure out what to do:                    */
  /*  1. if we have a non-periodic geometry, simply evaluate     */
  /*     the LDOS at each evaluation point at each frequency.    */
  /*                                                             */
  /*  2. if we have a periodic geometry and the user specified   */
  /*      a list of (Omega, kBloch) points, evaluate the LDOS    */
  /*      at each evaluation point at each (Omega, kBloch) point.*/
  /*                                                             */
  /*  3. otherwise (we have a periodic geometry and the user     */
  /*     did not request specific kBloch points), for each       */
  /*     frequency we perform a Brillouin-zone integral to get   */
  /*     the LDOS at each evaluation point.                      */
  /***************************************************************/

  /*--------------------------------------------------------------*/
  /*- non-PBC geometry: just do non-periodic DGF calculations at  */
  /*-                   each requested Omega value                */
  /*--------------------------------------------------------------*/
  if (LDim==0)
   {  
     for(int no=0; no<OmegaPoints->N; no++)
      GetLDOS( (void *)Data, OmegaPoints->GetEntry(no), 0, Result);
   }
  /*--------------------------------------------------------------*/
  /*- PBC structure with specified kBloch points: do a periodic   */
  /*- DGF calculation at each (Omega, kBloch) point               */
  /*--------------------------------------------------------------*/
  else if (OkBPoints)
   {  
     for(int nokb=0; nokb<OkBPoints->NR; nokb++)
      { 
        cdouble Omega=OkBPoints->GetEntry(nokb,0);

        double kBloch[3]={0.0, 0.0, 0.0};
        for(int d=0; d<LDim; d++)
         kBloch[d]=OkBPoints->GetEntryD(nokb,1+d);

        GetLDOS( (void *)Data, Omega, kBloch, Result);
      };
   }
  /*--------------------------------------------------------------*/
  /*- PBC structure without specified kBloch points: perform a    */
  /*- Brillouin-zone integration at each requested Omega value    */
  /*--------------------------------------------------------------*/
  else
   {
     /***************************************************************/
     /* complete the argument structure for Brillouin-zone          */
     /* that we started to initialize above                         */
     /***************************************************************/
     BZIArgs->BZIFunc     = GetLDOS;
     BZIArgs->UserData    = (void *)Data;
     BZIArgs->FDim        = FDim;
     UpdateBZIArgs(BZIArgs, Data->G->RLBasis, Data->G->RLVolume);

     /***************************************************************/
     /***************************************************************/
     /***************************************************************/
     for(int no=0; no<OmegaPoints->N; no++)
      { cdouble Omega=OmegaPoints->GetEntry(no);
        Log("Evaluating Brillouin-zone integral at omega=%s",z2s(Omega));
        GetBZIntegral(BZIArgs, Omega, Result);
        WriteData(Data, Omega, 0, FILETYPE_LDOS, Result, BZIArgs->BZIError);
      };
   };

}
