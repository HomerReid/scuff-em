#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <libhrutil.h>
#include <libscuff.h>
#include <libscuffInternals.h>
#include <libSGJC.h>
#include <BZIntegration.h>
#include "scuff-ldos.h"

#define II cdouble(0.0,1.0)
#define MAXEPFILES 100

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[]) 
{ 
  InstallHRSignalHandler();
  SetLogFileName("scuff-ldos.log");
  Log("scuff-ldos running on %s",GetHostName());

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
  cdouble Omega=0.0;   int nOmega=0;
  char *OmegaFile=0;
  char *OkBFile=0;
/**/
  char *EPFiles[MAXEPFILES];   int nEPFiles=0;
/**/
  bool GroundPlane=false;
  char *HalfSpace=0;
/**/
  double RelTol=1.0e-2;
  int MaxEvals=1000;
/**/
  char *FileBase=0;
  bool LDOSOnly=false;
/**/
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"geometry",    PA_STRING,  1, 1, (void *)&GeoFile,      0,  ".scuffgeo file"},
//
     {"GroundPlane", PA_BOOL,    1, 1, (void *)&GroundPlane,  0,  "simulate an infinite PEC ground plane at z=0"},
     {"HalfSpace",   PA_STRING,  1, 1, (void *)&HalfSpace,    0,  "simulate an infinite half-space for z<0 with the given material"},
//
     {"EPFile",      PA_STRING,  1, MAXEPFILES, (void *)EPFiles, &nEPFiles,  "list of evaluation points"},
//
     {"Omega",       PA_CDOUBLE, 1, 1, (void *)&Omega,  &nOmega,  "angular frequency"},
     {"OmegaFile",   PA_STRING,  1, 1, (void *)&OmegaFile,    0,  "list of omega points "},
     {"OmegakBlochFile",   PA_STRING,  1, 1, (void *)&OkBFile,    0,  "list of (omega,kBloch) points "},
//
     {"RelTol",      PA_DOUBLE,  1, 1, (void *)&RelTol,        0,  "relative tolerance"},
     {"MaxEvals",    PA_INT,     1, 1, (void *)&MaxEvals,      0,  "maximum number of integrand/summand evaluation"},
//
     {"FileBase",    PA_STRING,  1, 1, (void *)&FileBase,      0,  "base name for output files"},
     {"LDOSOnly",    PA_BOOL,    0, 1, (void *)&LDOSOnly,      0,  "omit DGF components from Brillouin-zone integration"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFile==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");
  if (nEPFiles==0)
   OSUsage(argv[0],OSArray,"you must specify at least one --EPFile");
  if (!FileBase)
   FileBase = vstrdup(GetFileBase(GeoFile));
  if (HalfSpace)
   { char *NewFileBase=vstrdup("%s.%s",FileBase,HalfSpace);
     FileBase=NewFileBase;
   };

  /***************************************************************/
  /* process --Omega, --OmegaFile, --OmegakBlochFile arguments   */
  /***************************************************************/
  HVector *OmegaPoints=0;
  HMatrix *OkBPoints=0;
  if ( nOmega!=0 )
   { if (OmegaFile)
      ErrExit("--Omega and --OmegaFile are incompatible");
     if (OkBFile)
      ErrExit("--Omega and --OmegakBlochFile are incompatible");
     OmegaPoints=new HVector(1, LHM_COMPLEX);
     OmegaPoints->SetEntry(0,Omega);
   }
  else if (OmegaFile)
   { if (OkBFile)
      ErrExit("--OmegaFile and --OmegakBlochFile are incompatible");
     OmegaPoints = new HVector(OmegaFile);
     if (OmegaPoints->ErrMsg)
      ErrExit(OmegaPoints->ErrMsg);
   }
  else if (OkBFile)
   { OkBPoints = new HMatrix(OkBFile);
     if (OkBPoints->ErrMsg)
      ErrExit(OkBPoints->ErrMsg);
   }
  else
   ErrExit("you must specify at least one frequency");

  /***************************************************************/
  /* create the SLDData structure that will be passed to all     */
  /* computational routines                                      */
  /***************************************************************/
  SLDData *Data     = CreateSLDData(GeoFile, EPFiles, nEPFiles);
  Data->RelTol      = RelTol;
  Data->MaxEvals    = MaxEvals;
  Data->FileBase    = FileBase;
  Data->LDOSOnly    = LDOSOnly;
  Data->GroundPlane = GroundPlane;

  if (HalfSpace)
   Data->HalfSpaceMP = new MatProp(HalfSpace);

  int LDim = Data->G->LDim;
  if ( HalfSpace && LDim!=2 )
   OSUsage(argv[0],OSArray,"--HalfSpace requires a 2D-periodic geometry");

  int NFun       = Data->LDOSOnly ? 2 : 38; // # outputs per eval pt
  int NX         = Data->TotalEvalPoints; 
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
        Omega=OkBPoints->GetEntry(nokb,0);

        double kBloch[2];
        kBloch[0]=OkBPoints->GetEntryD(nokb,1);
        if (LDim==2)
         kBloch[1]=OkBPoints->GetEntryD(nokb,2);

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
      { Omega=OmegaPoints->GetEntry(no);
        Log("Evaluating Brillouin-zone integral at omega=%s",z2s(Omega));
        GetBZIntegral(BZIArgs, Omega, Result);
        WriteData(Data, Omega, 0, FILETYPE_LDOS, Result, BZIArgs->BZIError);
      };
   };

}
