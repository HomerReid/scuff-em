#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#define HAVE_READLINE
#include <readline/readline.h>
#include <readline/history.h>

#include <libhrutil.h>
#include <libscuff.h>
#include <libscuffInternals.h>
#include <libSGJC.h>
#include "scuff-ldos.h"

#define II cdouble(0.0,1.0)

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
  /* process command-line arguments ******************************/
  /***************************************************************/
  char *GeoFile;
/**/
  cdouble Omega=0.0;   int nOmega=0;
  char *OmegaFile=0;
  char *OkBFile=0;
/**/
  char *EPFile=0;
/**/
  bool BZSymmetry=false;
  char *BZIString=0;
  double RelTol=1.0e-2;
  int MaxEvals=1000;
/**/
  char *HalfSpace=0;
/**/
  char *FileBase=0;
  bool LDOSPlus=false;
/**/
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"geometry",    PA_STRING,  1, 1, (void *)&GeoFile,      0,  ".scuffgeo file"},
//
     {"HalfSpace",   PA_STRING,  1, 1, (void *)&HalfSpace,    0,  "simulate an infinite half-space for z<0 with the given material"},
//
     {"EPFile",      PA_STRING,  1, 1, (void *)&EPFile,       0,  "list of evaluation points"},
//
     {"Omega",       PA_CDOUBLE, 1, 1, (void *)&Omega,  &nOmega,  "angular frequency"},
     {"OmegaFile",   PA_STRING,  1, 1, (void *)&OmegaFile,    0,  "list of omega points "},
//
     {"BZSymmetry",  PA_BOOL,    0, 1, (void *)&BZSymmetry,    0,  "assume BZ integrand is xy-symmetric"},
     {"BZIMethod",   PA_STRING,  1, 1, (void *)&BZIString,    0,  "Brillouin-zone integration method [DCUTRI | adaptive]"},
     {"RelTol",      PA_DOUBLE,  1, 1, (void *)&RelTol,       0,  "relative tolerance for Brillouin-zone integration"},
     {"MaxEvals",    PA_INT,     1, 1, (void *)&MaxEvals,     0,  "maximum number of Brillouin-zone samples"},
     {"OmegakBlochFile", PA_STRING,  1, 1, (void *)&OkBFile,  0,  "list of (omega, kx, ky) values"},
//
     {"FileBase",    PA_STRING,  1, 1, (void *)&FileBase,      0,  "base name for output files"},
     {"LDOSPlus",    PA_BOOL,    0, 1, (void *)&LDOSPlus,      0,  "LDOS plus"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFile==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");
  if (EPFile==0)
   OSUsage(argv[0],OSArray,"--EPFile option is mandatory");
  if (!FileBase)
   FileBase=vstrdup(GetFileBase(GeoFile));

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
     if (BZSymmetry)
      ErrExit("--BZSymmetry is incompatible with --OmegakBlochFile");
     if (BZIString)
      ErrExit("--BZIMethod is incompatible with --OmegakBlochFile");
   }
  else
   ErrExit("you must specify at least one frequency");

  /***************************************************************/
  /* create the SLDData structure that will be passed to all   */
  /* computational routines                                      */
  /***************************************************************/
  SLDData *Data     = CreateSLDData(GeoFile, EPFile);
  Data->RelTol      = RelTol;
  Data->MaxEvals    = MaxEvals;
  Data->FileBase    = FileBase;
  Data->LDOSOnly    = !LDOSPlus;
  Data->BZSymmetry  = BZSymmetry;

  if (HalfSpace)
   Data->HalfSpaceMP = new MatProp(HalfSpace);

  int LDim = Data->G->LDim;

  if ( HalfSpace && LDim!=2 )
   OSUsage(argv[0],OSArray,"--HalfSpace requires a 2D-periodic geometry");

  /***************************************************************/
  /* process arguments relating to brillouin-zone cubature *******/
  /***************************************************************/
  int BZIMethod=BZI_PCUBATURE;
  if (BZIString)
   { if (LDim==0)
      ErrExit("--BZIMethod may not be specified for non-PBC geometries");
     if (OkBPoints!=0)
      ErrExit("--BZIMethod and --OmegakBlochFile are mutually exclusive");
     if (!strcasecmp(BZIString, "DCUTRI"))
      { if (LDim==2)
         { Log("Using DCUTRI algorithm for Brillouin-zone integration.");
           BZIMethod=BZI_DCUTRI;
         }
        else
         { Warn("--BZIMethod DCUTRI not available for 1D periodicity "
                "(defaulting to --BZIMethod adaptive");
           BZIMethod=BZI_PCUBATURE;
         };
      }
     else if (!strncasecmp(BZIString, "FOTC",4) )
      {  
         BZIMethod=BZI_FOTC;
         Log("Using %s algorithm for Brillouin-zone integration.",BZIString);
      }
     else if (!strcasecmp(BZIString, "PCUBATURE") )
      {  Log("Using PCUBATURE algorithm for Brillouin-zone integration.");
         BZIMethod=BZI_PCUBATURE;
      }
     else 
      ErrExit("unknown Brillouin-zone integration scheme %s",BZIString);
   };

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
     Data->OutFileName=vstrdup("%s.LDOS",FileBase);
     WriteFilePreamble(Data->OutFileName, FILETYPE_LDOS, LDim);
     int NFun = 2*Data->XMatrix->NR;
     double *Result = new double[NFun];
     for(int no=0; no<OmegaPoints->N; no++)
      { Omega=OmegaPoints->GetEntry(no);
        GetLDOS(Data, Omega, 0, Result);
      };
     delete[] Result;
   }
  /*--------------------------------------------------------------*/
  /*- PBC structure with specified kBloch points: do a periodic   */
  /*- DGF calculation at each (Omega, kBloch) point               */
  /*--------------------------------------------------------------*/
  else if (OkBPoints)
   {  
     Data->ByKFileName=vstrdup("%s.byOmegakBloch",FileBase);
     WriteFilePreamble(Data->ByKFileName, FILETYPE_BYK, LDim);
     int NFun = 2*Data->XMatrix->NR;
     double *Result = new double[NFun];
     for(int nokb=0; nokb<OkBPoints->NR; nokb++)
      { 
        Omega=OkBPoints->GetEntry(nokb,0);

        double kBloch[2];
        kBloch[0]=OkBPoints->GetEntryD(nokb,1);
        if (LDim==2)
         kBloch[1]=OkBPoints->GetEntryD(nokb,2);

        GetLDOS(Data, Omega, kBloch, Result);
      };
     delete[] Result;
   }
  /*--------------------------------------------------------------*/
  /*- PBC structure without specified kBloch points: perform a    */
  /*- Brillouin-zone integration at each requested Omega value    */
  /*--------------------------------------------------------------*/
  else
   {
     Data->ByKFileName=vstrdup("%s.byOmegakBloch",FileBase);
     WriteFilePreamble(Data->ByKFileName, FILETYPE_BYK, LDim);
     Data->OutFileName=vstrdup("%s.LDOS",FileBase);
     WriteFilePreamble(Data->OutFileName, FILETYPE_LDOS, LDim);

     for(int no=0; no<OmegaPoints->N; no++)
      { Omega=OmegaPoints->GetEntry(no);
        switch(BZIMethod)
         { 
           case BZI_DCUTRI:
             GetBZIntegral_DCUTRI(Data, Omega);
             break;

           case BZI_FOTC:
             GetBZIntegral_FOTC(Data, Omega, BZIString);
             break;

           case BZI_PCUBATURE:
             GetBZIntegral_PCubature(Data, Omega);
             break;

           default:
             ErrExit("unknown Brillouin-zone integration scheme");
         }; 

      };

   };

}
