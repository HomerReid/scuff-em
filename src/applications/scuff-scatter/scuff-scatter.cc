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
 * a. options specifying the incident field: 
 * 
 *     --planewave
 * 
 *          incident field is a plane wave, traveing in the positive
 *          z direction, linearly polarized with E-field pointing in 
 *          the positive x-direction
 * 
 *     --gaussianbeam WW
 * 
 *          incident field is a gaussian beam, traveling in 
 *          the positive z direction, linearly polarized with
 *          E-field pointing in the positive x-direction, and with
 *          beam waist WW microns
 * 
 *     --pointsource xx yy zz nx ny nz 
 * 
 *          incident field is the field of a unit-strength point 
 *          electric dipole source at coordinates (xx,yy,zz) and
 *          pointing in the (nx,ny,nz) direction
 * 
 *       -------------------------------------------------
 * 
 * b. options specifying the output: 
 * 
 *     --console 
 * 
 *         if this option is specified, the code will enter a  
 *         console-query mode in which the user can enter field evaluation
 *         points and the code will print out the scattered E and H 
 *         fields at each point. 
 * 
 *     --EPFile MyEPFile
 * 
 *         a file containing a list of points at which the scattered 
 *         field is to be evaluated. (all field points must lie in the   
 *         EXTERIOR of all scattering objects in the geometry).
 * 
 *         each line of EPFile should contain 3 numbers, the cartesian
 *         components of the scattered field. (blank lines and comments,
 *         i.e. lines beginning with a '#', are skipped)
 * 
 * 
 *     --FluxMesh  MyFirstFluxMesh.msh
 *     --FluxMesh  MySecondFluxMesh.msh
 *        ...       
 *     --FluxMesh  MyNinthFluxMesh.msh
 * 
 *         GMSH .msh files containing 3D surface meshes (discretized
 *         into triangles) through which to compute poynting flux.
 * 
 *       -------------------------------------------------
 * 
 * c. options describing the frequency
 * 
 *     --frequency xx (value of real or imaginary frequency)
 *     --RF | --IF    (choice of real or imaginary frequency)
 *
 *       -------------------------------------------------
 * 
 * d. options describing the scatterer
 * 
 *     --geometry MyGeometry.rwggeo 
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

#include <libhrutil.h>
#include <libIncField.h>
#include <libhmat.h>

#include "libRWG.h"
#include "EMScatter.h"
 
/***************************************************************/
/***************************************************************/
/***************************************************************/
void usage(char *ProgramName, const char *format, ... )
{ 
  va_list ap;
  char buffer[1000];

  if (format)
   { va_start(ap,format);
     vsnprintf(buffer,1000,format,ap);
     va_end(ap);
     fprintf(stderr,"error: %s (aborting)\n\n",buffer);
   };

  fprintf(stderr,"usage: %s [incident field options] [scatterer options]\n",ProgramName);
  fprintf(stderr,"\n");
  fprintf(stderr," incident field options: \n\n");
  fprintf(stderr,"  --planewave \n");
  fprintf(stderr,"  --gaussianbeam WW \n");
  fprintf(stderr,"  --pointsource xx yy zz nx ny nz\n");
  fprintf(stderr,"  --MRICoil\n");
  fprintf(stderr,"\n");
  fprintf(stderr," output options: \n\n");
  fprintf(stderr,"  --console \n\n");
  fprintf(stderr,"  --OutsidePointsFile xx (list of field evaluation points in external region)\n");
  fprintf(stderr,"  --InsidePointsFile xx (list of field evaluation points inside first object)\n");
  fprintf(stderr,"  --FluxMesh FluxMesh1.msh [--FluxMesh FluxMesh2.msh ...]\n");
  fprintf(stderr,"\n");
  fprintf(stderr," frequency options: \n\n");
  fprintf(stderr,"  --frequency xx \n\n");
  fprintf(stderr,"  --rf (real frequency)\n");
  fprintf(stderr,"  --if (real frequency)\n");
  fprintf(stderr,"\n");
  fprintf(stderr," scatterer options: \n");
  fprintf(stderr,"  --geometry xx \n\n");
  fprintf(stderr,"\n");
  fprintf(stderr," miscellaneous options: \n");
  fprintf(stderr,"  --nThread xx \n");
  fprintf(stderr,"  --ExportMatrix \n");
  fprintf(stderr,"\n");
  fprintf(stderr,"\n");
  exit(1);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{

  int narg, nConv, DoConsole=0, SourceType=LIF_TYPE_PSEC, nThread;

  double Frequency, MuI, MuE;
  cdouble EpsI, EpsE;
  int lMax, RealFreq, ExportMatrix;

  PlaneWaveData MyPWD, *PWD=0;
  PointSourceData MyPSD, *PSD=0;
  GaussianBeamData MyGBD, *GBD=0;
  MagneticFrillData MyMFD, *MFD=0;
  double WW;
  double PSLocation[3], PSDirection[3];

  char *GeoFileName;

  #define MAXEPFILES 10
  char *OEPFiles[MAXEPFILES];
  int nOEPFiles=0;
  char *IEPFiles[MAXEPFILES];
  int nIEPFiles=0;
  
  #define MAXFLUXMESHES 10
  char *FluxMeshFiles[MAXFLUXMESHES];
  int nFluxMeshFiles=0;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble pwPol[MAXPWS * 3];
  double pwDir[MAXPWS * 3];
  double gbCenter[MAXGBS];
  ExtendedArgStruct ASArray[]=
   { {"pwPolarization", PA_CDOUBLE, 3, MAXPWS,  (void *)&pwPol,     &npwPol,        "plane wave polarization vector"},
     {"pwDirection",    PA_DOUBLE,  3, MAXPWS,  (void *)&pwDir,     &npwDir,        "plane wave direction vector"},
     {"gbCenter",       PA_DOUBLE,  3, MAXGBS,  (void *)&gbCenter,  &ngbCenter,     "gaussian beam center"},
     {"gbDirection",    PA_DOUBLE,  3, MAXGBS,  (void *)&gbDir,     &ngbDir,        "gaussian beam direction vector"},
     {"gbWaist",        PA_DOUBLE,  1, MAXGBS,  (void *)&gbWaist,   &ngbWaist,      "gaussian beam waist"},
     {"psLocation",     PA_DOUBLE,  3, MAXPSS,  (void *)&psLoc,     &npsLoc,       "point source location"},
     {"psOrientation",  PA_DOUBLE,  3, MAXPSS,  (void *)&psOrnt,    &npsOrnt,      "point source orientation vector"},
     {"Omega",          PA_DOUBLE,  1, MAXFREQ, (void *)&OmegaVals, &psOrientation, "point source orientation vector"},
     {0,0,0,0,0}

  /***************************************************************/
  /* process command line arguments not relating to scatterer    */
  /***************************************************************/
  nThread=-1;
  GeoFileName=0;
  RealFreq=1;
  ExportMatrix=0;
  for(narg=1; narg<argc; narg++)
   { 
     if ( !strcasecmp(argv[narg], "--geometry") )
      { 
        if ( narg+1 >= argc ) 
         usage(argv[0],"--geometry requires an argument");
        GeoFileName=argv[narg+1];
        narg++;
      }
     else if ( !strcasecmp(argv[narg], "--frequency") )
      { 
        if ( narg+1 >= argc ) 
         usage(argv[0],"--frequency requires an argument");
        nConv=sscanf(argv[narg+1],"%le",&Frequency);
        if (nConv!=1)
         usage(argv[0],"invalid value (%s) specified for frequency parameter",argv[narg]);
        printf("Working at frequency Omega=%e.\n",Frequency);
        narg++;
      }
     else if ( !strcasecmp(argv[narg], "--rf") )
      { 
        RealFreq=1;
        printf("Real-frequency computation selected.\n");
      }
     else if ( !strcasecmp(argv[narg], "--if") )
      { 
        RealFreq=0;
        printf("Imaginary-frequency computation selected.\n");
      }
     else if ( !strcasecmp(argv[narg], "--OutsidePointsFile") )
      { 
        if ( narg+1 >= argc ) 
         usage(argv[0],"--OutsidePointsFile requires an argument");
        if ( nOEPFiles == MAXEPFILES )
         ErrExit("too many EP files");
        OEPFiles[nOEPFiles++]=argv[narg+1];
        printf("Will evaluate fields at outside points in file %s.\n",argv[narg+1]);
        narg++;
      }
     else if ( !strcasecmp(argv[narg], "--InsidePointsFile") )
      { 
        if ( narg+1 >= argc ) 
         usage(argv[0],"--InsidePointsFile requires an argument");
        if ( nIEPFiles == MAXEPFILES )
         ErrExit("too many EP files");
        IEPFiles[nIEPFiles++]=argv[narg+1];
        printf("Will evaluate fields at inside points in file %s.\n",argv[narg+1]);
        narg++;
      }
    else if ( !strcasecmp(argv[narg], "--FluxMesh") )
      { 
        if ( narg+1 >= argc ) 
         usage(argv[0],"--FluxMesh requires an argument");
        if ( nFluxMeshFiles == MAXFLUXMESHES )
         ErrExit("too many flux meshes");
        FluxMeshFiles[nFluxMeshFiles++]=argv[narg+1];
        printf("Will generate flux plot for surface %s.\n",argv[narg+1]);
        narg++;
      }
     else if ( !strcasecmp(argv[narg], "--planewave") )
      { 
        if (PWD || GBD || PSD || MFD)
         usage(argv[0],"only one type of incident field may be specified");
        PWD=&MyPWD;
        printf("Incident field is z-traveling x-polarized plane wave.\n");
      }
     else if ( !strcasecmp(argv[narg], "--gaussianbeam") )
      { 
        if (PWD || GBD || PSD || MFD)
         usage(argv[0],"only one type of incident field may be specified");
        if ( narg+1 >= argc ) 
         usage(argv[0],"--gaussianbeam requires an argument");
        sscanf(argv[narg+1],"%le",&WW);
        printf("Incident field is z-traveling x-polarized gaussian beam\n");
        printf(" with beam waist %g um.\n",WW);
        GBD=&MyGBD;
        narg++;
      }
     else if ( !strcasecmp(argv[narg], "--pointsource") )
      { 
        if (PWD || GBD || PSD || MFD)
         usage(argv[0],"only one type of incident field may be specified");
        if ( narg+6 >= argc ) 
         usage(argv[0],"--pointsource requires 6 arguments");
        if (    1!=sscanf(argv[narg+1],"%le",PSLocation+0 )
             || 1!=sscanf(argv[narg+2],"%le",PSLocation+1 )
             || 1!=sscanf(argv[narg+3],"%le",PSLocation+2 )
             || 1!=sscanf(argv[narg+4],"%le",PSDirection+0 )
             || 1!=sscanf(argv[narg+5],"%le",PSDirection+1 )
             || 1!=sscanf(argv[narg+6],"%le",PSDirection+2 )
           )
         usage(argv[0],"invalid argument passed to --pointsource option");
 
        PSD=&MyPSD;
        printf("Incident field is field of point source \n");
        printf(" at (%g,%g,%g) pointing in (%g,%g,%g) direction.\n",
                PSLocation[0], PSLocation[1], PSLocation[2],
                PSDirection[0], PSDirection[1], PSDirection[2]);
        narg+=6;
      }
     else if ( !strcasecmp(argv[narg], "--MRICoil") )
      { 
        if (PWD || GBD || PSD || MFD)
         usage(argv[0],"only one type of incident field may be specified");
        printf("Incident field is field of MRI coil.\n");
        MFD=&MyMFD;
      }
     else if ( !strcasecmp(argv[narg], "--psmc") )
      { 
        SourceType=LIF_TYPE_PSMC;
        printf("Incident field is PSMC.\n");
      }
     else if ( !strcasecmp(argv[narg], "--nThread") )
      { 
        if ( narg+1 >= argc ) 
         usage(argv[0],"--nThread requires an argument");
        sscanf(argv[narg+1],"%i",&nThread);
        narg++;
      }
     else if ( !strcasecmp(argv[narg], "--ExportMatrix") )
      { 
        ExportMatrix=1;
        printf("Will export BEM matrix to binary file.\n");
      } 
     else
      { 
         usage(argv[0],"error: unknown option %s (aborting)",argv[narg]);
      };
   };

  if( GeoFileName==0 )
   usage(argv[0],"--geometry argument is mandatory ");

  if( PSD==0 && PWD==0 && GBD==0 && MFD==0 ) 
   usage(argv[0],"you must specify an incident field");

  if (nThread==-1)
   nThread=GetNumProcs();

  SetLogFileName("EMScatter.log");

  /*******************************************************************/
  /* create the RWGGeometry                                          */
  /*******************************************************************/
  RWGGeometry *G=new RWGGeometry(GeoFileName); 

  char buffer[1000];
  snprintf(buffer,1000,"%s.pp",GetFileBase(GeoFileName));
  G->WritePPMesh(buffer,GetFileBase(GeoFileName),0);

  Log("Precomputing...");
  printf("Precomputing tables for BEM matrix...\n");
  G->PreCompute(nThread);
  HMatrix *M=G->AllocateBEMMatrix(RealFreq);
  HVector *KN=G->AllocateRHSVector(RealFreq);

  /*******************************************************************/
  /* assemble the BEM matrix, export it to a binary data file if     */
  /* that was requested, then LU-factorize.                          */
  /*******************************************************************/
  Log("Assembling the BEM matrix...");
  printf("Assembling the BEM matrix...\n");
  G->AssembleBEMMatrix(Frequency, RealFreq, nThread, M);
  if (ExportMatrix)
   { void *pCC=HMatrix::OpenC2MLContext(GetFileBase(G->GeoFileName));
     M->ExportToMATLAB(pCC,"M");
     HMatrix::CloseC2MLContext(pCC);
   };

  Log("LU-factorizing BEM matrix...");
  printf("LU-factorizing BEM matrix...\n");
  M->LUFactorize();

  /***************************************************************/
  /* set up the incident field profile and create the RHS vector */
  /***************************************************************/
  Log("Assembling the RHS vector...");
  printf("Assembling the RHS vector...\n");
  if (PWD)
   { 
     PWD->E0[0]=1.0;        PWD->E0[1]=0.0;        PWD->E0[2]=0.0;
     PWD->nHat[0]=0.0;      PWD->nHat[1]=0.0;      PWD->nHat[2]=1.0;
     PWD->Frequency=Frequency;
     PWD->RealFreq=RealFreq;
     PWD->Eps=G->MP->GetEpsD(Frequency,RealFreq);
     PWD->Mu=G->MP->GetMu(Frequency,RealFreq);
     G->AssembleRHSVector(EHPlaneWave, (void *)PWD, RealFreq, nThread, KN);
   }
  else if (GBD)
   { 
     GBD->X0[0]=0.0;        GBD->X0[1]=0.0;        GBD->X0[2]=0.0;
     GBD->KProp[0]=0.0;     GBD->KProp[1]=0.0;     GBD->KProp[2]=1.0;
     GBD->E0[0]=1.0;        GBD->E0[1]=0.0;        GBD->E0[2]=0.0;
     GBD->Frequency=Frequency;
     GBD->RealFreq=RealFreq;
     GBD->Eps=G->MP->GetEpsD(Frequency,RealFreq);
     GBD->Mu=G->MP->GetMu(Frequency,RealFreq);
     GBD->W0=WW;
     G->AssembleRHSVector(EHGaussianBeam, (void *)GBD, RealFreq, nThread, KN);
   }
  else if (PSD) // PSD
   {
     memcpy(PSD->X0, PSLocation, 3*sizeof(double) );
     memcpy(PSD->nHat, PSDirection, 3*sizeof(double) );
     PSD->Frequency=Frequency;
     PSD->RealFreq=RealFreq;
     PSD->SourceType=SourceType;
     PSD->Eps=G->MP->GetEpsD(Frequency,RealFreq);
     PSD->Mu=G->MP->GetMu(Frequency,RealFreq);
     G->AssembleRHSVector(EHPointSource, (void *)PSD, RealFreq, nThread, KN);
   }
  else  // MFD
   {
     MFD->Frequency=Frequency;
     MFD->RealFreq=RealFreq;
     MFD->Eps=G->MP->GetEpsD(Frequency,RealFreq);
     MFD->Mu=G->MP->GetMu(Frequency,RealFreq);
     MFD->X0[0]=0.0;
     MFD->X0[1]=-1.0;
     MFD->X0[2]=0.0;
     MFD->Theta=0.0;
     MFD->Phi=0.0;
     MFD->RIn=0.01;
     MFD->ROut=0.02;
     G->AssembleRHSVector(EHMagneticFrill, (void *)MFD, RealFreq, nThread, KN);

     char buffer[1000];
     snprintf(buffer,1000,"%s.pp",GetFileBase(GeoFileName));
     DrawMagneticFrill( (void *)MFD, buffer);
   };


  /***************************************************************/
  /* solve the BEM system*****************************************/
  /***************************************************************/
  Log("Solving the BEM system...");
  printf("Solving the BEM system...\n");
  M->LUSolve(KN);

  /***************************************************************/
  /* and now we've solved the scattering problem, and we can     */
  /* proceed to compute the scattered field at arbitary points   */
  /* in space, which we do in one or more of three ways according*/
  /* to what the user requested.                                 */
  /***************************************************************/

  /***************************************************************/
  /* first create the data structure passed to the output modules*/
  /***************************************************************/
  EMSData MyEMSData, *EMSD=&MyEMSData;
  EMSD->G=G;
  EMSD->KN=KN;
  EMSD->Frequency=Frequency;
  EMSD->RealFreq=RealFreq;
  EMSD->nThread=nThread;
  EMSD->opPWD=(void *)PWD;
  EMSD->opGBD=(void *)GBD;
  EMSD->opPSD=(void *)PSD;
  EMSD->opMFD=(void *)MFD;

  /*--------------------------------------------------------------*/
  /* a. interactive console mode if user requested it             */
  /*--------------------------------------------------------------*/
  if (DoConsole)
   { Log("Entering console mode...");
     Console(EMSD);
   };

  /*--------------------------------------------------------------*/
  /*- b. fields at user-specified lists of evaluation points     -*/
  /*--------------------------------------------------------------*/
  int nepf;
  for(nepf=0; nepf<nOEPFiles; nepf++)
   ProcessEPFile(EMSD, OEPFiles[nepf], -1);
  for(nepf=0; nepf<nIEPFiles; nepf++)
   ProcessEPFile(EMSD, IEPFiles[nepf], 0);
  HMatrix *EPMatrix=0;

  /*--------------------------------------------------------------*/
  /*- c. flux plots on user-specified surface meshes -------------*/
  /*--------------------------------------------------------------*/
  int nfm;
  for(nfm=0; nfm<nFluxMeshFiles; nfm++)
   CreateFluxPlot(EMSD, FluxMeshFiles[nfm]);
   
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  printf("Thank you for your support.\n");

}
