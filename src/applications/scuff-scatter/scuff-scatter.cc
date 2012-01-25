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
 *     --pwPol Ex Ey Ez 
 *     --pwDir Dx Dy Dz 
 * 
 *          incident field is a plane wave with E-field polarization
 *          given by pwPol and propagating in the pwDir
 *          direction. (the components of the pwPol d
 *          
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
#define MAXPSS  10    // max number of point sources
#define MAXFREQ 10    // max number of frequencies 
#define MAXEPF  10    // max number of evaluation-point files
#define MAXFM   10    // max number of flux meshes

     {"psLocation",     PA_DOUBLE,  3, MAXPSS,  (void *)psLoc,      &npsLoc,      "point source location"},
     {"psOrientation",  PA_DOUBLE,  3, MAXPSS,  (void *)psOrnt,     &npsOrnt,     "point source orientation"},
     {"Omega",          PA_CDOUBLE, 1, MAXFREQ, (void *)OmegaVals,  &nOmegaVals,  "(angular) frequency"},
     {"OmegaFile",      PA_STRING,  1, 1,       (void *)&OmegaFile, &nOmegaFiles, "list of (angular) frequencies"},
     {"EPFile",         PA_STRING,  1, MAXEPF,  (void *)EPFiles,    &nEPFiles,    "list of evaluation points"},
     {"FluxMesh",       PA_STRING,  1, MAXFM,   (void *)FluxMeshes, &nFluxMeshes, "flux mesh"},
     {0,0,0,0,0,0,0}
   };
 
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
  fprintf(stderr," scatterer options: \n\n");
  fprintf(stderr,"  --geometry MyGeometry.scuffgeo\n");
  fprintf(stderr,"\n");
  fprintf(stderr," incident field options: \n\n");
  fprintf(stderr,"  --pwPolarization Ex Ey Ez \n");
  fprintf(stderr,"  --pwDirection Nx Ny Nz \n");
  fprintf(stderr,"  --gbCenter xx yy zz \n");
  fprintf(stderr,"  --gbDirection Nx Ny Nz\n");
  fprintf(stderr,"  --psLocation xx yy zz \n");
  fprintf(stderr,"  --psOrientation xx yy zz \n");
  fprintf(stderr,"\n");
  fprintf(stderr," output options: \n\n");
  fprintf(stderr,"  --EPFile xx \n");
  fprintf(stderr,"  --FluxMesh MyFluxMesh.msh \n");
  fprintf(stderr,"\n");
  fprintf(stderr," frequency options: \n\n");
  fprintf(stderr,"  --omega xx \n");
  fprintf(stderr,"  --omegaFile MyFile.dat\n");
  fprintf(stderr,"\n");
  fprintf(stderr," miscellaneous options: \n");
  fprintf(stderr,"  --nThread xx \n");
  fprintf(stderr,"  --ExportMatrix \n");
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
  /* process options *********************************************/
  /***************************************************************/
  cdouble pwPol[3];                  int npwPol;
  double pwDir[3];                   int npwDir;
  double gbCenter[3];                int ngbCenter;
  double gbDir[3];                   int ngbDir;
  double gbWaist[1];                 int ngbWaist;
  double psLoc[3*MAXPSS];            int npsLoc;
  double psOrnt[3*MAXPSS];           int npsOrnt;
  cdouble OmegaVals[MAXFREQ];        int nOmegaVals;
  char *OmegaFile;                   int nOmegaFiles;
  char *EPFile[MAXEPF];              int nEPFiles;
  char *FluxMeshes[MAXFM];           int nFluxMeshes;
  int nThread=0;
  int ExportMatrix=0;
  /* name               type     args  instances  storage           count         description*/
  OptStruct OSArray[]=
   { {"pwPolarization", PA_CDOUBLE, 3, 1,       (void *)pwPol,      &npwPol,      "plane wave polarization"},
     {"pwDirection",    PA_DOUBLE,  3, 1,       (void *)pwDir,      &npwDir,      "plane wave direction"},
     {"gbCenter",       PA_DOUBLE,  3, 1,       (void *)gbCenter,   &ngbCenter,   "gaussian beam center"},
     {"gbDirection",    PA_DOUBLE,  3, 1,       (void *)gbDir,      &ngbDir,      "gaussian beam direction"},
     {"gbWaist",        PA_DOUBLE,  1, 1,       (void *)gbWaist,    &ngbWaist,    "gaussian beam waist"},
     {"psLocation",     PA_DOUBLE,  3, MAXPSS,  (void *)psLoc,      &npsLoc,      "point source location"},
     {"psOrientation",  PA_DOUBLE,  3, MAXPSS,  (void *)psOrnt,     &npsOrnt,     "point source orientation"},
     {"Omega",          PA_CDOUBLE, 1, MAXFREQ, (void *)OmegaVals,  &nOmegaVals,  "(angular) frequency"},
     {"OmegaFile",      PA_STRING,  1, 1,       (void *)&OmegaFile, &nOmegaFiles, "list of (angular) frequencies"},
     {"EPFile",         PA_STRING,  1, MAXEPF,  (void *)EPFiles,    &nEPFiles,    "list of evaluation points"},
     {"FluxMesh",       PA_STRING,  1, MAXFM,   (void *)FluxMeshes, &nFluxMeshes, "flux mesh"},
     {"nThread",        PA_BOOL,    0, 1,       (void *)&nThread,   0,            "number of CPU threads to use"},
     {"ExportMatrix",   PA_BOOL,    0, 1,       (void *)&ExportMatrix, 0,         "export BEM matrix to file"},
     {0,0,0,0,0,0,0}
   };

  if (nThread==0)
   nThread=GetNumProcs();

  /*******************************************************************/
  /* process frequency-related options to construct a list of        */
  /* frequencies at which to run simulations                         */
  /*******************************************************************/
  HVector *OmegaList0=0, *OmegaList=0;
  int nFreq, NumFreqs=0;
  if (nOmegaFiles==1)
   { 
     OmegaList0=new HVector(OmegaFile,LHM_TXT);
     if (OmegaList0->ErrMsg)
      ErrExit(ErrMsg);

     NumFreqs=OmegaList0->N;

     // convert it to cdouble if it isn't already
     if (OmegaList0->RealComplex!=LHM_COMPLEX)
      { OmegaList=new HVector(NumFreqs, LHM_COMPLEX);
        for(nFreq=0; nFreq<NumFreqs; nFreq++)
         OmegaList->SetEntry(nFreq, OmegaList0->GetEntry(nFreq));
        delete OmegaList0;
        OmegaList0=OmegaList;
      };

   };

  if (nOmegaVals>0)
   { 
     NumFreqs += nOmegaVals;
     OmegaList=new HVector(NumFreqs, LHM_COMPLEX);
     nFreq=0;
     if (OmegaList0) 
      for(nFreq=0; nFreq<OmegaList0->N; nFreq++)
       OmegaList->SetEntry(nFreq, OmegaList0->GetEntry(nFreq));
     int nf;
     for(nf=0; nf<nOmegaVals; nf++)
      OmegaList->SetEntry(nFreq+nf, OmegaVals[nf]);
   };

  /*******************************************************************/
  /* process incident-field-related options to construct the data    */
  /* used to generate the incident field in our scattering problem   */
  /*******************************************************************/
  IncFieldData *IFDList=0, *IFD=0;
  int npw;
  if ( npwPol != npwDir )
   ErrExit("numbers of --pwPolarization and --pwDirection options must agree");
  if ( ngbCenter!=ngbDirection || ngbDirection!=ngbWaist )
   ErrExit("numbers of --gbCenter, --gbDirection, and --gbWaist options must agree ");
  if ( npsLoc!=npsOrnt )
   ErrExit("numbers of --psLocation and --psOrientation options must agree");
  if ( npwPol==1 )
   { IFD=new PlaneWaveData(npwPol, npwDir, OmegaList->GetEntry(0) );
   };

  /*******************************************************************/
  /* create the RWGGeometry                                          */
  /*******************************************************************/
  RWGGeometry *G=new RWGGeometry(GeoFileName); 
  HMatrix *M=G->AllocateBEMMatrix(RealFreq);
  HVector *KN=G->AllocateRHSVector(RealFreq);

  /*******************************************************************/
  /* loop over frequencies *******************************************/
  /*******************************************************************/
  char OmegaStr[MAXSTR];
  for(nFreq=0; nFreq<NumFreqs; nFreqs++)
   { 
     Omega = OmegaList->GetEntry(nFreq);

     if ( real(Omega)==0.0 )
      snprintf(OmegaStr,MAXSTR,"%gi",imag(Omega));
     else if ( imag(Omega)==0.0 )
      snprintf(OmegaStr,MAXSTR,"%g",real(Omega));
     else 
      snprintf(OmegaStr,MAXSTR,"%g+%gi",real(Omega),imag(Omega));

     /*******************************************************************/
     /* assemble the BEM matrix, export it to a binary data file if     */
     /* that was requested, then LU-factorize.                          */
     /*******************************************************************/
     Log("Working at frequency %s...",OmegaStr);
     Log("  Assembling the BEM matrix...");
     G->AssembleBEMMatrix(Omega, nThread, M);
     if (ExportMatrix)
      { void *pCC=HMatrix::OpenC2MLContext(GetFileBase(G->GeoFileName),"_%s",OmegaStr);
        M->ExportToMATLAB(pCC,"M");
        HMatrix::CloseC2MLContext(pCC);
      };

     Log("  LU-factorizing BEM matrix...");
     M->LUFactorize();

     /***************************************************************/
     /* set up the incident field profile and create the RHS vector */
     /***************************************************************/
     Log("  Assembling the RHS vector...");


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
