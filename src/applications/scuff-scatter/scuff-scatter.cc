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
 *     --FluxMesh  MyFluxMesh.msh
 * 
 *         (Note that --FluxMesh may be specified more than once.)
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

#include "scuff-scatter.h"

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define MAXPW   1     // max number of plane waves
#define MAXGB   1     // max number of gaussian beams
#define MAXPS   10    // max number of point sources
#define MAXFREQ 10    // max number of frequencies 
#define MAXEPF  10    // max number of evaluation-point files
#define MAXFM   10    // max number of flux meshes

#define MAXSTR  1000
 
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
  fprintf(stderr,"  --gbPolarization Nx Ny Nz\n");
  fprintf(stderr,"  --gbDirection Nx Ny Nz\n");
  fprintf(stderr,"  --gbCenter xx yy zz \n");
  fprintf(stderr,"  --gbWaist W \n");
  fprintf(stderr,"  --psLocation xx yy zz \n");
  fprintf(stderr,"  --psOrientation xx yy zz \n");
  fprintf(stderr,"\n");
  fprintf(stderr," output options: \n\n");
  fprintf(stderr,"  --EPFile xx \n");
  fprintf(stderr,"  --FluxMesh MyFluxMesh.msh \n");
  fprintf(stderr,"\n");
  fprintf(stderr," frequency options: \n\n");
  fprintf(stderr,"  --omega xx \n");
  fprintf(stderr,"  --omegaFile MyOmegaFile.dat\n");
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
void z2s(cdouble z, char *zStr)
{ 
  if ( real(z)==0.0 )
   snprintf(zStr,MAXSTR,"%gi",imag(z));
  else if ( imag(z)==0.0 )
   snprintf(zStr,MAXSTR,"%g",real(z));
  else 
   snprintf(zStr,MAXSTR,"%g+%gi",real(z),imag(z));
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  /***************************************************************/
  /* process options *********************************************/
  /***************************************************************/
  char *GeoFile=0;
  cdouble pwPol[3*MAXPW];            int npwPol;
  double pwDir[3*MAXPW];             int npwDir;
  cdouble gbPol[3*MAXGB];            int ngbPol;
  double gbDir[3*MAXGB];             int ngbDir;
  double gbCenter[3*MAXGB];          int ngbCenter;
  double gbWaist[MAXGB];             int ngbWaist;
  double psLoc[3*MAXPS];             int npsLoc;
  cdouble psStrength[3*MAXPS];       int npsStrength;
  cdouble OmegaVals[MAXFREQ];        int nOmegaVals;
  char *OmegaFile;                   int nOmegaFiles;
  char *EPFiles[MAXEPF];             int nEPFiles;
  char *FluxMeshes[MAXFM];           int nFluxMeshes;
  int nThread=0;
  int ExportMatrix=0;
  /* name               type     args  instances  storage           count         description*/
  OptStruct OSArray[]=
   { {"geometry",       PA_STRING,  1, 1,       (void *)&GeoFile,   0,            "geometry file"},
     {"pwPolarization", PA_CDOUBLE, 3, MAXPW,   (void *)pwPol,      &npwPol,      "plane wave polarization"},
     {"pwDirection",    PA_DOUBLE,  3, MAXPW,   (void *)pwDir,      &npwDir,      "plane wave direction"},
     {"gbPolarization", PA_CDOUBLE, 3, MAXGB,   (void *)gbPol,      &ngbPol,      "gaussian beam polarization"},
     {"gbDirection",    PA_DOUBLE,  3, MAXGB,   (void *)gbDir,      &ngbDir,      "gaussian beam direction"},
     {"gbCenter",       PA_DOUBLE,  3, MAXGB,   (void *)gbCenter,   &ngbCenter,   "gaussian beam center"},
     {"gbWaist",        PA_DOUBLE,  1, MAXGB,   (void *)gbWaist,    &ngbWaist,    "gaussian beam waist"},
     {"psLocation",     PA_DOUBLE,  3, MAXPS,   (void *)psLoc,      &npsLoc,      "point source location"},
     {"psStrength",     PA_CDOUBLE, 3, MAXPS,   (void *)psStrength, &npsStrength, "point source strength"},
     {"Omega",          PA_CDOUBLE, 1, MAXFREQ, (void *)OmegaVals,  &nOmegaVals,  "(angular) frequency"},
     {"OmegaFile",      PA_STRING,  1, 1,       (void *)&OmegaFile, &nOmegaFiles, "list of (angular) frequencies"},
     {"EPFile",         PA_STRING,  1, MAXEPF,  (void *)EPFiles,    &nEPFiles,    "list of evaluation points"},
     {"FluxMesh",       PA_STRING,  1, MAXFM,   (void *)FluxMeshes, &nFluxMeshes, "flux mesh"},
     {"nThread",        PA_INT,     1, 1,       (void *)&nThread,   0,            "number of CPU threads to use"},
     {"ExportMatrix",   PA_BOOL,    0, 1,       (void *)&ExportMatrix, 0,         "export BEM matrix to file"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  if (nThread==0)
   nThread=GetNumProcs();

  /*******************************************************************/
  /* process frequency-related options to construct a list of        */
  /* frequencies at which to run simulations                         */
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

  IncFieldData *IFDList=0, *IFD;
  int npw, ngb, nps;
  for(npw=0; npw<npwPol; npw++)
   { IFD=new PlaneWaveData(pwPol + 3*npw, pwDir + 3*npw);
     IFD->Next=IFDList;
     IFDList=IFD;
   };
#if 0
  for(ngb=0; ngb<ngbCenter; ngb++)
   { IFD=new GaussianBeamData(gbCenter + 3*ngb, gbDir + 3*ngb, gbPol + 3*ngb, gbWaist[ngb]);
     IFD->Next=IFDList;
     IFDList=IFD;
   };
#endif
  for(nps=0; nps<npsLoc; nps++)
   { IFD=new PointSourceData(psLoc + 3*nps, psStrength + 3*nps);
     IFD->Next=IFDList;
     IFDList=IFD;
   };

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  if ( (nEPFiles>0 || nFluxMeshes>0) && IFDList==0 )
   ErrExit("you must specify at least one incident field source");

  /*******************************************************************/
  /* create the RWGGeometry, allocate BEM matrix and RHS vector      */
  /*******************************************************************/
  RWGGeometry *G=new RWGGeometry(GeoFile);
#ifdef SCUFF
  HMatrix *M=G->AllocateBEMMatrix();
  HVector *KN=G->AllocateRHSVector();
#else
  HMatrix *M=G->AllocateBEMMatrix(REAL_FREQ);
  HVector *KN=G->AllocateRHSVector(REAL_FREQ);
#endif

  char GeoFileBase[MAXSTR];
  strncpy(GeoFileBase, GetFileBase(GeoFile), MAXSTR);

  /*******************************************************************/
  /* put all relevant information into a data structure to be passed */
  /* to the output modules                                           */
  /*******************************************************************/
  SSData MySSData, *SSD=&MySSData;
  SSD->G=G;
  SSD->KN=KN;
  SSD->opIFD=(void *)IFDList;

  /*******************************************************************/
  /* loop over frequencies *******************************************/
  /*******************************************************************/
  char OmegaStr[MAXSTR];
  cdouble Omega;
  cdouble Eps;
  double Mu;
  for(nFreq=0; nFreq<NumFreqs; nFreq++)
   { 
     Omega = OmegaList->GetEntry(nFreq);
     z2s(Omega, OmegaStr);

     /*******************************************************************/
     /* assemble the BEM matrix, export it to a binary data file if     */
     /* that was requested, then LU-factorize.                          */
     /*******************************************************************/
     Log("Working at frequency %s...",OmegaStr);
     Log("  Assembling the BEM matrix...");
#ifdef SCUFF
     G->AssembleBEMMatrix(Omega, nThread, M);
#else
     int RealFreq;
     if ( imag(Omega)==0.0 )
      RealFreq=1;
     else if ( real(Omega)==0.0 )
      RealFreq=0;
     else
      ErrExit("%s:%i:internal error",__FILE__,__LINE__);
     
     G->PreCompute(nThread);
     G->AssembleBEMMatrix( abs(Omega), RealFreq, nThread, M);
#endif
     if (ExportMatrix)
      { void *pCC=HMatrix::OpenC2MLContext(GeoFileBase,"_%s",OmegaStr);
        M->ExportToMATLAB(pCC,"M");
        HMatrix::CloseC2MLContext(pCC);
      };

     if ( nEPFiles==0 && nFluxMeshes==0 )
      continue;

     Log("  LU-factorizing BEM matrix...");
     M->LUFactorize();

     /***************************************************************/
     /* set up the incident field profile and assemble the RHS vector */
     /***************************************************************/
     Log("  Assembling the RHS vector..."); 
#ifdef SCUFF
     G->ExteriorMP->GetEpsMu(Omega,&Eps,&Mu);
#else
     G->MP->GetEpsMu(Omega,&Eps,&Mu);
#endif
     IFDList->SetFrequencyAndEpsMu(Omega,Eps,Mu);

#ifdef SCUFF
     G->AssembleRHSVector(EHIncField, (void *)IFDList, nThread, KN);
#else
     G->AssembleRHSVector(EHIncField, (void *)IFDList, RealFreq, nThread, KN);
#endif

     /***************************************************************/
     /* solve the BEM system*****************************************/
     /***************************************************************/
     Log("Solving the BEM system...");
     printf("Solving the BEM system...\n");
     M->LUSolve(KN);

     /***************************************************************/
     /* process requested outputs                                   */
     /***************************************************************/
     SSD->Omega=Omega;

     int nepf;
     for(nepf=0; nepf<nEPFiles; nepf++)
      ProcessEPFile(SSD, EPFiles[nepf], 0);

     int nfm;
     for(nfm=0; nfm<nFluxMeshes; nfm++)
      CreateFluxPlot(SSD, FluxMeshes[nfm]);

   }; //  for(nFreq=0; nFreq<NumFreqs; nFreqs++)

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  printf("Thank you for your support.\n");

}
