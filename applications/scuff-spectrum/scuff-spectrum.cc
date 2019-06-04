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
#include <libSGJC.h>
#include "libBeyn.h"

#define II cdouble(0.0,1.0)

using namespace scuff;

#define MAXEPF   10    // max number of evaluation-point files
#define MAXFVM   10    // max number of field visualization meshes

/***************************************************************/
/* prototypes for postprocessing routines in OutputModules.cc***/
/***************************************************************/
void ProcessEPFile(RWGGeometry *G, HVector *KN, cdouble Omega, double *kBloch,
                   char *EPFileName, char *OutFileBase);

void VisualizeFields(RWGGeometry *G, HVector *KN, cdouble Omega, double *kBloch,
                     char *OutFileBase, char *FVMesh, char *TransFile);

void VisualizeFields(RWGGeometry *G, HVector *KN, cdouble Omega, double *kBloch,
                     char *OutFileBase, double *Screen);

void VisualizeSurfaceFields(RWGGeometry *G, HVector *KN, cdouble Omega, double *kBloch,
                            char *OutFileBase);

void WriteCartesianMoments(RWGGeometry *G, HVector *KN, cdouble Omega, double *kBloch, char *CartesianMomentFile);

void WriteSphericalMoments(RWGGeometry *G, HVector *KN, cdouble Omega, int LMax, char *SphericalMomentFile);

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct BFData
 { 
   RWGGeometry *G;
   HMatrix *M;
   double kBloch[2];
   FILE *LogFile;
 } BFData;

void BeynFunc(cdouble Omega, void *UserData, HMatrix *VHat, HMatrix *MVHat)
{
  BFData *Data   = (BFData *)UserData;

  RWGGeometry *G = Data->G;
  HMatrix *M     = Data->M;
  double *kBloch = Data->kBloch;
  FILE *LogFile  = Data->LogFile;

  if (G->LDim==0)
   Log(" assembling BEM matrix at Omega=%s",CD2S(Omega));
  if (G->LDim==1)
   Log(" assembling BEM matrix at k={%e},Omega=%s", kBloch[0],CD2S(Omega));
  if (G->LDim==2)
   Log(" assembling BEM matrix at k={%e,%e},Omega=%s", kBloch[0],kBloch[1],CD2S(Omega));

  if (G->LDim==0)
   G->AssembleBEMMatrix(Omega, M);
  else
   G->AssembleBEMMatrix(Omega, kBloch, M);

  if (LogFile)
   fprintf(LogFile,"%e %e\n",real(Omega),imag(Omega));

  if (MVHat)
   M->Multiply(VHat, MVHat);
  else
   { Log(" LUFactorizing...");
     M->LUFactorize(); 
     Log(" LUSolving...");
     M->LUSolve(VHat);
   }
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[]) 
{ 
  InstallHRSignalHandler();
  InitializeLog(argv[0]);

  /***************************************************************/
  /* process command-line arguments ******************************/
  /***************************************************************/
  char *GeoFile=0;
//
  cdouble Omega0   = cdouble(1.0,-1.0);
  double Rx        = 0.1;
  double Ry        = 0.1;
  int    N         = 20;
  int    L         = 5;
  double kx        = 0.0;  int nkx;
  double ky        = 0.0;  int nky;
  char *ContourFile=0;
//
  char *FileBase=0;
//
  bool PlotContours=false;
  bool PlotSurfaceCurrents=false;
  bool PlotSurfaceFields=false;
  char *SphericalMomentFile=0;
  int LMax=3;
  char *CartesianMomentFile=0;
  char *EPFiles[MAXEPF];             int nEPFiles;
  double FVScreens[11*MAXFVM];       int nFVScreens;
  char *FVMeshes[MAXFVM];            int nFVMeshes;
  char *FVMeshTransFiles[MAXFVM];    int nFVMeshTransFiles;
  memset(FVMeshTransFiles, 0, MAXFVM*sizeof(char *));
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"geometry",           PA_STRING,  1, 1, (void *)&GeoFile,            0,    ".scuffgeo file"},
// 
     {"Omega0",             PA_CDOUBLE, 1, 1, (void *)&Omega0,             0,    "center of contour"},
     {"Rx",                 PA_DOUBLE,  1, 1, (void *)&Rx,                 0,    "horizonal radius"},
     {"Ry",                 PA_DOUBLE,  1, 1, (void *)&Ry,                 0,    "vertical radius"},
     {"N",                  PA_INT,     1, 1, (void *)&N,                  0,    "number of quadrature points"},
     {"L",                  PA_INT,     1, 1, (void *)&L,                  0,    "upper bound on expected number of eigenvalues in contour "},
//
     {"kx",                 PA_INT,     1, 1, (void *)&kx,                 &nkx, "x component of bloch vector"},
     {"ky",                 PA_INT,     1, 1, (void *)&ky,                 &nky, "y component of bloch vector"},
//
     {"ContourFile",        PA_STRING,  1, 1, (void *)&ContourFile,        0,   "list of contours"},
//
     {"PlotContours",       PA_BOOL,    0, 1, (void *)&PlotContours,       0,   "plot contours for visualization"},
//
     {"PlotSurfaceCurrents",PA_BOOL,    0, 1, (void *)&PlotSurfaceCurrents,0,   "generate visualization files for eigenmode currents"},
     {"PlotSurfaceFields",  PA_BOOL,    0, 1, (void *)&PlotSurfaceFields,  0,   "generate visualization files for eigenmode surface fields"},
//
     {"SphericalMomentFile",PA_STRING,  1, 1, (void *)&SphericalMomentFile,0,   "name of output file for spherical multipole moments"},
     {"LMax",               PA_INT,     1, 1, (void *)&LMax,               0,   "index of highest spherical multipole moment to retain"},
//
     {"CartesianMomentFile",PA_STRING,  1, 1, (void *)&CartesianMomentFile,0,   "name of output file for cartesian multipole moments"},
//
     {"EPFile",             PA_STRING,  1, MAXEPF, (void *)EPFiles,          &nEPFiles,           "list of evaluation points"},
//
     {"FVScreen",           PA_DOUBLE,  11, MAXFVM, (void *)FVScreens,       &nFVScreens,         "field visualization screen"},
     {"FVMesh",             PA_STRING,  1, MAXFVM, (void *)FVMeshes,         &nFVMeshes,          "field visualization mesh"},
     {"FVMeshTransFile",    PA_STRING,  1, MAXFVM, (void *)FVMeshTransFiles, &nFVMeshTransFiles,  "list of geometrical transformations applied to FVMesh"},
//
     {"FileBase",           PA_STRING,  1, 1, (void *)&FileBase,           0,   "base name for output files"},
//
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFile==0)
   OSUsage(argv[0], VERSION, OSArray,"--geometry option is mandatory");
  if (FileBase==0)
   FileBase=strdup(GetFileBase(GeoFile));

  /***************************************************************/
  /* initialize RWGGeometry, read list of contours               */
  /***************************************************************/
  RWGGeometry *G         = new RWGGeometry(GeoFile);
  HMatrix *M             = G->AllocateBEMMatrix();
  int D = G->TotalBFs;

  /***************************************************************/
  /* process contour specifications ******************************/
  /***************************************************************/
  HMatrix *ContourMatrix;
  if (ContourFile)
   ContourMatrix = new HMatrix(ContourFile);
  else
   { ContourMatrix = new HMatrix(1, 5 + G->LDim, LHM_COMPLEX);
     int nc=0;
     if (G->LDim>=1) ContourMatrix->SetEntry(0,nc++,kx);
     if (G->LDim==2) ContourMatrix->SetEntry(0,nc++,ky);
     ContourMatrix->SetEntry(0,nc++,Omega0);
     ContourMatrix->SetEntry(0,nc++,Rx);
     ContourMatrix->SetEntry(0,nc++,Ry);
     ContourMatrix->SetEntry(0,nc++,N);
     ContourMatrix->SetEntry(0,nc++,L);
   };
    
  /***************************************************************/
  /* open eigenfrequency output file and write file header       */
  /***************************************************************/
  FILE *f=vfopen("%s.ModeFrequencies","a",FileBase);
  fprintf(f,"#%s " VERSION " run on %s (%s)\n",argv[0], GetHostName(), GetTimeString());
  fclose(f);

  /***************************************************************/
  /* apply Beyn's method to each user-specified contour          */
  /***************************************************************/
  SetDefaultCD2SFormat("%.8e %.8e");
  for(int nr=0; nr<ContourMatrix->NR; nr++)
   {
     /***************************************************************/
     /* extract information on the contour **************************/
     /***************************************************************/
     struct BFData MyBFData = {G, M, {0,0}, 0};
     double *kBloch = (G->LDim > 0) ? MyBFData.kBloch : 0;
     for(int d=0; d<G->LDim; d++) 
      MyBFData.kBloch[d] = ContourMatrix->GetEntryD(nr, d);
     Omega0              = ContourMatrix->GetEntry (nr, G->LDim+0);
     Rx                  = ContourMatrix->GetEntryD(nr, G->LDim+1);
     Ry                  = ContourMatrix->GetEntryD(nr, G->LDim+2);
     N                   = (int)(ContourMatrix->GetEntryD(nr, G->LDim+3));
     L                   = (int)(ContourMatrix->GetEntryD(nr, G->LDim+4));

     // create a string identifier for this contour
     char kbStr[20]="";
     if (G->LDim==1)
      snprintf(kbStr,20,"kx%g_",kBloch[0]);
     if (G->LDim==2)
      snprintf(kbStr,20,"kx%g_ky%g_",kBloch[0],kBloch[1]);
     char ContourLabel[100];
     snprintf(ContourLabel, 100, "C%i_w%s_%sRX%g_RY%g_N%i",nr,z2s(Omega0),kbStr,Rx,Ry,N);
     Log("Looking for eigenvalues in contour %s",ContourLabel);

     // initialize visualization files showing the contour if that was requested
     if (PlotContours)
      { MyBFData.LogFile = vfopen("%s.contours","a",FileBase);
        fprintf(MyBFData.LogFile,"\n\n# (Omega0,Rx,Ry,N)=%s,%e,%e,%i\n",CD2S(Omega0),Rx,Ry,N);
      }

     /***************************************************************/
     /* run Beyn's algorithm for this contour                       */
     /***************************************************************/
     BeynSolver *Solver    = CreateBeynSolver(D, L);
     HVector *Eigenvalues  = Solver->Eigenvalues;
     HVector *EVErrors     = Solver->EVErrors;
     HMatrix *Eigenvectors = Solver->Eigenvectors;
     HVector *Residuals    = Solver->Residuals;
     int NumModes=BeynSolve(Solver, BeynFunc, (void *)&MyBFData, Omega0, Rx, Ry, N);

     if (PlotContours)
      fclose(MyBFData.LogFile);

     /***************************************************************/
     /* write eigenfrequency results to .ModeFrequencies file       */
     /***************************************************************/
     f=vfopen("%s.ModeFrequencies","a",FileBase);
     fprintf(f,"# For contour w0=%s, Rx=%e, Ry=%e, N=%i, L=%i",z2s(Omega0),Rx,Ry,N,L);
     if (kBloch)
      { fprintf(f,", kBloch=");
        fprintVec(f,kBloch,G->LDim);
      }
     fprintf(f,":\n");
     fprintf(f,"# re(w) im(w)   estimated error in re(w), im(w)    residual\n");
     for(int n=0; n<NumModes; n++)
      { fprintf(f,"%+12e %+12e  ",real(Eigenvalues->GetEntry(n)), imag(Eigenvalues->GetEntry(n)));
        if (n<EVErrors->N)
         fprintf(f,"%+12e %+12e  ",real(EVErrors->GetEntry(n)), imag(EVErrors->GetEntry(n)));
        fprintf(f,"%+e ",Residuals->GetEntryD(n));
        fprintf(f,"\n");
      }
     fclose(f);

     /***************************************************************/
     /* post-processing of eigenvector data                         */
     /***************************************************************/
     for(int nm=0; nm<NumModes; nm++)
      { 
        char OutFileBase[100];
        snprintf(OutFileBase,100,"%s_%s_Mode%i",FileBase,ContourLabel,nm);

        cdouble Omega = Eigenvalues->GetEntry(nm);
        HVector KN(D, LHM_COMPLEX, (cdouble *)Eigenvectors->GetColumnPointer(nm));

        // write cartesian multiple moments
        if (CartesianMomentFile)
         WriteCartesianMoments(G, &KN, Omega, kBloch, CartesianMomentFile);

        // write cartesian multiple moments
        if (SphericalMomentFile)
         WriteSphericalMoments(G, &KN, Omega, LMax, SphericalMomentFile);

        // write surface-current visualization files
        if (PlotSurfaceCurrents)
         G->PlotSurfaceCurrents(&KN, Omega, kBloch, "%s.pp", OutFileBase);

        if (PlotSurfaceFields)
         VisualizeSurfaceFields(G, &KN, Omega, kBloch, OutFileBase);

        // process user-specified lists of field-evaluation points
        for(int nepf=0; nepf<nEPFiles; nepf++)
         ProcessEPFile(G, &KN, Omega, kBloch, EPFiles[nepf], OutFileBase);

        // process user-specified field-visualization meshes
        for(int nfm=0; nfm<nFVMeshes; nfm++)
         VisualizeFields(G, &KN, Omega, kBloch, OutFileBase, FVMeshes[nfm], FVMeshTransFiles[nfm]);
        for(int ns=0; ns<nFVScreens; ns++)
         { char VFFileBase[100];
           snprintf(VFFileBase,100,"%s.Screen%i",OutFileBase,ns);
           VisualizeFields(G, &KN, Omega, kBloch, VFFileBase, FVScreens + 11*ns);
         }
      }

     DestroyBeynSolver(Solver);

   }

}
