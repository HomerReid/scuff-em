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
 * scuff-microstrip.cc -- scuff-EM module for RF modeling of microstrip
 *                     -- geometries
 * 
 * homer reid    -- 3/2018
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <libhrutil.h>
#include <libhmat.h>
#include <libscuff.h>
#include <libIncField.h>

#include "RWGPorts.h"

#define FREQ2OMEGA (2.0*M_PI/300.0)
#define OMEGA2FREQ (1/(FREQ2OMEGA))
#define II cdouble(0.0,1.0)

using namespace scuff;

#define MAXFREQ  10    // max number of frequencies
#define MAXSTR   1000
#define MAXEPF   10    // max number of field visualization files
#define MAXFVM   10    // max number of field visualization meshes

/***************************************************************/
/* output modules in OutputModules.cc **************************/
/***************************************************************/
void ProcessEPFile(RWGGeometry *G, HVector *KN, cdouble Omega,
                   RWGPortList *PortList, cdouble *PortCurrents,
                   char *EPFile, char *FileBase);

void ProcessFVMesh(RWGGeometry *G, HVector *KN, cdouble Omega,
                   RWGPortList *PortList, cdouble *PortCurrents,
                   char *FVMesh, char *FileBase);

void ComputeSZParms(RWGGeometry *G, RWGPortList *PortList, cdouble Omega,
                    HMatrix *M, HVector *KN, char *FileBase, bool SParms);
namespace scuff{
void GetReducedPotentials_Nearby(RWGSurface *S, const int ne,
                                 const double X0[3],  const cdouble k,
                                 cdouble *p, cdouble a[3],
                                 cdouble dp[3], cdouble da[3][3],
                                 cdouble ddp[3][3], cdouble dcurla[3][3],
                                 bool *IncludeTerm=0);
}

/***************************************************************/
/* main function   *********************************************/
/***************************************************************/  
int main(int argc, char *argv[])
{
  InstallHRSignalHandler();
  InitializeLog(argv[0]);

  /***************************************************************/
  /** process command-line arguments *****************************/
  /***************************************************************/
  char *GeoFile=0;
  char *PortFile=0;
  bool PlotPorts=false;
//
  char *FreqFile=0;
  double Freqs[MAXFREQ];    int nFreqs;
  double MinFreq;           int nMinFreq;
  double MaxFreq;	    int nMaxFreq;
  int NumFreqs;		    int nNumFreqs;
  bool LogFreq=false;
//
  char *PCFile=0;
//
  char *EPFiles[MAXEPF];    int nEPFiles;
  char *FVMeshes[MAXFVM];   int nFVMeshes;
  bool ZParms=false;
  bool SParms=false;
//
  char *FileBase=0;
  char *ContribOnly=0;
  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { {"geometry",       PA_STRING,  1, 1,       (void *)&GeoFile,    0,             "geometry file"},
//
     {"portfile",       PA_STRING,  1, 1,       (void *)&PortFile,   0,             "port file"},
     {"PlotPorts",      PA_BOOL,    0, 1,       (void *)&PlotPorts,  0,             "generate port visualization file"},
//
     {"freqfile",       PA_STRING,  1, 1,       (void *)&FreqFile,   0,             "list of frequencies"},
     {"frequency",      PA_DOUBLE,  1, MAXFREQ, (void *)Freqs,       &nFreqs,       "frequency (GHz)"},
     {"minfreq",        PA_DOUBLE,  1, 1,       (void *)&MinFreq,    &nMinFreq,     "starting frequency"},
     {"maxfreq",        PA_DOUBLE,  1, 1,       (void *)&MaxFreq,    &nMaxFreq,     "ending frequency"},
     {"numfreqs",       PA_INT,     1, 1,       (void *)&NumFreqs,   &nNumFreqs,    "number of frequencies"},
     {"logfreq",        PA_BOOL,    0, 1,       (void *)&LogFreq,    0,             "use logarithmic frequency steps"},
//
     {"portcurrentfile", PA_STRING,  1, 1,      (void *)&PCFile,     0,             "port current file"},
//
     {"ZParameters",    PA_BOOL,    0, 1,       (void *)&ZParms,     0,             "output z parameters"},
     {"SParameters",    PA_BOOL,    0, 1,       (void *)&SParms,     0,             "output s parameters"},
//
     {"EPFile",         PA_STRING,  1, MAXEPF,  (void *)EPFiles,     &nEPFiles,     "list of evaluation points"},
     {"FVMesh",         PA_STRING,  1, MAXFVM,  (void *)FVMeshes,    &nFVMeshes,    "field visualization mesh"},
//
     {"FileBase",       PA_STRING,  1, 1,       (void *)&FileBase,   0,             "base name for output files"},
//
     {"ContribOnly",    PA_STRING,  1, 1,       (void *)&ContribOnly,0,             "select port voltage contributors"},
//
//
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFile==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");
  if (PortFile==0)
   OSUsage(argv[0],OSArray,"--PortFileName option is mandatory");
  if (!FileBase)
   FileBase=strdup(GetFileBase(GeoFile));

  /***************************************************************/
  /* create the geometry                                         */
  /***************************************************************/
  RWGGeometry::UseHRWGFunctions=false;
  RWGGeometry *G=new RWGGeometry(GeoFile);
 
  HMatrix *M=G->AllocateBEMMatrix();
  HVector *KN=G->AllocateRHSVector();

  /***************************************************************/
  /* parse the port list and plot if requested *******************/
  /***************************************************************/
  RWGPortList *PortList=ParsePortFile(G, PortFile);
  int NumPorts = PortList->Ports.size();
  if (PlotPorts)
   { fprintf(stderr,"--PlotPorts option was specified; plotting ports ONLY.\n");
     PlotPortsInGMSH(G, PortList, "%s.ports.pp", GetFileBase(G->GeoFileName));
     fprintf(stderr,"RF ports plotted to file %s.ports.pp.\n",GetFileBase(G->GeoFileName));
     fprintf(stderr,"Thank you for your support.\n");
     exit(0);
   }

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
for(int npe=0; npe<PortList->PortEdges.size(); npe++)
{
RWGPortEdge *PE=PortList->PortEdges[npe];
RWGSurface *S=G->Surfaces[PE->ns];
double X[3]={0.1,0.2,0.0};
cdouble p[0]; cdouble a[3];
cdouble dp[3]; cdouble da[3][3];
cdouble ddp[3][3]; cdouble dcurla[3][3];
GetReducedPotentials_Nearby(S, PE->ne, X, 1.0, p, a, dp, da, ddp, dcurla);
printf("pe=%i: p=%e,%e foryaf.\n",npe,real(p[0]),imag(p[0]));
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  /***************************************************************/
  /* parse frequency specifications to create the FreqList vector*/
  /***************************************************************/
  HVector *FreqList=0;
  if (FreqFile) 
   FreqList = new HVector(FreqFile);
  else if (nFreqs)
   FreqList=new HVector(nFreqs, Freqs);
  else if (nMinFreq || nMaxFreq || nNumFreqs)
   { if ( !nMinFreq || !nMaxFreq || !nNumFreqs )
      ErrExit("--MinFreq, --MaxFreq, --NumFreqs must be all present or all absent");
     FreqList = LogFreq ? LogSpace(MinFreq, MaxFreq, NumFreqs) : LinSpace(MinFreq, MaxFreq, NumFreqs);
   }

  /***************************************************************/ 
  /* parse the port-current list if that was specified           */ 
  /***************************************************************/
  HMatrix *PCMatrix=0;
  cdouble *PortCurrents=0;
  if (PCFile)
   { if (FreqList)
      ErrExit("--PortCurrentFile may not be specified together with a frequency specification");
     PCMatrix=new HMatrix(PCFile);
     if (PCMatrix->NC != NumPorts+1)
      ErrExit("%s: expected %i columns (1 frequency, %i ports)\n",PCFile,NumPorts+1,NumPorts);
     FreqList=new HVector(PCMatrix->NR, LHM_REAL, PCMatrix->GetColumnPointer(0));
     PortCurrents = new cdouble[NumPorts];
   }

  /***************************************************************/
  /* sanity check input arguments ********************************/
  /***************************************************************/
  if (PCFile==0 && FreqList->N!=0 && (ZParms==0 && SParms==0) )
   OSUsage(argv[0],OSArray,"--ZParameters and/or --SParameters must be specified if a frequency specification is present");
  if (PCFile!=0 && (ZParms!=0 || SParms!=0) )
   OSUsage(argv[0],OSArray,"--ZParameters and --SParameters may not be used with --PortCurrentFile");
  if (PCFile!=0 && nEPFiles==0 && nFVMeshes==0)
   OSUsage(argv[0],OSArray,"--EPFile or --FVMesh must be specified if --PortCurrentFile is specified");
  if (PCFile==0 && (nEPFiles!=0 || nFVMeshes!=0) )
   OSUsage(argv[0],OSArray,"--EPFile and --FVMesh require --PortCurrentFile");

  /***************************************************************/
  /* loop over frequencies ***************************************/
  /***************************************************************/
  for (int nf=0; nf<FreqList->N; nf++)
   { 
     double Freq = FreqList->GetEntryD(nf);

     /*--------------------------------------------------------------*/
     /* assemble and factorize the BEM matrix at this frequency      */
     /*--------------------------------------------------------------*/
     cdouble Omega=FREQ2OMEGA * Freq;
     Log("Assembling BEM matrix at f=%g GHz...",Freq);
     G->AssembleBEMMatrix(Omega, M);
     Log("Factorizing...",Freq);
     M->LUFactorize();

     /*--------------------------------------------------------------*/
     /* switch off to output modules to handle various calculations -*/
     /*--------------------------------------------------------------*/
     if (ZParms || SParms)
      ComputeSZParms(G, PortList, Omega, M, KN, FileBase, SParms);

     /*--------------------------------------------------------------*/
     /*- if the user gave us driving port currents and asked for the */
     /*- radiated fields at a list of points/and or on a user-       */
     /*- specified flux-mesh surface, handle that                    */
     /*--------------------------------------------------------------*/
     if (nEPFiles!=0 || nFVMeshes!=0)
      { 
        Log(" Computing radiated fields..."); 

        /*--------------------------------------------------------------*/
        /* get the RHS vector corresponding to the user-specified port  */
        /* currents at this frequency, then solve the BEM sytem and     */
        /* handle post-processing field computations                    */
        /*--------------------------------------------------------------*/
        PCMatrix->GetEntries(nf,"1:end",PortCurrents);
        KN->Zero();
        Log("  assembling RHS vector");
        AddPortContributionsToRHS(G, PortList, PortCurrents, Omega, KN);
        Log("  solving the BEM system");
        M->LUSolve(KN);
         
        for(int n=0; n<nEPFiles; n++)
         ProcessEPFile(G, KN, Omega, PortList, PortCurrents, EPFiles[n], FileBase);

        for(int n=0; n<nFVMeshes; n++)
         ProcessFVMesh(G, KN, Omega, PortList, PortCurrents, FVMeshes[n], FileBase);
      } // if (nEPFiles!=0 || nFVMeshes!=0)

   }; // for (nf=0..)

  printf("Thank you for your support.\n");
}
