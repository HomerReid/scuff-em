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
#include <libSubstrate.h>
#include <RFSolver.h>

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
namespace scuff{
void GetReducedPotentials_Nearby(RWGSurface *S, const int ne,
                                 const double X0[3],  const cdouble k,
                                 cdouble *p, cdouble a[3],
                                 cdouble dp[3], cdouble da[3][3],
                                 cdouble ddp[3][3], cdouble dcurla[3][3],
                                 bool *IncludeTerm=0);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *GetFieldPSDMatrix(RWGGeometry *G, RWGPortList *PortList,
                           cdouble Omega, HVector *KN, cdouble *PortCurrents, HMatrix *FieldPSDMatrix)
{
  HMatrix *XMatrix = new HMatrix(3, G->TotalPanels);
  int NX = G->TotalPanels;
  for(int ns=0, nx=0; ns<G->NumSurfaces; ns++)
   for(int np=0; np<G->Surfaces[ns]->NumPanels; np++, nx++)
    XMatrix->SetEntriesD(":",nx,G->Surfaces[ns]->Panels[np]->Centroid);

  HMatrix *BFRPFMatrix=0, *PortRPFMatrix=0;
  GetMOIRPFMatrices(G, G->Substrate, PortList, Omega, XMatrix,
                    &BFRPFMatrix, &PortRPFMatrix);
  
  HMatrix *BFPFMatrix = new HMatrix(NPFC, NX, LHM_COMPLEX);
  HVector BFPFVector(NPFC*NX, BFPFMatrix->ZM);
  BFRPFMatrix->Apply(KN, &BFPFVector);

  HMatrix *PortPFMatrix = new HMatrix(NPFC, NX, LHM_COMPLEX);
  HVector PortPFVector(NPFC*NX, PortPFMatrix->ZM);
  HVector PCVector(PortList->Ports.size(), PortCurrents);
  PortRPFMatrix->Apply(&PCVector, &PortPFVector);

  if (FieldPSDMatrix && (FieldPSDMatrix->NR!=NX || FieldPSDMatrix->NC!=13) )
   { if (FieldPSDMatrix) delete FieldPSDMatrix;
     FieldPSDMatrix=0;
   }
  if (!FieldPSDMatrix)
   FieldPSDMatrix=new HMatrix(NX, 13, LHM_COMPLEX);
  FieldPSDMatrix->Zero();

  for(int nx=0; nx<NX; nx++)
   { FieldPSDMatrix->SetEntry(nx,4,BFPFMatrix->GetEntry(_PF_PHI,nx));
     FieldPSDMatrix->SetEntry(nx,5,BFPFMatrix->GetEntry(_PF_AX,nx));
     FieldPSDMatrix->SetEntry(nx,6,BFPFMatrix->GetEntry(_PF_AY,nx));
     FieldPSDMatrix->SetEntry(nx,7,BFPFMatrix->GetEntry(_PF_AZ,nx));
     FieldPSDMatrix->SetEntry(nx,8,PortPFMatrix->GetEntry(_PF_PHI,nx));
     FieldPSDMatrix->SetEntry(nx,9,PortPFMatrix->GetEntry(_PF_AX,nx));
     FieldPSDMatrix->SetEntry(nx,10,PortPFMatrix->GetEntry(_PF_AY,nx));
     FieldPSDMatrix->SetEntry(nx,11,PortPFMatrix->GetEntry(_PF_AZ,nx));
   }

  delete BFRPFMatrix;
  delete BFPFMatrix;
  delete PortRPFMatrix;
  delete PortPFMatrix;
  delete XMatrix;
  return FieldPSDMatrix;
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
  char *SubstrateFile=0;
  char *EpsStr = 0;
  double h     = 0.0;
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
     {"SubstrateFile",  PA_STRING,  1, 1,       (void *)&SubstrateFile,   0,        "substrate definition file"},
     {"Eps",            PA_STRING,  1, 1,       (void *)&EpsStr,     0,             "substrate permittivity"},
     {"h",              PA_DOUBLE,  1, 1,       (void *)&h,          0,             "substrate thickness"},
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
  /* process substrate-related options                           */
  /***************************************************************/
  if (SubstrateFile)
   G->Substrate = new LayeredSubstrate(SubstrateFile);
  else if (EpsStr!=0)
   { 
     char SubstrateDefinition[1000];
     if (h==0.0) // no ground plane
      snprintf(SubstrateDefinition,1000,"0.0 CONST_EPS_%s\n",EpsStr);
     else
      snprintf(SubstrateDefinition,1000,"0.0 CONST_EPS_%s\n%e GROUNDPLANE\n",EpsStr,-h);
     G->Substrate=CreateLayeredSubstrate(SubstrateDefinition);
   }

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
  //if (PCFile!=0 && nEPFiles==0 && nFVMeshes==0)
   //OSUsage(argv[0],OSArray,"--EPFile or --FVMesh must be specified if --PortCurrentFile is specified");
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
     AssembleMOIMatrix(G, G->Substrate, Omega, M);
     Log("Factorizing...");
     M->LUFactorize();

     /*--------------------------------------------------------------*/
     /* switch off to output modules to handle various calculations -*/
     /*--------------------------------------------------------------*/
     if (ZParms || SParms)
      ComputeSZParms(G, PortList, Omega, M, FileBase, SParms);

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
        Log("  assembling RHS vector");
        GetPortContributionToRHS(G, PortList, PortCurrents, Omega, KN);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
static HMatrix *PSDMatrix=0;
PSDMatrix=G->GetPanelSourceDensities(Omega, KN, PSDMatrix);
G->PlotSurfaceCurrents(PSDMatrix, Omega, "/tmp/BeforeWO.pp");

SetDefaultCD2SFormat("{%+.2e,%+.2e}");
cdouble Charge=0.0, DPM[3]={0.0, 0.0, 0.0};
for(int np=0; np<PSDMatrix->NR; np++)
 { cdouble Q = PSDMatrix->GetEntry(np,3) * PSDMatrix->GetEntry(np,4);
   Charge += Q;
   DPM[0] += PSDMatrix->GetEntry(np,0) * Q;
   DPM[1] += PSDMatrix->GetEntry(np,1) * Q;
   DPM[2] += PSDMatrix->GetEntry(np,2) * Q;
 }
printf("{Q,P} before (WO) = %s %s %s %s\n",CD2S(Charge),CD2S(DPM[0]),CD2S(DPM[1]),CD2S(DPM[2]));

AddPortContributionsToPSD(G, PortList, PortCurrents, Omega, PSDMatrix);
G->PlotSurfaceCurrents(PSDMatrix, Omega, "/tmp/BeforeW.pp");

Charge=DPM[0]=DPM[1]=DPM[2]=0.0;
for(int np=0; np<PSDMatrix->NR; np++)
 { cdouble Q = PSDMatrix->GetEntry(np,3) * PSDMatrix->GetEntry(np,4);
   Charge += Q;
   DPM[0] += PSDMatrix->GetEntry(np,0) * Q;
   DPM[1] += PSDMatrix->GetEntry(np,1) * Q;
   DPM[2] += PSDMatrix->GetEntry(np,2) * Q;
 }
printf("{Q,P} before (W) = %s %s %s %s\n",CD2S(Charge),CD2S(DPM[0]),CD2S(DPM[1]),CD2S(DPM[2]));

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
        Log("  solving the BEM system");
        M->LUSolve(KN);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

PSDMatrix=G->GetPanelSourceDensities(Omega, KN, PSDMatrix);
G->PlotSurfaceCurrents(PSDMatrix, Omega, "/tmp/AfterWO.pp");

Charge=DPM[0]=DPM[1]=DPM[2]=0.0;
for(int np=0; np<PSDMatrix->NR; np++)
 { cdouble Q = PSDMatrix->GetEntry(np,3) * PSDMatrix->GetEntry(np,4);
   Charge += Q;
   DPM[0] += PSDMatrix->GetEntry(np,0) * Q;
   DPM[1] += PSDMatrix->GetEntry(np,1) * Q;
   DPM[2] += PSDMatrix->GetEntry(np,2) * Q;
 }
printf("{Q,P} after (WO) = %s %s %s %s\n",CD2S(Charge),CD2S(DPM[0]),CD2S(DPM[1]),CD2S(DPM[2]));

AddPortContributionsToPSD(G, PortList, PortCurrents, Omega, PSDMatrix);
G->PlotSurfaceCurrents(PSDMatrix, Omega, "/tmp/AfterW.pp");

Charge=DPM[0]=DPM[1]=DPM[2]=0.0;
for(int np=0; np<PSDMatrix->NR; np++)
 { cdouble Q = PSDMatrix->GetEntry(np,3) * PSDMatrix->GetEntry(np,4);
   Charge += Q;
   DPM[0] += PSDMatrix->GetEntry(np,0) * Q;
   DPM[1] += PSDMatrix->GetEntry(np,1) * Q;
   DPM[2] += PSDMatrix->GetEntry(np,2) * Q;
 }
printf("{Q,P} after (W) = %s %s %s %s\n",CD2S(Charge),CD2S(DPM[0]),CD2S(DPM[1]),CD2S(DPM[2]));

static HMatrix *FieldPSDMatrix=0;
FieldPSDMatrix=GetFieldPSDMatrix(G, PortList, Omega, KN, PortCurrents, FieldPSDMatrix);
G->PlotSurfaceCurrents(FieldPSDMatrix, Omega, "/tmp/BFFields.pp");
cdouble iwAdJ[2]={0.0,0.0}, SigmaPhi[2]={0.0,0.0};
cdouble IW=II*Omega;
for(int np=0; np<FieldPSDMatrix->NC; np++)
 { double Area       = PSDMatrix->GetEntryD(np,3);
   cdouble Sigma     = PSDMatrix->GetEntry(np, 4);
   cdouble K[3];       PSDMatrix->GetEntries(np,"5:7",K);
   cdouble PhiBF     = FieldPSDMatrix->GetEntry(np, 4);
   cdouble ABF[3];     FieldPSDMatrix->GetEntries(np,"5:7",ABF);
   cdouble PhiPort   = FieldPSDMatrix->GetEntry(np, 8);
   cdouble APort[3];   FieldPSDMatrix->GetEntries(np,"9:11",APort);
   iwAdJ[0]    += Area*IW*(K[0]*ABF[0] + K[1]*ABF[1] + K[2]*ABF[2]);
   iwAdJ[1]    += Area*IW*(K[0]*APort[0] + K[1]*APort[1] + K[2]*APort[2]);
   SigmaPhi[0] += Area*Sigma*PhiBF;
   SigmaPhi[1] += Area*Sigma*PhiPort;
 }
FILE *f=fopen("/tmp/EJ.out","a");
fprintf(f,"%e ",real(Omega));
fprintf(f,"%e %e ",real(iwAdJ[0]),imag(iwAdJ[0]));
fprintf(f,"%e %e ",real(iwAdJ[1]),imag(iwAdJ[1]));
fprintf(f,"%e %e ",real(SigmaPhi[0]),imag(SigmaPhi[0]));
fprintf(f,"%e %e ",real(SigmaPhi[1]),imag(SigmaPhi[1]));
fprintf(f,"\n");
fclose(f);
for(int np=0; np<FieldPSDMatrix->NR; np++)
 { FieldPSDMatrix->SetEntry(np, 4, FieldPSDMatrix->GetEntry(np,8));
   FieldPSDMatrix->SetEntry(np, 5, FieldPSDMatrix->GetEntry(np,9));
   FieldPSDMatrix->SetEntry(np, 6, FieldPSDMatrix->GetEntry(np,10));
   FieldPSDMatrix->SetEntry(np, 7, FieldPSDMatrix->GetEntry(np,11));
 }
G->PlotSurfaceCurrents(FieldPSDMatrix, Omega, "/tmp/PortFields.pp");
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
         
        for(int n=0; n<nEPFiles; n++)
         ProcessEPFile(G, KN, Omega, PortList, PortCurrents, EPFiles[n], FileBase);

        for(int n=0; n<nFVMeshes; n++)
         ProcessFVMesh(G, KN, Omega, PortList, PortCurrents, FVMeshes[n], FileBase);
      } // if (nEPFiles!=0 || nFVMeshes!=0)

   }; // for (nf=0..)

  printf("Thank you for your support.\n");
}
