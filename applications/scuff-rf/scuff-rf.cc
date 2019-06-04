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
 * scuff-rf.cc   -- scuff-EM module for modeling of RF and microwave
 *               -- structures
 * 
 * homer reid    -- 9/2011
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
#define II cdouble(0.0,1.0)

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
extern int WriteLogFile;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

using namespace scuff;

#define MAXFREQ  10    // max number of frequencies
#define MAXCACHE 10    // max number of cache files for preload
#define MAXSTR   1000
#define MAXEPF   10    // max number of field visualization files
#define MAXFVM   10    // max number of field visualization meshes

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ProcessEPFile(RWGGeometry *G, HVector *KN, cdouble Omega, 
                   RWGPort **Ports, int NumPorts, cdouble *PortCurrents, 
                   char *EPFile, char *FileBase);

void ProcessFVMesh(RWGGeometry *G, HVector *KN, cdouble Omega,
                   RWGPort **Ports, int NumPorts, cdouble *PortCurrents, 
                   char *FVMesh, char *FileBase);

/***************************************************************/
/***************************************************************/
/***************************************************************/
void HVConcat(HVector *V1, HVector *V2)
{ 
  if ( V1==0 || V2==0 || V2->N==0 ) return;

  int NewN=V1->N + V2->N;
  V1->DV = (double *)realloc(V1->DV, NewN*sizeof(double));
  memcpy(V1->DV + V1->N, V2->DV, (V2->N)*sizeof(double));
  V1->N=NewN;

}

/***************************************************************/
/* main function   *********************************************/
/***************************************************************/  
int main(int argc, char *argv[])
{
  /***************************************************************/
  /** process command-line arguments *****************************/
  /***************************************************************/
  char *GeoFile=0;
  char *PortFile=0;
  char *PCFile=0;
  char *EPFiles[MAXEPF];                int nEPFiles;
  char *FVMeshes[MAXFVM];               int nFVMeshes;
  int PlotPorts=0;
  int ZParameters=0; 
  int SParameters=0;
  int Moments=0;
  double FrequencyValues[MAXFREQ];	int nFrequency;
  double MinFreq;			int nMinFreq;
  double MaxFreq;			int nMaxFreq;
  int NumFreqs;				int nNumFreqs;
  int LogFreq=0;
  char *FreqFile=0;
  char *Cache=0;
  char *ReadCache[MAXCACHE];         int nReadCache;
  char *WriteCache=0;
  char *ContribOnly=0;
  char *FileBase=0;
  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { {"geometry",       PA_STRING,  1, 1,       (void *)&GeoFile,    0,             "geometry file"},
//
     {"portfile",       PA_STRING,  1, 1,       (void *)&PortFile,   0,             "port file"},
     {"PlotPorts",      PA_BOOL,    0, 1,       (void *)&PlotPorts,  0,             "generate port visualization file"},
//
     {"frequency",      PA_DOUBLE,  1, MAXFREQ, (void *)FrequencyValues,  &nFrequency,   "frequency (GHz)"},
     {"minfreq",        PA_DOUBLE,  1, 1,       (void *)&MinFreq,    &nMinFreq,     "starting frequency"},
     {"maxfreq",        PA_DOUBLE,  1, 1,       (void *)&MaxFreq,    &nMaxFreq,     "ending frequency"},
     {"numfreqs",       PA_INT,     1, 1,       (void *)&NumFreqs,   &nNumFreqs,    "number of frequencies"},
     {"logfreq",        PA_BOOL,    0, 1,       (void *)&LogFreq,    0,             "use logarithmic frequency steps"},
     {"freqfile",       PA_STRING,  1, 1,       (void *)&FreqFile,   0,             "list of frequencies"},
//
     {"ZParameters",    PA_BOOL,    0, 1,       (void *)&ZParameters, 0,            "output z parameters"},
     {"SParameters",    PA_BOOL,    0, 1,       (void *)&SParameters, 0,            "output s parameters"},
     {"Moments",        PA_BOOL,    0, 1,       (void *)&Moments,     0,            "output dipole moments"},
//
     {"portcurrentfile", PA_STRING,  1, 1,      (void *)&PCFile,     0,             "port current file"},
     {"EPFile",         PA_STRING,  1, MAXEPF,  (void *)EPFiles,     &nEPFiles,             "list of evaluation points"},
//
     {"FVMesh",         PA_STRING,  1, MAXFVM,  (void *)FVMeshes,    &nFVMeshes,    "field visualization mesh"},
//
     {"FileBase",       PA_STRING,  1, 1,       (void *)&FileBase,     0,           "base name for output files"},
//
     {"Cache",          PA_STRING,  1, 1,       (void *)&Cache,      0,             "read/write cache"},
     {"ReadCache",      PA_STRING,  1, MAXCACHE,(void *)ReadCache,   &nReadCache,   "read cache"},
     {"WriteCache",     PA_STRING,  1, 1,       (void *)&WriteCache, 0,             "write cache"},
//
     {"WriteLogFile",   PA_BOOL,    0, 1,       (void *)&WriteLogFile, 0,           "write new log file"},
//
//
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFile==0)
   OSUsage(argv[0],VERSIONSTRING,OSArray,"--geometry option is mandatory");
  if (!FileBase)
   FileBase=strdup(GetFileBase(GeoFile));

  if (ContribOnly) SetContribOnly(ContribOnly);

  /***************************************************************/
  /* create the log file *****************************************/
  /***************************************************************/
  SetLogFileName("scuff-rf.log");
  int narg;
  Log("%s " VERSIONSTRING " running on %s with arguments ",argv[0],getenv("HOST"));
  Log("%s ",argv[0]);
  for(narg=1; narg<argc; narg++) 
   LogC("%s ",argv[narg]);

  /***************************************************************/
  /* create the geometry                                         */
  /***************************************************************/
  RWGGeometry::UseHRWGFunctions=false;
  RWGGeometry *G=new RWGGeometry(GeoFile);
 
  HMatrix *M=G->AllocateBEMMatrix();
  HVector *KN=G->AllocateRHSVector();

  /***************************************************************/ 
  /* parse the port list *****************************************/ 
  /***************************************************************/ 
  if (PortFile==0)
   OSUsage(argv[0],VERSIONSTRING,OSArray,"--PortFileName option is mandatory");

  int NumPorts;
  RWGPort **Ports=ParsePortFile(G, PortFile, &NumPorts);

  if (PlotPorts)
   { 
     fprintf(stderr,"--PlotPorts option was specified; plotting ports ONLY.\n");
     for(int np=0; np<NumPorts; np++)
      { fprintf(stderr," Port %2i, positive edge: %i edges (perimeter %e)\n", NumPorts, 
                         Ports[np]->NumPEdges, Ports[np]->PPerimeter);
        fprintf(stderr," Port %2i, negative edge: %i edges (perimeter %e)\n", NumPorts, 
                         Ports[np]->NumMEdges, Ports[np]->MPerimeter);
      };
     
     PlotPortsInGMSH(Ports, NumPorts, "%s.ports.pp", GetFileBase(G->GeoFileName));
     fprintf(stderr,"RF ports plotted to file %s.ports.pp.\n",GetFileBase(G->GeoFileName));
     fprintf(stderr,"Thank you for your support.\n");
     exit(1);
   };

  /***************************************************************/ 
  /* parse the port-current list if that was specified           */ 
  /***************************************************************/
  HMatrix *PCList=0;
  HVector *FreqList=0;
  int nf;
  if (PCFile)
   { 
     PCList=new HMatrix(PCFile,LHM_TEXT);
     if (PCList->ErrMsg)
      ErrExit("%s: %s",PCFile,PCList->ErrMsg);
     if (PCList->NC != 2*NumPorts+1)
      ErrExit("%s: expected %i columns (2*%i ports + 1)\n",PCFile,2*NumPorts+1,NumPorts);
     FreqList=new HVector(PCList->NR);
     for(nf=0; nf<FreqList->N; nf++)
      FreqList->SetEntry(nf,PCList->GetEntryD(nf,0));
   };

  /***************************************************************/
  /* parse frequency specifications to create the FreqList vector*/
  /***************************************************************/
  HVector *FV1=0, *FV2=0, *FV3=0;
  if (FreqFile)
   { FV1=new HVector(FreqFile);
     if (FV1->ErrMsg)
      ErrExit(FV1->ErrMsg);
   };
  if (nFrequency)
   { FV2=new HVector(nFrequency);
     memcpy(FV2->DV, FrequencyValues, nFrequency*sizeof(double));
   };
  if (nMinFreq || nMaxFreq || nNumFreqs)
   { if ( !nMinFreq || !nMaxFreq || !nNumFreqs )
      ErrExit("--MinFreq, --MaxFreq, --NumFreqs must be all present or all absent");
     FV3 = LogFreq ? LogSpace(MinFreq, MaxFreq, NumFreqs) : LinSpace(MinFreq, MaxFreq, NumFreqs);
   };
  if ( (FV1 || FV2 || FV3) )
   { if (PCFile!=0 )
      ErrExit("--portcurrentfile may not be combined with a frequency specification"); 
     FreqList = new HVector(0);
     HVConcat(FreqList, FV1);
     HVConcat(FreqList, FV2);
     HVConcat(FreqList, FV3);
     NumFreqs=FreqList->N;
   };

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  cdouble *PCVector = new cdouble[NumPorts];
  for(int np=0; np<NumPorts; np++)
   { KN->Zero();
     memset(PCVector, 0, NumPorts*sizeof(cdouble));
     PCVector[np]=1.0;
     AddPortContributionsToRHS(G, Ports, NumPorts, PCVector, FREQ2OMEGA * FreqList->GetEntry(0), KN);
     FILE *f=vfopen("/tmp/Port%i.Old","w",np);
     for(int ne=0; ne<KN->N; ne++)
      fprintf(f,"%5i   %+.12e   %+.12e\n",ne,real(KN->GetEntry(ne)),imag(KN->GetEntry(ne)));
     fclose(f);
    }
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  /***************************************************************/
  /* sanity check input arguments ********************************/
  /***************************************************************/
  if ( PCFile==0 && NumFreqs!=0 && (ZParameters==0 && SParameters==0) )
   OSUsage(argv[0],VERSIONSTRING,OSArray,"--zparameters and/or --sparameters must be specified if a frequency specification is present");
  if (PCFile!=0 && (ZParameters!=0 || SParameters!=0) )
   OSUsage(argv[0],VERSIONSTRING,OSArray,"--zparameters and --sparameters may not be used with --portcurrentfile");
  if (PCList!=0 && nEPFiles==0 && nFVMeshes==0)
   OSUsage(argv[0],VERSIONSTRING,OSArray,"--EPFile or --FVMesh must be specified if --portcurrentfile is specified");

  /***************************************************************/
  /* create output files *****************************************/
  /***************************************************************/
  char TimeString[200];
  time_t MyTime=time(0);
  struct tm *MyTm=localtime(&MyTime);
  strftime(TimeString,30,"%D::%T",MyTm);

  FILE *ZParFile=0;
  char ZParFileName[1000];
  HMatrix *ZMatrix=0;
  if (ZParameters)
   { 
     snprintf(ZParFileName,999,"%s.zparms",GetFileBase(G->GeoFileName));
     ZParFile=CreateUniqueFile(ZParFileName,1,ZParFileName);
     setbuf(ZParFile,0);
     fprintf(ZParFile,"# %s ran on %s (%s)\n",argv[0],getenv("HOST"),TimeString);
     fprintf(ZParFile,"# command-line arguments: \n#");
     for(narg=0; narg<argc; narg++)
      fprintf(ZParFile,"%s ",argv[narg]);
     fprintf(ZParFile,"\n#\n");
     fprintf(ZParFile,"# columns:\n");
     fprintf(ZParFile,"# 1 frequency (GHz)\n");
     fprintf(ZParFile,"# 2,3 real,imag Z_{11} \n");
     if (NumPorts>1)
      { fprintf(ZParFile,"# 4,5 real,imag Z_{12} \n");
        fprintf(ZParFile,"# ...\n");
        fprintf(ZParFile,"# %i, %i real, imag Z_{%i%i}\n",
                   2*NumPorts*NumPorts, 2*NumPorts*NumPorts+1, 
                   NumPorts+1,NumPorts+1);
      };

     ZMatrix=new HMatrix(NumPorts, NumPorts, LHM_COMPLEX);
   };

  FILE *SParFile=0;
  char SParFileName[1000];
  HMatrix *SMatrix=0;
  if (SParameters)
   { 
     snprintf(SParFileName,999,"%s.sparms",GetFileBase(G->GeoFileName));
     SParFile=CreateUniqueFile(SParFileName,1,SParFileName);
     setbuf(SParFile,0);
     fprintf(SParFile,"# %s ran on %s (%s)\n",argv[0],getenv("HOST"),TimeString);
     fprintf(SParFile,"# command-line arguments: \n#");
     for(narg=0; narg<argc; narg++)
      fprintf(SParFile,"%s ",argv[narg]);
     fprintf(SParFile,"\n#\n"); 
     fprintf(SParFile,"# columns:\n");
     fprintf(SParFile,"# 1 frequency (GHz)\n");
     fprintf(SParFile,"# 2,3 real,imag S_{11} \n");
     if (NumPorts>1)
      { fprintf(SParFile,"# 4,5 real,imag S_{12} \n");
        fprintf(SParFile,"# ...\n");
        fprintf(SParFile,"# %i, %i real, imag S_{%i%i}\n",
                   2*NumPorts*NumPorts,2*NumPorts*NumPorts+1, 
                   NumPorts+1,NumPorts+1);
      };

     SMatrix=new HMatrix(NumPorts, NumPorts, LHM_COMPLEX);
     if (ZMatrix==0) //ZMatrix is needed if s-parameters are computed
      ZMatrix=new HMatrix(NumPorts, NumPorts, LHM_COMPLEX); 
   };

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  FILE *MomentFile=0;
  //char PSDFileName[1000];
  HMatrix *PM=0, *PSD=0;
  if (Moments)
   { PM=new HMatrix(G->NumSurfaces, 6, LHM_COMPLEX);
     MomentFile=vfopen("%s.moments","w",GetFileBase(G->GeoFileName));
     PSD=new HMatrix(G->TotalPanels, 13, LHM_COMPLEX);
   }

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

  /***************************************************************/  
  /* sweep over frequencies                                      */
  /* note: Freq is measured in GHz, while Omega is measured in   */
  /* our natural unit for angular frequency, which is            */
  /* c / (1 mm) = 3e11 rad/sec                                   */
  /***************************************************************/
  int np, npp;
  double Freq; 
  cdouble Omega;
  cdouble *PortCurrents=new cdouble[NumPorts]; 
  cdouble *PortVoltages=new cdouble[NumPorts];
  for (nf=0, Freq=FreqList->GetEntryD(0); nf<FreqList->N; Freq=FreqList->GetEntryD(++nf) )
   { 
      /*--------------------------------------------------------------*/
      /* assemble and factorize the BEM matrix at this frequency      */
      /*--------------------------------------------------------------*/
      Omega=FREQ2OMEGA * Freq;
      Log("Assembling BEM matrix at f=%g GHz...",Freq);
      G->AssembleBEMMatrix(Omega, M);
      Log("Factorizing...",Freq);
      M->LUFactorize();

      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      if (WriteCache)
       { StoreCache( WriteCache );
         WriteCache=0;       
       };

      /*--------------------------------------------------------------*/
      /* if the user asked us to compute Z- or S- parameters, compute */
      /* the rows of the Z-matrix at this frequency                   */
      /*--------------------------------------------------------------*/
      if (ZParameters || SParameters)
       { for(np=0; np<NumPorts; np++)
          { 
            Log(" Computing row %i of the impedance matrix:",np+1);

            Log("  assembling RHS vector");
            memset(PortCurrents, 0, NumPorts*sizeof(cdouble));
            PortCurrents[np]=1.0;
            KN->Zero();
            AddPortContributionsToRHS(G, Ports, NumPorts, PortCurrents, Omega, KN);

            /*--------------------------------------------------------------*/
            /*--------------------------------------------------------------*/
            /*--------------------------------------------------------------*/
            if (Moments)
             { G->GetDipoleMoments(Omega, KN, PM);
               SetDefaultCD2SFormat("%+.6e %+.6e");
               for(int ns=0; ns<G->NumSurfaces; ns++)
                fprintf(MomentFile,"%e %i BEFORE %s %s %s %s %s %s %s \n",
                                    real(Omega),G->TotalPanels,G->Surfaces[ns]->Label,
                                    CD2S(PM->GetEntry(ns,0)), CD2S(PM->GetEntry(ns,1)), CD2S(PM->GetEntry(ns,2)),
                                    CD2S(PM->GetEntry(ns,3)), CD2S(PM->GetEntry(ns,4)), CD2S(PM->GetEntry(ns,5)));

               G->GetPanelSourceDensities(Omega, KN, PSD);
               G->PlotSurfaceCurrents(PSD, Omega, "%s_Before(WO).pp",GetFileBase(G->GeoFileName));

               SetDefaultCD2SFormat("{%+.2e,%+.2e}");
               cdouble Charge=0.0, DPM[3]={0.0, 0.0, 0.0};
               for(np=0; np<PSD->NR; np++)
                { cdouble Q = PSD->GetEntry(np,3) * PSD->GetEntry(np,4);
                  Charge += Q;
                  DPM[0] += PSD->GetEntry(np,0) * Q;
                  DPM[1] += PSD->GetEntry(np,1) * Q;
                  DPM[2] += PSD->GetEntry(np,2) * Q;
                }
               printf("{Q,P} before (WO) = %s %s %s %s\n",CD2S(Charge),CD2S(DPM[0]),CD2S(DPM[1]),CD2S(DPM[2]));

               AddPortContributionsToPSD(G, Ports, NumPorts, PortCurrents, Omega, PSD);
               G->PlotSurfaceCurrents(PSD, Omega, "%s_Before(W).pp",GetFileBase(G->GeoFileName));

               Charge=0.0; DPM[0]=DPM[1]=DPM[2]=0.0;
               for(np=0; np<PSD->NR; np++)
                { cdouble Q = PSD->GetEntry(np,3) * PSD->GetEntry(np,4);
                  Charge += Q;
                  DPM[0] += PSD->GetEntry(np,0) * Q;
                  DPM[1] += PSD->GetEntry(np,1) * Q;
                  DPM[2] += PSD->GetEntry(np,2) * Q;
                }
               printf("{Q,P} before (W) = %s %s %s %s\n",CD2S(Charge),CD2S(DPM[0]),CD2S(DPM[1]),CD2S(DPM[2]));
             };
                
            /*--------------------------------------------------------------*/
            /*--------------------------------------------------------------*/
            /*--------------------------------------------------------------*/
            Log("  solving the BEM system");
            M->LUSolve(KN);

            /*--------------------------------------------------------------*/
            /*--------------------------------------------------------------*/
            /*--------------------------------------------------------------*/
            if (Moments)
             { G->GetDipoleMoments(Omega, KN, PM);
               SetDefaultCD2SFormat("%+.6e %+.6e");
               for(int ns=0; ns<G->NumSurfaces; ns++)
                fprintf(MomentFile,"%e %i AFTER %s %s %s %s %s %s %s \n",
                                    real(Omega),G->TotalPanels,G->Surfaces[ns]->Label,
                                    CD2S(PM->GetEntry(ns,0)), CD2S(PM->GetEntry(ns,1)), CD2S(PM->GetEntry(ns,2)),
                                    CD2S(PM->GetEntry(ns,3)), CD2S(PM->GetEntry(ns,4)), CD2S(PM->GetEntry(ns,5)));

               G->GetPanelSourceDensities(Omega, KN, PSD);
               G->PlotSurfaceCurrents(PSD, Omega, "%s_AfterWO.pp",GetFileBase(G->GeoFileName));

               SetDefaultCD2SFormat("{%+.2e,%+.2e}");
               cdouble Charge=0.0, DPM[3]={0.0, 0.0, 0.0};
               for(np=0; np<PSD->NR; np++)
                { cdouble Q = PSD->GetEntry(np,3) * PSD->GetEntry(np,4);
                  Charge += Q;
                  DPM[0] += PSD->GetEntry(np,0) * Q;
                  DPM[1] += PSD->GetEntry(np,1) * Q;
                  DPM[2] += PSD->GetEntry(np,2) * Q;
                }
               printf("{Q,P} after (WO) = %s %s %s %s\n",CD2S(Charge),CD2S(DPM[0]),CD2S(DPM[1]),CD2S(DPM[2]));

               AddPortContributionsToPSD(G, Ports, NumPorts, PortCurrents, Omega, PSD);
               G->PlotSurfaceCurrents(PSD, Omega, "%s_AfterW.pp",GetFileBase(G->GeoFileName));

               Charge=DPM[0]=DPM[1]=DPM[2]=0.0;
               for(np=0; np<PSD->NR; np++)
                { cdouble Q = PSD->GetEntry(np,3) * PSD->GetEntry(np,4);
                  Charge += Q;
                  DPM[0] += PSD->GetEntry(np,0) * Q;
                  DPM[1] += PSD->GetEntry(np,1) * Q;
                  DPM[2] += PSD->GetEntry(np,2) * Q;
                }
               printf("{Q,P} after (W) = %s %s %s %s\n",CD2S(Charge),CD2S(DPM[0]),CD2S(DPM[1]),CD2S(DPM[2]));
             };

            Log("  computing port voltages");
            GetPortVoltages(G, KN, Ports, NumPorts, PortCurrents, Omega, PortVoltages);

            /* note: the entry in the Z-matrix is the complex conjugate */
            /* of the measured port voltage, because the Z-matrix is    */
            /* defined using the usual circuit theory convention in     */
            /* which all quantities have time dependence exp(+iwt),     */
            /* whereas scuff-EM uses the opposite sign convention.      */
            for(npp=0; npp<NumPorts; npp++)
             ZMatrix->SetEntry(np, npp, conj(PortVoltages[npp]));
          };

         /*--------------------------------------------------------------*/
         /*- write Z parameters to output file if that was requested    -*/
         /*--------------------------------------------------------------*/
         if (ZParameters)
          { fprintf(ZParFile,"%e ",Freq);
            for(np=0; np<NumPorts; np++)
             for(npp=0; npp<NumPorts; npp++)
              fprintf(ZParFile,"%e %e ",real(ZMatrix->GetEntry(np,npp)), 
                                        imag(ZMatrix->GetEntry(np,npp)));
            fprintf(ZParFile,"\n");
          };

         /*--------------------------------------------------------------*/
         /*- write S parameters to output file if that was requested    -*/
         /*--------------------------------------------------------------*/
         if (SParameters)
          { ZToS(ZMatrix, SMatrix);
            fprintf(SParFile,"%e ",Freq);
            for(np=0; np<NumPorts; np++)
             for(npp=0; npp<NumPorts; npp++)
              fprintf(SParFile,"%e %e ",real(SMatrix->GetEntry(np,npp)), imag(SMatrix->GetEntry(np,npp)));
            fprintf(SParFile,"\n");
          };
       }; // if (ZParameters || SParameters) 

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
         /* currents at this frequency                                   */
         /*--------------------------------------------------------------*/
         Log("  assembling RHS vector");
         for(np=0; np<NumPorts; np++)
          PortCurrents[np]=cdouble( PCList->GetEntryD(nf, 2*np+1), 
                                    PCList->GetEntryD(nf, 2*np+2));
         KN->Zero();
         AddPortContributionsToRHS(G, Ports, NumPorts, PortCurrents, Omega, KN);

         /*--------------------------------------------------------------*/
         /*- solve the system and evaluate the radiated fields   --------*/
         /*--------------------------------------------------------------*/
         Log("  solving the BEM system");
         M->LUSolve(KN);
         
         for(int n=0; n<nEPFiles; n++)
          ProcessEPFile(G, KN, Omega, Ports, NumPorts, PortCurrents, EPFiles[n], FileBase);

         for(int n=0; n<nFVMeshes; n++)
          ProcessFVMesh(G, KN, Omega, Ports, NumPorts, PortCurrents, FVMeshes[n], FileBase);

       };

   }; // for (nf=0, Freq=FreqList->GetEntryD(0); nf<FreqList->N; Freq=FreqList->GetEntryD(++nf) )

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (ZParameters)
   { fclose(ZParFile);
     printf("Z-parameters vs. frequency written to file %s\n",ZParFileName);
   }
  if (SParameters)
   { fclose(SParFile);
     printf("S-parameters vs. frequency written to file %s\n",SParFileName);
   }
  if (Moments)
   fclose(MomentFile);
  printf("Thank you for your support.\n");
}
