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
 * scuff-caspol.cc   -- a standalone code within the scuff-EM suite for
 *                   -- computing casimir-polder potentials for
 *                   -- polarizable molecules in the vicinity of material
 *                   -- bodies
 *
 * documentation at: 
 *  http://homerreid.com/scuff-em/scuff-caspol
 *
 * homer reid        -- 10/2006 -- 8/2013
 *
 */
#include "scuff-caspol.h"

#include <time.h>

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define MAXSTR  1000
#define MAXFREQ 10
#define MAXATOM 10

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define FILETYPE_BYXI 0
#define FILETYPE_OUT  1
void WriteFilePreamble(const char *FileName, int argc, char *argv[], 
                       int FileType, SCPData *SCPD)
{ 
  FILE *f = fopen(FileName, "a");
  if (!f) return;

  time_t MyTime=time(0);
  struct tm *MyTm=localtime(&MyTime);
  char TimeString[200];
  strftime(TimeString,30,"%D::%T",MyTm);

  fprintf(f,"# scuff-caspol running on %s (%s)\n",getenv("HOST"),TimeString);
  fprintf(f,"#\n");
  fprintf(f,"# command line:\n");
  fprintf(f,"#\n");
  fprintf(f,"# scuff-caspol ");
  for(int narg=0; narg<(argc-1); narg++)
   fprintf(f," %s%s",argv[narg+1],(narg%4)==3 ? "\n# " : " ");
  fprintf(f,"\n");
  fprintf(f,"#\n");
  fprintf(f,"# data columns:\n");
  fprintf(f,"#1: X (um) \n");
  fprintf(f,"#2: Y (um) \n");
  fprintf(f,"#3: Z (um) \n");
  
  int nc=4;
  if (FileType==FILETYPE_BYXI)
   { 
     fprintf(f,"#%i: Xi (imaginary frequency) (3e14 rad/sec)\n",nc++);

     for(int na=0; na<SCPD->NumAtoms; na++)
      { fprintf(f,"#%i: polarizability of %s (a_0^3)\n",
                   nc++,SCPD->PolModels[na]->Name);
        fprintf(f,"#%i: CP potential per unit frequency for %s (neV/3e14 rad/sec)\n",
                   nc++,SCPD->PolModels[na]->Name);
      };

   }
  else
   {
     for(int na=0; na<SCPD->NumAtoms; na++)
      fprintf(f,"#%i: Casimir-polder potential for %s (neV)\n",
                nc++,SCPD->PolModels[na]->Name);
   };

  fclose(f);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void Usage(const char *ProgramName, OptStruct *OSArray, const char *ErrMsg)
{
  (void) OSArray; // unused

  fprintf(stderr, "error: %s (aborting)\n",ErrMsg);
  fprintf(stderr, "\n");
  fprintf(stderr, "usage: %s [options] \n",ProgramName);
  fprintf(stderr, "\n");
  fprintf(stderr, "options: \n");
  fprintf(stderr, "\n");
  fprintf(stderr, " --geoFile     MyFile.scuffgeo   specify geometry file\n");
  fprintf(stderr, " --PECPlate                      use PEC plate\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " --EPFile      MyEvalPointFile   file specifying evaluation points\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " --Xi          2.17              imaginary angular frequency\n");
  fprintf(stderr, " --XiFile      MyXiFile          list of imaginary angular frequencies\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " --atom        Rubidium          atomic species (for built-in polarizability models)\n");
  fprintf(stderr, " --particle    MyParticle.pol    particle species (for a user-defined polarizability model)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " --temperature xx                compute at T=xx kelvin\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " --reltol      xx                relative error tolerance for frequency sums/integrals\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "built-in polarizability models: \n");
  fprintf(stderr, "\n");
  fprintf(stderr, " --atom Hydrogen   [--atom H ] \n");
  fprintf(stderr, " --atom Lithium    [--atom Li] \n");
  fprintf(stderr, " --atom Sodium     [--atom Na] \n");
  fprintf(stderr, " --atom Potassium  [--atom K ] \n");
  fprintf(stderr, " --atom Rubidium   [--atom Rb] \n");
  fprintf(stderr, " --atom Cesium     [--atom Cs] \n");
  fprintf(stderr, " --atom Francium   [--atom Fr] \n");
  fprintf(stderr, "\n");
  exit(1);

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
  bool PECPlate=false;
//
  char *Atoms[MAXATOM];                     int NumBIAtoms;
  char *Particles[MAXATOM];                 int NumParticles;
//
  char *EPFile=0;
//
  double Temperature = 0.0;                 int nTemperature;
  double XiVals[MAXFREQ];	            int nXiVals;
  char *XiFile=0;
//
  double RelTol = DEF_RELTOL;
  /* name           type  #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { 
     {"geometry",    PA_STRING,  1, 1,       (void *)&GeoFile,     0,             ".scuffgeo file"},
     {"PECPlate",    PA_BOOL,    0, 1,       (void *)&PECPlate,    0,             "use PEC plate"},
//
     {"Atom",        PA_STRING,  1, MAXATOM, (void *)Atoms,        &NumBIAtoms,   "atomic species (built-in polarizability models)"},
     {"Particle",    PA_STRING,  1, MAXATOM, (void *)Particles,    &NumParticles, "particle species (user-defined polarizability models)"},
//
     {"EPFile",      PA_STRING,  1, 1,       (void *)&EPFile,      0,             "list of evaluation points"},
//
     {"Temperature", PA_DOUBLE,  1, 1,       (void *)&Temperature, &nTemperature, "temperature in Kelvin"},
     {"Xi",          PA_DOUBLE,  1, MAXFREQ, (void *)XiVals,       &nXiVals,      "imaginary frequency"},
     {"XiFile",      PA_STRING,  1, 1,       (void *)&XiFile,      0,             "file containing Xi values"},
//
     {"RelTol",      PA_DOUBLE,  1, 1,       (void *)&RelTol,      0,             "relative error tolerance"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  /***************************************************************/
  /* some sanity checks on input arguments ***********************/
  /***************************************************************/
  if (GeoFile==0 && PECPlate==0)
   Usage(argv[0], OSArray,"either --geometry or --PECPlate must be specified");
  if (GeoFile!=0 && PECPlate!=0)
   ErrExit("geometry and --PECPlate are mutually exclusive");
  if (EPFile==0)
   Usage(argv[0], OSArray,"--EPFile option is mandatory");
  if ( XiFile && nXiVals )
   ErrExit("--XiFile and --Xi are mutually exclusive");
  if (nTemperature!=0 && (XiFile!=0 || nXiVals!=0) )
   ErrExit("--Xi and --XiFile options are incompatible with --Temperature");
  if (Temperature<0.0)
   ErrExit("--temperature may not be negative");
  int NumAtoms = NumBIAtoms + NumParticles;
  if (NumAtoms == 0)
   Usage(argv[0], OSArray,"you must specify at least one --atom or --particle");

  /*******************************************************************/
  /* create the RWGGeometry, allocate BEM matrix and RHS vector, and */
  /* set up the data structure passed to the computational routines  */
  /*******************************************************************/
  SCPData MySCPData, *SCPD=&MySCPData;
  if (GeoFile)
   { SCPD->G  = new RWGGeometry(GeoFile, SCUFF_TERSELOGGING);
     SCPD->M  = SCPD->G->AllocateBEMMatrix(SCUFF_PUREIMAGFREQ);
     SCPD->KN = SCPD->G->AllocateRHSVector(SCUFF_PUREIMAGFREQ);
   }
  else
   { SCPD->G  = 0; // in this case we take the 
     SCPD->M  = 0; // geometry to be a PEC plate in the xy plane
     SCPD->KN = 0;
     GeoFile = strdup("PECPlate");
   };
  SCPD->RelTol = RelTol;

  /*******************************************************************/
  /* create polarizability models for all specified atoms            */
  /*******************************************************************/
  SCPD->PolModels  = (PolModel **)malloc(NumAtoms * sizeof(PolModel *));
  SCPD->Alphas     = (HMatrix **) malloc(NumAtoms * sizeof(HMatrix *) );
  SCPD->NumAtoms   = NumAtoms;

  int na=0;
  for(int nbi=0; nbi<NumBIAtoms; nbi++)
   { SCPD->PolModels[na] = new PolModel(Atoms[nbi], PM_BUILTIN);
     if (SCPD->PolModels[na]->ErrMsg)
      Usage(argv[0], OSArray, SCPD->PolModels[na]->ErrMsg);
     SCPD->Alphas[na] = new HMatrix(3,3);
     na++;
   };
  for(int np=0; np<NumParticles; np++)
   { SCPD->PolModels[na] = new PolModel(Particles[np], PM_FILE);
     if (SCPD->PolModels[na]->ErrMsg)
      Usage(argv[0], OSArray, SCPD->PolModels[na]->ErrMsg);
     SCPD->Alphas[na] = new HMatrix(3,3);
     na++;
   };

  /*******************************************************************/
  /* process list of evaluation points *******************************/
  /*******************************************************************/
  HMatrix *EPMatrix = SCPD->EPMatrix = new HMatrix(EPFile,LHM_TEXT,"--ncol 3");
  if (EPMatrix->ErrMsg)
   ErrExit(EPMatrix->ErrMsg);

  int NumEvalPoints=EPMatrix->NR;
  // storage for values of CP potential at eval points
  double *U=(double *)mallocEC(NumEvalPoints * NumAtoms * sizeof(double));

  /*******************************************************************/
  /* process --Xi or --XiFile options if present                     */
  /*******************************************************************/
  HVector *XiList=0;
  if ( XiFile )
   { XiList = new HVector(XiFile,LHM_TEXT,"--nc 1 --strict");
     if (XiList->ErrMsg)
      ErrExit(XiList->ErrMsg);
     Log("Read %i Xi points from file %s.",XiList->N, XiFile);
   }
  else if ( nXiVals>0 )
   { XiList = new HVector(nXiVals);
     for(int n=0; n<nXiVals; n++)
      XiList->SetEntry(n,XiVals[n]);
   };

  /*******************************************************************/
  /* create the .byXi file *******************************************/
  /*******************************************************************/
  SCPD->ByXiFileName=vstrdup("%s.byXi",GetFileBase(GeoFile));
  WriteFilePreamble(SCPD->ByXiFileName, argc, argv, FILETYPE_BYXI, SCPD);

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  SetLogFileName("scuff-caspol.log");
  Log("scuff-caspol running on %s",GetHostName());
  if ( XiList )
   Log("Computing Casimir-Polder quantities at %i user-specified frequeny points.",XiList->N);
  else if ( Temperature!=0.0 )
   Log("Computing full Matsubara-summed Casimir quantities at T=%e Kelvin.",Temperature);
  else
   Log("Computing full zero-temperature Casimir quantities.");

  /*******************************************************************/
  /* main program body:                                              */
  /*                                                                 */
  /*  (a) if the user specified a set of frequencies, then we just   */
  /*      calculate the contributions to the casimir-polder potential*/
  /*      from those frequencies.                                    */
  /*                                                                 */
  /*  (b) if the user specified a temperature, then we evaluate the  */
  /*      matsubara sum for the full casimir-polder potential at     */
  /*      that temperature.                                          */
  /*                                                                 */
  /*  (c) if the user specified neither of the above, then we        */
  /*      evaluate the integral over the imaginary-frequency axis    */
  /*      to get the full zero-temperature casimir-polder potential. */
  /*                                                                 */
  /*******************************************************************/
  if (XiList)
   {
     for(int nXi=0; nXi<XiList->N; nXi++)
      GetCPIntegrand(SCPD, XiList->GetEntryD(nXi), U);
   }
  else if ( Temperature != 0.0 )
   { 
     EvaluateMatsubaraSum(SCPD, Temperature, U);
   }
  else
   { 
     EvaluateFrequencyIntegral(SCPD, U);
   };

  /***************************************************************/
  /* write frequency-integrated or Matsubara-summed results to   */
  /* output file                                                 */
  /***************************************************************/
  if (XiList==0)
   { 
     char OutFileName[MAXSTR];
     snprintf(OutFileName, MAXSTR, "%s.out",GetFileBase(GeoFile));
     WriteFilePreamble(OutFileName, argc, argv, FILETYPE_OUT, SCPD);

     FILE *f=fopen(OutFileName,"a");
     for(int nep=0; nep<EPMatrix->NR; nep++)
      { fprintf(f,"%e %e %e ", EPMatrix->GetEntryD(nep,0),
                               EPMatrix->GetEntryD(nep,1),
                               EPMatrix->GetEntryD(nep,2) );
        for(na=0; na<NumAtoms; na++)
         fprintf(f,"%e ", U[nep*NumAtoms + na]);

        fprintf(f,"\n");
      };
     fclose(f);

     if (Temperature>0.0)
      printf("Matsubara-summed data written to file %s.\n",OutFileName);
     else
      printf("Frequency-integrated data written to file %s.\n",OutFileName);
   };

  /***************************************************************/
  /***************************************************************/  
  /***************************************************************/
  printf("Frequency-resolved data written to file %s.\n",SCPD->ByXiFileName);
  printf("Thank you for your support.\n");

}
