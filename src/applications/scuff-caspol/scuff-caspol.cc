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
 *                   -- computing dyadic green's functions of material
 *                   -- bodies, as well as casimir-polder and van-der-Waals 
 *                   -- potentials for polarizable molecules in the 
 *                   -- vicinity of those bodies
 *
 * documentation at: 
 *  http://homerreid.com/scuff-em/scuff-caspol
 *
 * homer reid        -- 10/2006 -- 2/2012
 *
 */
#include "scuff-caspol.h"

#include <time.h>

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define MAXSTR  1000
#define MAXXI   10

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define FILETYPE_BYXI 0
#define FILETYPE_OUT  1
void WriteFilePreamble(FILE *f, int argc, char *argv[], int FileType)
{ 

  time_t MyTime=time(0);
  struct tm *MyTm=localtime(&MyTime);
  char TimeString[200];
  strftime(TimeString,30,"%D::%T",MyTm);

  fprintf(f,"# scuff-caspol running on %s (%s)\n",getenv("HOST"),TimeString);
  fprintf(f,"#");
  fprintf(f,"# command line:\n");
  fprintf(f,"# scuff-caspol ");
  for(int narg=0; narg<(argc-1); narg++)
   fprintf(f," %s%s",argv[narg+1],(narg%4)==3 ? "\n# " : " ");
  fprintf(f,"\n");
  fprintf(f,"#\n");
  fprintf(f,"# data columns:\n");
  fprintf(f,"#1: X \n");
  fprintf(f,"#2: Y \n");
  fprintf(f,"#3: Z \n");
  if (FileType==FILETYPE_BYXI)
   { fprintf(f,"#4: Xi \n");
     fprintf(f,"#5: casimir-polder potential per unit frequency at (X, Y, Z, Xi)\n");
   }
  else
   fprintf(f,"#4: casimir-polder potential at (X, Y, Z)\n");

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void Usage(char *ProgramName, OptStruct *OSArray, char *ErrMsg)
{
  fprintf(stderr, "error: %s (aborting)\n",ErrMsg);
  fprintf(stderr, "\n");
  fprintf(stderr, "usage: %s [options] \n");
  fprintf(stderr, "\n");
  fprintf(stderr, "options: \n");
  fprintf(stderr, " --geoFile  MyFile.scuffgeo         specify geometry file\n");
  fprintf(stderr, " --PECPlate                         use PEC plate\n");
  fprintf(stderr, " --EPFile   MyEvalPointFile         file specifying evaluation points\n");
  fprintf(stderr, " --atom     Rubidium                type of atom\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "atoms supported: \n");
  fprintf(stderr, " --atom Hydrogen   [--atom H]  \n");
  fprintf(stderr, " --atom Lithium    [--atom Li] \n");
  fprintf(stderr, " --atom Sodium     [--atom Na] \n");
  fprintf(stderr, " --atom Potassium  [--atom K]  \n");
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
  int PECPlate=0;
  char *Atom=0;
  char *EPFile=0;
  OptStruct OSArray[]=
   { {"geometry",       PA_STRING,  1, 1,       (void *)&GeoFile,  0,          ".scuffgeo file"},
     {"PECPlate",       PA_BOOL,    0, 1,       (void *)&PECPlate, 0,          "use PEC plate"},
     {"Atom",           PA_STRING,  1, 1,       (void *)&Atom,     0,          "type of atom"},
     {"EPFile",         PA_STRING,  1, 1,       (void *)&EPFile,   0,          "list of evaluation points"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  if (GeoFile==0 && PECPlate==0)
   Usage(argv[0], OSArray,"either --geometry or --PECPlate must be specified");
  if (GeoFile!=0 && PECPlate!=0)
   ErrExit("geometry and --PECPlate are mutually exclusive");
  if (Atom==0)
   Usage(argv[0], OSArray,"--Atom option is mandatory");
  if (EPFile==0)
   Usage(argv[0], OSArray,"--EPFile option is mandatory");

  /*******************************************************************/
  /* create the RWGGeometry, allocate BEM matrix and RHS vector, and */
  /* set up the data structure passed to the computational routines  */
  /*******************************************************************/
  SCPData MySCPData, *SCPD=&MySCPData;

  if (GeoFile)
   { SCPD->G  = new RWGGeometry(GeoFile); 
     SCPD->M  = SCPD->G->AllocateBEMMatrix(SCUFF_PUREIMAGFREQ);
     SCPD->KN = SCPD->G->AllocateRHSVector(SCUFF_PUREIMAGFREQ);
   }
  else
   { SCPD->G  = 0; // in this case we take the 
     SCPD->M  = 0; // geometry to be a PEC plate in the xy plane
     SCPD->KN = 0;
   };

  SetLogFileName("scuff-caspol.log");

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  SCPD->PM = new PolModel(Atom);
  if (SCPD->PM->ErrMsg)
   Usage(argv[0], SCPD->PM->ErrMsg);

  /*******************************************************************/
  /* process list of evaluation points *******************************/
  /*******************************************************************/
  HMatrix *EPMatrix = SCPD->EPMatrix = new HMatrix(EPFile,LHM_TEXT,"--ncol 3");
  if (EPMatrix->ErrMsg)
   ErrExit(EPMatrix->ErrMsg);

  int nEvalPoints=EPMatrix->NR;
  // storage for values of CP potential at eval points
  double *U=(double *)mallocEC(nEvalPoints * sizeof(double)); 

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  char GeoFileBase[MAXSTR], ByXiFileName[MAXSTR]; 
  if (GeoFile)
   strncpy(GeoFileBase, GetFileBase(GeoFile), MAXSTR);
  else
   sprintf(GeoFileBase,"PECPlate");
  snprintf(ByXiFileName, MAXSTR, "%s.%s.byXi",GeoFileBase,Atom);
  SCPD->ByXiFile=CreateUniqueFile(ByXiFileName,1,ByXiFileName);
  WriteFilePreamble(SCPD->ByXiFile, argc, argv, FILETYPE_BYXI);

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
     for(nXi=0; nXi<XiList->N; nXi++)
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
  fclose(SCPD->ByXiFile);

  /***************************************************************/
  /* write frequency-integrated or matsubara-summed results to   */
  /* output file if that was requested                           */
  /***************************************************************/
  if (!XiList)
   {
     char OutFileName[MAXSTR];
     snprintf(OutFileName, MAXSTR, "%s.out",GeoFileBase);
     FILE *OutFile=CreateUniqueFile(OutFileName,1,OutFileName);
     WriteFilePreamble(OutFile, argc, argv, FILETYPE_OUT);
     int nep;
     for(nep=0; nep<EPList->NR; nep++)
      fprintf(OutFile,"%e %e %e %e \n", EPList->GetEntryD(nep,0),
                                        EPList->GetEntryD(nep,1),
                                        EPList->GetEntryD(nep,2),
                                        U[nep]);
     fclose(OutFile);
     if (Temperature>0.0)
      printf("Matsubara-summed data written to file %s.\n",OutFileName);
     else
      printf("Frequency-integrated data written to file %s.\n",OutFileName);
   };

  /***************************************************************/
  /***************************************************************/  
  /***************************************************************/
  printf("Frequency-resolved data written to file %s.\n",ByXiFileName);
  printf("Thank you for your support.\n");

}
