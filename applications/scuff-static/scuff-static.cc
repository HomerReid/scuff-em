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
 * scuff-static.cc   -- a standalone code within the scuff-EM suite for 
 *                   -- solving electrostatics problems 
 *
 * homer reid        -- 10/2006 -- 5/2013
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "scuff-static.h"

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define MAXEPF   10    // max number of evaluation-point files
#define MAXFVM   10    // max number of field-visualization meshes
#define MAXCACHE 10    // max number of cache files for preload

#define MAX_MONOPOLES 10
#define MAX_DIPOLES 10

/***************************************************************/
/* Attempt to read or write the solution vector for a BEM      */
/* electrostatics problem to/from an HDF5 binary data file.    */
/***************************************************************/
#define FILEOP_READ  0
#define FILEOP_WRITE 1
bool FileOp(int Op, HVector *Sigma, char *SolutionFile, char *SolutionName, char *TransformLabel)
{  
  Log("Attempting to %s solution: file %s ",
       Op==FILEOP_READ ? "read" : "write", SolutionFile);

  char DataSetName[1000];
  if(TransformLabel)
   snprintf(DataSetName,1000,"%s_%s",SolutionName,TransformLabel);
  else
   snprintf(DataSetName,1000,"%s",SolutionName);
  
  switch(Op)
   { 
     case FILEOP_READ:
      Sigma->ImportFromHDF5(SolutionFile,DataSetName);
      break;

     case FILEOP_WRITE:
      Sigma->ExportToHDF5(SolutionFile,DataSetName);
      break;
   };

  if (Sigma->ErrMsg)
   { Log("failed! %s",Sigma->ErrMsg);
     free(Sigma->ErrMsg);
     return false;
   }
  Log("success!");
  return true;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  InstallHRSignalHandler();
  InitializeLog(argv[0]);

  /***************************************************************/
  /* process options *********************************************/
  /***************************************************************/
  InstallHRSignalHandler();
// geometry-definition stuff
  char *GeoFile     = 0;
  char *SubstrateFile = 0;
  char *TransFile   = 0;
// excitation-definition stuff
  char *PotFile     = 0;
  char *ConstField  = 0;
  double Monopoles[4*MAX_MONOPOLES]; int nMonopoles;
  double Dipoles[6*MAX_DIPOLES];     int nDipoles;
  char *PhiExt      = 0;
  char *ExcitationFile = 0;
// output-definition stuff
  char *EPFiles[MAXEPF];             int nEPFiles;
  char *FileBase    = 0;
  char *PolFile     = 0;
  char *CapFile     = 0;
  char *PlotFile    = 0;
  char *CMatrixFile = 0;
  char *CMatrixHDF5File = 0;
  int lMax          = 2;             int nlMax;
  char *FVMeshes[MAXFVM];            int nFVMeshes;
  char *FVMeshTransFiles[MAXFVM];    int nFVMeshTransFiles;
  memset(FVMeshTransFiles, 0, MAXFVM*sizeof(char *));
  bool SeparateOutputFiles = false;
// other misc stuff
  char *SolutionFile= 0;
  char *SolutionName= 0;
  bool DisableEquivalentPairs = false;
  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { 
     {"Geometry",       PA_STRING,  1, 1,       (void *)&GeoFile,    0,             "geometry file\n"},
     {"Substrate",      PA_STRING,  1, 1,       (void *)&SubstrateFile,  0,             "substrate file\n"},
     {"TransFile",      PA_STRING,  1, 1,       (void *)&TransFile,  0,             "list of geometry transformations\n"},
/**/
     {"PotentialFile",  PA_STRING,  1, 1,       (void *)&PotFile,    0,             "list of conductor potentials\n"},
     {"ConstField",     PA_STRING,  1, 1,       (void *)&ConstField, 0,             "direction of constant unit-strength E field (x,y,z)\n"},
     {"Monopole",       PA_DOUBLE,  4, MAX_MONOPOLES, (void *)Monopoles, &nMonopoles, "(x0,y0,z0,Q) for point-source excitation\n"},
     {"Dipole",         PA_DOUBLE,  6, MAX_DIPOLES,  (void *)Dipoles, &nDipoles, "(x0,y0,z0,px,py,pz) for point-dipole excitation\n"},
     {"PhiExt",         PA_STRING,  1, 1,       (void *)&PhiExt,     0,  "user-specified external potential\n"},
     {"ExcitationFile", PA_STRING,  1, 1,       (void *)&ExcitationFile,     0,  "list of excitations\n"},
/**/
     {"EPFile",         PA_STRING,  1, MAXEPF,  (void *)EPFiles,     &nEPFiles,     "list of evaluation points"},
     {"FileBase",       PA_STRING,  1, 1,       (void *)&FileBase,   0,             "base name for field output files"},
     {"PolFile",        PA_STRING,  1, 1,       (void *)&PolFile, 0,             "polarizability output file\n"},
     {"CapFile",        PA_STRING,  1, 1,       (void *)&CapFile,     0,         "capacitance-matrix output file\n"},
     {"PlotFile",       PA_STRING,  1, 1,       (void *)&PlotFile,   0,         "surface-charge visualization output file\n"},
     {"CMatrixFile",    PA_STRING,  1, 1,       (void *)&CMatrixFile, 0,            "C-matrix text output file"},
     {"CMatrixHDF5File", PA_STRING, 1, 1,       (void *)&CMatrixHDF5File, 0,        "C-matrix HDF5 output file"},
     {"lMax",           PA_INT,     1, 1,       (void *)&lMax,       &nlMax,        "maximum l-value of spherical harmonic in C-matrix\n"},
     {"FVMesh",         PA_STRING,  1, MAXFVM,  (void *)FVMeshes,    &nFVMeshes,    "field visualization mesh"},
     {"FVMeshTransFile", PA_STRING,  1, MAXFVM,  (void *)FVMeshTransFiles,    &nFVMeshTransFiles,    "list of geometrical transformations for FVMesh\n"},
     {"SeparateOutputFiles", PA_BOOL, 1, 0,     (void *)&SeparateOutputFiles, 0,    "write separate output files for each transformation and excitation"},
/**/
     {"SolutionFile",   PA_STRING,  1, 1,       (void *)&SolutionFile, 0,           "name of HDF5 file for solution input/output"},
     {"SolutionName",   PA_STRING,  1, 1,       (void *)&SolutionName, 0,           "name of dataset within HDF5 file\n"},
     {"DisableEquivalentPairs", PA_BOOL, 1, 0,  (void *)&DisableEquivalentPairs, 0, "do not identify equivalent surface pairs\n"},
/**/
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  if (GeoFile==0)
   OSUsage(argv[0], VERSION, OSArray, "--geometry option is mandatory");

  if (nlMax && (CMatrixFile==0 || CMatrixHDF5File) )
   ErrExit("--lMax option can only be used with --CMatrixFile or --CMatrixHDF5File");

  /*******************************************************************/
  /* sanity checks on input arguments.                               */
  /* Type-1 calculations (polarizability, capacitance, C-matrix) do  */
  /*  not accept excitation specifications.                          */
  /* Type-2 calculations (fields, field visualization, surface-charge*/
  /*  visualization) require excitation specifications.              */
  /*******************************************************************/
  bool HaveType1Outputs = (PolFile || CapFile || CMatrixFile || CMatrixHDF5File);
  bool HaveType2Outputs = (nEPFiles>0 || nFVMeshes>0 || PlotFile!=0);
  bool HaveType2Inputs  = (    ExcitationFile!=0 || PotFile!=0  || ConstField!=0 
                            || nMonopoles!=0     || nDipoles!=0 || PhiExt!=0
                          );

  if ( (!HaveType1Outputs) && (!HaveType2Outputs) )
   OSUsage(argv[0], VERSION, OSArray, "you have not selected any type of calculation");
  if (HaveType1Outputs && HaveType2Outputs)
   ErrExit("{--EPFile,--FVMesh{ may not be used with {--polfilebase, --capfilebase}");
  if (HaveType1Outputs && HaveType2Inputs)
   ErrExit("potential/field specifications may not be used with {--polfilebase, --capfilebase}");
  if (HaveType2Outputs && !HaveType2Inputs)
   ErrExit("you have not specified any external potential or field");
  if (HaveType1Outputs && SolutionFile) 
   ErrExit("--SolutionFile may not be used with {--polfilebase, --capfilebase}");

  if (SolutionFile && SolutionName==0)
   SolutionName = strdup("DEFAULT");

  if (!FileBase) FileBase=strdup(GetFileBase(GeoFile));

  /*******************************************************************/
  /* create the ScuffStaticGeometry **********************************/
  /*******************************************************************/
  StaticSolver *SS   = new StaticSolver(GeoFile, SubstrateFile);
  SS->SeparateOutputFiles = SeparateOutputFiles;

  HMatrix *M      = 0;
  HVector *Sigma  = SS->AllocateRHSVector();
  RWGGeometry *G  = SS->G;

  /*******************************************************************/
  /* read in information about excitations ***************************/
  /*******************************************************************/
  StaticExcitation **SEList=0;
  int NumExcitations=0;
  if (ExcitationFile)
   { 
     SEList = ReadExcitationFile(SS, ExcitationFile, &NumExcitations);
   };
  if (PotFile || ConstField || nMonopoles>0 || nDipoles>0 || PhiExt )
   { 
     if (ExcitationFile)
      ErrExit("--ExcitationFile is incompatible with [--PotFile | --ConstField | --Monopole | --Dipole | --PhiExt]");
     SEList = CreateSimpleSEList(G, PotFile, ConstField,
                                 nMonopoles, Monopoles,
                                 nDipoles,   Dipoles, PhiExt);
     NumExcitations=1;
   };

  /****************************************************************/
  /*- read the transformation file if one was specified and check */
  /*- that it plays well with the specified geometry file.        */
  /*- note if TransFile==0 then this code snippet still works; in */
  /*- this case the list of GTComplices is initialized to contain */
  /*- a single empty GTComplex and the check automatically passes.*/
  /****************************************************************/
  GTCList GTCs=ReadTransFile(TransFile);
  int NT = GTCs.size();
  char *ErrMsg=G->CheckGTCList(GTCs);
  if (ErrMsg)
   ErrExit("file %s: %s",TransFile,ErrMsg);

  /*******************************************************************/
  /* if we have more than one geometrical transformation,            */
  /* allocate storage for BEM matrix blocks                          */
  /* 20170619 just do this by default now that we have equivalent-   */
  /*  pair detection                                                 */
  /*******************************************************************/
  //BMAccelerator *BMA = (NT==1) ? 0 : CreateBMAccelerator(SS);
  // TODO figure out interaction between transfiles and equivalent pairs
  if (NT>1) DisableEquivalentPairs=true; 
  BMAccelerator *BMA = CreateBMAccelerator(SS, GTCs, !DisableEquivalentPairs);

  /*******************************************************************/
  /* loop over transformations.                                      */
  /* if no --TransFile was specified, this loop iterates just once,  */
  /* for the 'default' (identity) transformation.                    */
  /*******************************************************************/
  for(int nt=0; nt<NT; nt++)
  { 
     if (TransFile)
      { G->Transform(GTCs[nt]);
        SS->TransformLabel=strdup(GTCs[nt]->Tag);
        Log("Working at transformation %s...",GTCs[nt]->Tag);
      }

     /*******************************************************************/
     /*******************************************************************/
     /*******************************************************************/
     bool HaveSolution=false;
     if (SolutionFile)
      HaveSolution=FileOp(FILEOP_READ, Sigma, SolutionFile, SolutionName, SS->TransformLabel);

     if (!HaveSolution)
      ReassembleBEMMatrix(SS, &M, BMA, nt);

     /*******************************************************************/
     /* now switch off depending on the type of calculation the user    */
     /* requested. The first three options below do not require         */
     /* user-specified external fields.                                 */
     /*******************************************************************/
     if (PolFile)
      WritePolarizabilities(SS, M, Sigma, PolFile);
     if (CapFile)
      WriteCapacitanceMatrix(SS, M, Sigma, CapFile);
     if (CMatrixFile)
      WriteCMatrix(SS, M, Sigma, lMax, CMatrixFile, CMatrixHDF5File);

     /******************************************************************/
     /* The remaining options do require user-specified external       */
     /* fields unless the user supplied a precomputed solution vector. */
     /******************************************************************/
     for(int n=0; n<NumExcitations; n++)
      { 
        StaticExcitation *SE = SEList[n];
        SS->ExcitationLabel = SE->Label;
        Log("Handling excitation %s",SE->Label ? SE->Label : "(default)");
        if (!HaveSolution)
         { 
           SS->AssembleRHSVector(SE, Sigma);
           M->LUSolve(Sigma);
           if (SolutionFile)
            FileOp(FILEOP_WRITE, Sigma, SolutionFile, SolutionName, SS->TransformLabel);
         };

        if (PlotFile)
         SS->PlotChargeDensity(Sigma, PlotFile);

        if (nEPFiles>0 || PlotFile )
         WriteFields(SS, Sigma, SE, EPFiles, nEPFiles, FileBase);

        for(int nfm=0; nfm<nFVMeshes; nfm++)
         SS->VisualizeFields(Sigma, FVMeshes[nfm], FileBase, SE,
                             FVMeshTransFiles[nfm]);
      };

     /*******************************************************************/
     /*******************************************************************/
     /*******************************************************************/
     G->UnTransform();

   }; // for(int nt=0; nt<NT; nt++)

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  printf("Thank you for your support.\n");

}
