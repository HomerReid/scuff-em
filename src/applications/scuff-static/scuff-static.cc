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

#include <libscuff.h>
#include <SSSolver.h>
#include <cmatheval.h>

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define MAXEPF   10    // max number of evaluation-point files
#define MAXCACHE 10    // max number of cache files for preload

#define MAXSTR   1000

/***************************************************************/
/* routines in OutputModules.cc ********************************/
/***************************************************************/
void WritePolarizabilities(SSSolver *SSS, HMatrix *M,
                           HVector *Sigma, char *FileName);

void WriteCapacitanceMatrix(SSSolver *SSS, HMatrix *M,
                            HVector *Sigma, char *CapFile);

void WriteCMatrix(SSSolver *SSS, HMatrix *M,
                  HVector *Sigma, int lMax,
                  char *TextFileName, char *HDF5FileName);

void WriteFields(SSSolver *SSS, HMatrix *M, HVector *Sigma,
                 char *PotFile, char *PhiExt, int ConstFieldDirection,
                 char *PlotFile, char **EPFiles, int nEPFiles,
                 char *FileBase);

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
  char *GeoFile     = 0;
  char *TransFile   = 0;
  char *PolFile     = 0;
  char *CapFile     = 0;
  char *CMatrixFile = 0;
  char *CMatrixHDF5File = 0;
  int lMax          = 2;             int nlMax;
  char *PotFile     = 0;
  char *PhiExt      = 0;
  char *EPFiles[MAXEPF];             int nEPFiles;
  char *FileBase    = 0;
  char *PlotFile    = 0;
  char *Cache       = 0;
  char *ReadCache[MAXCACHE];         int nReadCache;
  char *WriteCache  = 0;
  char *ConstField  = 0;
  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { 
     {"Geometry",       PA_STRING,  1, 1,       (void *)&GeoFile,    0,             "geometry file"},
/**/
     {"TransFile",      PA_STRING,  1, 1,       (void *)&TransFile,  0,             "list of geometry transformations"},
/**/
     {"PolFile",        PA_STRING,  1, 1,       (void *)&PolFile,    0,             "polarizability output file"},
/**/
     {"CapFile",        PA_STRING,  1, 1,       (void *)&CapFile,    0,             "capacitance matrix output file"},
/**/
     {"CMatrixFile",    PA_STRING,  1, 1,       (void *)&CMatrixFile, 0,            "C-matrix text output file"},
     {"CMatrixHDF5File", PA_STRING, 1, 1,       (void *)&CMatrixHDF5File, 0,        "C-matrix HDF5 output file"},
     {"lMax",           PA_INT,     1, 1,       (void *)&lMax,       &nlMax,        "maximum l-value of spherical harmonic in C-matrix "},
/**/
     {"PotentialFile",  PA_STRING,  1, 1,       (void *)&PotFile,    0,             "list of conductor potentials"},
/**/
     {"PhiExternal",    PA_STRING,  1, 1,       (void *)&PhiExt,     0,             "functional form of external potential"},
     {"ConstField",     PA_STRING,  1, 1,       (void *)&ConstField, 0,             "direction of constant unit-strength E field (x,y,z) "},
/**/
     {"EPFile",         PA_STRING,  1, MAXEPF,  (void *)EPFiles,     &nEPFiles,     "list of evaluation points"},
     {"FileBase",       PA_STRING,  1, 1,       (void *)&FileBase,   0,             "base filename for EP file output"},
     {"PlotFile",       PA_STRING,  1, 1,       (void *)&PlotFile,   0,             "surface charge visualization output file"},
/**/
     {"Cache",          PA_STRING,  1, 1,       (void *)&Cache,      0,             "read/write cache"},
     {"ReadCache",      PA_STRING,  1, MAXCACHE,(void *)ReadCache,   &nReadCache,   "read cache"},
     {"WriteCache",     PA_STRING,  1, 1,       (void *)&WriteCache, 0,             "write cache"},
/**/
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  if (GeoFile==0)
   OSUsage(argv[0], OSArray, "--geometry option is mandatory");

  if (FileBase==0)
   FileBase=vstrdup(GetFileBase(GeoFile));

  if (nlMax && (CMatrixFile==0 || CMatrixHDF5File) )
   ErrExit("--lMax option can only be used with --CMatrixFile or --CMatrixHDF5File");

  /*******************************************************************/
  /* sanity check on input arguments *********************************/
  /*******************************************************************/
  if ( (PolFile || CapFile || CMatrixFile || CMatrixHDF5File) && (nEPFiles>0 || PotFile!=0 || PhiExt!=0 ) )
   ErrExit("(--EPFile,--PotFile,--PhiExternal) may not be used with (--polfile, --capfile)");
  if ( nEPFiles==0 && !PolFile && !CapFile && !CMatrixFile && !CMatrixHDF5File && !PlotFile )
   OSUsage(argv[0], OSArray, "you have not selected any type of calculation");

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  int ConstFieldDirection=-1;
  if (ConstField)
   { switch(tolower(ConstField[0]))
      { case 'x': ConstFieldDirection=0; break;
        case 'y': ConstFieldDirection=1; break;
        case 'z': ConstFieldDirection=2; break;
        default:  ErrExit("invalid --ConstField specification");
      };
   };

  /*******************************************************************/
  /* create the ScuffStaticGeometry **********************************/
  /*******************************************************************/
  SSSolver *SSS   = new SSSolver(GeoFile);
  HMatrix *M      = SSS->AllocateBEMMatrix();
  HVector *Sigma  = SSS->AllocateRHSVector();

  RWGGeometry *G  = SSS->G;

  /****************************************************************/
  /*- read the transformation file if one was specified and check */
  /*- that it plays well with the specified geometry file.        */
  /*- note if TransFile==0 then this code snippet still works; in */
  /*- this case the list of GTComplices is initialized to contain */
  /*- a single empty GTComplex and the check automatically passes.*/
  /****************************************************************/
  int NT;
  GTComplex **GTCList=ReadTransFile(TransFile, &NT);
  char *ErrMsg=G->CheckGTCList(GTCList, NT);
  if (ErrMsg)
   ErrExit("file %s: %s",TransFile,ErrMsg);

  /*******************************************************************/
  /* if we have more than one geometrical transformation,            */
  /* allocate storage for BEM matrix blocks                          */
  /*******************************************************************/
  HMatrix **TBlocks=0, **UBlocks=0;
  int NS=G->NumSurfaces;
  if (NT>1)
   { int NADB = NS*(NS-1)/2; // number of above-diagonal blocks
     TBlocks  = (HMatrix **)mallocEC(NS*sizeof(HMatrix *));
     UBlocks  = (HMatrix **)mallocEC(NADB*sizeof(HMatrix *));
     for(int ns=0, nb=0; ns<NS; ns++)
      { 
        int NBF=G->Surfaces[ns]->NumPanels;

        // allocate diagonal block for surface #ns, unless
        // it has a mate in which case we reuse the mate's block
        int nsMate = G->Mate[ns];
        if ( nsMate!=-1 )
         TBlocks[ns] = TBlocks[nsMate];
        else
         TBlocks[ns] = new HMatrix(NBF, NBF, M->RealComplex);

        // allocate off-diagonal blocks for surfaces ns,nsp>ns
        for(int nsp=ns+1; nsp<NS; nsp++, nb++)
         { int NBFp = G->Surfaces[nsp]->NumPanels;
           UBlocks[nb] = new HMatrix(NBF, NBFp, M->RealComplex);
         };
      };
   };

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

  /*******************************************************************/
  /* loop over transformations ***************************************/
  /* (if no --TransFile was specified, this loop iterates just once, */
  /* for the 'default' (identity) transformation                     */
  /*******************************************************************/
  for(int nt=0; nt<NT; nt++)
   { 
     char TransformStr[100]="";
     if (TransFile)
      { 
        G->Transform(GTCList[nt]);
        SSS->TransformLabel=GTCList[nt]->Tag;
        Log("Working at transformation %s...",SSS->TransformLabel);
        snprintf(TransformStr,100,"_%s",SSS->TransformLabel);
      };

     /*******************************************************************/
     /* assemble BEM matrix, either (a) all at once if we have only one */
     /* geometric transformation, or (b) with the diagonal and off-     */
     /* diagonal blocks computed separately so that the former can be   */
     /* reused for multiple geometric transformations                   */
     /*******************************************************************/
     if (NT==1)
      SSS->AssembleBEMMatrix(M);
     else
      { 
        if (nt==0)
         for(int ns=0, nb=0; ns<NS; ns++)
          SSS->AssembleBEMMatrixBlock(ns, ns, TBlocks[nb]);

        for(int ns=0, nb=0; ns<NS; ns++)
         for(int nsp=ns+1; nsp<NS; nsp++, nb++)
          if (nt==0 || G->SurfaceMoved[ns] || G->SurfaceMoved[nsp] )
           SSS->AssembleBEMMatrixBlock(ns, nsp, UBlocks[nb]);

        /*******************************************************************/
        /* stamp blocks into BEM matrix and LU-factorize *******************/
        /*******************************************************************/
        for(int ns=0, nb=0; ns<NS; ns++)
         { 
           int RowOffset=G->PanelIndexOffset[ns];
           M->InsertBlock(TBlocks[ns], RowOffset, RowOffset);
           for(int nsp=ns+1; nsp<NS; nsp++, nb++)
            { int ColOffset=G->PanelIndexOffset[nsp];
              M->InsertBlock(UBlocks[nb], RowOffset, ColOffset);
              M->InsertBlockTranspose(UBlocks[nb], ColOffset, RowOffset);
            };
         };
      };
     M->LUFactorize();

     /*******************************************************************/
     /* now switch off depending on the type of calculation the user    */
     /* requested                                                       */
     /*******************************************************************/
     if (PolFile)
      WritePolarizabilities(SSS, M, Sigma, PolFile);
     if (CapFile)
      WriteCapacitanceMatrix(SSS, M, Sigma, CapFile);
     if (CMatrixFile)
      WriteCMatrix(SSS, M, Sigma, lMax, CMatrixFile, CMatrixHDF5File);
     if (nEPFiles>0 || PlotFile )
      WriteFields(SSS, M, Sigma,
                  PotFile, PhiExt, ConstFieldDirection,
                  PlotFile, EPFiles, nEPFiles, FileBase);

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
