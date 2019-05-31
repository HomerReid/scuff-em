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
 * scuff-tmatrix.cc 
 *
 * This program computes the T-matrix, in a spherical-wave basis,  
 * for an arbitary scatterer. (What this means is the following: 
 * For each spherical wave up to a user-specified maximum L value, 
 * the code scatters that wave off of the geometry and decomposes
 * the resulting scattered field into a basis of spherical waves.)
 * The implementation of this code demonstrates, among other things, 
 * how to define your own type of incident field routine.
 *
 * Homer Reid    -- 6/2012
 *
 * Update, 11/2017: 
 *
 *   -- fixed incorrect normalization conventions
 *   -- updated spherical-moment calculation to use modern 'Panel Cubature' framework
 *   -- updated documentation at
 *      http://homerreid.github.io/scuff-em-documentation/applications/scuff-tmatrix/scuff-tmatrix
 */

#ifdef VERSION
#define VERSIONSTRING VERSION " "
#else
#define VERSIONSTRING ""
#endif

#include <stdio.h>
#include <stdlib.h>

#include "libscuff.h"
#include "libhmat.h"
#include "libSpherical.h"

using namespace scuff;
#define II cdouble (0.0,1.0)

/***************************************************************/
/***************************************************************/
/***************************************************************/
HVector *GetSphericalMoments(RWGGeometry *G, cdouble k, int LMax,
                             HVector *KN, HVector *MomentVector);

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  InstallHRSignalHandler();
  InitializeLog(argv[0]);

  /*--------------------------------------------------------------*/
  /*- process command-line options -------------------------------*/
  /*--------------------------------------------------------------*/
  char *GeoFileName=0;       // .scuffgeo file
  int LMax=3;                // maximum L-value of spherical wave computed
  cdouble Omega=0;           // angular frequency at which to run the computation
  char *OmegaFile=0;         // list of angular frequencies
  char *Cache=0;             // scuff cache file 
  char *FileBase=0;          // base filename for output file
  bool WriteHDF5Files=false; // write T-matrix data to HDF5 files
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"geometry",       PA_STRING,  1, 1, (void *)&GeoFileName,    0,  ".scuffgeo file"},
     {"LMax",           PA_INT,     1, 1, (void *)&LMax,           0,  "maximum l-value"},
     {"Omega",          PA_CDOUBLE, 1, 1, (void *)&Omega,          0,  "angular frequency"},
     {"OmegaFile",      PA_STRING,  1, 1, (void *)&OmegaFile,      0,  "list of angular frequencies"},
     {"Cache",          PA_STRING,  1, 1, (void *)&Cache,          0,  "scuff cache file"},
     {"FileBase",       PA_STRING,  1, 1, (void *)&FileBase,       0,  "base filename for output files"},
     {"WriteHDF5Files", PA_BOOL,    0, 1, (void *)&WriteHDF5Files, 0,  "write HDF5 output files"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFileName==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");

  /*--------------------------------------------------------------*/
  /*- process frequency options ----------------------------------*/
  /*--------------------------------------------------------------*/
  HVector *OmegaVector=0;
  if (OmegaFile)
   OmegaVector=new HVector(OmegaFile);
  else if ( Omega!=0.0 )
   { OmegaVector=new HVector(1, LHM_COMPLEX);
     OmegaVector->SetEntry(0,Omega);
   }
  else
   OSUsage(argv[0],OSArray,"either --omega or --omegafile must be specified");

  /*--------------------------------------------------------------*/
  /* create the RWGGeometry from the .scuffgeo file               */
  /*--------------------------------------------------------------*/
  RWGGeometry *G=new RWGGeometry(GeoFileName);
  G->SetLogLevel(SCUFF_VERBOSELOGGING);
  if (Cache)
   PreloadCache(Cache);

  /*--------------------------------------------------------------*/
  /* preallocate BEM matrix and RHS vector                        */
  /*--------------------------------------------------------------*/
  HMatrix *M  = G->AllocateBEMMatrix();
  HVector *KN = G->AllocateRHSVector();

  /*--------------------------------------------------------------*/
  /*- preallocate an HMatrix to store the T-matrix data           */
  /*--------------------------------------------------------------*/
  int NumLMs = (LMax+1)*(LMax+1) - 1;
  int NumMoments= 2*NumLMs;
  HMatrix *TMatrix = new HMatrix(NumMoments, NumMoments, LHM_COMPLEX);

  /*--------------------------------------------------------------*/
  /* instantiate a SphericalWave structure (we will set the L, M, */
  /* and P fields later)                                          */
  /*--------------------------------------------------------------*/
  SphericalWave SW;

  /*--------------------------------------------------------------*/
  /* open output file and write preamble    ----------------------*/
  /*--------------------------------------------------------------*/
  if (!FileBase)
   FileBase = vstrdup(GetFileBase(GeoFileName));
  FILE *f=vfopen("%s.TMatrix","a",FileBase);
  fprintf(f,"# scuff-tmatrix " VERSIONSTRING "run on %s (%s)\n",GetHostName(),GetTimeString());
  fprintf(f,"# columns:\n");
  fprintf(f,"# 1 omega\n");
  fprintf(f,"# 2,3,4,5 (alpha, {L,M,P}_alpha)   (T-matrix row index)\n");
  fprintf(f,"# 6,7,8,9 ( beta, {L,M,P}_beta)    (T-matrix column index)\n");
  fprintf(f,"# 10, 11  real, imag T_{alpha, beta}\n");

  /*--------------------------------------------------------------*/
  /*- outer loop over frequencies. -------------------------------*/
  /*--------------------------------------------------------------*/
  for(int nOmega=0; nOmega<OmegaVector->N; nOmega++)
   {
     Log("Computing T matrix at frequency %s...",z2s(Omega));

     /*--------------------------------------------------------------*/
     /* assemble and factorize the BEM matrix at this frequency      */
     /*--------------------------------------------------------------*/
     Omega=OmegaVector->GetEntry(nOmega);
     G->AssembleBEMMatrix(Omega, M);
     M->LUFactorize();

     /*--------------------------------------------------------------*/
     /*- inner loop over incident spherical waves (i.e. over columns */
     /*- of the T matrix; note Betais a running column index)        */
     /*--------------------------------------------------------------*/
     TMatrix->Zero();
     for(int LBeta=1, Beta=0; LBeta<=LMax; LBeta++)
      for(int MBeta=-LBeta; MBeta<=LBeta; MBeta++)
       for(int PBeta=0; PBeta<2; PBeta++, Beta++)
        { 
           SW.SetL(LBeta);
           SW.SetM(MBeta);
           SW.SetP(PBeta);

           // solve the scattering problem for this incident spherical wave
           Log("Solving scattering problem with incident spherical wave #%i: (L,M,P)=(%i,%i,%i)",Beta,LBeta,MBeta,PBeta);
           G->AssembleRHSVector(Omega, &SW, KN);
           M->LUSolve(KN);

           // compute the full vector of spherical multipole moments induced by the
           // incident wave on the object and store it as the Betath column of the T-matrix
           HVector TColumn(NumMoments, LHM_COMPLEX, (cdouble *)TMatrix->GetColumnPointer(Beta));
           GetSphericalMoments(G, Omega, LMax, KN, &TColumn);
        }; // for (nc=l=0...)

     /*--------------------------------------------------------------*/
     /*- write the full content of the T-matrix at this frequency to */
     /*- the text output file                                        */
     /*--------------------------------------------------------------*/
     for(int LAlpha=1, Alpha=0; LAlpha<=LMax; LAlpha++)
      for(int MAlpha=-LAlpha; MAlpha<=LAlpha; MAlpha++)
       for(int PAlpha=0; PAlpha<2; PAlpha++, Alpha++)
        for(int LBeta=1, Beta=0; LBeta<=LMax; LBeta++)
         for(int MBeta=-LBeta; MBeta<=LBeta; MBeta++)
          for(int PBeta=0; PBeta<2; PBeta++, Beta++)
           fprintf(f,"%s  %i %i %+i %i   %i %i %+i %i  %+15.8e %+15.8e\n",
                      z2s(Omega), 
                      Alpha,LAlpha,MAlpha,PAlpha,Beta,LBeta,MBeta,PBeta,
                      real(TMatrix->GetEntry(Alpha,Beta)),
                      imag(TMatrix->GetEntry(Alpha,Beta)));

     if (WriteHDF5Files)
      { char FileName[100];
        snprintf(FileName,100,"%s_w%s.HDF5",FileBase,z2s(Omega));
        TMatrix->ExportToHDF5(FileName,"T");
      };

    }; // for( nOmega= ... )

  fclose(f);
      
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (Cache)
   StoreCache(Cache);
  printf("Thank you for your support.\n");
  
}
