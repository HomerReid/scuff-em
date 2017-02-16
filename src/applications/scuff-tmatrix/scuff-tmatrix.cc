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
 */

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
HVector *GetSphericalMoments(RWGGeometry *G, cdouble k, int lMax,
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
  char *GeoFileName=0;  // .scuffgeo file
  int lMax=5;           // maximum L-value of spherical wave computed
  cdouble Omega=0;      // angular frequency at which to run the computation
  char *OmegaFile=0;    // list of angular frequencies
  char *Cache=0;        // scuff cache file 
  char *FileBase=0;     // base filename for output file
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"geometry",  PA_STRING,  1, 1, (void *)&GeoFileName,  0,  ".scuffgeo file"},
     {"lMax",      PA_INT,     1, 1, (void *)&lMax,         0,  "maximum l-value"},
     {"Omega",     PA_CDOUBLE, 1, 1, (void *)&Omega,        0,  "angular frequency"},
     {"OmegaFile", PA_STRING,  1, 1, (void *)&OmegaFile,    0,  "list of angular frequencies"},
     {"Cache",     PA_STRING,  1, 1, (void *)&Cache,        0,  "scuff cache file"},
     {"FileBase",  PA_STRING,  1, 1, (void *)&FileBase,     0,  "base filename for output file"},
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
  int NumMoments= 2*(lMax+1)*(lMax+1);
  HMatrix *TMatrix = new HMatrix(NumMoments, NumMoments, LHM_COMPLEX);
  HVector *AVector = new HVector(NumMoments, LHM_COMPLEX);

  /*--------------------------------------------------------------*/
  /* instantiate a SphericalWave structure (we will set the l, m, */
  /* and Type fields later)                                       */
  /*--------------------------------------------------------------*/
  SphericalWave SW;

  /*--------------------------------------------------------------*/
  /*- outer loop over frequencies --------------------------------*/
  /*--------------------------------------------------------------*/
  int Type, l, m;
  int TypeP, lP, mP;
  const char *TypeChar="ME";
  int nr, nc;
  FILE *TextOutputFile;
  if (FileBase)
   TextOutputFile=vfopen("%s.TMatrix","w",FileBase);
  else
   TextOutputFile=vfopen("%s.TMatrix","w",GetFileBase(GeoFileName));
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
     /*- inner loop over incident spherical waves (i.e. over columns-*/
     /*- of the T matrix; note nc is a running column index)        -*/
     /*--------------------------------------------------------------*/
     TMatrix->Zero();
     for(nc=l=0; l<=lMax; l++)
      for(m=-l; m<=l; m++)
       for(Type=SW_MAGNETIC; Type<=SW_ELECTRIC; Type++, nc++)
        { 
           if (l==0) 
            continue;

           SW.SetType(Type);
           SW.SetL(l);
           SW.SetM(m);

           // solve the scattering problem for this incident spherical wave
           Log("Solving scattering problem with incident spherical wave P(l,m)=%c(%i,%i)",TypeChar[Type],l,m);
           G->AssembleRHSVector(Omega, &SW, KN);
           M->LUSolve(KN);

           // compute the spherical multipole moments induced by the 
           // incident wave on the object 
           GetSphericalMoments(G, Omega, lMax, KN, AVector);

           // stamp in the vector of moments as the ncth row of the T-matrix
           // NOTE: i don't know here the missing factor of -1.0 is coming
           // from here...
           for(nr=0; nr<NumMoments; nr++)
            TMatrix->SetEntry(nr, nc, -1.0*Omega*AVector->GetEntry(nr));

        }; // for (nc=l=0...)

     /*--------------------------------------------------------------*/
     /*- write the full content of the T-matrix at this frequency to */
     /*- the text output file                                        */
     /*--------------------------------------------------------------*/
     for(nr=l=0; l<=lMax; l++)
      for(m=-l; m<=l; m++)
       for(Type=SW_MAGNETIC; Type<=SW_ELECTRIC; Type++, nr++)
        for(nc=lP=0; lP<=lMax; lP++)
         for(mP=-lP; mP<=lP; mP++)
          for(TypeP=SW_MAGNETIC; TypeP<=SW_ELECTRIC; TypeP++, nc++)
           { 
             if (l==0 || lP==0) continue;
             fprintf(TextOutputFile,"%s %c %i %i %c %i %i  %+15.8e %+15.8e\n",
                                     z2s(Omega),
                                     TypeChar[Type],l,m,
                                     TypeChar[TypeP],lP,mP,
                                     real(TMatrix->GetEntry(nr,nc)),
                                     imag(TMatrix->GetEntry(nr,nc)));
           };

    }; // for( nOmega= ... )

  fclose(TextOutputFile);
      
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (Cache)
   StoreCache(Cache);
  printf("Thank you for your support.\n");
  
}
