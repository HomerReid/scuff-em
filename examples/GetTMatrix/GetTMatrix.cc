/*
 * GetTMatrix.cc -- a slightly more complicated example of how to use the
 *               -- scuff-EM C++ API
 *
 *               -- This program computes the T-matrix, in a spherical-wave basis,  
 *               -- for an arbitary scatterer. (What this means is the following: 
 *               -- For each spherical wave up to a user-specified maximum L value, 
 *               -- the code scatters that wave off of the geometry and decomposes
 *               -- the resulting scattered field into a basis of spherical waves.)
 *               -- The implementation of this code demonstrates, among other things, 
 *               -- how to define your own type of incident field routine.
 *
 * Homer Reid    -- 6/2012
 */

#include <stdio.h>
#include <stdlib.h>

#include <libscuff.h>
using namespace scuff;

#include "GetTMatrix.h"

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  /*--------------------------------------------------------------*/
  /*- process command-line options -------------------------------*/
  /*--------------------------------------------------------------*/
  char *GeoFileName=0;  // .scuffgeo file
  int lMax=5;           // maximum L-value of spherical wave computed
  cdouble Omega=0;      // angular frequency at which to run the computation
  char *OmegaFile=0;    // list of angular frequencies
  char *Cache;          // scuff cache file 
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"geometry",  PA_STRING,  1, 1, (void *)&GeoFileName,  0,  ".scuffgeo file"},
     {"lMax",      PA_INT,     1, 1, (void *)&lMax,         0,  "maximum l-value"},
     {"Omega",     PA_CDOUBLE, 1, 1, (void *)&Omega,        0,  "angular frequency"},
     {"OmegaFile", PA_STRING,  1, 1, (void *)&OmegaFile,    0,  "list of angular frequencies"},
     {"Cache",     PA_STRING,  1, 1, (void *)&Cache,        0,  "scuff cache file"},
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
     OmegaVector->SetEntry(1,Omega);
   }
  else
   OSUsage(argv[0],OSArray,"either --omega or --omegafile must be specified");

  /*--------------------------------------------------------------*/
  /* create the RWGGeometry from the .scuffgeo file               */
  /*--------------------------------------------------------------*/
  SetLogFileName("%s.log",GetFileBase(GeoFileName));
  RWGGeometry *G=new RWGGeometry(GeoFileName);
  G->SetLogLevel(SCUFF_VERBOSELOGGING);
  PreloadCache(Cache);

  /*--------------------------------------------------------------*/
  /* preallocate BEM matrix and RHS vector                        */
  /*--------------------------------------------------------------*/
  HMatrix *M  = G->AllocateBEMMatrix();
  HVector *KN = G->AllocateRHSVector();

  /*--------------------------------------------------------------*/
  /*- preallocate an HMatrix to store the T-matrix data           */
  /*--------------------------------------------------------------*/
  int Dimension = 2*(lMax+1)*(lMax+1);
  HMatrix *TMatrix = new HMatrix(Dimension, Dimension, LHM_COMPLEX);

  /*--------------------------------------------------------------*/
  /* instantiate a SphericalWave structure (we will set the l, m, */
  /* and Type fields later)                                       */
  /*--------------------------------------------------------------*/
  SphericalWave SW(0, 0, SW_ELECTRIC);

  /*--------------------------------------------------------------*/
  /*- outer loop over frequencies --------------------------------*/
  /*--------------------------------------------------------------*/
  double Omega;
  int Type, l, m;
  int TypeP, lP, mP;
  int nr, nc;
  for(int nOmega=0; nOmega<OmegaVector->N; nOmega++)
   {
     /*--------------------------------------------------------------*/
     /* assemble and factorize the BEM matrix at this frequency      */
     /*--------------------------------------------------------------*/
     Omega=OmegaVector->GetEntry(nOmega);
     G->AssembleBEMMatrix(Omega, M);
     M->LUFactorize();

     /*--------------------------------------------------------------*/
     /*- inner loop over incident spherical waves (i.e. over rows   -*/
     /*- of the T matrix; note nr is a running row index)           -*/
     /*--------------------------------------------------------------*/
     for(nr=Type=0; Type<=1; Type++)
      for(l=0; l<=lMax; l++)
       for(m=-l; m<=l; m++, nr++)
        { if (l==0) 
           continue;

           SW.SetType(Type);
           SW.SetL(l);
           SW.SetM(m);

           // solve the scattering problem for this incident spherical wave
           G->AssembleRHSVector(Omega, &SW, KN);
           M->LUSolve(KN);

           // compute the projection of the scattered field onto each
           // spherical wave

        }; // for (nr=Type= ... )

    }; // for( nOmega= ... )
      
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  fclose(f);
  StoreCache(Cache);
  printf("Thank you for your support.\n");
  
}
