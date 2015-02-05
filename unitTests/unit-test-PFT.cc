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
 * unit-test-PFT.cc -- SCUFF-EM unit tests for PFT calculations
 * 
 * homer reid       -- 11/2005 -- 10/2011
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>
#include "libscuff.h"
#include "PFTOptions.h"
#include "libscuffInternals.h"
#include "libIncField.h"

using namespace scuff;

#define II cdouble (0.0,1.0)
#define TENTHIRDS 3.33333333333333333334
#define RT1_2     0.70710678118654752440

/***************************************************************/
/* tables of power and force data for spheres of radius R=1um  */
/* and various material compositions.                          */
/*                                                             */
/* in each case, powers and forces are for the case of         */
/* illumination by a unit-magnitude z-traveling plane wave     */
/* with left-circular polarization, ie the incident E-field is */
/*  E(x,y,z) = (1/sqrt{2}) (\vec{x} + i \vec{y}) e^{ikz}       */
/*                                                             */
/* data table entries:                                         */
/*  (1) omega (angular frequency, units 3e14 rad/sec)          */
/*  (2) absorbed  power (watts)                                */
/*  (2) scattered power (watts)                                */
/*  (3) Z-force (nanoNewtons)                                  */
/*  (4) Z-torque (nanoNewtons*microns)                         */
/***************************************************************/
#if 0
#define NUMOMEGAS 2
double OmegaList[] = {1.0e-1, 1.0};
#endif
#define NUMOMEGAS 11
double OmegaList[] = 
 { 0.01000000, 0.01584893, 0.02511886, 0.03981072, 0.06309573,
   0.10000000, 0.15848932, 0.25118864, 0.39810717, 0.63095734,
   1.00000000};

double PECPFT[] =
   { 1.0e-2, 3.320868e-08, 0.0, 4.641536e-08,
     1.0e-1, 5.980295e+00, 0.0, 6.020131e+00,
     1.0e-0, 6.012197e+00, 0.0, 5.880804e+00
   };

double GoldPFT[] =
   { 1.0e-2, 3.111826e-08, 3.293537e-03, 3.293579e-03,
     1.0e-1, 3.244591e-04, 8.972653e-03, 9.419389e-03,
     1.0e-0, 2.125119e+00, 1.983620e-02, 2.509075e+00
   };

double SiO2PFT[] =
   { 1.0e-2, 1.504039e-08, 3.348068e-07, 3.498465e-07,
     1.0e-1, 1.527705e-04, 3.648611e-05, 1.885546e-04,
     1.0e-0, 1.465431e+00, 4.180791e-03, 8.951046e-01 
   };

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{ 
  SetLogFileName("scuff-unit-tests.log");
  Log("SCUFF-EM PFT unit tests running on %s",GetHostName());

  /***************************************************************/
  /* set up incident field ***************************************/
  /***************************************************************/
  const cdouble E0[3]  = { RT1_2, II*RT1_2, 0.0 };
  const double nHat[3] = { 0.0, 0.0, 1.0 };
  PlaneWave *PW = new PlaneWave(E0, nHat);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  PFTOptions *Options=InitPFTOptions();

  /***************************************************************/
  /* loop over all three material geometries *********************/
  /***************************************************************/
  #define NUMCASES 3
  const char *GeoFileNames[NUMCASES] = { "PECSphere_501.scuffgeo",
                                         "SiO2Sphere_501.scuffgeo",
                                         "GoldSphere_501.scuffgeo"
                                       };
  double *ExactData[NUMCASES] = { PECPFT, SiO2PFT, GoldPFT };
  int NumTests = NUMCASES * NUMOMEGAS;
  int PassedTests = 0 ;
  FILE *DataLogFile=fopen("unit-test-PFT.data","w");
  for(int nCase=0; nCase<NUMCASES; nCase++)
   { 
     RWGGeometry *G = new RWGGeometry(GeoFileNames[nCase]);
     G->SetLogLevel(SCUFF_VERBOSELOGGING);
     RWGSurface  *S = G->Surfaces[0];
     HMatrix *TIn   = S->IsPEC ? 0 : G->AllocateBEMMatrix();
     HMatrix *TOut  = G->AllocateBEMMatrix();
     HMatrix *M     = G->AllocateBEMMatrix();
     HVector *KN    = G->AllocateRHSVector();
     HVector *RHS   = G->AllocateRHSVector();

     if (DataLogFile)
      fprintf(DataLogFile,"# geometry %s: \n",G->GeoFileName);

     for (int nOmega=0; nOmega<NUMOMEGAS; nOmega++)
      { 
        // solve the scattering problem
        double Omega=OmegaList[nOmega];

        if (S->IsPEC)
         {
           G->AssembleBEMMatrix(Omega, TOut);
           M->Copy(TOut);
         }
        else
         { for(int nr=0; nr<G->NumRegions; nr++)
            G->RegionMPs[nr]->Zero();

           G->RegionMPs[1]->UnZero();
           G->AssembleBEMMatrix(Omega, TIn);
           G->RegionMPs[1]->Zero();

           G->RegionMPs[0]->UnZero();
           G->AssembleBEMMatrix(Omega, TOut);
           G->RegionMPs[0]->Zero();
 
           for(int nr=0; nr<G->NumRegions; nr++)
            G->RegionMPs[nr]->UnZero();

           M->Copy(TOut);
           M->AddBlock(TIn,0,0);
         };

        M->LUFactorize();
        G->AssembleRHSVector(Omega, PW, RHS);
        KN->Copy(RHS);
        M->LUSolve(KN);

        // compute overlap PFT
        InitPFTOptions(Options);
        Options->PFTMethod=SCUFF_PFT_OVERLAP;
        Options->RHSVector=RHS;
        double OPFT[8];
        G->GetPFT(0, PW, KN, Omega, OPFT, Options);
  
        // compute DSIPFT
        InitPFTOptions(Options);
        Options->PFTMethod=SCUFF_PFT_DSI;
        double DSIPFT[8];
        G->GetPFT(0, PW, KN, Omega, DSIPFT, Options);

        // compute EP without using precomputed matrices
        InitPFTOptions(Options);
        Options->PFTMethod=SCUFF_PFT_EP;
        double EPPFTWithout[8];
        G->GetPFT(0, PW, KN, Omega, EPPFTWithout, Options);

        // compute EPPFT using precomputed matrices
        InitPFTOptions(Options);
        Options->PFTMethod=SCUFF_PFT_EP;
        Options->TInterior=TIn;
        Options->TExterior=TOut;
        double EPPFTWith[8];
        G->GetPFT(0, PW, KN, Omega, EPPFTWith, Options);

        fprintf(DataLogFile,"%e %e %e %e %e %e %e %e %e %e %e %e %e \n",Omega,
                             OPFT[0],OPFT[1],OPFT[4],OPFT[7],
                             DSIPFT[0],DSIPFT[1],DSIPFT[4],DSIPFT[7],
                             EPPFTWith[0],EPPFTWith[1],
                             EPPFTWithout[0],EPPFTWithout[1]);
        fflush(DataLogFile);

#if 0
        if (    RD( OPFT[0], ExactPFT[0]     ) > 0.05
             || RD( OPFT[1], ExactPFT[1]     ) > 0.05
             || RD( OPFT[4], ExactPFT[4]     ) > 0.05
             || RD( OPFT[7], ExactPFT[7]     ) > 0.05
           )
          { Log("OPFT test failed for %s at Omega=%e: ",
                G->GeoFileName,Omega);
            Log(" (%e,%e) (%e,%e) (%e,%e) \n", 
                  QScatScuff,  QScatExact, 
                  QAbsScuff,   QAbsExact,
                  QForceScuff, QForceExact);
          } 
         else
          PassedTests++;

        if (DataLogFile)
         { fprintf(DataLogFile,"%e  %e %e (%e) %e %e (%e) %e %e (%e)\n",
                   Omega,
                   QScatScuff,  QScatExact,  RD(QScatScuff,  QScatExact), 
                   QAbsScuff,   QAbsExact,   RD(QAbsScuff,   QAbsExact),
                   QForceScuff, QForceExact, RD(QForceScuff, QForceExact));
           fflush(DataLogFile);
         };
#endif

      }; // for (int n=0; n<NUMOMEGAS; n++)

     if(DataLogFile)
      fprintf(DataLogFile,"\n\n");

     delete G; 
     delete TIn;
     delete TOut;
     delete M; 
     delete KN;
     delete RHS;

   }; // for(int nCase=0; nCase<NUMCASES; nCase++)

  if(DataLogFile)
   fclose(DataLogFile);
  
  /***************************************************************/ 
  /***************************************************************/ 
  /***************************************************************/ 
  printf("%i/%i tests successfully passed.\n",NumTests,PassedTests);

  int FailedCases=NumTests - PassedTests;
  if (FailedCases>0)
   abort();

  return 0;

}
