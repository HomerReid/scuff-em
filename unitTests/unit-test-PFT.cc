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
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>
#include "libscuff.h"
#include "libscuffInternals.h"
#include "libIncField.h"

using namespace scuff;

#define II cdouble (0.0,1.0)
#define TENTHIRDS 3.33333333333333333334

/***************************************************************/
/* tables of power and force data for spheres of radius R=1um  */
/* and various material composition.                           */
/*                                                             */
/* in each case, powers and forces are for the case of         */
/* illumination by a z-traveling plane wave with E-field       */
/* polarized in the +x direction.                              */
/*                                                             */
/* data table entries:                                         */
/*  (1) omega (angular frequency, units 3e14 rad/sec)          */
/*  (2) efficiency for power scattering                        */
/*  (3) efficiency for power absorption                        */
/*  (4) efficiency for z-force                                 */
/*                                                             */
/* where efficiencies are cross-sections divided by \pi R^2.   */
/*                                                             */
/* cross-sections for {absorbed, scattered} power are the      */
/* total power {absorbed by, scattered from} the object        */
/* divided by incident power flux. the cross-section for the   */
/* force is the z-force on the object divided by the incident  */
/* momentum flux.                                              */
/***************************************************************/
#define NUMOMEGAS 4
double OmegaList[] = { 1.0e-3, 1.0e-2, 1.0e-1, 1.0 };

double GoldPFT[] =
   { 1.0e-3, 2.841416e-12, 7.066437e-04, 7.066437e-04,
     1.0e-2, 3.111826e-08, 3.293537e-03, 3.293579e-03,
     1.0e-1, 3.244591e-04, 8.972653e-03, 9.419389e-03,
     1.0e-0, 2.125119e+00, 1.983620e-02, 2.509075e+00
   };

double PECPFT[] =
   { 1.0e-3, 3.216470e-12, 0.0, 4.427315e-12,
     1.0e-2, 3.320868e-08, 0.0, 4.641536e-08,
     1.0e-1, 5.980295e+00, 0.0, 6.020131e+00,
     1.0e-0, 6.012197e+00, 0.0, 5.880804e+00
   };

double SiCPFT[] =
   { 1.0e-3, 1.503808e-12, 3.345234e-09, 3.346738e-09,
     1.0e-2, 1.504039e-08, 3.348068e-07, 3.498465e-07,
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
  const cdouble E0[3]  = { 1.0, 0.0, 0.0 };
  const double nHat[3] = { 0.0, 0.0, 1.0 };
  PlaneWave *PW = new PlaneWave(E0, nHat);
  double IncFlux = 1.0/(2.0*ZVAC);

  /***************************************************************/ 
  /* loop over all three material geometries *********************/ 
  /***************************************************************/ 
  #define NUMCASES 3
  const char *GeoFileNames[NUMCASES] = { "PECSphere_474.scuffgeo", 
                                         "GoldSphere_474.scuffgeo",
                                         "SiCSphere_474.scuffgeo"
                                       };
  double *ExactData[NUMCASES] = { PECPFT, GoldPFT, SiCPFT };
  int NumTests = NUMCASES * NUMOMEGAS;
  int PassedTests = 0 ;
  FILE *DataLogFile=fopen("unit-test-PFT.data","w");
  for(int nCase=0; nCase<NUMCASES; nCase++)
   { 
     RWGGeometry *G = new RWGGeometry(GeoFileNames[nCase]);
     G->SetLogLevel(SCUFF_VERBOSELOGGING);
     HMatrix *M   = G->AllocateBEMMatrix();
     HVector *KN  = G->AllocateRHSVector();
     HVector *RHS = G->AllocateRHSVector();

     if (DataLogFile)
      fprintf(DataLogFile,"# geometry %s: \n",G->GeoFileName);

     for (int nOmega=0; nOmega<NUMOMEGAS; nOmega++)
      { 
        double Omega=OmegaList[nOmega];
        G->AssembleBEMMatrix(Omega, M);
        G->AssembleRHSVector(Omega, PW, RHS);
        KN->Copy(RHS);
        M->LUFactorize();
        M->LUSolve(KN);

        double PFT[8];
        G->GetPFT(KN, RHS, Omega, 0, PFT);

        double Denom = M_PI*IncFlux;
        double QScatScuff  = G->GetScatteredPower(KN, Omega, 0) / Denom;
        double QAbsScuff   = PFT[0] / Denom;
        double QForceScuff = PFT[5] / (TENTHIRDS*Denom);

        double QScatExact  = ExactData[nCase][4*nOmega + 1];
        double QAbsExact   = ExactData[nCase][4*nOmega + 2];
        double QForceExact = ExactData[nCase][4*nOmega + 3];
  
        if (    RD( QScatScuff , QScatExact ) > 0.05
             || RD( QAbsScuff  , QScatExact ) > 0.05
             || RD( QForceScuff, QForceExact) > 0.05
           )
          { Log("PFT test failed for %s at Omega=%e: ",
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
                   QScatScuff,  QScatExact, RD(QScatScuff,  QScatExact), 
                   QAbsScuff,   QAbsExact,  RD(QAbsScuff,  QAbsExact),
                   QForceScuff,   QForceExact,  RD(QForceScuff,  QForceExact));
           fflush(DataLogFile);
         };

      }; // for (int n=0; n<NUMOMEGAS; n++)

     if(DataLogFile)
      fprintf(DataLogFile,"\n\n");

     delete G; 
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
