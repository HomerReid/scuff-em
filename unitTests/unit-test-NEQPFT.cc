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
 * unit-test-NEQPFT.cc -- SCUFF-EM unit tests for non-equilibrium
 *                     -- power, force, and torque calculations
 * 
 * homer reid          -- 11/2014
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
/* tables of non-equilibrium PFT data for SiC spheres of       */
/* radius R=1 um separated by a (center--center) distance of   */
/* 10 um in the Z direction.                                   */
/***************************************************************/

#define NUMOMEGAS 3
double OmegaList[] = { 1.0e-2, 1.0e-1, 1.0 };

// PowerFlux[no][0] == 0->0 power at frequency #no
// PowerFlux[no][1] == 0->1 power at frequency #no
// PowerFlux[no][2] == 1->0 power at frequency #no
// PowerFlux[no][3] == 1->1 power at frequency #no
// 
// ForceFlux[no][0] == 0->0 force at frequency #no
// ForceFlux[no][1] == 0->1 force at frequency #no
// ForceFlux[no][2] == 1->0 force at frequency #no
// ForceFlux[no][3] == 1->1 force at frequency #no

double PowerFlux[NUMOMEGAS][4]
 { };

double ForceFlux[NUMOMEGAS][4]
 { };

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{ 
  (void )argc;
  (void )argv; 
  SetLogFileName("scuff-unit-tests.log");
  Log("SCUFF-EM NEQPFT unit tests running on %s",GetHostName());

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

        double OPFT[7], PTot;
        G->GetOPFT(0, Omega, KN, RHS, 0, OPFT, &PTot);

        double PFT[8];
        PFT[0]=OPFT[0];
        PFT[1]=PTot-PFT[0];
        memcpy(PFT+2, OPFT+1, 6*sizeof(double));

        double Denom = M_PI*IncFlux;
        //double QScatScuff  = G->GetScatteredPower(KN, Omega, 0) / Denom;
        double QScatScuff  = 0.0; // FIXME
        double QAbsScuff   = PFT[0] / Denom;
        double QForceScuff = PFT[5] / (TENTHIRDS*Denom);

        double QScatExact  = ExactData[nCase][4*nOmega + 1];
        double QAbsExact   = ExactData[nCase][4*nOmega + 2];
        double QForceExact = ExactData[nCase][4*nOmega + 3];
  
        if (    RD( QScatScuff , QScatExact ) > 0.05
             || RD( QAbsScuff  , QAbsExact )  > 0.05
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
                   QScatScuff,  QScatExact,  RD(QScatScuff,  QScatExact), 
                   QAbsScuff,   QAbsExact,   RD(QAbsScuff,   QAbsExact),
                   QForceScuff, QForceExact, RD(QForceScuff, QForceExact));
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
