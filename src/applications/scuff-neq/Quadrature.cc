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
 * FrequencyIntegral.cc -- scuff-neq module for numerical quadrature 
 *                      -- over frequencies
 *
 * homer reid           -- 5/2012
 *
 */

#include <libSGJC.h>
#include "scuff-neq.h"

// how this works: 
//  a. temperature in eV = kT = 8.6173e-5 * (T in kelvin)
//  b. temperature in our internal energy units
//     = (kT in eV) / (0.1973 eV)
//     = (8.6173e-5 / 0.1973 ) * (T in Kelvin)
#define BOLTZMANNK 4.36763e-4

// for example: suppose in the real world we have 
// omega=3e14 rad/sec, T = 300 kelvin. then
// \hbar \omega / (kT) 
//   = (6.6e-16 ev s )(3e14 s^{-1}) / (0.026 ev) 
//   = 7.6
// whereas in this code we would have Omega==1, T=300, and hence
// Omega/(BOLTZMANNK*T) = (1/(4.36763e-4*300)) = 7.6. 

/***************************************************************/
/***************************************************************/
/***************************************************************/
double Theta(double Omega, double T)
{ 
  return Omega / ( exp( (Omega/BOLTZMANNK*T) ) - 1.0 );
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct FIData 
 {
   SNEQData *SNEQD;
   double OmegaMin;
   int Infinite;
   double *TObjects;
   double TEnvironment;

 } FIData;

void SGJCIntegrand(unsigned ndim, const double *x, void *params,
                   unsigned fdim, double *fval)
{
  (void) ndim; // unused
  (void) fdim;

  FIData *FID         = (FIData *)params;

  SNEQData *SNEQD     = FID->SNEQD; 
  double OmegaMin     = FID->OmegaMin;
  int Infinite        = FID->Infinite;
  double *TObjects    = FID->TObjects;
  double TEnvironment = FID->TEnvironment;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double Omega, Jacobian; 
  if(Infinite)
   { Omega    = OmegaMin + x[0] / (1.0-x[0]);
     Jacobian = 1.0 / ( (1.0-x[0]) * (1.0-x[0]) );
   }
  else
   { Omega    = x[0];
     Jacobian = 1.0;
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  GetFrequencyIntegrand(SNEQD, Omega, fval);

  /*--------------------------------------------------------------*/
  /*- quantities arising from sources inside object nop are       */
  /*- weighted by a factor [Theta(T) - Theta(TEnv)]                */
  /*--------------------------------------------------------------*/
  int nt, nq, no, nop;
  int NO = SNEQD->G->NumObjects;
  int NT = SNEQD->NumTransformations;
  int NQ = SNEQD->NQ;
  int NONQ = NO*NQ;
  int NO2NQ = NO*NO*NQ;
  double DeltaTheta;
  for(nop=0; nop<NO; nop++)
   { DeltaTheta = Theta(Omega, TObjects[nop]) - Theta(Omega, TEnvironment);
     for(nt=0; nt<NT; nt++)
      for(no=0; no<NO; no++)
       for(nq=0; nq<NQ; nq++)
        fval[ nt*NO2NQ + no*NONQ + nop*NQ + nq ] *= Jacobian*DeltaTheta/M_PI;
   };
     
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void EvaluateFrequencyIntegral(SNEQData *SNEQD, 
                               double OmegaMin, double OmegaMax,
                               double *TObjects, double TEnvironment, 
                               double AbsTol, double RelTol,
                               double *I, double *E)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  FIData MyFIData, *FID=&MyFIData;
  FID->SNEQD = SNEQD;
  FID->OmegaMin=OmegaMin;

  if (OmegaMax==-1.0)
   { FID->Infinite=1; 
     OmegaMin=0.0;
     OmegaMax=1.0;
   }
  else
   FID->Infinite=0;

  FID->TObjects=TObjects;
  FID->TEnvironment=TEnvironment;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  RWGGeometry *G = SNEQD -> G;
  int NO = G->NumObjects;
  int NT = SNEQD->NumTransformations;
  int NQ = SNEQD->NQ;
  int fdim = NT*NO*NO*NQ;
  adapt_integrate_log(fdim, SGJCIntegrand, (void *)FID, 1, &OmegaMin, &OmegaMax,
                      1000, AbsTol, RelTol, I, E, "scuff-neq.SGJClog",15);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int nt, nq, no, nop;
  int NONQ = NO*NQ;
  int NO2NQ = NO*NO*NQ;
  FILE *f;
  char FileName[1000];
  for(no=0; no<NO; no++)
   for(nop=0; nop<NO; nop++)
    { 
      snprintf(FileName,1000,"From%sTo%s.out",G->Objects[no]->Label,G->Objects[nop]->Label);
      f=CreateUniqueFile(FileName,1);
      for(nt=0; nt<NT; nt++)
       { fprintf(f,"%s ",SNEQD->GTCList[nt]->Tag);
         for(nq=0; nq<NQ; nq++)
          fprintf(f,"%e %e ", I[ nt*NO2NQ + no*NONQ + nop*NQ + nq ],
                              E[ nt*NO2NQ + no*NONQ + nop*NQ + nq ] );
         fprintf(f,"\n");
       };
      fclose(f);
    };

}
