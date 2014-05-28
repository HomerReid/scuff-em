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
 * SumsIntegrals.cc -- routines for evaluating matsubara sums,
 *                  -- imaginary-frequency integrals, and 
 *                  -- brillouin-zone integrations
 *
 * homer reid       -- 2/2012
 *
 */

#include "scuff-cas3D.h"
#include "libscuffInternals.h"
#include <libSGJC.h>

#define XIMIN  0.001
#define XIMAX 10.000

/*****************************************************************/
/* 3, 6, 10, 15-point monkhorst-pack grids for square-lattice BZ */
/*****************************************************************/
double MP3[]=
 { 7.853982e-01, 7.853982e-01, 0.25,
   7.853982e-01, 2.356194e+00, 0.5,
   2.356194e-01, 2.356194e+00, 0.25
 };

double MP6[]=
 { 5.235988e-01, 5.235988e-01, 0.111111111111,
   5.235988e-01, 1.570796e+00, 0.222222222222,
   5.235988e-01, 2.617994e+00, 0.222222222222,
   1.570796e+00, 1.570796e+00, 0.111111111111,
   1.570796e+00, 2.617994e+00, 0.222222222222,
   2.617994e+00, 2.617994e+00, 0.111111111111
 };

double MP10[]=
 { 3.926991e-01, 3.926991e-01, 0.0625,
   3.926991e-01, 1.178097e+00, 0.1250, 
   3.926991e-01, 1.963495e+00, 0.1250,
   3.926991e-01, 2.748894e+00, 0.1250,
   1.178097e+00, 1.178097e+00, 0.0625,
   1.178097e+00, 1.963495e+00, 0.1250,
   1.178097e+00, 2.748894e+00, 0.1250,
   1.963495e+00, 1.963495e+00, 0.0625,
   1.963495e+00, 2.748894e+00, 0.1250,
   2.748894e+00, 2.748894e+00, 0.0625
 };

double MP15[]=
 { 3.141593e-01, 3.141593e-01, 0.04,
   3.141593e-01, 9.424778e-01, 0.08,
   3.141593e-01, 1.570796e+00, 0.08,
   3.141593e-01, 2.199115e+00, 0.08,
   3.141593e-01, 2.827433e+00, 0.08,
   9.424778e-01, 9.424778e-01, 0.04,
   9.424778e-01, 1.570796e+00, 0.08,
   9.424778e-01, 2.199115e+00, 0.08,
   9.424778e-01, 2.827433e+00, 0.08,
   1.570796e+00, 1.570796e+00, 0.04,
   1.570796e+00, 2.199115e+00, 0.08,
   1.570796e+00, 2.827433e+00, 0.08,
   2.199115e+00, 2.199115e+00, 0.04,
   2.199115e+00, 2.827433e+00, 0.08,
   2.827433e+00, 2.827433e+00, 0.04
 };

/***************************************************************/
/* CacheRead: attempt to bypass an entire GetXiIntegrand       */
/* calculation by reading results from the .byXi file.         */
/* Returns 1 if successful (which means the values of the      */
/* energy/force/torque integrand for ALL transformations at    */
/* this value of Xi were successfully read from the file) or 0 */
/* on failure).                                                */
/***************************************************************/
int CacheRead(const char *ByXiFileName, SC3Data *SC3D, double Xi, double *EFT)
{ 
  if (SC3D->UseExistingData==false)
   return 0;

  FILE *f;
  double Q[4];
  char Line[1000], fTag[1000];
  double fXi;
  int nt, ntnq, nRead, FoundFirst;

  /*----------------------------------------------------------*/
  /* 0. try to open the cache file. --------------------------*/
  /*----------------------------------------------------------*/
  if ( !(f=fopen(ByXiFileName,"r")) )
   return 0;

  /*----------------------------------------------------------*/
  /*----------------------------------------------------------*/
  /*----------------------------------------------------------*/
  for(;;)
   {
     /*----------------------------------------------------------*/
     /* 1. skip down through the cache file until we find a line */
     /*    whose Xi and Tag values equal Xi and the first        */
     /*    tag in the workspace structure.                       */
     /*----------------------------------------------------------*/
     FoundFirst=0;
     while( !FoundFirst && fgets(Line,1000,f) )
      { sscanf(Line,"%s %le",fTag,&fXi);
        if ( fabs(fXi-Xi) < 1.0e-8*Xi && !strcmp(fTag,SC3D->GTCList[0]->Tag) )
         FoundFirst=1;
      };
     if ( !FoundFirst ) 
      { fclose(f); 
        return 0;
      };
   
     Log(" found (Tag,Xi)=(%s,%e) in cache file...",fTag,Xi);
   
     /*----------------------------------------------------------*/
     /* 2. verify that the line we just read from the cache file */
     /*    contains data for all the quantities we need          */
     /*----------------------------------------------------------*/ 
     nRead=sscanf(Line,"%s %le %le %le %le %le",fTag,&fXi,Q,Q+1,Q+2,Q+3); 
     if ( nRead != SC3D->NumQuantities+2 )
      { Log(" ...but number of quantities is wrong (skipping)");
        continue;
      };
     memcpy(EFT,Q,SC3D->NumQuantities*sizeof(double));
     ntnq=SC3D->NumQuantities;
   
     /*----------------------------------------------------------*/
     /* 3. ok, since that worked, now keep going ----------------*/
     /*----------------------------------------------------------*/
     for(nt=1; nt<SC3D->NumTransformations; nt++)
      { 
        /* check for premature end of file */
        if ( !fgets(Line,1000,f) )
         { Log(" ...but data for some transforms were missing (skipping)");
           break;
         };
   
        nRead=sscanf(Line,"%s %le %le %le %le %le",fTag,&fXi,Q,Q+1,Q+2,Q+3);
   
        /* check for incorrect number of quantities */
        if ( nRead != SC3D->NumQuantities+2 )
         { Log(" ...but number of quantities is wrong (skipping)");
           break;
         };
   
        /* check for tag and/or Xi mismatch */
        if ( fabs(fXi-Xi)>1.0e-8*Xi || strcmp(fTag,SC3D->GTCList[nt]->Tag) )
         { Log(" ...but tag #%i did not match (%s != %s) (skipping)",
               nt,fTag,SC3D->GTCList[nt]->Tag);
           break;
         };
   
        memcpy(EFT+ntnq,Q,SC3D->NumQuantities*sizeof(double));
        ntnq+=SC3D->NumQuantities;
      };
   
     if (ntnq==SC3D->NTNQ)
      { Log(" ...and successfully read data for all quantities at all transforms");
        fclose(f);
        return 1;
      };

   };

}

/***************************************************************/
/* wrapper around GetCasimirIntegrand with correct prototype   */
/* prototype for adapt_integrate()                             */
/***************************************************************/
int GetCasimirIntegrand2(unsigned ndim, const double *x, void *params,
                         unsigned fdim, double *fval)
{
  (void) ndim; // unused
  (void) fdim; // unused

  SC3Data *SC3D = (SC3Data *)params;
  double *EFT = fval;

  double kBloch[2];

  kBloch[0] = x[0]*SC3D->RLBasisVectors[0][0] + x[1]*SC3D->RLBasisVectors[1][0];
  kBloch[1] = x[0]*SC3D->RLBasisVectors[0][1] + x[1]*SC3D->RLBasisVectors[1][1];

  GetCasimirIntegrand(SC3D, SC3D->Xi, kBloch, EFT);

  return 0;

}


/***************************************************************/
/* get the contribution of a single imaginary angular frequency*/
/* to the Casimir quantities. For periodic geometries, this    */
/* means integrating over the Brillouin zone at the Xi value   */
/* in question.                                                */
/***************************************************************/
void GetXiIntegrand(SC3Data *SC3D, double Xi, double *EFT)
{
  /***************************************************************/
  /* attempt to bypass calculation by reading data from .byXi file */
  /***************************************************************/
  if ( CacheRead(SC3D->ByXiFileName, SC3D, Xi, EFT) )
   return; 

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (SC3D->G->NumLatticeBasisVectors==0)
   { 
     GetCasimirIntegrand(SC3D, Xi, 0, EFT);
   }
  else if (SC3D->BZIMethod == BZIMETHOD_ADAPTIVE)
   {
     double Lower[2] = {0.0, 0.0};
     double Upper[2] = {1.0, 1.0};
     double *Error = new double[SC3D->NTNQ];

     SC3D->Xi = Xi;
     pcubature(SC3D->NTNQ, GetCasimirIntegrand2, (void *)SC3D, 2, Lower, Upper,
               SC3D->MaxkBlochPoints, SC3D->AbsTol, SC3D->RelTol, ERROR_INDIVIDUAL,
               EFT, Error);

     for(int ntnq=0; ntnq<SC3D->NTNQ; ntnq++)
      { EFT[ntnq] *= SC3D->BZVolume;
        Error[ntnq] *= SC3D->BZVolume;

        if (    (Error[ntnq] > 10.0*SC3D->AbsTol) 
             || (Error[ntnq] > 10.0*SC3D->RelTol*fabs(EFT[ntnq]))
           )
         Warn("potentially large errors (Q%i: %.1e %%) in BZ integration",
               ntnq,Error[ntnq]/fabs(EFT[ntnq]));
      };
     
     delete[] Error;
   }
  else if (SC3D->BZIMethod == BZIMETHOD_CC5917 )
   { 
      ErrExit("CC5917 not yet implemented");
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  FILE *f=fopen(SC3D->ByXiFileName,"a");
  for(int ntnq=0, nt=0; nt<SC3D->NumTransformations; nt++)
   { fprintf(f,"%s %.6e ",SC3D->GTCList[nt]->Tag,Xi);
     for(int nq=0; nq<SC3D->NumQuantities; nq++, ntnq++) 
      fprintf(f,"%.8e ",EFT[ntnq]);
     fprintf(f,"\n");
     fflush(f);
  };
  fclose(f);

}

/***************************************************************/
/* wrapper with correct prototype for adapt_integrate          */
/***************************************************************/
int GetXiIntegrand2(unsigned ndim, const double *x, void *params, 
                    unsigned fdim, double *fval)
{
  (void) ndim; // unused
  (void) fdim; // unused

  double Xi = x[0]/(1.0-x[0]); 
  double Jacobian = 1.0/( (1.0-x[0])*(1.0-x[0]) );
  SC3Data *SC3D = (SC3Data *)params;
  double *EFT = fval;

  GetXiIntegrand(SC3D, Xi, EFT);

  for(int ntnq=0; ntnq<SC3D->NTNQ; ntnq++)
   EFT[ntnq]*=Jacobian;

  return 0;

}

/***************************************************************/
/* Evaluate the Matsubara sum to get total Casimir quantities  */
/* at temperature Temperature degrees Kelvin.                  */
/*                                                             */
/* how the temperature conversion works:                       */
/*  a. temperature in eV = kT = 8.6173e-5 * (T in kelvin)      */
/*  b. temperature in our internal energy units                */
/*     = (kT in eV) / (0.1973 eV)                              */
/*     = (8.6173e-5 / 0.1973 ) * (T in Kelvin)                 */
/***************************************************************/
#define BOLTZMANNK 4.36763e-4
void GetMatsubaraSum(SC3Data *SC3D, double Temperature, double *EFT, double *Error)
{ 
  int n, ntnq, NTNQ=SC3D->NTNQ;
  double Xi, Weight;

  double *dEFT = new double [NTNQ]; 
  double *LastEFT = new double[NTNQ]; 
  double RelDelta;
  int AllConverged=0;
  int *ConvergedIters = new int [NTNQ];

  double kT = BOLTZMANNK * Temperature;

  memset(EFT,0,NTNQ*sizeof(double));
  memset(SC3D->Converged,0,NTNQ*sizeof(int));
  memset(ConvergedIters,0,NTNQ*sizeof(int));

  Log("Beginning Matsubara sum at T=%g kelvin...",Temperature);

  for(n=0; n<SC3D->MaxXiPoints; n++)
   { 
     /***************************************************************/
     /* compute the next matsubara frequency ************************/
     /***************************************************************/
     if (n==0)
      { Weight=0.5;
        // NOTE: we assume that the integrand is constant for Xi < XIMIN
        Xi=XIMIN;
      }
     else
      { Weight=1.0;
        Xi=2.0*M_PI*kT*((double)n);
      };

     /***************************************************************/
     /* evaluate the frequency integrand at this matsubara frequency*/
     /***************************************************************/
     GetXiIntegrand(SC3D, Xi, dEFT);

     /***************************************************************/
     /* accumulate contributions to the sum.                        */
     /*                                                             */
     /* how it works: the matsubara sum is                          */
     /*  2\pi kT *  \sum_n^\prime F(\xi_n)                          */
     /* where \xi_n is the nth matsubara frequency and F(\xi) is    */
     /* the casimir integrand (and the primed sum means that the    */
     /* n==0 term is weighted with a factor of 1/2).                */
     /*                                                             */
     /* however, my GetXiIntegrand() routine returns the quantity   */
     /* FI = F(\xi_n) / (2\pi). (this is so that the integral of FI */
     /* over all \xi returns the correct casimir quantity with no   */
     /* additional multiplicative prefactors.)                      */
     /*                                                             */
     /* thus the matsubara sum is                                   */
     /*  4\pi^2 kT *  \sum_n^\prime FI(\xi_n)                       */
     /*                                                             */
     /* where FI is what is returned by GetXiIntegrand.             */
     /***************************************************************/
     memcpy(LastEFT,EFT,NTNQ*sizeof(double));
     for(ntnq=0; ntnq<NTNQ; ntnq++)
      EFT[ntnq] += Weight * 4.0*M_PI*M_PI* kT * dEFT[ntnq];

     /*********************************************************************/
     /* convergence analysis.                                             */
     /* how it works: if the relative change in output quantity #ntnq is  */
     /* less than EPSREL, we increment ConvergedIters[ntnq]; otherwise we */
     /* set ConvergedIters[ntnq] to 0.                                    */
     /* when ConvergedIters[ntnq] hits 2, we mark quantity #ntnq as       */
     /* having converged.                                                 */
     /* when all quantities have converged, we are done.                  */
     /*********************************************************************/
     for(AllConverged=1, ntnq=0; ntnq<NTNQ; ntnq++)
      { 
        if ( SC3D->Converged[ntnq] ) 
         continue;

        Error[ntnq] = fabs(EFT[ntnq] - LastEFT[ntnq]);
        RelDelta = Error[ntnq] / fabs(EFT[ntnq]);
        if ( RelDelta < 1.0e-6 )
         ConvergedIters[ntnq]++;
        else
         ConvergedIters[ntnq]=0;

        if ( ConvergedIters[ntnq]>=2 )
         SC3D->Converged[ntnq]=1;
        else  
         AllConverged=0;
      }; 

     if (AllConverged==1)
      break;

   }; /* for (n=0 ... */

  delete[] dEFT;   
  delete[] LastEFT;
  delete[] ConvergedIters;
  
  if (AllConverged==0)
   { 
     fprintf(stderr,"\n*\n* WARNING: Matsubara sum unconverged after %i frequency samples.\n*\n",n);
     Log("Matsubara sum UNCONVERGED at n=%i samples",n); 
   } 
  else 
   Log("Matsubara sum converged after summing n=%i frequency points.",n);
    
} 

/***************************************************************/
/* Integrate over the positive imaginary frequency axis to get */
/* the total Casimir quantities at zero temperature.           */
/***************************************************************/
void GetXiIntegral(SC3Data *SC3D, double *EFT, double *Error)
{
  double Lower[1] = {0.0}; 
  double Upper[1] = {1.0};

  pcubature(SC3D->NTNQ, GetXiIntegrand2, (void *)SC3D, 1, Lower, Upper,
            SC3D->MaxXiPoints, SC3D->AbsTol, SC3D->RelTol,
            ERROR_INDIVIDUAL, EFT, Error);
   
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetXiIntegral2(SC3Data *SC3D, int NumIntervals, double *I, double *E)
{ 
  int fdim = SC3D->NTNQ;
  double *fLeft  = new double[fdim];
  double *fMid   = new double[fdim];
  double *fRight = new double[fdim];
  double Xi;

  /*--------------------------------------------------------------*/
  /*- evaluate integrand at leftmost frequency point and estimate */
  /*- the integral from 0 to XIMIN by assuming that the integrand */
  /*- is constant in that range                                   */
  /*--------------------------------------------------------------*/
  Xi=XIMIN;
  GetXiIntegrand(SC3D, Xi, fLeft);
  for(int nf=0; nf<fdim; nf++)
   I[nf] = fLeft[nf] * XIMIN;
  memset(E,0,SC3D->NTNQ*sizeof(double));

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double Delta = (XIMAX - XIMIN ) / NumIntervals;
  double *ISimp  = new double[fdim];
  double *ITrap  = new double[fdim];
  for(int nIntervals=0; nIntervals<NumIntervals; nIntervals++)
   { 
     // evaluate integrand at midpoint of interval 
     Xi += 0.5*Delta;
     GetXiIntegrand(SC3D, Xi, fMid);

     // evaluate integrand at right end of interval 
     Xi += 0.5*Delta;
     GetXiIntegrand(SC3D, Xi, fRight);

     // compute the simpson's rule and trapezoidal rule
     // estimates of the integral over this interval  
     // and take their difference as the error
     for(int nf=0; nf<fdim; nf++)
      { ISimp[nf] = (fLeft[nf] + 4.0*fMid[nf] + fRight[nf])*Delta/6.0;
        ITrap[nf] = (fLeft[nf] + 2.0*fMid[nf] + fRight[nf])*Delta/4.0;
        I[nf] += ISimp[nf];
        E[nf] += fabs(ISimp[nf] - ITrap[nf]);
      };

     // prepare for next iteration
     memcpy(fLeft, fRight, fdim*sizeof(double));

   };
  delete[] fLeft;
  delete[] fMid;
  delete[] fRight;
  delete[] ISimp;
  delete[] ITrap;

}
