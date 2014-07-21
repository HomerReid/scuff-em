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
#include <libTriInt.h>

#define XIMIN  0.001
#define XIMAX 10.000

/***************************************************************/
/* helper routine for CacheRead that, given a line of text     */
/* from a cache file, does the following:                      */
/*  1. checks that the first token on the line matches Tag     */
/*     and that the next NumbersToMatch tokens match the       */
/*     corresponding entries of Numbers                        */ 
/*  2. if (1) failed, returns false                            */
/*  3. otherwise, fills in Values with the remaining numbers   */
/*     on the line (up to 4) and returns the number of those   */
/*     numbers that were successfully read as *NumValues       */
/***************************************************************/
bool LineMatches(char *Line, char *Tag, double *Keys, int NumKeys,
                 double *Values, int *NumValues)
{
  char MyTag[1000];
  double Numbers[7];
  int nRead=sscanf(Line,"%s %le %le %le %le %le %le %le",
                         MyTag, Numbers+0, Numbers+1,
                         Numbers+2, Numbers+3, Numbers+4,
                         Numbers+5, Numbers+6);
  
  if ( nRead < (NumKeys+2) )
   return false;

  if ( strcasecmp(Tag, MyTag) ) 
   return false;

  for(int nk=0; nk<NumKeys; nk++)
   if ( fabs(Numbers[nk]-Keys[nk]) > 1.0e-6*fabs(Keys[nk]) )
    return false;

  *NumValues = nRead - 1 - NumKeys;
  for(int nv=0; nv < (*NumValues); nv++)
   Values[nv] = Numbers[NumKeys + nv];

  return true;

}

/***************************************************************/
/* CacheRead: attempt to bypass an entire GetXiIntegrand       */
/* calculation by reading results from the .byXi or .byXikbloch*/
/* file. Returns true if successful (which means the values of */
/* the energy/force/torque integrand for ALL transformations   */
/* at this value of Xi were successfully read from the file)   */
/* or false on failure.                                        */
/***************************************************************/
bool CacheRead(SC3Data *SC3D, double Xi, double *kBloch, double *EFT)
{ 
  if (SC3D->UseExistingData==false)
   return false;

  if (    (kBloch==0 && SC3D->G->NumLatticeBasisVectors>0)
       || (kBloch!=0 && SC3D->G->NumLatticeBasisVectors==0)
     ) ErrExit("%s:%i: internal error",__FILE__,__LINE__);

  /*----------------------------------------------------------*/
  /*----------------------------------------------------------*/
  /*----------------------------------------------------------*/
  int NumKeys = 1 + SC3D->G->NumLatticeBasisVectors;
  double Keys[3];
  Keys[0]=Xi;
  if(NumKeys>=2) Keys[1]=kBloch[0];
  if(NumKeys>=3) Keys[2]=kBloch[1];

  /*----------------------------------------------------------*/
  /* 0. try to open the cache file. --------------------------*/
  /*----------------------------------------------------------*/
  FILE *f;
  if (kBloch)
   f=fopen(SC3D->ByXiKFileName,"r");
  else
   f=fopen(SC3D->ByXiFileName,"r");
  if (f==0) return false;

  /*----------------------------------------------------------*/
  /* 1. skip down through the cache file until we find a line */
  /*    whose...                                              */
  /*----------------------------------------------------------*/
  double Values[7];
  int NumValues;
  int LineNum=0;
  char Line[1000];
  char *Tag = SC3D->GTCList[0]->Tag;
  bool FoundFirst=false;
  int NQ = SC3D->NumQuantities;
  while( !FoundFirst && fgets(Line,1000,f) )
   { LineNum++;
     if ( !LineMatches(Line,Tag,Keys,NumKeys,Values,&NumValues) )
      continue;
     if ( NumValues != NQ )
      { Log(" found matching (Tag,freqs) on line %i of cache file"
            " but number of quantities is wrong (%i, %i)",
            LineNum,NumValues,NQ);
        continue;
      };
     FoundFirst=true;
   };
  if (!FoundFirst)
   { fclose(f);
     return false;
   };

  Log("...found matching dataset on line %i...",LineNum);
  memcpy(EFT,Values,NQ*sizeof(double));
   
  /*----------------------------------------------------------*/
  /* 3. ok, since that worked, now keep going ----------------*/
  /*----------------------------------------------------------*/
  for(int nt=1; nt<SC3D->NumTransformations; nt++)
   { 
     LineNum++;
     Tag = SC3D->GTCList[nt]->Tag;

     if (    !fgets(Line,1000,f)  
          || !LineMatches(Line,Tag,Keys,NumKeys,Values,&NumValues) 
          || (NumValues!=NQ)
        )
      { Log(" ...but data for some transforms were missing (line %i) (aborting)",LineNum);
        fclose(f);
        return false;
      };

     memcpy(EFT + nt*NQ,Values,NQ*sizeof(double));
   };
   
  fclose(f);
  Log(" ...and successfully read data for all quantities at all transforms");
  return true;

}

/***************************************************************/
/* wrapper around GetCasimirIntegrand with correct prototype   */
/* prototype for adapt_integrate()                             */
/***************************************************************/
int kBlochIntegrand(unsigned ndim, const double *x, void *params,
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
/* to the Casimir quantities. For non-periodic geometries, this*/
/* involves making a single call to GetCasimirIntegrand. For   */
/* periodic geometries, this means integrating over the        */
/* Brillouin zone at the Xi value in question.                 */
/***************************************************************/
void GetXiIntegrand(SC3Data *SC3D, double Xi, double *EFT) 
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (SC3D->G->NumLatticeBasisVectors==0)
   { 
     GetCasimirIntegrand(SC3D, Xi, 0, EFT);
   }
  else
   { 
     /************************************************************/
     /* perform brillouin-zone integration.                      */
     /*                                                          */
     /* \int_{BZ} f( \vec k ) d^2 k                              */
     /*  = J*\int_0^1 du \int_0^1 dv f(\vec k(u,v)) du dv        */
     /*  =4J*\int_0^{1/2} du \int_0^{1/2} dv f(\vec k(u,v)) du dv*/
     /*                                                          */
     /* where \vec k(u,v) = u*\Gamma_1 + v*\Gamma_2              */
     /* Jacobian: d^2 k = J du dv                                */
     /* where J = (det Gamma) = V_{BZ}.                          */
     /* and thus the integral reads                              */
     /***************************************************************/
     double Lower[2] = {0.0, 0.0};
     double Upper[2] = {0.5, 0.5};
     double *Error = new double[SC3D->NTNQ];
     if (SC3D->BZICutoff!=0.0)
      { Upper[0]*=SC3D->BZICutoff;
        Upper[1]*=SC3D->BZICutoff;
      };
     
     SC3D->Xi = Xi;
     if (SC3D->BZIOrder == 0)
      { 
        pcubature(SC3D->NTNQ, kBlochIntegrand, (void *)SC3D, 2, Lower, Upper,
                  SC3D->MaxkBlochPoints, SC3D->AbsTol, SC3D->RelTol, ERROR_INDIVIDUAL,
                  EFT, Error);
      }
     else
      { 
        ECC2D(SC3D->BZIOrder, Lower, Upper, kBlochIntegrand, (void *)SC3D,
              SC3D->NTNQ, SC3D->BZSymmetry, SC3D->BZIValues, 0,
              EFT, Error);
      };

     for(int ntnq=0; ntnq<SC3D->NTNQ; ntnq++)
      { 
        EFT[ntnq] *= 4.0*SC3D->BZVolume;
        Error[ntnq] *= 4.0*SC3D->BZVolume;

        if ( Error[ntnq] > 10.0*SC3D->RelTol*fabs(EFT[ntnq]) )
         Warn("potentially large errors (Q%i: %i %%) in BZ integration",
               ntnq,ceil(100.0*Error[ntnq]/fabs(EFT[ntnq])));
      };
     
     /***************************************************************/
     /* write data to .byXi file                                    */
     /***************************************************************/
     FILE *f=fopen(SC3D->ByXiFileName,"a");
     for(int ntnq=0, nt=0; nt<SC3D->NumTransformations; nt++)
      { fprintf(f,"%s %.6e ",SC3D->GTCList[nt]->Tag,Xi);
        for(int nq=0; nq<SC3D->NumQuantities; nq++, ntnq++) 
         fprintf(f,"%.8e %.8e ",EFT[ntnq],Error[ntnq]);
        fprintf(f,"\n");
        fflush(f);
     };
     fclose(f);

     delete[] Error;
   };

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
     /*  2\pi kT *  \sum_n^\prime FI(\xi_n)                       */
     /*                                                             */
     /* where FI is what is returned by GetXiIntegrand.             */
     /***************************************************************/
     memcpy(LastEFT,EFT,NTNQ*sizeof(double));
     for(ntnq=0; ntnq<NTNQ; ntnq++)
      EFT[ntnq] += Weight * 2.0*M_PI* kT * dEFT[ntnq];

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
