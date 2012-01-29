/*
 * GetCPIntegrand.cc -- routine to compute the contribution of a single 
 *                   -- imaginary frequency to the casimir-polder potential 
 *
 * homer reid        -- 10/2006 -- 2/2012
 *
 */

#include <stdio.h>
#include <math.h>
#include <stdarg.h>

#include <scuff-caspol.h>

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetCPIntegrand(SCPData *SCPD, double Xi, double *U)
{
  RWGGeometry *G     = SCPD->G;
  HMatrix *M         = SCPD->M;
  HVector *KN        = SCPD->KN;
  MatProp **AlphaMP  = SCPD->AlphaMP;
  HMatrix *EPList    = SCPD->EPList;
  int nThread        = SCPD->nThread;
  FILE *ByXiFile     = SCPD->ByXiFile;

  /***************************************************************/ 
  /* assemble and factorize the BEM matrix at this frequency     */ 
  /***************************************************************/
  Log("Assembling BEM matrix at Xi=%g\n",Xi);
  G->AssembleBEMMatrix(Xi, IMAG_FREQ, nThread, M);
  M->LUFactorize();

  /***************************************************************/ 
  /* get polarizability values at this frequency                 */ 
  /***************************************************************/ 
  double Alpha[9];
  memset(Alpha, 0, 9*sizeof(double));

  if (SCPD->AlphaFunc) 
   { 
     SCPD->AlphaFunc(Xi, Alpha);
   }
  else
   { int i, j;
     for(i=0; i<3; i++)
      for(j=0; j<3; j++)
       if (AlphaMP[i+3*j])
        Alpha[ i + 3*j] = AlphaMP[ i + 3*j]->GetEpsD(Xi, IMAG_FREQ);
   };

  /***************************************************************/ 
  /* explain me **************************************************/ 
  /***************************************************************/ 
  double PreFac = -Xi / (2.0*M_PI);

  /***************************************************************/ 
  /* loop over all evaluation points to get the contributione of */ 
  /* this frequency to the CP potential at each point            */
  /***************************************************************/ 
  int nep;
  double R[3];
  cdouble GE[3][3], GM[3][3];
  for(nep=0; nep<EPList->NR; nep++)
   { 
      /* get the dyadic GF at this eval point */
      R[0]=M->GetEntryD(nep, 0);
      R[1]=M->GetEntryD(nep, 1);
      R[2]=M->GetEntryD(nep, 2);
      G->GetGij(R, M, 0, KN, Xi, IMAG_FREQ, nThread, GE, GM);

      /* compute the trace of Alpha \cdot G */
      for(U[nep]=0.0, i=0; i<3; i++)
       for(j=0; j<3; j++)
        U[nep] += PreFac * Alpha[i][j] * real(GE[j][i]);
      
      fprintf(SCPD,"%e %e %e %e %e\n",R[0],R[1],R[2],Xi,U[nep]);
 
   }; 
}

/***************************************************************/
/* evaluate the matsubara sum at Temperature kelvin.           */
/***************************************************************/
// how this works: 
//  a. temperature in eV = kT = 8.6173e-5 * (T in kelvin)
//  b. temperature in our internal energy units
//     = (kT in eV) / (0.1973 eV)
//     = (8.6173e-5 / 0.1973 ) * (T in Kelvin)
#define BOLTZMANNK 4.36763e-4
void EvaluateMatsubaraSum(SCPData *SCPD, double Temperature, double *U)
{ 
  int nXi, NEP=SCPD->EPList->NR;  // 'number of evaluation points' 
  double Xi, Weight, Delta;
  double dU[NEP], LastU[NEP]; 
  int ConvergedIters[NEP], AllConverged;
  double kT=Temperature*BOLTZMANNK;

  memset(U,0,NEP*sizeof(double));
  memset(ConvergedIters,0,NEP*sizeof(int));

  Log("Beginning Matsubara sum at T=%g kelvin...",T);

  for(nXi=0; nXi<100000; nXi++)
   { 
     /***************************************************************/
     /* compute the next matsubara frequency ************************/
     /***************************************************************/
     if (nXi==0)
      { Weight=0.5;
        // NOTE: we assume that the integrand is constant for Xi < XIMIN
        Xi=SCPD->XiMin;
      }
     else
      { Weight=1.0;
        Xi=2.0*M_PI*kT*((double)nXi);
      };

     /***************************************************************/
     /* evaluate the frequency integrand at this matsubara frequency*/
     /***************************************************************/
     GetCPIntegrand(SCPD, Xi, dU);

     /***************************************************************/
     /* accumulate contributions to the sum.                        */
     /* how it works: the matsubara sum is                          */
     /*  2\pi kT *  \sum_n^\prime F(\xi_n)                          */
     /* where \xi_n is the nth matsubara frequency and F(\xi) is    */
     /* the frequency integrand.                                    */
     /* however, my FrequencyIntegrand() routine returns the        */
     /* quantity FI = F(\xi_n) / (2\pi). (this is so that the       */
     /* integral of FI over all \xi returns the correct casimir     */
     /* quantity with no additional multiplicative prefactors.)     */
     /* thus the matsubara sum is                                   */
     /*  4\pi^2 kT *  \sum_n^\prime F(\xi_n)                        */
     /* where the primed sum means that the n==0 term is weighted   */
     /* with a factor of 1/2.                                       */
     /***************************************************************/
     memcpy(LastU,U,NEP*sizeof(double));
     for(nep=0; nep<NEP; nep++)
      U[nep] += Weight * 4.0*M_PI*M_PI* kT * dU[nep];

     /*********************************************************************/
     /* convergence analysis.                                             */
     /* how it works: if the absolute or relative change in the potential */
     /* at point #nep is within our tolerances, we increment              */
     /* ConvergedIters[nep]; otherwise we set ConvergedIters[nep] to 0.   */
     /* when ConvergedIters[nep] hits 3, we mark quantity #nep as         */
     /* having converged.                                                 */
     /* when all quantities have converged, we are done.                  */
     /*********************************************************************/
     for(AllConverged=1, nep=0; nep<NEP; nep++)
      { 
        Delta = fabs( (U[nep]-LastU[nep]) );
        if ( Delta < AbsTol || Delta < RelTol*fabs(U[nep]) )
         ConvergedIters[nep]++;
        else
         ConvergedIters[nep]=0;
   
        if (ConvergedIters[nep]<3) AllConverged=0;
      }; 

     if (AllConverged==1)
      break;

   }; /* for (n=0 ... */
  
  if (AllConverged==0)
   { 
     fprintf(stderr,"\n*\n* WARNING: Matsubara sum unconverged after %i frequency samples.\n*\n",nXi);
     Log("Matsubara sum UNCONVERGED at n=%i samples",nXi); 
   } 
  else 
   Log("Matsubara sum converged after summing n=%i frequency points.",nXi);
    
} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SGJCIntegrand(unsigned ndim, const double *x, void *params,
                   unsigned fdim, double *fval)
{
  double Xi = x[0] / (1.0-x[0]);

  SCPData *SCPD = (SCPData *)params;
  GetCPIntegrand(SCPD, Xi, fval);

  int nf;
  double J  = 1.0 / ((1.0-x[0])*(1.0-x[0])); // jacobian
  for(nf=0; nf<fdim; nf++)
   fval[nf]*=J;
}

void EvaluateFrequencyIntegral(SCPData *SCPD, double *U)
{
  double Lower=0.0;
  double Upper=1.0;

  int fdim=SCPD->EPList->NR; // dimension of integrand vector 
  double Error[fdim];

  adapt_integrate_log(fdim, SGJCIntegrand, (void *)SCPD, 1, 
                      &Lower, &Upper, 0, SCPD->AbsTol, SCPD->RelTol,
                      U, Error, "scuff-caspol.cubaturelog", 15);
}
