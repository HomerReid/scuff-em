/*
 * GetCPIntegrand.cc -- routine to compute the contribution of a single 
 *                   -- imaginary frequency to the casimir-polder potential 
 *
 * homer reid        -- 10/2006 -- 2/2012
 *
 */

#include <complex>

#include <libSGJC.h>
#include <libIncField.h>

#include "scuff-caspol.h"

// prefactor that enters into the calculation of the 
// casimir-polder potential.
// below, we obtain the casimir-polder potential as the quantity
//  PF * \int d\xi * \xi^2 * Tr(\alpha * G) 
// where 
//  -- \xi has units of 3e14 rad/sec 
//  -- \alpha has units of a0^3 (bohr radius)
//  -- G has units of 1/um 
// thus if we put 
//  PF = (\hbar c / 1um) * (a0 / 1um)^3
//     = (0.1973 ev) * (0.529177e-4)^3
// then the result of our calculation will be an 
// energy in units of ev.                         
#define PREFAC 2.9237e-14


// factor that converts temperature from degrees
// kelvin to my internal energy units 
// (in which '1' == 0.1973 ev.)
// how this works: 
//  a. temperature in eV = kT = 8.6173e-5 * (T in kelvin)
//  b. temperature in my internal energy units
//     = (kT in eV) / (0.1973 eV)
//     = (8.6173e-5 / 0.1973 ) * (T in Kelvin)
#define BOLTZMANNK 4.36763e-4

/***************************************************************/
/* compute the dyadic green's function at a distance Z above a */
/* PEC plate using the method of images.                       */
/***************************************************************/
void GetPECPlateDGF(double Z, double Xi, cdouble GE[3][3])
{
  // construct a point source at the image location
  double X0[3]={0.0, 0.0, -Z};
  cdouble P0[3] = {0,0,0}; 
  PointSourceData MyPSD(X0, P0);
  MyPSD.Omega=cdouble(0,Xi);
  MyPSD.Eps=1.0;
  MyPSD.Mu=1.0;

  // get each column of the DGF
  cdouble EH[6];
  double X[3]={0.0, 0.0, Z}; 

  MyPSD.P[0] = -1.0;   MyPSD.P[1] =  0.0; MyPSD.P[2] = 0.0;
  MyPSD.GetFields(X, EH);
  GE[0][0]=EH[0]; GE[1][0]=EH[1]; GE[2][0]=EH[2];

  MyPSD.P[0] =  0.0;   MyPSD.P[1] = -1.0; MyPSD.P[2] = 0.0;
  MyPSD.GetFields(X, EH);
  GE[0][1]=EH[0]; GE[1][1]=EH[1]; GE[2][1]=EH[2];
  
  MyPSD.P[0] =  0.0;   MyPSD.P[1] =  0.0; MyPSD.P[2] = 1.0;
  MyPSD.GetFields(X, EH);
  GE[0][2]=EH[0]; GE[1][2]=EH[1]; GE[2][2]=EH[2];

  // normalization factor needed to convert from 
  // my normalization of the DGF, which has units
  // of electric field / surface current, to the 
  // usual normalization in which the DGF has units 
  // of inverse length 
  double Factor = -ZVAC*Xi;
  int i, j;
  for(i=0; i<3; i++)
   for(j=0; j<3; j++)
    GE[i][j] /= Factor;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetCPIntegrand(SCPData *SCPD, double Xi, double *U)
{
  RWGGeometry *G     = SCPD->G;
  HMatrix *M         = SCPD->M;
  HVector *KN        = SCPD->KN;
  PolModel *PM       = SCPD->PM;
  HMatrix *EPList    = SCPD->EPList;
  int nThread        = SCPD->nThread;
  FILE *ByXiFile     = SCPD->ByXiFile;

  /***************************************************************/ 
  /* assemble and factorize the BEM matrix at this frequency     */ 
  /* (unless we are doing the PEC plate case...)                 */ 
  /***************************************************************/
  if (G)
   { Log("Assembling BEM matrix at Xi=%g\n",Xi);
     G->AssembleBEMMatrix(cdouble(0,Xi), M, nThread);
     M->LUFactorize();
   };

  /***************************************************************/ 
  /* get polarizability values at this frequency                 */ 
  /***************************************************************/ 
  double Alpha[9];
  memset(Alpha, 0, 9*sizeof(double));
  PM->GetPolarizability(Xi, Alpha);

  /***************************************************************/ 
  /* loop over all evaluation points to get the contribution  of */ 
  /* this frequency to the CP potential at each point            */
  /***************************************************************/ 
  int i, j, nep;
  double R[3];
  cdouble GE[3][3], GM[3][3];
  for(nep=0; nep<EPList->NR; nep++)
   { 
      /* get the dyadic GF at this eval point */
      R[0]=EPList->GetEntryD(nep, 0);
      R[1]=EPList->GetEntryD(nep, 1);
      R[2]=EPList->GetEntryD(nep, 2);
      if (G) {
	fprintf(stderr, "error: GetGij not implemented\n");
	exit(1);
	// G->GetGij(R, M, 0, KN, Xi, SCUFF_PUREIMAGFREQ, nThread, GE, GM);
      }
      else
       GetPECPlateDGF(R[2], Xi, GE);

      /* compute the trace of Alpha \cdot G */
      for(U[nep]=0.0, i=0; i<3; i++)
       for(j=0; j<3; j++)
        U[nep] += PREFAC * Xi * Xi * Alpha[i+3*j] * real(GE[j][i]);
      
      fprintf(SCPD->ByXiFile,"%e %e %e %e %e\n",R[0],R[1],R[2],Xi,U[nep]);
 
   }; 
}

/***************************************************************/
/* evaluate the matsubara sum at Temperature kelvin.           */
/***************************************************************/
void EvaluateMatsubaraSum(SCPData *SCPD, double Temperature, double *U)
{ 
  int nXi, NEP=SCPD->EPList->NR;  // 'number of evaluation points' 
  double Xi, Weight, Delta;
  double *dU = new double[NEP], *LastU = new double[NEP]; 
  int *ConvergedIters = new int[NEP], AllConverged;
  double kT=Temperature*BOLTZMANNK;

  memset(U,0,NEP*sizeof(double));
  memset(ConvergedIters,0,NEP*sizeof(int));

  Log("Beginning Matsubara sum at T=%g kelvin...",Temperature);

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
     /* accumulate contributions to the matsubara sum,              */
     /* how it works: the matsubara sum is                          */
     /*  2\pi kT *  \sum_n^\prime F(\xi_n)                          */
     /* where \xi_n is the nth matsubara frequency and F(\xi) is    */
     /* the frequency integrand, and where the primed sum means     */
     /* he n==0 term is weighted                                    */
     /* with a factor of 1/2.                                       */
     /***************************************************************/
     memcpy(LastU,U,NEP*sizeof(double));
     int nep;
     for(nep=0; nep<NEP; nep++)
      U[nep] += Weight * 2.0*M_PI*kT * dU[nep];

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
        if ( Delta < SCPD->AbsTol || Delta < SCPD->RelTol*fabs(U[nep]) )
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

  delete[] ConvergedIters;
  delete[] LastU;
  delete[] dU;
} 

/***************************************************************/
/* integrand routine used to evaluate imaginary frequency      */
/* integral by adaptive quadrature                             */
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
  double *Error = new double[fdim];

  adapt_integrate_log(fdim, SGJCIntegrand, (void *)SCPD, 1, 
                      &Lower, &Upper, 0, SCPD->AbsTol, SCPD->RelTol,
                      U, Error, "scuff-caspol.cubaturelog", 15);

  delete[] Error;
}
