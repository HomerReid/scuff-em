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
 * GetCPIntegrand.cc -- routine to compute the contribution of a single 
 *                   -- imaginary frequency to the casimir-polder potential 
 *
 * homer reid        -- 10/2006 -- 2/2012
 *
 */

#include <complex>

#include <libSGJC.h>
#include <libIncField.h>
#include <libTriInt.h>

#include "scuff-caspol.h"

// prefactor that enters into the calculation of the 
// casimir-polder potential.
// The CP potential is 
//  (2\hbar) \int d\xi \xi^2 Tr(\alpha*G)
// where 
//  -- \xi has units of 3e14 rad/sec 
//  -- \alpha has units of a0^3 (bohr radius)
//  -- G has units of 1/um.
// thus if we put 
//  PF = (2\hbar * 3e14 s^{-1}) * (a0 / 1um)^3
//     = (2*6.582119e-16 ev s) * (3e14 s^{-1} * (0.529177e-4)^3
//     = 
// then the result of our calculation will be an 
// energy in units of ev.
// #define PREFAC 5.852209e-14
// update 20130825: let's report potentials in nanoelectronvolts,
// so we'll up that prefactor by 1e-9
#define PREFAC 5.852209e-5

// factor that converts temperature from degrees
// kelvin to my internal energy units 
// (in which '1' == 0.1973 ev.)
// how this works: 
//  a. temperature in eV = kT = 8.6173e-5 * (T in kelvin)
//  b. temperature in my internal energy units
//     = (kT in eV) / (0.1973 eV)
//     = (8.6173e-5 / 0.1973 ) * (T in Kelvin)
#define BOLTZMANNK 4.36763e-4

#define ABSTOL 0.0
#define XIMIN  1.0e-6

/***************************************************************/
/* compute the dyadic green's function at a distance Z above a */
/* PEC plate using the method of images.                       */
/***************************************************************/
void GetPECPlateDGF(double Z, double Xi, cdouble GE[3][3])
{
  // construct a point source at the image location
  double X0[3]={0.0, 0.0, -Z};
  cdouble P0[3] = {0,0,0}; 
  PointSource MyPSD(X0, P0);
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
  double Factor = -Xi*Xi;
  int i, j;
  for(i=0; i<3; i++)
   for(j=0; j<3; j++)
    GE[i][j] /= Factor;
}

/***************************************************************/
/* the first parameter here needs to be void * so the routine  */
/* can be passed as the integrand function to GetBZIntegral()  */
/***************************************************************/
void GetCPIntegrand(void *pSCPD, cdouble Omega, double *kBloch, double *U)
{
  SCPData *SCPD        = (SCPData *)pSCPD;

  RWGGeometry *G       = SCPD->G;
  HMatrix *M           = SCPD->M;
  HVector *KN          = SCPD->KN;
  int NumAtoms         = SCPD->NumAtoms;
  HMatrix **Alphas     = SCPD->Alphas;   
  HMatrix *EPMatrix    = SCPD->EPMatrix;
  void **ABMBCache     = SCPD->ABMBCache;

  if (imag(Omega) < XIMIN)
   Omega = cdouble( real(Omega), XIMIN );
  double Xi=imag(Omega);

  /***************************************************************/
  /* assemble and factorize the BEM matrix at this               */
  /* (frequency, kBloch) point, unless we are doing the PEC      */
  /* plate case.                                                 */
  /***************************************************************/
  if (G)
   { 
     if (G->LDim==0)
      { Log("Assembling BEM matrix at Xi=%g",Xi);
        G->AssembleBEMMatrix(Omega, M);
      }
     else
      { if (G->LDim==1)
         Log("Assembling BEM matrix at (Xi,kx)=(%g,%g)",Xi,kBloch[0]);
        else
         Log("Assembling BEM matrix at (Xi,kx,ky)=(%g,%g,%g)",Xi,kBloch[0],kBloch[1]);
        int NS = G->NumSurfaces;
        for(int ns=0, nb=0; ns<NS; ns++)
         for(int nsp=ns; nsp<NS; nsp++, nb++)
          { 
            int RowOffset = G->BFIndexOffset[ns];
            int ColOffset = G->BFIndexOffset[nsp];
            G->AssembleBEMMatrixBlock(ns, nsp, Omega, kBloch,
                                      M, 0, RowOffset, ColOffset,
                                      ABMBCache[nb], false);

            if (nsp>ns)
             G->AssembleBEMMatrixBlock(nsp, ns, Omega, kBloch,
                                       M, 0, ColOffset, RowOffset,
                                       ABMBCache[nb], false);
          };
      };
     Log("LU-factorizing...");
     M->LUFactorize();
   };

  /***************************************************************/ 
  /* loop over all evaluation points to get the contribution of  */ 
  /* this frequency to the CP potential at each point            */
  /***************************************************************/ 
  double R[3];
  cdouble GE[3][3], GM[3][3];
  Log("Computing CP potential at %i eval points...",EPMatrix->NR);
  FILE *f = kBloch ? fopen(SCPD->ByXikFileName, "a") : 0;
  for(int nep=0; nep<EPMatrix->NR; nep++)
   { 
      /* get the dyadic GF at this eval point */
      R[0]=EPMatrix->GetEntryD(nep, 0);
      R[1]=EPMatrix->GetEntryD(nep, 1);
      R[2]=EPMatrix->GetEntryD(nep, 2);

      if (G)
       G->GetDyadicGFs(R, Omega, kBloch, M, KN, GE, GM);
      else
       GetPECPlateDGF(R[2], Xi, GE);

      if (f) fprintf(f,"%e %e %e %e %e ",R[0],R[1],R[2],Xi,kBloch[0]);
      if (G && G->LDim==2) fprintf(f,"%e ",kBloch[1]);
      for(int na=0; na<NumAtoms; na++)
       { 
         HMatrix *Alpha = Alphas[na];
         double UValue=0.0;
         for(int i=0; i<3; i++)
          for(int j=0; j<3; j++)
           UValue += PREFAC * Xi * Xi * Alpha->GetEntryD(i,j) * real(GE[j][i]);

         U[nep*NumAtoms + na] = UValue;

         if (f) fprintf(f,"%e %e ",Alpha->GetEntryD(0,0),UValue);
       };
      if (f) fprintf(f,"\n");
   }; 
  if (f) fclose(f);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetXiIntegrand(SCPData *SCPD, double Xi, double *U)
{ 
  RWGGeometry *G       = SCPD->G;
  int NumAtoms         = SCPD->NumAtoms;
  PolModel **PolModels = SCPD->PolModels;
  HMatrix **Alphas     = SCPD->Alphas;   
  HMatrix *EPMatrix    = SCPD->EPMatrix;

  if (Xi<=XIMIN)
   Xi=XIMIN;

  /***************************************************************/
  /* now the polarizability tensor for each atomic species must  */
  /* have been looked up by LoadPolarizability at this frequency */
  /***************************************************************/

  /***************************************************************/ 
  /***************************************************************/ 
  /***************************************************************/ 
  cdouble Omega=cdouble(0.0,Xi);
  int LDim = G ? G->LDim : 0;
  if (LDim==0)
   GetCPIntegrand((void *)SCPD, Omega, 0, U);
  else 
   GetBZIntegral(SCPD->BZIArgs, Omega, U);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  FILE *f=fopen(SCPD->ByXiFileName,"a");
  for(int nep=0; nep<EPMatrix->NR; nep++)
   { 
      double R[3];
      EPMatrix->GetEntriesD(nep, ":", R);
      fprintf(f,"%e %e %e %e ",R[0],R[1],R[2],Xi);
      for(int na=0; na<NumAtoms; na++)
       { fprintf(f,"%e ",Alphas[na]->GetEntryD(0,0));
         fprintf(f,"%e ",U[nep*NumAtoms + na]);
       };
      fprintf(f,"\n");
   }; 
  fclose(f);

}

/***************************************************************/
/* evaluate the matsubara sum at Temperature kelvin.           */
/***************************************************************/
void EvaluateMatsubaraSum(SCPData *SCPD, double Temperature, double *U)
{ 
  double kT=Temperature*BOLTZMANNK;

  int NEP=SCPD->EPMatrix->NR;  // 'number of evaluation points'
  int NA=SCPD->NumAtoms;
  int NU=NA*NEP;

  double *dU = new double[NU];
  double *LastU = new double[NU];
  double *MaxdU= new double[NU];
  memset(U,0,NU*sizeof(double));
  memset(MaxdU,0,NU*sizeof(double));

  int *ConvergedIters = new int[NU], AllConverged;
  memset(ConvergedIters,0,NU*sizeof(int));

  Log("Beginning Matsubara sum at T=%g kelvin...",Temperature);

  int nXi;
  for(nXi=0; nXi<100000; nXi++)
   { 
     /***************************************************************/
     /* compute the next matsubara frequency ************************/
     /***************************************************************/
     double Xi, Weight, Delta;
     if (nXi==0)
      { Weight=0.5;
        // NOTE: we assume that the integrand is constant for Xi < XIMIN
        Xi=XIMIN;
      }
     else
      { Weight=1.0;
        Xi=2.0*M_PI*kT*((double)nXi);
      };

     /***************************************************************/
     /* evaluate the frequency integrand at this matsubara frequency*/
     /***************************************************************/
     LoadPolarizability(SCPD, Xi);
     GetXiIntegrand(SCPD, Xi, dU);

     /***************************************************************/
     /* accumulate contributions to the matsubara sum,              */
     /* how it works: the matsubara sum is                          */
     /*  2\pi kT * \sum_n^\prime F(\xi_n)                           */
     /* where \xi_n is the nth matsubara frequency and F(\xi) is    */
     /* the frequency integrand, and where the primed sum means     */
     /* the n==0 term is weighted with a factor of 1/2.             */
     /***************************************************************/
     memcpy(LastU,U,NU*sizeof(double));
     for(int nu=0; nu<NU; nu++)
      U[nu] += Weight * 2.0*M_PI*kT * dU[nu];

     /*********************************************************************/
     /* convergence analysis.                                             */
     /* how it works: if the absolute or relative change in the potential */
     /* at point #nu is within our tolerances, we increment               */
     /* ConvergedIters[nu]; otherwise we set ConvergedIters[nu] to 0.    */
     /* when ConvergedIters[nep] hits 3, we mark quantity #nep as         */
     /* having converged.                                                 */
     /* when all quantities have converged, we are done.                  */
     /*********************************************************************/
     AllConverged=1;
     for(int nu=0; nu<NU; nu++)
      { 
        Delta = fabs( (U[nu]-LastU[nu]) );
        if ( Delta < ABSTOL || Delta < (SCPD->RelTol)*fabs(U[nu]) )
         ConvergedIters[nu]++;
        else
         ConvergedIters[nu]=0;
        
        // 20130827 second convergence test. for small temperatures, 
        // each incremental addition to the integral can be small
        // even though the integral is far from converged, whereupon 
        // the preceding test can detect false convergence.
        // to correct for this, we will additionally require that 
        // the increment to the integral be less than RELTOL*its
        // maximum value at any previous frequency.
        // This prevents false detection of convergence in 
        // stretches where we are making a large number of 
        // small constant contributions to the integral.
        if ( fabs(dU[nu]) > MaxdU[nu] )
         MaxdU[nu] = fabs(dU[nu]);
        if ( fabs(dU[nu]) > (SCPD->RelTol)*MaxdU[nu] )
         ConvergedIters[nu]=0;
   
        if (ConvergedIters[nu]<3) AllConverged=0;
      }; 
     if (AllConverged==1)
      break;

   }; /* for (nXi=0 ... */
  
  if (AllConverged==0)
   { 
     fprintf(stderr,"\n*\n* WARNING: Matsubara sum unconverged after %i frequency samples.\n*\n",nXi);
     Log("Matsubara sum UNCONVERGED at n=%i samples",nXi); 
   } 
  else 
   Log("Matsubara sum converged after summing n=%i frequency points.",nXi);

  delete[] ConvergedIters;
  delete[] MaxdU;
  delete[] LastU;
  delete[] dU;
} 

/***************************************************************/
/* integrand routine used to evaluate imaginary frequency      */
/* integral by adaptive quadrature                             */
/***************************************************************/
int SGJCIntegrand(unsigned ndim, const double *x, void *params,
                  unsigned fdim, double *fval)
{
  (void) ndim; // unused

  if (x[0]==1.0)
   { memset(fval, 0, fdim*sizeof(double));  
     return 0;
   };

  double Xi = x[0] / (1.0-x[0]);

  SCPData *SCPD = (SCPData *)params;
  LoadPolarizability(SCPD, Xi);
  GetXiIntegrand(SCPD, Xi, fval);

  unsigned nf;
  double J  = 1.0 / ((1.0-x[0])*(1.0-x[0])); // jacobian
  for(nf=0; nf<fdim; nf++)
   fval[nf]*=J;

  return 0;

}
 
/***************************************************************/
/***************************************************************/
/***************************************************************/
void EvaluateFrequencyIntegral(SCPData *SCPD, double *U)
{
  double Lower=0.0;
  double Upper=1.0;
  int fdim = (SCPD->EPMatrix->NR * SCPD->NumAtoms);
  double *Error = new double[fdim];

  pcubature_log(fdim, SGJCIntegrand, (void *)SCPD, 1, 
                &Lower, &Upper, 0, ABSTOL, SCPD->RelTol,
                ERROR_INDIVIDUAL, U, Error, "pcubature.log");

  delete[] Error;
}

void LoadPolarizability(SCPData *SCPD, double Xi)
{
  if (Xi<=XIMIN)
   Xi=XIMIN;
  PolModel **PolModels = SCPD->PolModels;
  HMatrix **Alphas     = SCPD->Alphas;
  for(int na=0; na<(SCPD->NumAtoms); na++)
  {
   PolModels[na]->GetPolarizability(Xi, Alphas[na]);
   printf("@Xi=%e, Alphas[%d][0,0] = %e\n", Xi, na, Alphas[na]->GetEntryD(0,0));
  }
}
