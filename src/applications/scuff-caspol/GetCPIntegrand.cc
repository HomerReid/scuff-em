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
void GetCPIntegrand(SCPData *SCPD, double Xi, HMatrix *EPList, double *U)
{
  RWGGeometry *G = SCPD->G;
  HMatrix *M     = SCPD->M;
  HVector *KN    = SCPD->KN;
  int nThread    = SCPD->nThread;

  /***************************************************************/ 
  /***************************************************************/ 
  /***************************************************************/ 
  G->AssembleBEMMatrix(Xi, IMAG_FREQ, nThread, M);
  M->LUFactorize();

  /***************************************************************/ 
  /* get polarizability values at this frequency                 */ 
  /***************************************************************/ 
  double Alpha[9];
  memset(Alpha, 0, 9*sizeof(double));

  int i, j;
  for(i=0; i<3; i++)
   for(j=0; j<3; j++)
    if (AlphaMP[i+3*j])
     Alpha[ i + 3*j] = AlphaMP[ i + 3*j]->GetEpsD(Xi, IMAG_FREQ);

  /***************************************************************/ 
  /***************************************************************/ 
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
 
   };

}
