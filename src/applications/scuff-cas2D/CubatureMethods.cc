#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <libhrutil.h>
#include <libSGJC.h>

#include "libTDRT.h"
#include "scuff-cas2D.h"

/***************************************************************/
/***************************************************************/
/***************************************************************/
int XiIntegrand(unsigned ndim, const double *x, void *params, 
                unsigned fdim, double *fval)
{ 
  (void) ndim; // unused 

  C2DWorkspace *W=(C2DWorkspace *)params;
  double Xi, Jacobian;
  unsigned ntnq;

  if (x[0]==1.0)
   { memset(fval, 0, fdim*sizeof(double));
     return 0; 
   };
   
  Xi = W->XQMin + x[0] / (1.0-x[0]);
  Jacobian=1.0 / ((1.0-x[0])*(1.0-x[0]));
  
  XQIntegrand(W, Xi, W->FixedQ, fval);
 
  for(ntnq=0; ntnq<fdim; ntnq++)
   fval[ntnq]*=Jacobian;

  return 0;

}  

/***************************************************************/
/***************************************************************/
/***************************************************************/
void EvaluateXiIntegral(C2DWorkspace *W, double Q, double *I, double *E)
{ 
  int ntnq;
  double *I1 = new double [W->NTNQ]; 
  double *I2 = new double [W->NTNQ];

  /*--------------------------------------------------------------*/
  /*- estimate the integral over the range Xi=[0,XQMIN] by       -*/
  /*- assuming that the integrand is constant over that range    -*/
  /*--------------------------------------------------------------*/
  XQIntegrand(W, W->XQMin, Q, I1);
  for(ntnq=0; ntnq<W->NTNQ; ntnq++)
   I1[ntnq]*=(W->XQMin);

  /*--------------------------------------------------------------*/
  /*- now use quadrature to integrate over the range             -*/
  /*- Xi=[XQMIN, Infinity].                                      -*/
  /*--------------------------------------------------------------*/
  double Lower=0.0;
  double Upper=1.0;
  W->FixedQ=Q;
  pcubature(W->NTNQ, XiIntegrand, (void *)W, 1, &Lower, &Upper,
            0, W->AbsTol, W->RelTol, ERROR_INDIVIDUAL, I2, E);

  for(ntnq=0; ntnq<W->NTNQ; ntnq++)
   I[ntnq] = I1[ntnq] + I2[ntnq];

  delete[] I1;
  delete[] I2;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int QIntegrand(unsigned ndim, const double *x, void *params, 
                unsigned fdim, double *fval)
{ 
  (void) ndim; // unsigned 
  C2DWorkspace *W=(C2DWorkspace *)params;
  double Q, Jacobian;
  unsigned ntnq;

  if (x[0]==1.0)
   { memset(fval, 0, fdim*sizeof(double));
     return 0;
   };
   
  Q=x[0] / (1.0-x[0]);
  Jacobian=1.0 / ((1.0-x[0])*(1.0-x[0]));
  
  XQIntegrand(W, W->FixedXi, Q, fval);
 
  for(ntnq=0; ntnq<fdim; ntnq++)
   fval[ntnq]*=Jacobian;

  return 0;

}  

/***************************************************************/
/***************************************************************/
/***************************************************************/
void EvaluateQIntegral(C2DWorkspace *W, double Xi, double *EF, double *Error)
{ 
  /***************************************************************/
  /* first look through the .byXi file to see if the Q integral  */
  /* has already been evaluated at this value of Xi.             */
  /***************************************************************/
  if ( CacheRead(W, Xi, -1.0, EF) )
   return;

  /***************************************************************/
  /* otherwise, evaluate the Q integral **************************/
  /***************************************************************/
  if (TDRTGeometry::LogLevel>=1)
   Log("Evaluating Q integral at Xi=%g...",Xi);

  double Lower=0.0;
  double Upper=1.0;
  W->FixedXi=Xi;
  pcubature_log(W->NTNQ, QIntegrand, (void *)W, 1, &Lower, &Upper,
                0, W->AbsTol, W->RelTol, ERROR_INDIVIDUAL, 
                EF, Error,"QIntegral.log");

  /***************************************************************/
  /* and write the results to the .byXi file                     */
  /***************************************************************/
  char DateString[200];
  time_t MyTime;
  struct tm *MyTm;
  MyTime=time(0);
  MyTm=localtime(&MyTime);
  strftime(DateString,30,"%D::%T",MyTm);

  FILE *f=fopen(W->ByXiFileName,"a");
  fprintf(f,"# %s\n",DateString);
  if (!f) return;
  int nt, nq, ntnq;
  for(ntnq=nt=0; nt<W->NumTransforms; nt++)
   { fprintf(f,"%s %.15e  ",W->Tags[nt],Xi);
     for(nq=0; nq<W->NumQuantities; nq++, ntnq++)
      fprintf(f,"%.15e ",EF[ntnq]);
     ntnq-=W->NumQuantities;
     for(nq=0; nq<W->NumQuantities; nq++, ntnq++)
      fprintf(f,"%.2e ",Error[ntnq]);
     fprintf(f,"\n");
   };
  fclose(f);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int RXQIntegrand(unsigned ndim, const double *x, void *params, 
                 unsigned fdim, double *fval)
{
  (void )ndim; //unused 

  if (x[0]==1.0)
   { memset(fval, 0, fdim*sizeof(double));
     return 0;
   };
   
  C2DWorkspace *W = (C2DWorkspace *)params;

  double RXQ = W->XQMin + x[0] / (1.0 - x[0]);
  double Jacobian= 0.5 * M_PI * RXQ / ((1.0-x[0])*(1.0-x[0]));

  XQIntegrand(W, RXQ, 0, fval);

  for(unsigned ntnq=0; ntnq<fdim; ntnq++)
   fval[ntnq]*=Jacobian;

  return 0;
}

int XQIntegrand(unsigned ndim, const double *x, void *params, 
                unsigned fdim, double *fval)
{
  (void )ndim; //unused 
  double Xi, q, Jacobian;
  unsigned ntnq;

  if ( x[0]==1.0 || x[1]==1.0 )
   { memset(fval, 0, fdim*sizeof(double));
     return 0;
   };
   
  C2DWorkspace *W = (C2DWorkspace *)params;
  Xi = W->XQMin + x[0] / (1.0 - x[0]);
   q = x[1] / (1.0 - x[1]);
  Jacobian= 1.0 / ((1.0-x[0])*(1.0-x[0])*(1.0-x[1])*(1.0-x[1]));

  XQIntegrand(W, Xi, q, fval);

  for(ntnq=0; ntnq<fdim; ntnq++)
   fval[ntnq]*=Jacobian;

  return 0;
}

void EvaluateXQIntegral(C2DWorkspace *W, double *I, double *E)
{
  if (W->G->AllPEC)
   { 
     /***************************************************************/
     /* assume that the integrand is constant for Xi < XQMIN        */
     /***************************************************************/
     double Xi = W->XQMin;
     RXQIntegrand(1, &Xi, (void *)W, W->NTNQ, I);
     for(int ntnq=0; ntnq<W->NTNQ; ntnq++)
      I[ntnq] *= W->XQMin;

     /***************************************************************/
     /* now get the integral from XQMIN to infinity *****************/
     /***************************************************************/
     double Lower[1] = {0.0};
     double Upper[1] = {1.0};
     double *DeltaI = new double [W->NTNQ];
     pcubature_log(W->NTNQ, RXQIntegrand, (void *)W, 1, Lower, Upper, 
                   0,W->AbsTol,W->RelTol,ERROR_INDIVIDUAL,
                   DeltaI,E,"RXQIntegral.log");

     for(int ntnq=0; ntnq<W->NTNQ; ntnq++)
      I[ntnq] += DeltaI[ntnq];

     free(DeltaI);
   }
  else
   { 
     double Lower[2] = {0.0, 0.0}; 
     double Upper[2] = {1.0, 1.0};
     pcubature_log(W->NTNQ, XQIntegrand, (void *)W, 2, Lower, Upper, 
                   0,W->AbsTol,W->RelTol,ERROR_INDIVIDUAL,
                   I,E,"XQIntegral.log");
   };

}

/***************************************************************/
/* evaluate the matsubara sum at a given temperature           */
/***************************************************************/
void EvaluateMatsubaraSum(C2DWorkspace *W, double T, double *I, double *E)
{
  int n, ntnq, NTNQ=W->NTNQ;
  double Xi, Weight;
  double *dI = new double[NTNQ];
  double *LastI = new double[NTNQ];
  int *ConvergedIters = new int[NTNQ];
  double AbsDelta, RelDelta;
  int NumConverged;

  // how this works: 
  //  a. temperature in eV = kT = 8.6173e-5 * (T in kelvin)
  //  b. temperature in our internal energy units
  //     = (kT in eV) / (0.1973 eV)
  //     = (8.6173e-5 / 0.1973 ) * (T in Kelvin)
  double kT = 4.36763e-4 * T;

  /***************************************************************/
  /* sum over matsubara frequencies ******************************/
  /***************************************************************/
  if (TDRTGeometry::LogLevel>=1)
   Log("Beginning Matsubara sum at T=%g kelvin...",T);
  memset(I,0,NTNQ*sizeof(double));
  memset(W->Converged,0,NTNQ*sizeof(int));
  memset(ConvergedIters,0,NTNQ*sizeof(int));
  NumConverged=0;
  for(n=0; n<10000; n++)
   {
     /***************************************************************/
     /* compute the next matsubara frequency ************************/
     /***************************************************************/
     if (n==0)
      { Weight=0.5;
        // NOTE: we assume that the integrand is constant for Xi < XQMIN
        Xi = W->XQMin;
      }
     else
      { Weight=1.0;
        Xi=2.0*M_PI*kT*((double)n);
      };

     /***************************************************************/
     /* evaluate the Q integral with Xi set equal to this matsubara */
     /* frequency                                                   */
     /***************************************************************/
     EvaluateQIntegral(W,Xi,dI,E);

     /***************************************************************/
     /* accumulate contributions to the sums. ***********************/
     /* how the normalization works:                                */
     /*  (a) the sum we want to compute is                          */
     /*       2*pi*kT * \sum^\prime_\xi (1/pi) * \int_0 I(xi,q)     */
     /***************************************************************/
     memcpy(LastI,I,NTNQ*sizeof(double));
     for(ntnq=0; ntnq<NTNQ; ntnq++)
      I[ntnq] += Weight * 2.0*M_PI * kT * dI[ntnq];

     /*********************************************************************/
     /* convergence analysis.                                             */
     /* how it works: if the relative change in output quantity #ntnq is  */
     /* less than EPSREL, we increment ConvergedIters[ntnq]; otherwise we */
     /* set ConvergedIters[ntnq] to 0.                                    */
     /* when ConvergedIters[ntnq] hits 3, we mark quantity #ntnq as       */
     /* having converged.                                                 */
     /* when all quantities have converged, we are done.                  */
     /*********************************************************************/
     if (n>2)
      { 
        for(ntnq=0; ntnq<NTNQ; ntnq++)
         { 
           if (W->Converged[ntnq]) 
            continue;

           AbsDelta = fabs( I[ntnq]-LastI[ntnq] );
           RelDelta = AbsDelta / fabs(I[ntnq]);
           if ( RelDelta < 1.0e-4*W->RelTol ) // || AbsDelta < W->AbsTol )
            { ConvergedIters[ntnq]++;
              if (ConvergedIters[ntnq]>2)
               { W->Converged[ntnq]=1;
                 if (TDRTGeometry::LogLevel>=1)
                  Log("Quantity %i converged...",ntnq); 
                 NumConverged++;
               };
            }
           else
            ConvergedIters[ntnq]=0;
         }; 
      };
     if (NumConverged==W->NTNQ)
      break;

   }; // for(n=0; n<10000; n++)
  
  /***************************************************************/ 
  /***************************************************************/ 
  /***************************************************************/ 
  if (NumConverged<W->NTNQ)
   { 
     fprintf(stderr,"\n*\n* WARNING: Matsubara sum unconverged after %i frequency samples.\n*\n",n);
     if (TDRTGeometry::LogLevel>=1)
      Log("Matsubara sum UNCONVERGED at n=%i samples",n); 
   } 
  else 
   if (TDRTGeometry::LogLevel>=1)
    Log("Matsubara sum converged after summing n=%i frequency points.",n);

  /***************************************************************/ 
  /* return the error as the difference between the most recent  */ 
  /* two iterations                                              */ 
  /***************************************************************/ 
  for(ntnq=0; ntnq<NTNQ; ntnq++)
   E[ntnq]=fabs(I[ntnq] - LastI[ntnq]);

  delete[] dI;
  delete[] LastI;
  delete[] ConvergedIters;
    
} 

/***************************************************************/
/* read (xi,q) points from the specified file and              */
/* evaluate the integrand at each point                        */
/***************************************************************/
void ProcessXQList(C2DWorkspace *W, char *XQListName)
{
  double *I = new double[W->NTNQ];
  char buffer[200], *p;
  double Xi, Q;
  int LineNum, nConv, nRead;
  FILE *f;
  
  f=fopen(XQListName,"r"); /* we have already checked that this works */

  LineNum=0;
  while( fgets(buffer,1000,f) )
   { 
      LineNum++;

      /* skip blank lines and comments */
      p=buffer;
      while( isspace(*p) ) p++;
      if ( *p==0 || *p=='#' )
       continue;

      /* try to interpret line as two numbers */
      nConv=sscanf(p,"%le %le %n\n",&Xi,&Q,&nRead);
      if ( nConv!=2 )
       { fprintf(stderr,"%s:%i: syntax error (skipping)\n",XQListName,LineNum);
         continue;
       };
      p+=nRead;
      while( *p && isspace(*p) ) p++;
      if ( *p!=0 )
       { fprintf(stderr,"%s:%i: extra characters on line (skipping)\n",XQListName,LineNum);
         continue;
       };

      /* evaluate the integrand at the given point and print results on console */
      XQIntegrand( W, Xi, Q, I);

      printf("\n ** At (Xi,Q)=(%15.12e,%15.12e):\n",Xi,Q);
      PrintConsoleOutput(W,I);
      printf("\n");

   }; // while( fgets(buffer,1000,f) )

  fclose(f);
  delete[] I;

}
