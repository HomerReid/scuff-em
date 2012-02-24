/*****************************************************************/
/* the old (libRWG) versions of the In_EIKROverR, etc. functions */
/*****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex>

#include <gsl/gsl_sf.h>

#include <libhrutil.h>

#define II cdouble(0.0,1.0)

/***************************************************************/
/* complex version of the 'relative exponential' function,     */
/* which is e^x, minus the first n terms in the taylor series  */
/* for e^x, divided by the n+1th term in the taylor series for */
/* e^x (so that cExpRel2(n,0)==1 for all n.)                   */
/***************************************************************/
#define EXPRELTOL  1.0e-8
#define EXPRELTOL2 EXPRELTOL*EXPRELTOL
cdouble cExpRel2(int n, cdouble x)
{
  int m;
  cdouble Term, Sum;
  double mag2Term, mag2Sum;

  /*--------------------------------------------------------------*/
  /*- small-x expansion                                          -*/
  /*--------------------------------------------------------------*/
  Sum=1.0;
  for(Term=1.0, m=1; m<100; m++)
   { Term*=x/((double)(m+n));
     Sum+=Term;
     mag2Term=norm(Term);
     mag2Sum=norm(Sum);
     if ( mag2Term < EXPRELTOL2*mag2Sum )
      break;
    };
  return Sum;
} 


/***************************************************************/
/* quick factorial function needed below                       */
/***************************************************************/
static inline double factorial(int n) 
{ 
  static double FactTable[7]={1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 620.0}; 

  if (n<=6)
   return FactTable[n];
  else
   return ((double)n)*factorial(n-1);
}

// In function for g(r) = e^{-Kappa*r}/r 
cdouble In_EMKROverR_libRWG(int n, cdouble GParam, double X)
{ 
  double RetVal, Kappa = imag(GParam);
  double KX=Kappa*X;

  /* 20100610 note: if KX>40.0, then 
      exp(-KX)*gsl_sf_exprel_n(n,KX) 
     agrees with
      n! / (KX)^n
     to 12 digits (tested this for n up to 5).
     this modification essentially allows us to 
     go up to arbitrarily large values of kappa.
  */
  if ( KX > 40.0  )
   RetVal = factorial(n-1) / (X * pow(KX,n));
  else
   RetVal = exp(-KX)*gsl_sf_exprel_n(n,KX) / (X*(double)(n));
  
  return RetVal;
}

// g(r) = -(kr+1.0)*exp(-kr)/r^3
cdouble In_GradEMKROverR_libRWG(int n, cdouble GParam, double R)
{ 
  double RetVal, Kappa = imag(GParam);
  double KR=Kappa*R;
  RetVal=-Kappa*exp(-KR)/(R*R)
         *(  gsl_sf_exprel_n(n-1,KR)/(    (double)(n-1) )
           + gsl_sf_exprel_n(n-2,KR)/( KR*(double)(n-2) )
          );
  return RetVal;
}


// g(r) = exp(I*K*r)/r
cdouble In_EIKROverR_libRWG(int n, cdouble K, double R)
{ 
  cdouble KR=K*R;

  /* 20100610 note: if imag(KR) > 40.0, then
      \int_0^1 dw w^n e^{w*I*K*R} / (w*R)
     agrees with
      (n-1)! / ( R* (-I*K*R)^(n) )
     to high accuracy.
  */

  if ( (imag(KR)) > 40.0  )
   return factorial(n-1) / (R*pow(-II*K*R,n));
  else
   return exp(II*K*R)*cExpRel2(n,-II*K*R) / ((double)(n)*R);
}

// g(r) = (ikr-1.0)*exp(ikr)/r^3
cdouble In_GradEIKROverR_libRWG(int n, cdouble K, double R)
{ 
  cdouble KR=K*R;

  return II*K*exp(II*KR)/(R*R)
          *(      cExpRel2(n-1,-II*KR)/(       (double)(n-1) )
             -1.0*cExpRel2(n-2,-II*KR)/( II*KR*(double)(n-2) )
           );
}
