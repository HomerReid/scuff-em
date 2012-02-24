/*
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex>

#include <libhrutil.h>

#define EXPRELTOL  1.0e-8
#define EXPRELTOL2 EXPRELTOL*EXPRELTOL

/***************************************************************/
/***************************************************************/
/***************************************************************/
cdouble cExpRel_SmallZ(int n, cdouble Z)
{
  int m;
  cdouble Term, Sum;
  double mag2Term, mag2Sum;

  Sum=1.0;
  for(Term=1.0, m=1; m<100; m++)
   { Term*=Z/((double)(m+n));
     Sum+=Term;
     mag2Term=norm(Term);
     mag2Sum=norm(Sum);
     if ( mag2Term < EXPRELTOL2*mag2Sum )
      break;
   };
  return Sum;
} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
cdouble cExpRel3(int n, cdouble Z)
{
  cdouble Term, Sum;
  int m;

  Sum=exp(Z);
  for(Term=1.0, m=0; m<n; m++)
   { 
     Sum-=Term;
     Term*=Z/((double)(m+1));
   };
  return Sum / Term;
} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
cdouble cEZpRel4(int n, cdouble Z)
{
  int m;
  cdouble Term, Sum;
  double mag2Term, mag2Sum;

  /*--------------------------------------------------------------*/
  /*- small-Z expansion                                          -*/
  /*--------------------------------------------------------------*/
  if ( abs(Z) < 0.1 )
   { Sum=1.0;
     for(Term=1.0, m=1; m<100; m++)
      { Term*=Z/((double)(m+n));
        Sum+=Term;
        mag2Term=norm(Term);
        mag2Sum=norm(Sum);
        if ( mag2Term < EXPRELTOL2*mag2Sum )
         break;
      };
     return Sum;
   }
  else
   {
     Sum=exp(Z);
     for(Term=1.0, m=0; m<n; m++)
      { 
        Sum-=Term;
        Term*=Z/((double)(m+1));
      };
     return Sum / Term;
   };
} 

/***************************************************************/
/* *************************************************************/
/***************************************************************/
double ExpRel4(int n, double Z)
{
  int m;
  double Term, Sum;

  /*--------------------------------------------------------------*/
  /*- small-Z expansion                                          -*/
  /*--------------------------------------------------------------*/
  if ( fabs(Z) < 0.1 )
   { Sum=1.0;
     for(Term=1.0, m=1; m<100; m++)
      { Term *= Z/((double)(m+n));
        Sum+=Term;
        if ( fabs(Term) < EXPRELTOL*fabs(Sum))
         break;
      };
     return Sum;
   }
  else
   {
     Sum=exp(Z);
     for(Term=1.0, m=0; m<n; m++)
      { 
        Sum-=Term;
        Term*=Z/((double)(m+1));
      };
     return Sum / Term;
   };
} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{ 
  /*--------------------------------------------------------------*/
  /*- process options  -------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble Z;
  int n;
  int NTimes=1000;
  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { {"Z",      PA_CDOUBLE, 1, 1, (void *)&Z,       0,  "Z"},
     {"n",      PA_INT,     1, 1, (void *)&n,       0,  "n"},
     {"NTimes", PA_INT,     1, 1, (void *)&NTimes,  0,  "number of times"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble ERSmallZ, ER
  double Elapsed2, Elapsed3;
  int nTimes;

  Tic();
  for(nTimes=0; nTimes<NTimes; nTimes++)
   ER2=cExpRel4(n, Z);
  Elapsed2=Toc();

  Tic();
  for(nTimes=0; nTimes<NTimes; nTimes++)
   ER3=cdouble(ExpRel4(n, real(Z)),0);
  Elapsed3=Toc();

  printf("Complex Z: %.3f us\n",Elapsed2*1.0e6);
  printf("Real    Z: %.3f us\n",Elapsed3*1.0e6);

  SetDefaultCD2SFormat("(%+16.8e,%+16.8e)");
  printf("%35s %35s %s\n","                HR2                ","                HR3                "," RD "); 
  printf("%35s-%35s-%s\n","-----------------------------------","-----------------------------------","-----"); 
  printf("%s %s %e\n",CD2S(ER2),CD2S(ER3),RD(ER2,ER3));

}
