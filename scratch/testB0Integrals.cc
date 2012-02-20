/*
 * 
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <libhrutil.h>
#include <libSGJC.h>

#define NFUN 8

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct GIBFData 
 { 
   double a;
   double vp0;
   double b2;

 } GIBFData;

void GetIntegralsBF_Integrand(unsigned ndim, const double *x, void *params,
			      unsigned fdim, double *fval)
{
  double vp=x[0];

  GIBFData *GIBFD=(GIBFData *)params;

  double a=GIBFD->a;
  double vp0=GIBFD->vp0;
  double b2=GIBFD->b2;

  double R2 = a*a*((vp+vp0)*(vp+vp0) + b2);
  double R  = sqrt(R2);
  double R3 = R2*R;

  fval[0] = 1.0 / R3;
  fval[1] = vp  / R3;
  fval[2] = 1.0 / R;
  fval[3] = vp  / R;
  fval[4] = R;      
  fval[5] = vp*R;      
  fval[6] = R2;      
  fval[7] = vp*R2;      
} 

void GetIntegrals_BF(double a, double u, double up, double vp0, double b, double Result[NFUN])
{ 
  GIBFData MyGIBFData, *GIBFD=&MyGIBFData;

  GIBFD->a=a;
  GIBFD->vp0=vp0;
  GIBFD->b2=b*b;

  double Lower[1], Upper[1];
  double Error[NFUN];

  Lower[0] = 0.0;
  Upper[0] = up;

  adapt_integrate(NFUN, GetIntegralsBF_Integrand, (void *)GIBFD, 1,
		  Lower, Upper, 0, 0.0, 1.0e-8, Result, Error);

  int nf;
  for(nf=0; nf<NFUN; nf++)
   Result[nf]*=u;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetIntegrals_Analytic(double a, double u, double up, double vp0, double b, double Result[NFUN])
{
  double OneRM3Int, vpRM3Int; 
  double OneRM1Int, vpRM1Int; 
  double OneR1Int, vpR1Int; 
  double OneR2Int, vpR2Int; 

  double vp02=vp0*vp0;
  double vp03=vp02*vp0;
  double vp04=vp03*vp0;
  double b2=b*b;
  double a2=a*a, a3=a2*a;

  if ( b2 > 1.0e-10 )
   { 
     double S1, S2, LogFac, Sum, Sum3, Sum4;
     S1=sqrt( b2 + vp02 );
     S2=sqrt( b2 + (vp0+up)*(vp0+up) );
     LogFac=log( (S2 + (up+vp0)) / (S1 + vp0) );
     Sum=vp0+up;
     Sum3=Sum*Sum*Sum;
     Sum4=Sum3*Sum;
   
     OneRM3Int = u*( (up+vp0)/S2 - vp0/S1 ) / (a3*b2);
     vpRM3Int  = -vp0*OneRM3Int + u*( 1.0/S1 - 1.0/S2 ) / a3;
     OneRM1Int = u*LogFac/a;
     vpRM1Int  = -vp0*OneRM1Int + u*(S2-S1)/a;
     OneR1Int  = u*a*0.5*(b2*LogFac + (up+vp0)*S2 - vp0*S1);
     vpR1Int   = -vp0*OneR1Int + u*a*(S2*S2*S2 - S1*S1*S1) / 3.0;
     OneR2Int  = u*a2*( (Sum3-vp03)/3.0 + up*b2 );
     vpR2Int   = -vp0*OneR2Int + u*a2*( (Sum4-vp04)/4.0 + up*(vp0 + 0.5*up)*b2 );
   }
  else // b is close to zero
   {
printf("**b==0 case\n");
     double vp0pup = vp0 + up;
     double vp0pup2 = vp0pup*vp0pup;
     double vp0pup3 = vp0pup2*vp0pup;
     double vp0pup4 = vp0pup3*vp0pup;
     OneRM3Int = u*fabs( 1.0/vp02 - 1.0/vp0pup2 ) / (2.0*a3);
     vpRM3Int  = -vp0*OneRM3Int + u*fabs(1.0/vp0 - 1.0/vp0pup) / a3;
     OneRM1Int = u*fabs(log(vp0pup/vp0)) / a;
     vpRM1Int  = -vp0*OneRM1Int + u*up/a;
     OneR1Int  =  u*a*fabs(vp0pup2 - vp02) / 2.0;
     OneR2Int  = u*a2*fabs(vp0pup3 - vp03) / 3.0;
     vpR1Int   = -vp0*OneR1Int + OneR2Int/a;
     vpR2Int   = -vp0*OneR2Int + u*a2*fabs(vp0pup4 - vp04) / 4.0;
   };
  
  Result[0] = OneRM3Int;
  Result[1] =  vpRM3Int;
  Result[2] = OneRM1Int;
  Result[3] =  vpRM1Int;
  Result[4] = OneR1Int;
  Result[5] =  vpR1Int;
  Result[6] = OneR2Int;
  Result[7] =  vpR2Int;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  /*--------------------------------------------------------------*/
  /*- process command-line arguments -----------------------------*/
  /*--------------------------------------------------------------*/
  double a, u, up, vp0, b;
  ArgStruct ASArray[]=
   { {"a",         PA_DOUBLE, (void *)&a,        "-1.0",  "a"},
     {"u",         PA_DOUBLE, (void *)&u,         "1.0",  "u"},
     {"up",        PA_DOUBLE, (void *)&up,       "-1.0",  "up"},
     {"vp0",       PA_DOUBLE, (void *)&vp0,      "-1.0",  "vp0"},
     {"b",         PA_DOUBLE, (void *)&b,        "-1.0",  "b"},
     {0,0,0,0,0}
   };
  ProcessArguments(argc, argv, ASArray);

  srand48(time(0));

  /*--------------------------------------------------------------*/
  /* choose random values for any parameters that were left unset */
  /*--------------------------------------------------------------*/

  // choose a between 0.1 and 3 
  if (a==-1.0)
   a = 0.1 + (3.0 - 0.1) * drand48();

  // choose up between 0 and 1 
  if (up==-1.0)
   up = drand48();

  // choose b between -3 and 3 
  if (b==-1.0)
   b = -3.0 + 6.0*drand48();

  // choose vp0 between -3 and 2 but excluding the interval [-up:0]
  if (vp0==-1.0)
   if ( lrand48()%2 == 1)
    vp0=2.0*drand48();
   else
    vp0=-up-2.0*drand48();

  printf("--a %e --up %e --vp0 %e --b %e\n",a,up,vp0,b);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double IntegralsBF[NFUN], IntegralsAnalytic[NFUN];

  GetIntegrals_BF(a, u, up, vp0, b, IntegralsBF);
  GetIntegrals_Analytic(a, u, up, vp0, b, IntegralsAnalytic);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int ni;
  printf("i | %15s | %15s | %5s\n", "    Analytic   ","       BF      "," RD  ");
  for(ni=0; ni<NFUN; ni++)
   printf("%i | %+15.8e | %+15.8e | %5.2e\n",
           ni,IntegralsAnalytic[ni],IntegralsBF[ni],
           RD(IntegralsAnalytic[ni],IntegralsBF[ni]) );

}
