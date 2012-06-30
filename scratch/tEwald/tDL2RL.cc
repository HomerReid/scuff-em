/*
 * DL2RL.cc 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex>

#include <libhrutil.h>

#define NFIRSTROUND 10
#define NSUM 1
#define NMAX 100000

typedef std::complex<double> cdouble;

#define II cdouble(0.0,1.0)

/***************************************************************/
/* global variables ********************************************/
/***************************************************************/

#ifdef SQUARELATTICE 
// direct lattice vectors 
double L1[2]  = { 2.0, 0.0 };
double L2[2]  = { 0.0, 3.0 };

// reciprocal lattice vectors 
double Gamma1[2] = { 2.0*M_PI/2.0, 0.0      };
double Gamma2[2] = {          0.0,  2.0*M_PI/3.0 };

// volume of 2D Brillouin zone 
double AGamma=4.0*M_PI*M_PI/6.0;
#endif

#define HEXAGONALLATTICE
#ifdef HEXAGONALLATTICE

#define LL 1.0
#define RT3 1.73205080756888

// direct lattice vectors 
double L1[2]  = { LL,     0.0        };
double L2[2]  = { 0.5*LL, 0.5*RT3*LL };

// reciprocal lattice vectors 
double Gamma1[2] = { 2.0*M_PI/LL, 2.0*M_PI/(LL*RT3) };
double Gamma2[2] = {         0.0, 4.0*M_PI/(LL*RT3) };

double AGamma=8.0*M_PI*M_PI/(LL*LL*RT3);
#endif

int nCells;

double ABSTOL=1.0e-12;
double RELTOL=1.0e-8;

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddDLTerm(double s, cdouble k, double *P, double *X,  
               int n1, int n2, cdouble *Sum)
{
  double L[2];

  L[0] = n1*L1[0] + n2*L2[0];
  L[1] = n1*L1[1] + n2*L2[1];
  
  double XpL2 = (X[0]+L[0])*(X[0]+L[0]) + (X[1]+L[1])*(X[1]+L[1]) + X[2]*X[2];
  double PdL  = P[0]*L[0] + P[1]*L[1];

  double s2=s*s;
  Sum[0] += exp( II*PdL ) * exp( -s2*XpL2 + k*k/(4.0*s2) );
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddRLTerm(double s, cdouble k, double *P, double *X,  
               int n1, int n2, cdouble *Sum)
{
  double Gamma[2];

  Gamma[0] = n1*Gamma1[0] + n2*Gamma2[0];
  Gamma[1] = n1*Gamma1[1] + n2*Gamma2[1];
  
  cdouble GammapP2pK2 =  (Gamma[0]-P[0])*(Gamma[0]-P[0])
                       + (Gamma[1]-P[1])*(Gamma[1]-P[1])
                       - k*k;
                   
  double GammadX = Gamma[0]*X[0] + Gamma[1]*X[1];

  double s2=s*s;
  Sum[0] += exp( II*GammadX ) * exp( -0.25*GammapP2pK2/s2 - s2*X[2]*X[2]);
  
}
              
/***************************************************************/
/* DL = 1 --> direct lattice sum *******************************/
/* DL = 0 --> reciprocal  lattice sum **************************/
/***************************************************************/
cdouble GetSum(int DL, double s, cdouble k, double *P, double *X)
{

  nCells=0;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int n1, n2;
  cdouble Sum[NSUM];
  memset(Sum,0,NSUM*sizeof(cdouble));
  for (n1=-NFIRSTROUND; n1<=NFIRSTROUND; n1++)
   for (n2=-NFIRSTROUND; n2<=NFIRSTROUND; n2++, nCells++)
    DL ? AddDLTerm(s, k, P, X, n1, n2, Sum) : AddRLTerm(s, k, P, X, n1, n2, Sum); 
         
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble LastSum[NSUM];
  memcpy(LastSum,Sum,NSUM*sizeof(cdouble));
  int ConvergedIters=0;
  int NN;
  double Delta, AbsDelta, RelDelta, MaxAbsDelta, MaxRelDelta, AbsSum;
  for(NN=NFIRSTROUND+1; ConvergedIters<3 && NN<=NMAX; NN++)
   {  
     for(n1=-NN; n1<=NN; n1++)
      for(n2=-NN; n2<=NN; n2++)
        { 
          if ( (abs(n1)<NN) && (abs(n2)<NN) )
           continue;

          DL ? AddDLTerm(s, k, P, X, n1, n2, Sum) : AddRLTerm(s, k, P, X, n1, n2, Sum); 

          nCells++;

        };

     /*--------------------------------------------------------------*/
     /* convergence analysis ----------------------------------------*/
     /*--------------------------------------------------------------*/
     MaxAbsDelta=MaxRelDelta=0.0;
     for(int i=0; i<NSUM; i++)
      { Delta=abs(Sum[i]-LastSum[i]);
        if ( Delta>MaxAbsDelta )
         MaxAbsDelta=Delta;
        AbsSum=abs(Sum[i]);
        if ( AbsSum>0.0 && (Delta > MaxRelDelta*AbsSum) )
         MaxRelDelta=Delta/AbsSum;
      };
     if ( MaxAbsDelta<ABSTOL || MaxRelDelta<RELTOL )
      ConvergedIters++;
     else
      ConvergedIters=0;

     memcpy(LastSum,Sum,NSUM*sizeof(cdouble));
   };

  if (DL)
   return 2.0*Sum[0]/(4.0*M_PI*sqrt(M_PI));
  else
   return AGamma*exp(II*(-P[0]*X[0]-P[1]*X[1]))*Sum[0]/(s*s*8.0*M_PI*M_PI*sqrt(M_PI));

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  /*--------------------------------------------------------------*/
  /*- process options  -------------------------------------------*/
  /*--------------------------------------------------------------*/
  double s=0.0;                   int ns;
  cdouble k=0.0;                  int nk;
  double P[2]={0.0, 0.0};         int npx, npy;
  double X[3]={0.0, 0.0, 0.0};    int nX, nY, nZ;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"s",       PA_DOUBLE,  1, 1, (void *)&s,      &ns,       "s"},
     {"k",       PA_CDOUBLE, 1, 1, (void *)&k,      &nk,       "k"},
     {"px",      PA_DOUBLE,  1, 1, (void *)(P+0),  &npx,       "px"},
     {"py",      PA_DOUBLE,  1, 1, (void *)(P+1),  &npy,       "py"},
     {"x",       PA_DOUBLE,  1, 1, (void *)(X+0),   &nX,       "x",},
     {"y",       PA_DOUBLE,  1, 1, (void *)(X+1),   &nY,       "y",},
     {"z",       PA_DOUBLE,  1, 1, (void *)(X+2),   &nZ,       "z",},
     {"ABSTOL",  PA_DOUBLE,  1, 1, (void *)&ABSTOL,   0,       "abstol",},
     {"RELTOL",  PA_DOUBLE,  1, 1, (void *)&RELTOL,   0,       "reltol",},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  srand48(time(0));
  if (ns==0)
   { s=drand48(); s=s/(1.0-s); }  
  if (nk==0)
   { real(k)=drand48(); imag(k)=drand48(); }
  if (npx==0)
   { P[0]=2.0*(drand48()-0.5); }
  if (npy==0)
   { P[1]=2.0*(drand48()-0.5); }
  if (nX==0)
   { X[0]=5.0*(drand48()-0.5); }
  if (nY==0)
   { X[1]=5.0*(drand48()-0.5); }
  if (nZ==0)
   { X[2]=5.0*(drand48()-0.5); }
 
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  nCells=0;
  cdouble DLSum = GetSum(1, s, k, P, X);
  printf("Direct:     %i cells\n",nCells);

  nCells=0;
  cdouble RLSum = GetSum(0, s, k, P, X);
  printf("Reciprocal: %i cells\n",nCells);

  printf("--s %g --k %g+%gi --px %g --py %g --x %g --y %g --z %g\n",s,real(k),imag(k),P[0],P[1],X[0],X[1],X[2]);
  printf("DL: (%+10.3e , %+10.3e )\n",real(DLSum),imag(DLSum));
  printf("RL: (%+10.3e , %+10.3e )\n",real(RLSum),imag(RLSum));
  printf("RD: %e\n",abs(DLSum-RLSum)/abs(DLSum));
  printf("Ratio: (%+10.3e , %+10.3e )\n",real(DLSum/RLSum),imag(DLSum/RLSum));


}
