#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define X10 0.7
#define X20 1.1
#define A1 2.0
#define A2 3.0

double Phi0(double X1, double X2)
{ 
  return exp( -A1*(X1-X10)*(X1-X10)- A2*(X2-X20)*(X2-X20));
} 


void Phi1(double X1, double X2, double *PhiVD)
{
  double Delta=1.0e-3;
  PhiVD[0] = Phi0(X1, X2);
  PhiVD[1] = (Phi0(X1+Delta, X2) - Phi0(X1-Delta,X2))/(2.0*Delta);
  PhiVD[2] = (Phi0(X1, X2+Delta) - Phi0(X1,X2-Delta))/(2.0*Delta);
  PhiVD[3] = (   Phi0(X1+Delta, X2+Delta) 
               - Phi0(X1+Delta, X2-Delta)
               - Phi0(X1-Delta, X2+Delta)
               + Phi0(X1-Delta, X2-Delta)
             ) /(4.0*Delta*Delta);
}

void Phi2(double X1, double X2, double *PhiVD)
{
  double X1Bar, X2Bar;
  X1Bar=X1-X10;
  X2Bar=X2-X20;
  PhiVD[0] = exp( -A1*X1Bar*X1Bar - A2*X2Bar*X2Bar );
  PhiVD[1] = -2.0*A1*X1Bar*PhiVD[0];
  PhiVD[2] = -2.0*A2*X2Bar*PhiVD[0];
  PhiVD[3] = 4.0*A1*A2*X1Bar*X2Bar*PhiVD[0];
}

int main()
{ 
  double X1, X2, PhiVD1[4], PhiVD2[4]; 
  char buffer[100];

  srand48(time(0));
  while( fgets(buffer,100,stdin) )
   { 
     X1=2.0*drand48();
     X2=2.0*drand48();
     printf("\n Enter: ");
     sscanf(buffer,"%le %le",&X1,&X2);
     Phi1(X1,X2,PhiVD1);
     Phi2(X1,X2,PhiVD2);
     printf("(X1,X2)=(%e %e): \n",X1,X2);
     printf("Phi:   %+15.12e \n",PhiVD1[0]);
     printf("       %+15.12e \n",PhiVD2[0]);
     printf("PhiX:  %+15.12e \n",PhiVD1[1]);
     printf("       %+15.12e \n",PhiVD2[1]);
     printf("PhiY:  %+15.12e \n",PhiVD1[2]);
     printf("       %+15.12e \n",PhiVD2[2]);
     printf("PhiZ:  %+15.12e \n",PhiVD1[3]);
     printf("       %+15.12e \n",PhiVD2[3]);
     printf("       %+15.12e \n",PhiVD2[3]);
   };

}
