/*
 * tPointSource.cc 
 *
 * This program prints the E-field of a 
 */
#include <stdio.h>
#include <stdlib.h>

#include <libhrutil.h>
#include <libIncField.h>

int main(int argc, char *argv[])
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble Omega=0.001;         // angular frequency
  double X[3]={1.0, 0.0, 0.0}; // field evaluation point
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"Omega",    PA_CDOUBLE, 1, 1, (void *)&Omega,  0, "angular frequency"},
     {"X",        PA_DOUBLE,  3, 1, (void *)X,       0, "evaluation point "},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double X0[3] = {0.0, 0.0, 0.0};
  cdouble P[3] = {0.0, 0.0, 1.0};
  PointSource *PS = new PointSource(X0, P);

  cdouble EH[6]; 
  PS->SetFrequencyAndEpsMu(Omega, 1.0, 1.0, false);
  PS->GetFields(X, EH);

  printf("Ex=(%+10.2e,%+10.2e)\n",real(EH[0]),imag(EH[0]));
  printf("Ey=(%+10.2e,%+10.2e)\n",real(EH[1]),imag(EH[1]));
  printf("Ez=(%+10.2e,%+10.2e)\n",real(EH[2]),imag(EH[2]));
  printf("Hx=(%+10.2e,%+10.2e)\n",real(EH[3]),imag(EH[3]));
  printf("Hy=(%+10.2e,%+10.2e)\n",real(EH[4]),imag(EH[4]));
  printf("Hz=(%+10.2e,%+10.2e)\n",real(EH[5]),imag(EH[5]));

}
