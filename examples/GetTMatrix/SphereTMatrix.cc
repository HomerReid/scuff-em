/*
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libhrutil.h>
#include "libSpherical.h"
#include "libscuff.h"

using namespace scuff;

void GetSphereTMatrix(int l, cdouble Omega, cdouble Eps, cdouble Mu, double R,  
                      cdouble *TMM, cdouble *TEE)
{

 cdouble n  = sqrt(Eps*Mu);
 cdouble ka = Omega*R;
 cdouble nka = n*Omega*R;
 cdouble jlka, jlnka, hlka;
 cdouble jlkaSlash, jlnkaSlash, hlkaSlash; 
 GetRadialFunction(l,   Omega, R, LS_REGULAR,  &jlka,  0, &jlkaSlash);
 GetRadialFunction(l, n*Omega, R, LS_REGULAR,  &jlnka, 0, &jlnkaSlash);
 GetRadialFunction(l,   Omega, R, LS_OUTGOING, &hlka,  0,  &hlkaSlash);

 *TMM =  - (jlka * (nka*jlnkaSlash)  -  Mu*jlnka*(ka*jlkaSlash) )
         / (hlka * (nka*jlnkaSlash)  -  Mu*jlnka*(ka*hlkaSlash) );

 *TEE =  - (jlka * (nka*jlnkaSlash)  -  Eps*jlnka*(ka*jlkaSlash) )
         / (hlka * (nka*jlnkaSlash)  -  Eps*jlnka*(ka*hlkaSlash) );

}

int main()
{
  double Omega=0.01;

  FILE *f=fopen("E10Sphere_TMatrix.dat","w");
  cdouble TMM, TEE;
  for(Omega=0.01; Omega<10.0; Omega+=0.01)
   { 
     fprintf(f,"%e ",Omega);
     for(int l=1; l<=5; l++)
      { GetSphereTMatrix(l, Omega, 10.0, 1.0, 1.0, &TMM, &TEE);
        fprintf(f,"%e %e %e %e ",real(TMM),imag(TMM),real(TEE),imag(TEE));
      };
     fprintf(f,"\n");
   };
  fclose(f);

}
