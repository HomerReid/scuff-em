/*
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libhrutil.h>
#include "libscuff.h"
#include "libSpherical.h"
#include "libMatProp.h"

using namespace scuff;

/***************************************************************/
/* get the diagonal T-matrix elements for a sphere at a given  */
/* frequency                                                   */
/***************************************************************/
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

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  char *Material=0;
  int lMax=3;
  int R=1.0;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"Material", PA_STRING,  1, 1, (void *)&Material,     0,  "material designation"},
     {"lMax",     PA_INT,     1, 1, (void *)&lMax,         0,  "maximum L-value"},
     {"R",        PA_DOUBLE,  1, 1, (void *)&R,            0,  "sphere radius"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (Material==0)
   OSUsage(argv[0],OSArray,"--Material option is mandatory");

  MatProp *MP = new MatProp(Material);
  if (MP->ErrMsg)
   ErrExit(MP->ErrMsg);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  FILE *f=vfopen("%s.SphereTMatrix","w",Material);
  cdouble TMM, TEE;
  double Omega=0.01;
  cdouble Eps, Mu;
  for(Omega=0.01; Omega<10.0; Omega+=0.01)
   { 
     MP->GetEpsMu(Omega, &Eps, &Mu);
     fprintf(f,"%e ",Omega);
  
     for(int l=1; l<=lMax; l++)
      { GetSphereTMatrix(l, Omega, Eps, Mu, 1.0, &TMM, &TEE);
        fprintf(f,"%e %e %e %e ",real(TMM),imag(TMM),real(TEE),imag(TEE));
      };
     fprintf(f,"\n");
   };
  fclose(f);

}
