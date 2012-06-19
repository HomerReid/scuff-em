/*
 * PlateWithHole.cc -- a simple example of how to use the
 *                  -- scuff-EM C++ API
 *
 * This program computes the diagonal elements of the
 * scattering portion of the dyadic Green's function
 * for a plate-with-hole geometry.
 *
 * Homer Reid       -- 6/2012
 */

#include <stdio.h>
#include <stdlib.h>

#include <libscuff.h>
using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  /*--------------------------------------------------------------*/
  /* create the RWGGeometry from the .scuffgeo file               */
  /*--------------------------------------------------------------*/
  RWGGeometry *G=new RWGGeometry("PlateWithHole.scuffgeo");
  SetLogFileName("PlateWithHole.log");
  G->SetLogLevel(SCUFF_VERBOSELOGGING);
  PreloadCache("PlateWithHole.scuffcache");

  /*--------------------------------------------------------------*/
  /* preallocate BEM matrix and RHS vector                        */
  /*--------------------------------------------------------------*/
  HMatrix *M  = G->AllocateBEMMatrix();
  HVector *KN = G->AllocateRHSVector();

  /*--------------------------------------------------------------*/
  /* preallocate point-source incident field (the location and    */
  /* direction are arbitrary as we will change them later)        */
  /*--------------------------------------------------------------*/
  double X0[3] = {0.0, 0.0, 0.0};  // point source location 
  cdouble P[3] = {1.0, 0.0, 0.0};  // point source direction
  PointSource PS(X0, P);

  /*--------------------------------------------------------------*/
  /*- outer loop over frequencies --------------------------------*/
  /*--------------------------------------------------------------*/
  double Omega;
  cdouble EH[6];
  FILE *f=fopen("PlateWithHole.out","w");
  for(Omega=0.01; Omega<=1.1; Omega*=sqrt(10.0))
   {
     /*--------------------------------------------------------------*/
     /* assemble and factorize the BEM matrix at this frequency      */
     /*--------------------------------------------------------------*/
     G->AssembleBEMMatrix(Omega, M);
     M->LUFactorize();

    /*--------------------------------------------------------------*/
    /*- inner loop over Z coordinate of point source ---------------*/
    /*--------------------------------------------------------------*/
    for ( X0[2] = 0.0; X0[2] < 4.0; X0[2]+=0.1 )
     { 
       fprintf(f,"%e %e ",Omega,X0[2]);

       // set the location of the point source 
       PS.SetX0(X0);

       // solve three separate scattering problems (one for each 
       // possible direction of the point source) and in each case 
       // extract the corresponding electric field component
       for(int i=0; i<3; i++)
        { 
           // configure the point source to point in the i direction
           memset(P,0,3*sizeof(cdouble));
           P[i]=1.0;
           PS.SetP(P);  

           // solve the scattering problem for this point source
           G->AssembleRHSVector(Omega, &PS, KN);
           M->LUSolve(KN);

           // compute the scattered E and H fields at X0 
           // and extract the ith component of the E field
           G->GetFields(0, KN, Omega, X0, EH);
           fprintf(f,"%e %e ",real(EH[i]),imag(EH[i]));

        }; // for (i=0; ... 
       fprintf(f,"\n");

     }; // for( X0[0]= ... 
    fprintf(f,"\n\n");  // put 2 newlines between frequencies in the data file
      
   }; // for (Omega= ... );

 fclose(f);
 StoreCache("PlateWithHole.scuffcache");
 printf("Thank you for your support.\n");
  
}
