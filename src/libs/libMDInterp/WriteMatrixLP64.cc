#include <stdio.h>
#include <stdlib.h>

#include <libhmat.h>
#include <libhrutil.h>

/***************************************************************/
/* legendre polynomials ****************************************/
/***************************************************************/
double P(int Order, double u)
{ 
  double u2=u*u;

  switch(Order)
   { case 0: return 1.0;
     case 1: return u;
     case 2: return 0.5*(3.0*u2-1.0);
     case 3: return 0.5*(5.0*u2-3.0)*u;
     case 4: return 0.125*(  (  35.0*u2- 30.0)*u2 +   3.0);
     case 5: return 0.125*(  (  63.0*u2- 70.0)*u2 +  15.0)*u;
     case 6: return 0.0625*( ((231.0*u2-315.0)*u2 + 105.0)*u2 - 5.0);
  };
  return 0.0;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
double dP(int Order, double u)
{ 
  double u2=u*u, u3=u2*u, u4=u2*u2;

  switch(Order)
   { case 0: return 0.0;
     case 1: return 1.0;
     case 2: return 3.0*u;
     case 3: return 0.5*(3.0*5.0*u2-3.0);
     case 4: return 0.125*( 4.0*35.0*u3- 2.0*30.0*u);
     case 5: return 0.125*( 5.0*63.0*u4 - 3.0*70.0*u2 +  15.0);
  };
  return 0.0;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{ 
   int ncp;
   double X1Bar, X2Bar, X3Bar;
   double X1Powers[4], X2Powers[4], X3Powers[4];
   double X1Derivs[4], X2Derivs[4], X3Derivs[4];
   int nm;
   int p, q, r;
   HMatrix *M;

   /*--------------------------------------------------------------*/
   /*- allocate the matrix ----------------------------------------*/
   /*--------------------------------------------------------------*/
   M=new HMatrix(64, 64);
   M->Zero();

   X1Powers[0]=X1Powers[2]=1.0;
   X2Powers[0]=X2Powers[2]=1.0;
   X3Powers[0]=X3Powers[2]=1.0;

   /*--------------------------------------------------------------*/
   /*- construct the matrix.                                      -*/
   /*- the outer loop (actually three loops) loops over the 8     -*/
   /*- corners of the grid cell.                                  -*/
   /*--------------------------------------------------------------*/
   for(ncp=0, X1Bar=-1.0; X1Bar<=1.0; X1Bar+=2.0)
    for(X2Bar=-1.0; X2Bar<=1.0; X2Bar+=2.0)
     for(X3Bar=-1.0; X3Bar<=1.0; X3Bar+=2.0, ncp++)
      { 

/*
         X1Powers[1]=X1Powers[3]=X1Bar;
         X2Powers[1]=X2Powers[3]=X2Bar;
         X3Powers[1]=X3Powers[3]=X3Bar;
*/

        for(p=0; p<4; p++)
         { X1Powers[p]=P(p,X1Bar);
           X1Derivs[p]=dP(p,X1Bar);
           X2Powers[p]=P(p,X2Bar);
           X2Derivs[p]=dP(p,X2Bar);
           X3Powers[p]=P(p,X3Bar);
           X3Derivs[p]=dP(p,X3Bar);
         };

        nm=0;
        for(p=0; p<4; p++)
         for(q=0; q<4; q++)
          for(r=0; r<4; r++)
           { 
             /* the nmth monomial is x1^p x2^q x3^r            */
             /* its x1 derivative is p x1^{p-1} x2^q x3^r, etc */
             M->SetEntry(8*ncp + 0, nm, X1Powers[p]*X2Powers[q]*X3Powers[r]);
             M->SetEntry(8*ncp + 1, nm, X1Derivs[p]*X2Powers[q]*X3Powers[r]);
             M->SetEntry(8*ncp + 2, nm, X1Powers[p]*X2Derivs[q]*X3Powers[r]);
             M->SetEntry(8*ncp + 3, nm, X1Powers[p]*X2Powers[q]*X3Derivs[r]);
             M->SetEntry(8*ncp + 4, nm, X1Derivs[p]*X2Derivs[q]*X3Powers[r]);
             M->SetEntry(8*ncp + 5, nm, X1Derivs[p]*X2Powers[q]*X3Derivs[r]);
             M->SetEntry(8*ncp + 6, nm, X1Powers[p]*X2Derivs[q]*X3Derivs[r]);
             M->SetEntry(8*ncp + 7, nm, X1Derivs[p]*X2Derivs[q]*X3Derivs[r]);

             nm++;
           };
      }; 

   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   M->ExportToText("MLP64");
   printf("Thank you for your support.\n");

}
