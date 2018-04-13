#include <stdio.h>
#include <stdlib.h>

#include <libhmat.h>
#include <libhrutil.h>

int main(int argc, char *argv[])
{ 
   int ncp;
   double X1Bar, X2Bar;
   double X1Powers[5], X2Powers[5];
   int nm, nMonomials;
   int p, q;
   HMatrix *M;

   /*--------------------------------------------------------------*/
   /*- allocate the matrix ----------------------------------------*/
   /*--------------------------------------------------------------*/
   nMonomials=15;
   M=new HMatrix(12, nMonomials);
   M->Zero();

   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   X1Powers[0]=X1Powers[2]=X1Powers[4]=1.0;
   X2Powers[0]=X2Powers[2]=X2Powers[4]=1.0;

   for(ncp=0, X1Bar=-1.0; X1Bar<=1.0; X1Bar+=2.0)
    for(X2Bar=-1.0; X2Bar<=1.0; X2Bar+=2.0, ncp++)
     { 
       X1Powers[1]=X1Powers[3]=X1Bar;
       X2Powers[1]=X2Powers[3]=X2Bar;

       nm=0;
       for(p=0; p<5; p++)
        for(q=0; q<5; q++)
         { 
           if ( p+q > 4 ) continue;

           /* the nmth monomial is x1^p x2^q x3^r            */
           /* its x1 derivative is p x1^{p-1} x2^q x3^r, etc */
           M->SetEntry(3*ncp + 0, nm, X1Powers[p]*X2Powers[q]);
           if (p>0)
            M->SetEntry(3*ncp + 1, nm, p*X1Powers[p-1]*X2Powers[q]);
           if (q>0)
            M->SetEntry(3*ncp + 2, nm, q*X1Powers[p]*X2Powers[q-1]);

           nm++;
         };
if (nm!=15) printf("bawonkatage! nm=%i\n",nm);
      }; 

   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   M->ExportToText("M2D");
   printf("Thank you for your support.\n");

}
