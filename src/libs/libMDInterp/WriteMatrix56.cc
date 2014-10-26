#include <stdio.h>
#include <stdlib.h>

#include <libhmat.h>
#include <libhrutil.h>

int main(int argc, char *argv[])
{ 
   int ncp;
   double X1Bar, X2Bar, X3Bar;
   double X1Powers[6], X2Powers[6], X3Powers[6];
   int nm;
   int p, q, r;
   HMatrix *M;

   /*--------------------------------------------------------------*/
   /*- allocate the matrix ----------------------------------------*/
   /*--------------------------------------------------------------*/
   M=new HMatrix(56, 56);
   M->Zero();

   /*--------------------------------------------------------------*/
   /*- construct the matrix.                                      -*/
   /*- the outer loop (actually three loops) loops over the 8     -*/
   /*- corners of the grid cell.                                  -*/
   /*- the inner loop (actually three loops) loops over the 32    -*/
   /*- or 35 monomials in the interpolating polynomial.           -*/
   /*--------------------------------------------------------------*/
   X1Powers[0]=X1Powers[2]=X1Powers[4]=1.0;
   X2Powers[0]=X2Powers[2]=X2Powers[4]=1.0;
   X3Powers[0]=X3Powers[2]=X3Powers[4]=1.0;

   for(ncp=0, X1Bar=-1.0; X1Bar<=1.0; X1Bar+=2.0)
    for(X2Bar=-1.0; X2Bar<=1.0; X2Bar+=2.0)
     for(X3Bar=-1.0; X3Bar<=1.0; X3Bar+=2.0, ncp++)
      { 
        X1Powers[1]=X1Powers[3]=X1Powers[5]=X1Bar;
        X2Powers[1]=X2Powers[3]=X2Powers[5]=X2Bar;
        X3Powers[1]=X3Powers[3]=X3Powers[5]=X3Bar;

        nm=0;
        for(p=0; p<6; p++)
         for(q=0; q<6; q++)
          for(r=0; r<6; r++)
           { 
             if ( p+q+r > 5 ) continue;

             /* the nmth monomial is x1^p x2^q x3^r            */
             /* its x1 derivative is p x1^{p-1} x2^q x3^r, etc */
             M->SetEntry(7*ncp + 0, nm, X1Powers[p]*X2Powers[q]*X3Powers[r]);
             if (p>0)
              M->SetEntry(7*ncp + 1, nm, p*X1Powers[p-1]*X2Powers[q]*X3Powers[r]);
             if (q>0)
              M->SetEntry(7*ncp + 2, nm, q*X1Powers[p]*X2Powers[q-1]*X3Powers[r]);
             if (r>0)
              M->SetEntry(7*ncp + 3, nm, r*X1Powers[p]*X2Powers[q]*X3Powers[r-1]);
             if (p>0 && q>0)
              M->SetEntry(7*ncp + 4, nm, p*q*X1Powers[p-1]*X2Powers[q-1]*X3Powers[r]);
             if (p>0 && r>0)
              M->SetEntry(7*ncp + 5, nm, p*r*X1Powers[p-1]*X2Powers[q]*X3Powers[r-1]);
             if (q>0 && r>0)
              M->SetEntry(7*ncp + 6, nm, q*r*X1Powers[p]*X2Powers[q-1]*X3Powers[r-1]);

             nm++;
           };
if (nm!=56) printf("bawonkatage!! nm=%i\n",nm);
      }; 

   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   M->ExportToText("WriteMatrix56.dat");
   printf("Thank you for your support.\n");

}
