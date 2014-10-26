
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
   double X1Powers[6], X2Powers[6], X3Powers[6];
   double X1Derivs[6], X2Derivs[6], X3Derivs[6];
   int nm, nMonomials;
   int p, q, r, pqrMax;
   HMatrix *M;

   /*--------------------------------------------------------------*/
   /*- process command-line arguments -----------------------------*/
   /*--------------------------------------------------------------*/
   ArgStruct ASArray[]=
    { {"nMonomials", PA_INT, (void *)&nMonomials,  "32",   "number of monomials"},
      {0,0,0,0,0}
    };
   ProcessArguments(argc, argv, ASArray);

   /*--------------------------------------------------------------*/
   /*- allocate the matrix ----------------------------------------*/
   /*--------------------------------------------------------------*/
   M=new HMatrix(32, nMonomials);
   M->Zero();

   /*--------------------------------------------------------------*/
   /*- construct the matrix.                                      -*/
   /*- the outer loop (actually three loops) loops over the 8     -*/
   /*- corners of the grid cell.                                  -*/
   /*- the inner loop (actually three loops) loops over the 32    -*/
   /*- or 35 monomials in the interpolating polynomial.           -*/
   /*--------------------------------------------------------------*/

   pqrMax = (nMonomials == 32) ? 4 : 5;

   for(ncp=0, X1Bar=-1.0; X1Bar<=1.0; X1Bar+=2.0)
    for(X2Bar=-1.0; X2Bar<=1.0; X2Bar+=2.0)
     for(X3Bar=-1.0; X3Bar<=1.0; X3Bar+=2.0, ncp++)
      { 
        for(p=0; p<6; p++)
         { X1Powers[p]=P(p,X1Bar);
           X1Derivs[p]=dP(p,X1Bar);
           X2Powers[p]=P(p,X2Bar);
           X2Derivs[p]=dP(p,X2Bar);
           X3Powers[p]=P(p,X3Bar);
           X3Derivs[p]=dP(p,X3Bar);
         };

        nm=0;
        for(p=0; p<pqrMax; p++)
         for(q=0; q<pqrMax; q++)
          for(r=0; r<pqrMax; r++)
           { 
             if ( p+q+r > 4 ) continue;

             /* the nmth monomial is x1^p x2^q x3^r            */
             /* its x1 derivative is p x1^{p-1} x2^q x3^r, etc */
             M->SetEntry(4*ncp + 0, nm, X1Powers[p]*X2Powers[q]*X3Powers[r]);
             M->SetEntry(4*ncp + 1, nm, X1Derivs[p]*X2Powers[q]*X3Powers[r]);
             M->SetEntry(4*ncp + 2, nm, X1Powers[p]*X2Derivs[q]*X3Powers[r]);
             M->SetEntry(4*ncp + 3, nm, X1Powers[p]*X2Powers[q]*X3Derivs[r]);

             nm++;
           };
      }; 

   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   M->ExportToText("WriteMatrixLP.dat");
   printf("Thank you for your support.\n");

}
