/*
 * libAmosBessel.cc -- a C++ wrapper around the original fortran routines 
 *                  -- of T. E. Amos for evaluating bessel-type functions
 *                  -- at complex arguments 
 *
 * homer reid       -- 2/2011
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex>
#include <cmath>

#include "libAmosBessel.h"

/***************************************************************/
/* prototypes for fortran functions ****************************/
/***************************************************************/
extern "C" {

 void zbesj_(double *zr, double *zi, double *fnu, int *kode, int *n, 
             double *cyr, double *cyi, int *nz, int *ierr);

 void zbesy_(double *zr, double *zi, double *fnu, int *kode, int *n, 
             double *cyr, double *cyi, int *nz, double *cwrkr, 
             double *cwrki, int *ierr);

 void zbesi_(double *zr, double *zi, double *fnu, int *kode, int *n, 
             double *cyr, double *cyi, int *nz, int *ierr);

 void zbesk_(double *zr, double *zi, double *fnu, int *kode, int *n, 
             double *cyr, double *cyi, int *nz, int *ierr);

 void zbesh_(double *zr, double *zi, double *fnu, int *kode, int *m,int *n, 
             double *cyr, double *cyi, int *nz, int *ierr);

 void zairy_(double *zr, double *zi, int *id, int *kode, 
             double *air, double *aii, int *nz, int *ierr);

 void zbiry_(double *zr, double *zi, int *id, int *kode, 
             double *air, double *aii, int *nz, int *ierr);
}

/***************************************************************/
/* WhichFunction:                                              */
/*                                                             */
/*  'J': regular cylindrical bessel J_n                        */
/*  'Y': irregular cylindrical bessel Y_n                      */
/*  'I': modified regular cylindrical bessel I_n               */
/*  'K': modified irregular cylindrical bessel K_n             */
/*  'O': type-1 cylindrical hankel function H^{(1)}_n          */
/*  'T': type-2 cylindricak hankel function H^{(2)}_n          */
/*                                                             */
/*  'j': regular cylindrical bessel j_n                        */
/*  'y': irregular cylindrical bessel y_n                      */
/*  'i': modified regular cylindrical bessel i_n               */
/*  'k': modified irregular cylindrical bessel k_n             */
/*  'o': type-1 spherical hankel function h^{(1)}_n            */
/*  't': type-2 spherical hankel function h^{(2)}_n            */
/***************************************************************/
int AmosBessel(char WhichFunction, cdouble z, double MinOrder, int NumOrders, 
               int Scale, cdouble *f)
{ 
  int no, nz, ierr, Type;
  double zr=real(z), zi=imag(z);
  double fr[NumOrders], fi[NumOrders];
  double cwrkr[NumOrders], cwrki[NumOrders];

  int kode = Scale ? 2 : 1;

  int Spherical=0;

  /***************************************************************/
  /* special treatment for j_n(0) and i_n(0)                     */
  /***************************************************************/
  if ( (WhichFunction=='j' || WhichFunction=='i') && abs(z)==0.0 )
   { memset(f,0,NumOrders*sizeof(cdouble));
     if (MinOrder==0)
      f[0]=cdouble(1.0,0.0);
     return 0;
   };
  
  switch(WhichFunction)
   { 
     case 'J':
        zbesj_( &zr, &zi, &MinOrder, &kode, &NumOrders, fr, fi, &nz, &ierr);
        break;

     case 'Y':
        zbesy_( &zr, &zi, &MinOrder, &kode, &NumOrders, fr, fi, &nz, cwrkr, cwrki, &ierr);
        break;

     case 'I':
        zbesi_( &zr, &zi, &MinOrder, &kode, &NumOrders, fr, fi, &nz, &ierr);
        break;

     case 'K':
        zbesk_( &(real(z)), &(imag(z)), &MinOrder, &kode, &NumOrders,
                fr, fi, &nz, &ierr);
        break;

     case 'O':
        Type=1;
        zbesh_( &(real(z)), &(imag(z)), &MinOrder, &kode, &Type, &NumOrders,
                fr, fi, &nz, &ierr);
        break;

     case 'T':
        Type=2;
        zbesh_( &(real(z)), &(imag(z)), &MinOrder, &kode, &Type, &NumOrders,
                fr, fi, &nz, &ierr);
        break;

     case 'j':
        Spherical=1;
        MinOrder+=0.5;
        zbesj_( &zr, &zi, &MinOrder, &kode, &NumOrders, fr, fi, &nz, &ierr);
        break;

     case 'y':
        Spherical=1;
        MinOrder+=0.5;
        zbesy_( &zr, &zi, &MinOrder, &kode, &NumOrders, fr, fi, &nz, cwrkr, cwrki, &ierr);
        break;

     case 'i':
        Spherical=1;
        MinOrder+=0.5;
        zbesi_( &zr, &zi, &MinOrder, &kode, &NumOrders, fr, fi, &nz, &ierr);
        break;

     case 'k':
        Spherical=1;
        MinOrder+=0.5;
        zbesk_( &(real(z)), &(imag(z)), &MinOrder, &kode, &NumOrders,
                fr, fi, &nz, &ierr);
        break;

     case 'o':
        Spherical=1;
        MinOrder+=0.5;
        Type=1;
        zbesh_( &(real(z)), &(imag(z)), &MinOrder, &kode, &Type, &NumOrders,
                fr, fi, &nz, &ierr);
        break;

     case 't':
        Spherical=1;
        MinOrder+=0.5;
        Type=2;
        zbesh_( &(real(z)), &(imag(z)), &MinOrder, &kode, &Type, &NumOrders,
                fr, fi, &nz, &ierr);
        break;

     default:
       fprintf(stderr,"%s:%i: internal error (WhichFunction=%c)",
                       __FILE__,__LINE__,WhichFunction);
       exit(1);
   };

  if (Spherical)
   { cdouble PreFactor;
     if (WhichFunction=='k') 
      PreFactor=sqrt( 2.0 / (M_PI*z) );
     else
      PreFactor=sqrt( M_PI/(2.0*z) );
     for(no=0; no<NumOrders; no++)
      f[no] = PreFactor*cdouble(fr[no],fi[no]);
   }
  else
   { for(no=0; no<NumOrders; no++)
      f[no] = cdouble( fr[no],fi[no] );
   };

  return ierr;

}

/***************************************************************/
/* WhichFunction:                                              */
/*  'A': value of Airy function Ai(z)                          */
/*  'a': derivative of Airy function dAi(z) /dz                */
/*  'B': value of Airy function Bi(z)                          */
/*  'b': derivative of Airy function dBi(z) /dz                */
/***************************************************************/
int AmosAiry(char WhichFunction, cdouble z, int Scale, cdouble *f)
{
  double zr=real(z), zi=imag(z);
  double fr, fi;
  int id, nz, ierr;

  int kode = Scale ? 2 : 1;

  switch(WhichFunction)
   { 
     case 'A':
        id=0;
        zairy_( &zr, &zi, &id, &kode, &fr, &fi, &nz, &ierr);
        break;

     case 'a':
        id=1;
        zairy_( &zr, &zi, &id, &kode, &fr, &fi, &nz, &ierr);
        break;

     case 'B':
        id=0;
        zbiry_( &zr, &zi, &id, &kode, &fr, &fi, &nz, &ierr);
        break;

     case 'b':
        id=1;
        zbiry_( &zr, &zi, &id, &kode, &fr, &fi, &nz, &ierr);
        break;

     default:
       fprintf(stderr,"%s:%i: internal error (WhichFunction=%c)",
                       __FILE__,__LINE__,WhichFunction);
       exit(1);
   };

  f[0]=cdouble(fr,fi);

 
 return ierr;

}
