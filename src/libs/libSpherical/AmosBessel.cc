/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This file is part of SCUFF-EM.
 *
 * SCUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * SCUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * AmosBessel.cc -- a C++ wrapper around the original fortran routines
 *               -- of T. E. Amos for evaluating bessel-type functions
 *               -- at complex arguments
 *
 * homer reid    -- 2/2011
 *
 * --------------------------------------------------------------
 *
 * there are two routines here:
 *  
 * (1) int AmosBessel(char WhichFunction, cdouble z,
 *                    double MinOrder, int NumOrders,
 *                    bool Scale, cdouble *f, double *Workspace=0);
 *
 * (2) int AmosAiry(char WhichFunction, cdouble z,
 *                  bool Scale, cdouble *f);
 * 
 * --------------------------------------------------------------
 * 
 * the function parameters are:
 *
 *  WhichFunction: a character that selects which function you 
 *                 want to compute (see below) 
 *
 *              z: the complex argument of the bessel/airy function
 *
 *       MinOrder: 
 *      NumOrders: bessel functions of orders 
 *                 MinOrder, MinOrder+1, ... MinOrder+NumOrders-1
 *                 will be computed. (Note that MinOrder need not be 
 *                 an integer or even a half-integer, but may be  
 *                 any arbitrary real number >=0 ).
 *                 
 *          Scale: if Scale==true, scaled versions of the bessel functions 
 *                 are computed. 
 *                 the definition of 'scaled' depends on the function 
 *                 in question (see below). 
 *
 *              f: on entry, must point to a buffer with room for 
 *                 NumOrders cdoubles (for AmosBessel) or 1 cdouble
 *                 (for AmosAiry).
 *            
 *                 on return, f is filled in with the values of the 
 *                 function requested. for example, if WhichFunction='J',
 *                 MinOrder=3.25, NumOrders=3, then on return we have         
 *            
 *                 f[0]=J_{3.25}(z), f[1]=J_{4.25}(z), f[2]=J_{5.25}(z)
 *            
 *                 where J_{\nu}(z) is the regular cylindrical bessel 
 *                 function of order \nu.
 *            
 *     Workspace:  an optional user-allocated buffer with enough
 *                 room to store 4*NumOrders doubles. (Note that
 *                 Workspace has type double*, not cdouble*). If
 *                 this parameter is omitted, the workspace is
 *                 internally allocated and deleted.
 * --------------------------------------------------------------
 *            
 *  the return value from both functions is the IERR parameter in
 *  the fortran routine (=0 for a successful computation) 
 *            
 * --------------------------------------------------------------
 *            
 * Values of the 'WhichFunction' field:
 *
 * == for AmosBessel:
 *
 *  'J': regular cylindrical bessel J_n
 *  'Y': irregular cylindrical bessel Y_n
 *  'I': modified regular cylindrical bessel I_n
 *  'K': modified irregular cylindrical bessel K_n
 *  'O': type-1 cylindrical hankel function H^(1)_n
 *  'T': type-2 cylindricak hankel function H^(2)_n
 *
 *  'j': regular spherical bessel j_n
 *  'y': irregular spherical bessel y_n
 *  'i': modified regular spherical bessel i_n
 *  'k': modified irregular spherical bessel k_n
 *  'o': type-1 spherical hankel function h^(1)_n
 *  't': type-2 spherical hankel function h^(2)_n
 *
 * == for AmosAiry:
 *
 *  'A': value of Airy function Ai(z)
 *  'a': derivative of Airy function dAi(z) /dz
 *  'B': value of Airy function Bi(z)
 *  'b': derivative of Airy function dBi(z) /dz
 *
 * --------------------------------------------------------------
 *
 * as noted above, the meaning of 'scaling' differs for different
 * functions. if you set Scale=true on entry to AmosBessel or AmosAiry,  
 * the function returned will be: 
 *
 * WhichFunction | Quantity returned 
 * ---------------------------------
 *
 *     'J'       | e^{ -|im z| } * J_\nu(z) 
 *     'j'       | e^{ -|im z| } * j_\nu(z) 
 *
 *     'Y'       | e^{ -|im z| } * Y_\nu(z) 
 *     'y'       | e^{ -|im z| } * y_\nu(z) 
 *
 *     'I'       | e^{ -|re z| } * I_\nu(z) 
 *     'i'       | e^{ -|re z| } * i_\nu(z) 
 *
 *     'K'       | e^{ +z } * K_\nu(z) 
 *     'k'       | e^{ +z } * k_\nu(z) 
 *
 *     'O'       | e^{ -iz } * H^(1)_\nu(z) 
 *     'o'       | e^{ -iz } * h^(1)_\nu(z) 
 *
 *     'T'       | e^{ +iz } * H^(2)_\nu(z) 
 *     't'       | e^{ +iz } * h^(2)_\nu(z) 
 *
 *
 *     'A'       | e^{ 2/3 * z^(3/2) }* Ai(z)
 *     'a'       | e^{ 2/3 * z^(3/2) }* dAi(z)/dz
 *
 *     'B'       | e^{ -2/3 * z^(3/2) }* Bi(z)
 *     'b'       | e^{ -2/3 * z^(3/2) }* dBi(z)/dz
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex>
#include <cmath>
#include "libhrutil.h" // for Warn()

#ifndef cdouble
  typedef std::complex<double> cdouble;
#endif

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
/*                                                             */
/* Workspace is an an optionally user-supplied workspace       */
/* buffer. It may be NULL, but if it is non-null it must       */
/* contain enough space to hold 4*NumOrders doubles.           */
/* Note that Workspace is of type double*, not cdouble*.       */
/***************************************************************/
int AmosBessel(char WhichFunction, cdouble z,
               double MinOrder, int NumOrders,
               bool Scale, cdouble *f,
               double *Workspace)
{ 
  /***************************************************************/
  /* special treatment for j_n(0) and i_n(0)                     */
  /***************************************************************/
  if ( (WhichFunction=='j' || WhichFunction=='i') && abs(z)==0.0 )
   { memset(f,0,NumOrders*sizeof(cdouble));
     if (MinOrder==0)
      f[0]=cdouble(1.0,0.0);
     return 0;
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double *fr, *fi, *cwrkr, *cwrki;
  if (Workspace)
   { fr    = Workspace + 0*NumOrders;
     fi    = Workspace + 1*NumOrders;
     cwrkr = Workspace + 2*NumOrders;
     cwrki = Workspace + 3*NumOrders;
   }
  else
   { 
     fr    = new double[4*NumOrders];
     fi    = fr + 1*NumOrders;
     cwrkr = fr + 2*NumOrders;
     cwrki = fr + 3*NumOrders;
   };
  memset(fr,    0, (NumOrders)*sizeof(double));
  memset(fi,    0, (NumOrders)*sizeof(double));
  memset(cwrkr, 0, (NumOrders)*sizeof(double));
  memset(cwrki, 0, (NumOrders)*sizeof(double));

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int no, nz, ierr, Type;
  double zr=real(z), zi=imag(z);
  int kode = Scale ? 2 : 1;
  int Spherical=0;
  
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
        zbesk_( &zr, &zi, &MinOrder, &kode, &NumOrders,
                fr, fi, &nz, &ierr);
        break;

     case 'O':
        Type=1;
        zbesh_( &zr, &zi, &MinOrder, &kode, &Type, &NumOrders,
                fr, fi, &nz, &ierr);
        break;

     case 'T':
        Type=2;
        zbesh_( &zr, &zi, &MinOrder, &kode, &Type, &NumOrders,
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
        zbesk_( &zr, &zi, &MinOrder, &kode, &NumOrders,
                fr, fi, &nz, &ierr);
        break;

     case 'o':
        Spherical=1;
        MinOrder+=0.5;
        Type=1;
        zbesh_( &zr, &zi, &MinOrder, &kode, &Type, &NumOrders,
                fr, fi, &nz, &ierr);
        break;

     case 't':
        Spherical=1;
        MinOrder+=0.5;
        Type=2;
        zbesh_( &zr, &zi, &MinOrder, &kode, &Type, &NumOrders,
                fr, fi, &nz, &ierr);
        break;

     default:
       fprintf(stderr,"%s:%i: internal error (WhichFunction=%c)\n",
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

  if (!Workspace)
   delete[] fr;

  if (ierr!=0)
   Warn("error in AmosBessel(%c,|z|=%g) (ierr=%i)",WhichFunction,abs(z),ierr);

  return ierr;

}

/***************************************************************/
/* WhichFunction:                                              */
/*  'A': value of Airy function Ai(z)                          */
/*  'a': derivative of Airy function dAi(z) /dz                */
/*  'B': value of Airy function Bi(z)                          */
/*  'b': derivative of Airy function dBi(z) /dz                */
/***************************************************************/
int AmosAiry(char WhichFunction, cdouble z, bool Scale, cdouble *f)
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
       fprintf(stderr,"%s:%i: internal error (WhichFunction=%c)\n",
                       __FILE__,__LINE__,WhichFunction);
       exit(1);
   };

  f[0]=cdouble(fr,fi);

 
 return ierr;

}
