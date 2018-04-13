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
 * tlibVVQAG.cc 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#include <libvvqag.h>

int count = 0;
const double radius = 0.50124145262344534123412; /* random */

void f_test(unsigned dim, const double *x, void *data, 
            unsigned nfun, double *fvec)
{
     unsigned i;
     ++count;

     /* first component of integrand: volume of hypersphere */
     fvec[0] = 0;
     for (i = 0; i < dim; ++i)
      fvec[0] += x[i] * x[i];
     fvec[0] = fvec[0] < radius * radius;

     /* second component of integrand: */
     /* simple smooth (separable) objective: prod. cos(x[i]). */
     fvec[1]=1.0;
     for (i = 0; i < dim; ++i)
      fvec[1] *= cos(x[i]);

     /* third component of integrand: */
     /* integral of exp(-x^2), rescaled to (0,infinity) limits */
     double scale = 1.0;
     fvec[2]=0.0;
     for (i = 0; i < dim; ++i) 
      {
       double z = (1 - x[i]) / x[i];
       fvec[2] += z * z;
       scale *= M_2_SQRTPI / (x[i] * x[i]);
      }
     fvec[2]= exp(-fvec[2]) * scale;

     /* fourth component of integrand: */
     fvec[3]=0.0;
     for (i = 0; i < dim; ++i) 
      fvec[3] += (x[i]-0.5)*(x[i]-0.5);
     fvec[3]= 1.0 / (fvec[3] + 1.0e-3);

}

void myvvfun(int nfun, int nquad, double *x, void *params, double *f)
{ 
  int nq;
  for (nq=0; nq<nquad; nq++)
   f_test(1,x+nq,0,4,f+4*nq);
} 

/* surface area of n-dimensional unit hypersphere */
static double S(unsigned n)
{
     double val;
     int fact = 1;
     if (n % 2 == 0) { /* n even */
	  val = 2 * pow(M_PI, n * 0.5);
	  n = n / 2;
	  while (n > 1) fact *= (n -= 1);
	  val /= fact;
     }
     else { /* n odd */
	  val = (1 << (n/2 + 1)) * pow(M_PI, n/2);
	  while (n > 2) fact *= (n -= 2);
	  val /= fact;
     }
     return val;
}

static double exact_integral(unsigned dim, const double *xmax, int which) {
     unsigned i;
     double val;
     switch(which) {
	 case 0:
	      val = dim == 0 ? 1 : S(dim) * pow(radius * 0.5, dim) / dim;
	      break;
	 case 1:
	      val = 1;
	      for (i = 0; i < dim; ++i)
		   val *= sin(xmax[i]);
	      break;
	 case 2:
	      val = 1;
	      break;
	 case 3:
              val=18.0478;
	      break;
	 default:
	      fprintf(stderr, "unknown integrand %d\n", which);
              exit(EXIT_FAILURE);
     }
     return val;
}

int main(int argc, char **argv)
{
     double *xmin, *xmax;
     double tol, val[4], err[4];
     unsigned i, nf, dim, maxEval;

     /* HR 9/07 just hard-code these parameters */
     dim=1;
     tol=1.0e-8;
     maxEval=0;

     xmin = (double *) malloc(dim * sizeof(double));
     xmax = (double *) malloc(dim * sizeof(double));
     for (i = 0; i < dim; ++i) {
	  xmin[i] = 0;
	  xmax[i] = 1.0;
     }

     printf("%u-dim integral, tolerance = %g\n", dim, tol);
int error_type[4];
vvgsl_integration_workspace *w=vvgsl_integration_workspace_alloc (10000, 4);
vvqag(myvvfun, 0, 4, 0.0, 1.0,
      0.0, tol, 10000, w, val, err, error_type);

// adapt_integrate(4, f_test, 0, dim, xmin, xmax, maxEval, 0, tol, val, err);

     for(nf=0; nf<4; nf++)
      { printf("integrand %i: \n",nf);
        printf(" integration val = %.15g, exact=%.15g \n",
                 val[nf], exact_integral(dim,xmax,nf));
      };

     printf("#evals = %d\n", count);

     free(xmax);
     free(xmin);

     return 0;
}
