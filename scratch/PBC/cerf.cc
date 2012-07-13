#include "libhrutil.h"

extern "C" 
 { 
   void wofz_(double *XI, double *YI, double *U, double *V, int *FLAG);
 };

cdouble cerfc(cdouble z)
{ 
  double zr, zi, wr, wi;
  int flag;

  /* handle some special cases that can cause the wofz routine to barf */
  if ( real(z)<-4.0 && fabs(imag(z))<1.0 )
   return cdouble(2.0,0.0);

  /* note the wofz_ function interprets its input as i*Z */
  zr = -imag(z);
  zi = real(z);

  wofz_(&zr, &zi, &wr, &wi, &flag);

  if (flag)
{ fprintf(stderr,"** warning: wofz(%e,%e) barfed\n",zr,zi);
   return 0.0;
};

  return exp(-z*z) * cdouble(wr, wi);
 
} 

cdouble cerf(cdouble z)
{ 
  return cdouble(1.0,0.0) - cerfc(z);
} 
