#include <stdio.h>
#include <math.h>
#include <complex>

typedef std::complex<double> cdouble;
#define II cdouble(0.0,1.0)

cdouble SphericalHankelH1(int l, cdouble z)
{ 
  if (l<0) return 0.0;

  cdouble ooz=1.0/z;
  cdouble ExpFac = ooz*exp(II*z);

  switch(l)
   {
     case 0:  return  -II*ExpFac;
     case 1:  return -1.0*ExpFac*(1.0 + II*ooz);
     case 2:  return   II*ExpFac*(1.0 + ooz*(3.0*II - 3.0*ooz));
     case 3:  return      ExpFac*(1.0 + ooz*(6.0*II + ooz*(-15.0 -15.0*II*ooz)));
     default: return (2.0*((double)l)-1.0)*SphericalHankelH1(l-1,z)/z - SphericalHankelH1(l-2,z);
   };
}

void Gethl(int l, cdouble z, cdouble *phl, cdouble *phlPrime)
{
  cdouble hl   = SphericalHankelH1(l, z);
  cdouble hlM1 = SphericalHankelH1(l-1, z);
  *phl = hl;
  *phlPrime = hlM1 - ( double(l) + 1.0 )*hl/z;
}

int main(int argc, char *argv[])
{ 
  cdouble z;
  sscanf(argv[1], "%le", &( real(z) ));
  sscanf(argv[2], "%le", &( imag(z) ));

  cdouble hl, hlPrime;
  for(int l=1; l<6; l++)
   { Gethl(l, z, &hl, &hlPrime);
     printf("%i: (%+9.2e,%+9.2e), (%+9.2e,%+9.2e)\n",
             l,real(hl),imag(hl),real(hlPrime),imag(hlPrime));
   };

}
