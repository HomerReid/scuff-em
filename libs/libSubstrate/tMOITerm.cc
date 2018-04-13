#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <complex>
#include <cmath>

typedef std::complex<double> cdouble;

#define II cdouble(0.0,1.0)

bool NeeddRho=true, Needdz=true;
bool PPIsOnly=false;
int NumSGFs=4;

#define _SGF_APAR  0
#define _SGF_PHI   1
#define _SGF_DRPHI 2
#define _SGF_DZPHI 3

cdouble AddSGFTerm_MOI(cdouble k, double Rho, double z, double h,
                       cdouble Eta, int n, cdouble EtaFac, cdouble *VVector)
{ 

  double zn = z + (n==0 ? 0.0 : 2.0*n*h);
  double r2 = Rho*Rho + zn*zn;
  if (r2==0.0)
   return 0.0;
  double r=sqrt(r2);
  cdouble ikr=II*k*r, ikr2=ikr*ikr, ikr3=ikr2*ikr;
  double Sign = (n%2) ? -1.0 : 1.0;
  cdouble ExpFac[4];
  ExpFac[0] = Sign*exp(ikr) / (4.0*M_PI*r);
  ExpFac[1] = ExpFac[0] * (ikr - 1.0)/r2;
  ExpFac[2] = ExpFac[0] * (ikr2 - 3.0*ikr + 3.0)/(r2*r2);
  ExpFac[3] = ExpFac[0] * (ikr3 - 6.0*ikr2 + 15.0*ikr - 15.0)/(r2*r2*r2);

  for(int dz=0, nvd=0; dz<=(Needdz?1:0); dz++)
   for(int dRho=0; dRho<=(NeeddRho?1:0); dRho++, nvd++)
    { 
      cdouble *V = VVector + nvd*NumSGFs;
      double zFac = (dz ? zn : 1.0), RhoFac = (dRho ? Rho : 1.0), zRho=zFac*RhoFac;
      int Index = dRho+dz;

      // contributions to V^A_parallel
      if (n<=1) V[_SGF_APAR] += zRho*ExpFac[Index];

      // contributions to V^Phi
      V[_SGF_PHI] += EtaFac*zRho*ExpFac[Index];

      if (PPIsOnly) continue;

      // contributions to V^Phi derivatives
      V[_SGF_DRPHI] += EtaFac*(zRho*Rho*ExpFac[Index+1] + dRho*zFac*ExpFac[Index]);
      V[_SGF_DZPHI] += EtaFac*(zRho*zn *ExpFac[Index+1] + dz*RhoFac*ExpFac[Index]);
    }

  return Eta*EtaFac;
}

int main(int argc, char *argv[])
{
  cdouble k=0.4;
  double Rho=0.3;
  double z=1.3;
  double h=0.9;

  cdouble VArray[16];

  for(int n=0; n<3; n++)
   { memset(VArray, 0, 16*sizeof(cdouble));
     AddSGFTerm_MOI(k, Rho, z, h, 1.0, 0, 1.0, VArray);
     printf("**n=%i: ",n);
     for(int dz=0, nvd=0; dz<=(Needdz?1:0); dz++)
      for(int dRho=0; dRho<=(NeeddRho?1:0); dRho++, nvd++)
       { cdouble *V = VArray + nvd*4;
         char Prefix[10]="";
         if (dRho) strcat(Prefix,"dR");
         if (dz) strcat(Prefix,"dz");
         for(int p=0; p<4; p++)
          printf("%sV[%i]={%e,%e}\n",Prefix,p,real(V[p]),imag(V[p])); 
         printf("\n\n");
       }
   }
}
