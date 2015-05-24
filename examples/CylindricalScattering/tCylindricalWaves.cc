/* Copyright (C) 2005-2011 M. T. Homer ReidElectricOnly=true;
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

#include <stdio.h>
#include <stdlib.h>

#include <libscuff.h>
#include "CylindricalWave.h"

using namespace scuff;
using namespace std;

#define II cdouble(0.0,1.0) 

void GetCylMN(cdouble k0, int Nu, double kz, int WaveType,
              double RPZ[3], cdouble MVec[3], cdouble NVec[3]);

void CoordinateCar2Cyl(double XYZ[3], double RPZ[3]);

void VectorCyl2Car(double Phi, cdouble VCyl[3], cdouble VCar[3]);

/***************************************************************/
/***************************************************************/
/***************************************************************/
void Compare(cdouble *V1, cdouble *V2, int N,
             const char *str1, const char *str2)
{ 
  //char FStr[10];
  //snprintf(FStr,10,"%+.%ie",Precision);

  printf(" n | %-25s | %-25s | RD      | Ratio\n",str1,str2);
  for(int n=0; n<N; n++)
   printf("%2i | (%+.4e,%+.4e) | (%+.4e,%+.4e) | %.1e | %.3e\n",n,
    real(V1[n]),imag(V1[n]), real(V2[n]),imag(V2[n]),
    RD(V1[n],V2[n]), abs(V1[n]/V2[n]));
  printf("\n");
}

typedef void (*vvfun)(double X[3], cdouble F[3], void *UserData);

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetFAndCurlF(vvfun Func, double X[3], void *UserData,
                  cdouble F[3], cdouble CurlF[3])
{ 
  cdouble Fp[3], Fm[3], dF[3][3];
  
  for(int Mu=0; Mu<3; Mu++)
   { 
     double XX[3];
     XX[0]=X[0]; 
     XX[1]=X[1];
     XX[2]=X[2];

     double Delta = (X[Mu]==0.0) ? 1.0e-4 : ( 1.0e-4 * fabs(X[Mu]) );
     XX[Mu] += Delta;
     Func(XX, Fp, UserData);
     XX[Mu] -= 2.0*Delta;
     Func(XX, Fm, UserData);

     dF[Mu][0] = (Fp[0] - Fm[0]) / (2.0*Delta);
     dF[Mu][1] = (Fp[1] - Fm[1]) / (2.0*Delta);
     dF[Mu][2] = (Fp[2] - Fm[2]) / (2.0*Delta);
     
   };

  Func(X, F, UserData);

  CurlF[0] = dF[1][2] - dF[2][1];
  CurlF[1] = dF[2][0] - dF[0][2];
  CurlF[2] = dF[0][1] - dF[1][0];

}

typedef struct MyData
 { cdouble k0; 
   int Nu;
   double kz; 
   int WaveType;
 } MyData;

void MFunc(double X[3], cdouble F[3], void *UserData)
{ 
  cdouble M[3], N[3];
  MyData *Data = (MyData *) UserData;
  GetCylMN(Data->k0, Data->Nu, Data->kz, Data->WaveType, X, M, N);

  double RPZ[3];
  CoordinateCar2Cyl(X, RPZ);
  VectorCyl2Car(RPZ[1], M, F);
}

void NFunc(double X[3], cdouble F[3], void *UserData)
{
  cdouble M[3], N[3];
  MyData *Data=(MyData *)UserData;
  GetCylMN(Data->k0, Data->Nu, Data->kz, Data->WaveType, X, M, N);

  double RPZ[3];
  CoordinateCar2Cyl(X, RPZ);
  VectorCyl2Car(RPZ[1], N, F);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble Omega=1.0;
  double kz=0.0;
  double XYZ[3];

  srand48(time(0));
  kz = drand48();
  XYZ[0] = 2.0*( drand48() - 0.5 );
  XYZ[1] = 2.0*( drand48() - 0.5 );
  XYZ[2] = 2.0*( drand48() - 0.5 );

  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"Omega",          PA_CDOUBLE, 1, 1, (void *)&Omega, 0, "angular frequency"},
     {"kz",             PA_DOUBLE,  1, 1, (void *)&kz,    0, "z wavenumber"},
     {"xyz",            PA_DOUBLE,  3, 1, (void *)XYZ,    0, "eval point"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  cdouble IK0=II*Omega;
  
  printf("--Omega %s ",z2s(Omega));
  printf("--kz    %g ",kz);
  printf("--xyz   %g %g %g ",XYZ[0],XYZ[1],XYZ[2]);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  MyData DataBuffer, *Data=&DataBuffer;
  SetDefaultCD2SFormat("{%+.2e,%+.2e}");
  for(int Nu=0; Nu<3; Nu++)
   for(int WaveType=CW_REGULAR; WaveType<=2; WaveType++)
    {  
      Data->k0=Omega;
      Data->kz=kz;
      Data->Nu=Nu;
      Data->WaveType=WaveType;

      cdouble M[3], CurlM[3], N[3], CurlN[3];
      GetFAndCurlF(MFunc, XYZ, (void *)Data, M, CurlM);
      GetFAndCurlF(NFunc, XYZ, (void *)Data, N, CurlN);

      cdouble k0=Omega;
      double IKMMag   = norm(k0*M[0]) + norm(k0*M[1]) + norm(k0*M[2]);
      double IKNMag   = norm(k0*N[0]) + norm(k0*N[1]) + norm(k0*N[2]);
      double CurlMMag = norm(CurlM[0]) + norm(CurlM[1]) + norm(CurlM[2]);
      double CurlNMag = norm(CurlN[0]) + norm(CurlN[1]) + norm(CurlN[2]);

      printf("\n**\n** (Nu,Wave)=(%i,%i):\n**\n\n",Nu,WaveType);

      N[0] *= -1.0*IK0;
      N[1] *= -1.0*IK0;
      N[2] *= -1.0*IK0;
//      Compare(CurlM, N, 3, "CurlM","-ikN");
//      printf("\n");
      printf("(Nu,Wave)=(%i,%i): Mag2(CurlM, -IKN)=(%.2e,%.2e) (ratio %.2e)\n",
              Nu,WaveType,CurlMMag, IKNMag, CurlMMag / IKNMag);
//      printf("\n");

      M[0] *= 1.0*IK0;
      M[1] *= 1.0*IK0;
      M[2] *= 1.0*IK0;
//      Compare(CurlN, M, 3, "CurlN","+ikM");
//      printf("\n");
      printf("(Nu,Wave)=(%i,%i): Mag2(CurlN, +IKM)=(%.2e,%.2e) (ratio %.2e)\n",
              Nu,WaveType,CurlNMag, IKMMag, CurlNMag / IKMMag);
//      printf("\n");

    };

}
