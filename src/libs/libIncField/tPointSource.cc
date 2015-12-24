/*
 * tPointSource.cc 
 *
 * This program prints the E-field of a 
 */
#include <stdio.h>
#include <stdlib.h>

#include <libhrutil.h>
#include <libIncField.h>

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
typedef struct SummandData
 {
   double *X, *X0;
   double kBloch;
   PointSource *PS;
 } SummandData

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
void MySummand(double *L, void *UserData, double *Sum)
{
  SummandData *Data = (SummandData *UserData);
  double *X         = Data->X;
  double *X0        = Data->X0;
  double *kBloch    = Data->kBloch;
  int LDim          = Data->LDim;
 
  double KDotL = 0.0;
  for(int nd=0; nd<LDim; nd++)
   KDotL += kBloch[nd] * L[nd];
  cdouble BlochPhase = exp(II*KDotL);

  double X0pL[3];
  X0pL[0]=X0[0] + L[0];
  X0pL[1]=X0[1] + L[1];
  X0pL[2]=X0[2] + L[2];
  PS->SetX(X0pL);

  cdouble EH[6];
  PS->GetFields(X0, EH);

  cdouble *zSum = (cdouble *)Sum;
  VecPlusEquals(zSum, BlochPhase, EH);

}

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble Omega=0.01;              // angular frequency
  cdouble EpsRel=1.0;
  cdouble MuRel=1.0;
  double X[3]={1.0, 0.0, 0.0};     // field evaluation point
  double X0[3]={0.0, 0.0, 0.0};    // point source location
  cdouble P[3]={0.0, 0.0, 1.0};    // point source location
  char *BasisFile=0;
  double kx=0.0;   int nkx=0;
  double ky=0.0;   int nky=0;
  double kz=0.0;   int nkz=0;
  double RelTol=1.0e-3;
  int MaxCells=1000;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"Omega",    PA_CDOUBLE, 1, 1, (void *)&Omega,   0, "angular frequency"},
     {"Eps",      PA_CDOUBLE, 1, 1, (void *)&Eps,     0, "material permittivitty"},
     {"Mu",       PA_CDOUBLE, 1, 1, (void *)&Mu,      0, "material permeability"},
     {"X",        PA_DOUBLE,  3, 1, (void *)X,        0, "evaluation point "},
     {"X0",       PA_DOUBLE,  3, 1, (void *)X0,       0, "source point"},
     {"kx",       PA_DOUBLE,  1, 1, (void *)&kx,    &nkx, "x component of Bloch vector"},
     {"ky",       PA_DOUBLE,  1, 1, (void *)&ky,    &nky, "y component of Bloch vector"},
     {"kz",       PA_DOUBLE,  1, 1, (void *)&kz,    &nkz, "z component of Bloch vector"},
     {"RelTol",   PA_DOUBLE,  1, 1, (void *)&RelTol,   0, ""},
     {"MaxCells", PA_INT,     1, 1, (void *)&MaxCells, 0, ""},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  HMatrix *LBasis=0;
  double kBloch[3];
  if (BasisFile)
   LBasis=new HMatrix(BasisFile);
  if (nkx || nky || nkz)
   { kBloch[0]=kx;
     kBloch[1]=ky;
     kBloch[2]=kz;
     int LDim = (nkz>0 ? 3 : nky>0 ? 2 : 1);
     LBasis=new HMatrix(3,LDim);
     for(int nd=0; nd<LDim; nd++)
      LBasis->SetEntry(nd,nd,1.0);
   };
  
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  PointSource *PS = new PointSource(X0, P);

  cdouble EH[6];
  PS->SetFrequencyAndEpsMu(Omega, EpsRel, MuRel, false);
  PS->GetFields(X, EH);
  if (LBasis)
   { PS->SetBasis(LBasis);
     PS->SetkBloch(kBloch);
   };

  printf("Ex=(%+10.2e,%+10.2e)\n",real(EH[0]),imag(EH[0]));
  printf("Ey=(%+10.2e,%+10.2e)\n",real(EH[1]),imag(EH[1]));
  printf("Ez=(%+10.2e,%+10.2e)\n",real(EH[2]),imag(EH[2]));
  printf("Hx=(%+10.2e,%+10.2e)\n",real(EH[3]),imag(EH[3]));
  printf("Hy=(%+10.2e,%+10.2e)\n",real(EH[4]),imag(EH[4]));
  printf("Hz=(%+10.2e,%+10.2e)\n",real(EH[5]),imag(EH[5]));

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (LBasis && Test)
   { 
      cdouble EHBF[6];
      struct SummandData MyData, *Data=&MyData;
      Data->X=X;
      Data->X0=X0;
      Data->kBloch=kBloch;
      Data->PS = new PointSource(X0, P);
      Data->PS->SetFrequencyAndEpsMu(Omega, EpsRel, MuRel, false);
      GetLatticeSum(MySummand, (void *)Data, 12, LBasis, 
                    (double *)EHBF, 0.0, RelTol, MaxCells);
      Compare(EH, EHBF, 6, "Ewald", "BF");
   };

}
