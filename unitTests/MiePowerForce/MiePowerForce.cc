/*
 * MiePowerForce.cc
 *
 * Homer Reid    -- 1/2013
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <libscuff.h>
#include <libhrutil.h>
#include <libSGJC.h>
#include <libSpherical.h>

using namespace scuff;

#define II cdouble(0.0,1.0)

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetMieSurfaceCurrents(double Theta, double Phi,
                           cdouble Omega, cdouble Epsilon, cdouble Mu, 
                           int nMax, cdouble K[3], cdouble N[3]);

void GetMieCrossSections(double Omega, cdouble Epsilon, cdouble Mu, int nMax, 
                         double *SigmaScat, double *SigmaTot, double *SigmaForce);

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WriteST(RWGSurface *S, RWGPanel *P, double Value, FILE *f)
{
  double *VV[3];
  VV[0] = S->Vertices + (3*P->VI[0]);
  VV[1] = S->Vertices + (3*P->VI[1]);
  VV[2] = S->Vertices + (3*P->VI[2]);
 
  fprintf(f,"ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n",
             VV[0][0], VV[0][1], VV[0][2],
             VV[1][0], VV[1][1], VV[1][2],
             VV[2][0], VV[2][1], VV[2][2], 
             Value,Value,Value);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int nCalls=0;
void GetATPowerDensity(double Theta, double Phi, 
                       cdouble Omega, cdouble Epsilon, cdouble Mu, 
                       int nMax, double *pAPower, double *pTPower)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble KS[3], NS[3];
  GetMieSurfaceCurrents(Theta,Phi,Omega,Epsilon,Mu,nMax,KS,NS);

  /***************************************************************/
  /* incident fields *********************************************/
  /***************************************************************/
  cdouble EC[3], HC[3];
  cdouble ExpFac = exp(II*Omega*cos(Theta));
  EC[0]=ExpFac; EC[1]=0.0;         EC[2]=0.0;
  HC[0]=0.0;    HC[1]=ExpFac/ZVAC; HC[2]=0.0;
  cdouble ES[3], HS[3];
  VectorC2S(Theta, Phi, EC, ES);
  VectorC2S(Theta, Phi, HC, HS);

  /***************************************************************/
  /* Absorbed power = (1/2) K^* \dot nHat \times N               */
  /***************************************************************/
  *pAPower = 0.5*real( -conj(KS[1])*NS[2] + conj(KS[2])*NS[1] );

  /***************************************************************/
  /* Total power                                                 */
  /***************************************************************/
  *pTPower = 0.5*real( conj(ES[1])*KS[1] + conj(ES[2])*KS[2]
                      +conj(HS[1])*NS[1] + conj(HS[2])*NS[2] );

  nCalls++;

} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct IntegrandData 
 { 
   cdouble Omega;
   cdouble Epsilon, Mu;
   int nMax;
 } IntegrandData; 
 
void Integrand(unsigned ndim, const double *x, void *params,
               unsigned fdim, double *fval)
{
  double Theta = acos(x[0]);
  double Phi   = x[1];
 
  IntegrandData *ID = (IntegrandData *)params;
  cdouble Omega   = ID->Omega;
  cdouble Epsilon = ID->Epsilon;
  cdouble Mu      = ID->Mu; 
  int nMax        = ID->nMax;

  GetATPowerDensity(Theta, Phi, Omega, Epsilon, Mu, nMax, fval+0, fval+1);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetATPower(cdouble Omega, cdouble Epsilon, cdouble Mu, int nMax,
                double ATPower[2], double Error[2])
{
  IntegrandData MyID, *ID = &MyID;
  ID->Omega               = Omega;
  ID->Epsilon             = Epsilon;
  ID->Mu                  = Mu;
  ID->nMax                = nMax;

  double Lower[2] = {-1.0, 0.0      };
  double Upper[2] = { 1.0, 2.0*M_PI };

  adapt_integrate(2, Integrand, (void *)ID, 2, Lower, Upper, 
                  100000, 0.0, 1.0e-2, ATPower, Error);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  /*--------------------------------------------------------------*/
  /*- process command-line options -------------------------------*/
  /*--------------------------------------------------------------*/
  double Omega=1.0;
  int nMax=10;
  char *GeoFileName=0;
  char *Material=0;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"Material", PA_STRING,  1, 1, (void *)&Material,    0,  "Material designation"},
     {"nMax",     PA_INT,     1, 1, (void *)&nMax,        0,  "nMax"},
     {"Omega",    PA_DOUBLE,  1, 1, (void *)&Omega,       0,  "angular frequency"},
     {"geometry", PA_STRING,  1, 1, (void *)&GeoFileName, 0,  "geometry"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  MatProp *MP;
  if (Material)
   { MP=new MatProp(Material);
     if (MP->ErrMsg)
      ErrExit(MP->ErrMsg);
   }
  else
   MP = new MatProp(MP_PEC);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if(GeoFileName)
   { RWGGeometry *G=new RWGGeometry(GeoFileName);
     RWGSurface  *S=G->Surfaces[0];
     cdouble Epsilon, Mu; 
     MP->GetEpsMu(Omega,&Epsilon,&Mu);
     FILE *f=vfopen("MiePower.E%s.w%g.pp","w",z2s(Epsilon),Omega);
     fprintf(f,"View \"%s\" {\n","Absorbed power density");
     for (int np=0; np<S->NumPanels; np++)
      { 
        double r, Theta, Phi;
        RWGPanel *P=S->Panels[np];
        CoordinateC2S(P->Centroid, &r, &Theta, &Phi);

        double APower, TPower;
        GetATPowerDensity(Theta, Phi, Omega, Epsilon, Mu, nMax, &APower, &TPower);

        WriteST(S,P,APower,f);
 
      };
     fprintf(f,"};\n");
     fclose(f);
   }
  else
   { 
     MP->SetFreqUnit(3.0e14);
     FILE *f=vfopen("%s.MPF","w",MP->Name);
     fprintf(f,"# data file columns: \n");
     fprintf(f,"#  1  Omega \n");
     fprintf(f,"#  2  re(Epsilon)\n");
     fprintf(f,"#  3  im(Epsilon)\n");
     fprintf(f,"#  4  nCalls\n");
     fprintf(f,"#  5  QScat  (Mie series)\n");
     fprintf(f,"#  6  QTot   (Mie series) \n");
     fprintf(f,"#  7  QAbs   (Mie series) \n");
     fprintf(f,"#  8  QForce (Mie series)\n");
     fprintf(f,"#  9  QAbs   (numerical integration)\n");
     fprintf(f,"#  10 QAbs   (numerical integration error)\n");
     fprintf(f,"#  11 QTot   (numerical integration)\n");
     fprintf(f,"#  12 QTot   (numerical integration error)\n");
     for(Omega=0.001; Omega<=10.1; Omega*=exp(0.05*log(10.0)) )
      { 
        cdouble Epsilon, Mu; 
        MP->GetEpsMu(Omega,&Epsilon,&Mu);
        double ATPower[2], Error[2], SigmaScat, SigmaTot, SigmaForce;
        nCalls=0;
        GetATPower(Omega, Epsilon, Mu, nMax, ATPower, Error);
        GetMieCrossSections(Omega, Epsilon, Mu, nMax, &SigmaScat, &SigmaTot, &SigmaForce);

        double IncFlux=1.0/(2.0*376.73031346177);
        double QScatSeries     = SigmaScat / M_PI;
        double QTotSeries      = SigmaTot  / M_PI;
        double QAbsSeries      = (SigmaTot-SigmaScat)/ M_PI;
        double QForceSeries    = SigmaForce/ M_PI;
        double QAbsNumInt      = ATPower[0] / (M_PI*IncFlux);
        double QAbsNumIntError = Error[0] / (M_PI*IncFlux);
        double QTotNumInt      = ATPower[1] / (M_PI*IncFlux);
        double QTotNumIntError = Error[1] / (M_PI*IncFlux);

        fprintf(f,"%e %e %e %i ",Omega,real(Epsilon),imag(Epsilon),nCalls);
        fprintf(f,"%e %e %e %e ", QScatSeries, QTotSeries, QAbsSeries, QForceSeries);
        fprintf(f,"%e %e %e %e \n",QAbsNumInt, QAbsNumIntError, QTotNumInt, QTotNumIntError);
        fflush(f);

      };
     fclose(f);
   };

}
