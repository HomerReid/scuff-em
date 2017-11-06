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

#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>
#include <unistd.h>

#include "libhrutil.h"
#include "libSubstrate.h"
#include "libscuff.h"

// infinite silicon half-space
const char SISubstrateFile[]=
 "0.0 CONST_EPS_11.7\n";

using namespace scuff;

#define II cdouble(0.0,1.0)

/***************************************************************/
namespace scuff{
void CalcGC(double R[3], cdouble Omega,
            cdouble EpsR, cdouble MuR,
            cdouble GMuNu[3][3], cdouble CMuNu[3][3],
            cdouble GMuNuRho[3][3][3], cdouble CMuNuRho[3][3][3]);
}

void Getg0112(LayeredSubstrate *Substrate, cdouble Omega, double qMag,
              double zDest, double zSource, HMatrix *RTwiddle,
              HMatrix *WMatrix, HMatrix *STwiddle, HMatrix *g0112[4]);

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddScriptG0(LayeredSubstrate *S, cdouble Omega,
                 HMatrix *XMatrix, HMatrix *GMatrix)
{
  for(int nx=0; nx<XMatrix->NR; nx++)
   { 
     double XD[6], *XS=XD+3;
     XMatrix->GetEntriesD(nx, "0:5", XD);
     int nrDest   = S->GetRegionIndex(XD[2]);
     int nrSource = S->GetRegionIndex(XS[2]);
     if (nrDest!=nrSource) continue;

     double R[3];
     VecSub(XD, XS, R);
     cdouble EpsRel, MuRel;
     S->MPLayer[nrDest]->GetEpsMu(Omega, &EpsRel, &MuRel);
     cdouble k=sqrt(EpsRel*MuRel)*Omega;

     if (S->StaticLimit)
      { double r2=VecNorm2(R), r=sqrt(r2);
        cdouble G0 = II*exp(II*k*r)/(4.0*M_PI*r*r*r);
        cdouble EFactor = ZVAC*MuRel*G0/k;
        cdouble HFactor = G0/(k*ZVAC*EpsRel);
        for(int Mu=0; Mu<3; Mu++)
         for(int Nu=0; Nu<3; Nu++)
          { double Bracket = 3.0*R[Mu]*R[Nu]/r2 - ((Mu==Nu) ? 1.0 : 0.0);
            GMatrix->AddEntry(6*(0+Mu) + 0*Nu, nx, EFactor*Bracket);
            GMatrix->AddEntry(6*(3+Mu) + 3*Nu, nx, HFactor*Bracket);
          };

        if (nrDest==S->NumInterfaces && !isinf(S->zGP))
         { R[2] = XD[2] + XS[2] - 2.0*S->zGP;
           r2=VecNorm2(R), r=sqrt(r2);
           G0 = II*exp(II*k*r)/(4.0*M_PI*r*r*r);
           EFactor = ZVAC*MuRel*G0/k;
           HFactor = G0/(k*ZVAC*EpsRel);
           for(int Mu=0; Mu<3; Mu++)
            for(int Nu=0; Nu<3; Nu++)
             { double Sign = (Nu==2) ? -1.0 : 1.0;
               double Bracket = 3.0*R[Mu]*R[Nu]/r2 - ((Mu==Nu) ? 1.0 : 0.0);
               GMatrix->AddEntry(6*(0+Mu) + 0*Nu, nx, +1.0*Sign*EFactor*Bracket);
               GMatrix->AddEntry(6*(3+Mu) + 3*Nu, nx, -1.0*Sign*HFactor*Bracket);
             };
         };
        continue;
      };

     cdouble G[3][3], C[3][3], dG[3][3][3], dC[3][3][3];
     scuff::CalcGC(R, Omega, EpsRel, MuRel, G, C, dG, dC);
     cdouble EEPreFac = II*Omega*ZVAC*MuRel;
     cdouble EMPreFac = II*k;
     cdouble MEPreFac = -1.0*II*k;
     cdouble MMPreFac = II*Omega*EpsRel/ZVAC;
     for(int Mu=0; Mu<3; Mu++)
      for(int Nu=0; Nu<3; Nu++)
       { GMatrix->AddEntry(6*(0+Mu) + 0+Nu, nx, EEPreFac*G[Mu][Nu]);
         GMatrix->AddEntry(6*(0+Mu) + 3+Nu, nx, EMPreFac*C[Mu][Nu]);
         GMatrix->AddEntry(6*(3+Mu) + 0+Nu, nx, MEPreFac*C[Mu][Nu]);
         GMatrix->AddEntry(6*(3+Mu) + 3+Nu, nx, MMPreFac*G[Mu][Nu]);
       };

     if (nrDest==S->NumInterfaces && !isinf(S->zGP))
      { R[2] = XD[2] + XS[2] - 2.0*S->zGP;
        scuff::CalcGC(R, Omega, EpsRel, MuRel, G, C, dG, dC);
        const double ImageSign[3]={-1.0, -1.0, +1.0};
        for(int Mu=0; Mu<3; Mu++)
         for(int Nu=0; Nu<3; Nu++)
          { GMatrix->AddEntry(6*(0+Mu) + 0+Nu, nx, +1.0*ImageSign[Nu]*EEPreFac*G[Mu][Nu]);
            GMatrix->AddEntry(6*(0+Mu) + 3+Nu, nx, -1.0*ImageSign[Nu]*EMPreFac*C[Mu][Nu]);
            GMatrix->AddEntry(6*(3+Mu) + 0+Nu, nx, +1.0*ImageSign[Nu]*MEPreFac*C[Mu][Nu]);
            GMatrix->AddEntry(6*(3+Mu) + 3+Nu, nx, -1.0*ImageSign[Nu]*MMPreFac*G[Mu][Nu]);
          };
      };

   };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct qFunctionSCData
 {
   HMatrix *XMatrix;
   HMatrix *RTwiddle;
   HMatrix *WMatrix;
   HMatrix *STwiddle;
   FILE *byqFile;
 } qFunctionSCData;

void qFunctionSC1D(LayeredSubstrate *Substrate,
                   double q2D[2], cdouble Omega,
                   void *UserData, cdouble *Integrand)
{
  qFunctionSCData *Data = (qFunctionSCData *)UserData;

  HMatrix *XMatrix    = Data->XMatrix;
  HMatrix *RTwiddle   = Data->RTwiddle;
  HMatrix *WMatrix    = Data->WMatrix;
  HMatrix *STwiddle   = Data->STwiddle;
  FILE *byqFile       = Data->byqFile;

  int EntryOnly       = Substrate->EntryOnly;
  bool EEOnly         = Substrate->EEOnly;
  bool XYOnly         = Substrate->XYOnly;
  int NumInterfaces   = Substrate->NumInterfaces;
  double *zInterface  = Substrate->zInterface;
  cdouble *EpsLayer   = Substrate->EpsLayer;

  double qMag = q2D[0];

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int nx=0, ni=0; nx<XMatrix->NR; nx++)
   { 
     double Rho[2];
     Rho[0]          = XMatrix->GetEntryD(nx,0) - XMatrix->GetEntryD(nx,3);
     Rho[1]          = XMatrix->GetEntryD(nx,1) - XMatrix->GetEntryD(nx,4);
     double RhoMag   = sqrt( Rho[0]*Rho[0] + Rho[1]*Rho[1] );
     double CosTheta = (RhoMag==0.0) ? 0.0 : Rho[0] / RhoMag;
     double SinTheta = (RhoMag==0.0) ? 0.0 : Rho[1] / RhoMag;
     double zDest    = XMatrix->GetEntryD(nx,2);
     double zSource  = XMatrix->GetEntryD(nx,5);

     // fetch scalar quantities
     cdouble gBuf[4*36];
     HMatrix g0( 6, 6, LHM_COMPLEX, gBuf +  0  );
     HMatrix g1x(6, 6, LHM_COMPLEX, gBuf + 36  );
     HMatrix g1y(6, 6, LHM_COMPLEX, gBuf + 72  );
     HMatrix g2( 6, 6, LHM_COMPLEX, gBuf + 108 );
     HMatrix *g0112[4];
     g0112[0]=&g0;
     g0112[1]=&g1x;
     g0112[2]=&g1y;
     g0112[3]=&g2;
     Getg0112(Substrate, Omega, qMag, zDest, zSource, RTwiddle, WMatrix, STwiddle, g0112);

     // fetch bessel functions
     double TT=Secs();
bool NeedRhoDerivatives=false;
     cdouble J[4], dJ[4];
     GetJFactors(qMag, RhoMag, J, NeedRhoDerivatives ? dJ : 0);
     Substrate->Times[BESSELTIME] += (Secs()-TT);
  
cdouble J0, J1oqRho, J1[2], J2[2][2]; 
  J0       = J[0];
  J1[0]    = II*CosTheta*J[1];
  J1[1]    = II*SinTheta*J[1];
  J1oqRho  = J[2];
  J2[0][0] = -1.0*CosTheta*CosTheta*J[3] + J1oqRho;
  J2[0][1] = -1.0*CosTheta*SinTheta*J[3];
  J2[1][0] = -1.0*SinTheta*CosTheta*J[3];
  J2[1][1] = -1.0*SinTheta*SinTheta*J[3] + J1oqRho;
#if 0
dJ0, dJ1oqRho, dJ1[2], dJ2[2][2];
  dJData->J0       = dJ[0];
  dJData->J1oqRho  = dJ1oqRho;
  dJData->J1[0]    = II*CosTheta*dJ[1];
  dJData->J1[1]    = II*SinTheta*dJ[1];
  dJData->J2[0][0] = -1.0*CosTheta*CosTheta*dJ[2] + dJ1oqRho;
  dJData->J2[0][1] = -1.0*CosTheta*SinTheta*dJ[2];
  dJData->J2[1][0] = -1.0*SinTheta*CosTheta*dJ[2];
  dJData->J2[1][1] = -1.0*SinTheta*SinTheta*dJ[2] + dJ1oqRho;
#endif

     // check for zDest, zSource both on a layer boundary
     // note that this calculation assumes Rho=(x,0))
     cdouble Gxx0=0.0, Gyy0=0.0;
     bool OnInterface=false;
     if ( EqualFloat(zDest, zSource) )
      for(int nz=0; nz<NumInterfaces && !OnInterface; nz++)
       if ( fabs(zDest-zInterface[nz]) < 1.0e-6 )
        { cdouble EpsA = EpsLayer[nz], EpsB=EpsLayer[nz+1];
          cdouble Factor = (-0.5*II*ZVAC/Omega)*(EpsA-EpsB)/(EpsA+EpsB);
          Gxx0 = Factor*qMag*(J0 - J1oqRho);
          Gyy0 = Factor*qMag*J1oqRho;
          OnInterface=true;
        };
     
     // assemble GTwiddle
     cdouble GTwiddle[6][6];
     memset( (cdouble *)GTwiddle, 0, 36*sizeof(cdouble));
     int pqMax = (EEOnly ? 1 : 2);
     int ijMax = (XYOnly ? 2 : 3);
     for(int p=0; p<pqMax; p++)
      for(int q=0; q<pqMax; q++)
       for(int i=0; i<ijMax; i++)
        for(int j=0; j<ijMax; j++)
         { 
           int Mu=3*p+i, Nu=3*q+j;
           if (EntryOnly!=-1 && (EntryOnly!=(6*Mu+Nu)) )
            continue;
            
           cdouble g0MuNu  = g0. GetEntry(Mu,Nu);
           cdouble g1xMuNu = g1x.GetEntry(Mu,Nu);
           cdouble g1yMuNu = g1y.GetEntry(Mu,Nu);
           cdouble g2MuNu  = g2. GetEntry(Mu,Nu);

           GTwiddle[Mu][Nu]   = g0MuNu * J0;

           if (i<2 && j<2)
            GTwiddle[Mu][Nu] += g2MuNu * J2[i][j];
           else if ( (i<2 && j==2) || (i==2 && j<2) )
            GTwiddle[Mu][Nu] += (g1xMuNu*J1[0] + g1yMuNu*J1[1]);
         };

     if (OnInterface)
      { GTwiddle[0][0] -= Gxx0;
        GTwiddle[1][1] -= Gyy0;
      };

     // stamp into integrand vector
     for(int Mu=0; Mu<6; Mu++)
      for(int Nu=0; Nu<6; Nu++)
       Integrand[ni++] = GTwiddle[Mu][Nu];

     Substrate->Times[STAMPTIME]+=Secs() - TT;

     if (byqFile)
      { fprintf(byqFile,"%e ",qMag);
        fprintf(byqFile,"%e %e %e %e ",Rho[0],Rho[1],zDest,zSource);
 //      fprintf(byqFile,"%e %e ",real(J[0]),imag(J[0]));
 //       fprintf(byqFile,"%e %e ",real(J1oqRho),imag(J1oqRho));
 //       fprintf(byqFile,"%e %e ",real(J[2]),imag(J[2]));
        fprintVec(byqFile,(cdouble *)GTwiddle, 36);
        fprintf(byqFile,"%e %e %e %e ",real(Gxx0), imag(Gxx0), real(Gyy0), imag(Gyy0));
        fprintf(byqFile,"\n");
        fflush(byqFile);
      };

   }; // for(int nx=0, ni=0; nx<XMatrix->NR; nx++)

}

void qFunctionSC2D(LayeredSubstrate *Substrate,
                   double q2D[2], cdouble Omega,
                   void *UserData, cdouble *Integrand)
{
  qFunctionSCData *Data = (qFunctionSCData *)UserData;

  HMatrix *XMatrix    = Data->XMatrix;
  HMatrix *RTwiddle   = Data->RTwiddle;
  HMatrix *WMatrix    = Data->WMatrix;
  HMatrix *STwiddle   = Data->STwiddle;
  FILE *byqFile       = Data->byqFile;

  int EntryOnly       = Substrate->EntryOnly;
  bool EEOnly         = Substrate->EEOnly;
  bool XYOnly         = Substrate->XYOnly;
  int NumInterfaces   = Substrate->NumInterfaces;
  double *zInterface  = Substrate->zInterface;
  cdouble *EpsLayer   = Substrate->EpsLayer;

  cdouble gBuf[36];
  HMatrix GMatrix(6,6,LHM_COMPLEX,gBuf);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int nx=0, ni=0; nx<XMatrix->NR; nx++)
   { 
     double Rho[2];
     Rho[0]         = XMatrix->GetEntryD(nx,0) - XMatrix->GetEntryD(nx,3);
     Rho[1]         = XMatrix->GetEntryD(nx,1) - XMatrix->GetEntryD(nx,4);
     double zDest   = XMatrix->GetEntryD(nx,2);
     double zSource = XMatrix->GetEntryD(nx,5);
     Substrate->GetScriptGTwiddle(Omega, q2D, zDest, zSource,
                                  RTwiddle, WMatrix, STwiddle,
                                  &GMatrix); 

     cdouble ExpFac=exp(II*(q2D[0]*Rho[0] + q2D[1]*Rho[1]));
     for(int Mu=0; Mu<6; Mu++)
      for(int Nu=0; Nu<6; Nu++)
       Integrand[ni++] = ExpFac*GMatrix.GetEntry(Mu,Nu);

     if (byqFile)
      { fprintf(byqFile,"%e %e ",q2D[0],q2D[1]);
        fprintf(byqFile,"%e %e %e %e ",Rho[0],Rho[1],zDest,zSource);
 //      fprintf(byqFile,"%e %e ",real(J[0]),imag(J[0]));
 //       fprintf(byqFile,"%e %e ",real(J1oqRho),imag(J1oqRho));
 //       fprintf(byqFile,"%e %e ",real(J[2]),imag(J[2]));
        fprintVec(byqFile,GMatrix.ZM, 36);
        fprintf(byqFile,"\n");
        fflush(byqFile);
      };

   }; // for(int nx=0, ni=0; nx<XMatrix->NR; nx++)

}


/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetSubstrateDGF_SC(LayeredSubstrate *Substrate,
                        cdouble Omega,
                        HMatrix *XMatrix,
                        HMatrix *GMatrix,
                        bool Force2D=false)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int NX = XMatrix->NR;
  int FDim = 36*NX;

  UpdateCachedEpsMu(Omega);

  HMatrix *RTwiddle  = new HMatrix(6,               4*NumInterfaces, LHM_COMPLEX);
  HMatrix *WMatrix   = new HMatrix(4*NumInterfaces, 4*NumInterfaces, LHM_COMPLEX);
  HMatrix *STwiddle  = new HMatrix(4*NumInterfaces, 6,               LHM_COMPLEX);

  qFunctionSCData MyData, *Data=&MyData;
  Data->XMatrix    = XMatrix;
  Data->RTwiddle   = RTwiddle;
  Data->WMatrix    = WMatrix;
  Data->STwiddle   = STwiddle;
  Data->byqFile    = fopen("/tmp/qIntegral.log","w");

  if (Force2D)
   qIntegrate(Omega, qFunctionSC2D, (void *)Data, GMatrix->ZM, FDim, false);
  else
   qIntegrate(Omega, qFunctionSC1D, (void *)Data, GMatrix->ZM, FDim);
  
  delete RTwiddle;
  delete WMatrix;
  delete STwiddle;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  InstallHRSignalHandler();
  InitializeLog(argv[0]);

  /*--------------------------------------------------------------*/
  /*- process command-line arguments -----------------------------*/
  /*--------------------------------------------------------------*/
  char *GeoFile=0;
  char *SubstrateFile=0;
  cdouble Omega=0.1;
  char *EPFile=0;
  double XDS[6]={1.0, 0.0, 1.0, 0.0, 0.0, 0.5}; int nXDS=1;
  bool FreeSpace=false;
  bool OmitFreeSpace=false;
  int EntryOnly=-1;
  bool EEOnly=false; 
  bool XYOnly=false;
  bool Force2D=false;
  /* name, type, #args, max_instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"Geometry",      PA_STRING,  1, 1, (void *)&GeoFile,    0, ".scuffgeo file"},
     {"SubstrateFile", PA_STRING,  1, 1, (void *)&SubstrateFile, 0, ".substrate file"},
     {"Omega",         PA_CDOUBLE, 1, 1, (void *)&Omega,      0, "angular frequency"},
     {"EPFile",        PA_STRING,  1, 1, (void *)&EPFile,     0, "list of evaluation points"},
     {"XDS",           PA_DOUBLE,  6, 1, (void *)XDS,         &nXDS, ""},
     {"FreeSpace",     PA_BOOL,    1, 1, (void *)&FreeSpace,  0, ""},
     {"OmitFreeSpace", PA_BOOL,    0, 1, (void *)&OmitFreeSpace, 0, ""},
     {"EntryOnly",     PA_INT,     1, 1, (void *)&EntryOnly,  0, ""},
     {"EEOnly",        PA_BOOL,    0, 1, (void *)&EEOnly,     0, ""},
     {"XYOnly",        PA_BOOL,    0, 1, (void *)&XYOnly,     0, ""},
     {"Force2D",       PA_BOOL,    0, 1, (void *)&Force2D,    0, ""},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  /***************************************************************/
  /* read substrate                                              */
  /***************************************************************/
  bool OwnsSubstrateFile=false;
  LayeredSubstrate *S=0;
  RWGGeometry *G=0;
  if (GeoFile)
   { G=new RWGGeometry(GeoFile);
     S=G->Substrate;
     printf("Read substrate definition from file %s\n",GeoFile);
   }
  else 
   { 
     if (!SubstrateFile)
      { SubstrateFile=strdup("XXXXXX");
        if ( mkstemp(SubstrateFile) == -1 )
         ErrExit("could not create temporary file");
        FILE *f=fopen(SubstrateFile,"w");
        fprintf(f,SISubstrateFile);
        fclose(f);
        OwnsSubstrateFile=true;
        printf("Using build-in substrate definition.\n");
      };
     S=new LayeredSubstrate(SubstrateFile);
     if (S->ErrMsg)
      ErrExit(S->ErrMsg);
     printf("Read substrate definition from file %s\n",SubstrateFile);
     if (OwnsSubstrateFile)unlink(SubstrateFile);
   };
  S->Describe();

  if (FreeSpace) 
   { S->ForceFreeSpace=true;
     printf("Doing the free-space case foryaf.\n");
   };

  S->EntryOnly = EntryOnly;
  S->EEOnly    = EEOnly;
  S->XYOnly    = XYOnly;
  
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  HMatrix *XMatrix;
  if (EPFile)
   XMatrix = new HMatrix(EPFile);
  else if (nXDS)
   { 
     XMatrix = new HMatrix(1,6);
     XMatrix->SetEntriesD(0,":",XDS);
     printf("Rho={%e,%e}, ZDest=%e, ZSource=%e\n",XDS[0]-XDS[3],XDS[1]-XDS[4],XDS[2],XDS[5]);
   }
  else
   { 
     #define NUMPTS  150
     #define XMIN   -1.5
     #define XMAX    1.5
     #define DX      (XMAX - XMIN)/(NUMPTS-1)
     XMatrix = new HMatrix(NUMPTS, 6);
     double XXP[6] = {0.0, 0.0, 0.0, 0.25, 0.5, 1.0};
     for(int nx=0; nx<NUMPTS; nx++)
      { XXP[2] = XMIN + ((double)nx)*DX;
        XMatrix->SetEntriesD(nx, "0:5", XXP);
      };
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  FILE *ff=fopen("tFullWaveSubstrate.gp","w");
  fprintf(ff,"XD(x)=(column(1))\n");
  fprintf(ff,"YD(x)=(column(2))\n");
  fprintf(ff,"ZD(x)=(column(3))\n");
  fprintf(ff,"XS(x)=(column(4))\n");
  fprintf(ff,"YS(x)=(column(5))\n");
  fprintf(ff,"ZS(x)=(column(6))\n");
  const char *EM="EM";
  const char *xyz="xyz";
  for(int pMu=0, nc=7; pMu<6; pMu++)
   for(int qNu=0; qNu<6; qNu++, nc+=2)
    { 
      int p  = pMu/3;
      int Mu = pMu%3;
      int q  = qNu/3;
      int Nu = qNu%3;

      int Offset=0;
      fprintf(ff,"rG%c%c%c%cSL(x)=($%i)\n",EM[p],EM[q],xyz[Mu],xyz[Nu],Offset+nc);
      fprintf(ff,"iG%c%c%c%cSL(x)=($%i)\n",EM[p],EM[q],xyz[Mu],xyz[Nu],Offset+nc+1);
      fprintf(ff,"mG%c%c%c%cSL(x)=(D2($%i,$%i))\n",EM[p],EM[q],xyz[Mu],xyz[Nu],Offset+nc,Offset+nc+1);

      Offset+=72;
      fprintf(ff,"rG%c%c%c%cSC(x)=($%i)\n",EM[p],EM[q],xyz[Mu],xyz[Nu],Offset+nc);
      fprintf(ff,"iG%c%c%c%cSC(x)=($%i)\n",EM[p],EM[q],xyz[Mu],xyz[Nu],Offset+nc+1);
      fprintf(ff,"mG%c%c%c%cSC(x)=(D2($%i,$%i))\n",EM[p],EM[q],xyz[Mu],xyz[Nu],Offset+nc,Offset+nc+1);

      Offset+=72;
      fprintf(ff,"rG%c%c%c%cFS(x)=($%i)\n",EM[p],EM[q],xyz[Mu],xyz[Nu],Offset+nc);
      fprintf(ff,"iG%c%c%c%cFS(x)=($%i)\n",EM[p],EM[q],xyz[Mu],xyz[Nu],Offset+nc+1);
      fprintf(ff,"mG%c%c%c%cFS(x)=(D2($%i,$%i))\n",EM[p],EM[q],xyz[Mu],xyz[Nu],Offset+nc,Offset+nc+1);

    };
  fclose(ff);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SetConsoleLogging();

  HMatrix *GSL = new HMatrix(36, XMatrix->NR, LHM_COMPLEX);
  if (FreeSpace)
   AddScriptG0(S, Omega, XMatrix, GSL);
  else 
   S->GetSubstrateDGF(Omega, XMatrix, GSL, STATIC_LIMIT);

  HMatrix *GFS = new HMatrix(GSL->NR, GSL->NC, LHM_COMPLEX);
  GFS->Zero();
  AddScriptG0(S, Omega, XMatrix, GFS);

  HMatrix *GSC = S->GetSubstrateDGF(Omega, XMatrix, SURFACE_CURRENT);

  HMatrix *GBF = new HMatrix(36, XMatrix->NR, LHM_COMPLEX);
  GetSubstrateDGF_SC(Substrate, Omega, XMatrix, GBF, Force2D);
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/

  FILE *f=fopen("tFullWaveSubstrate.out","w");
  for(int nx=0; nx<XMatrix->NR; nx++)
   { 
     double XXP[6];
     XMatrix->GetEntriesD(nx,"0:5",XXP);
     fprintVec(f,XXP,6);

     cdouble GVector[36];

     GSL->GetEntries(":",nx,GVector);
     fprintVec(f,GVector,36);

     GSC->GetEntries(":",nx,GVector);
     fprintVec(f,GVector,36);

     GFS->GetEntries(":",nx,GVector);
     fprintVec(f,GVector,36);
     fprintf(f,"\n");

     if (!FreeSpace)
      GFS->Zero();

     if (XMatrix->NR == 1)
      { 
        SetDefaultCD2SFormat("{%+.2e,%+.2e}");
        cdouble *GSCij = (cdouble *)GSC->GetColumnPointer(0);
        cdouble *GFSij = (cdouble *)GFS->GetColumnPointer(0);
        cdouble *GSLij = (cdouble *)GSL->GetColumnPointer(0);

        for(int P=0; P<2; P++)
         for(int Q=0; Q<2; Q++)
          { double Norm=0.0, Diff=0.0;
            for(int i=0; i<3; i++)
             for(int j=0; j<3; j++)
              { cdouble gFS = GFSij[ (3*P+i)*6 + (3*Q+j) ];
                cdouble gSC = GSCij[ (3*P+i)*6 + (3*Q+j) ];
                cdouble gSL = GSLij[ (3*P+i)*6 + (3*Q+j) ];
		Norm += norm(gSL);
		Diff += norm(gSC+gFS-gSL);
                GSCij[ (3*P+i)*6 + (3*Q+j)] += gFS;
              };

            printf("\n ** Quadrant %c%c: RD=%e \n\n",EM[P],EM[Q],Diff/Norm);
            for(int i=0; i<3; i++)
             printf("%s %s %s || %s %s %s \n",
                    CD2S( GSC->GetEntry( (P*3+i)*6 + (Q*3+0) , 0)),
                    CD2S( GSC->GetEntry( (P*3+i)*6 + (Q*3+1) , 0)),
                    CD2S( GSC->GetEntry( (P*3+i)*6 + (Q*3+2) , 0)),
                    CD2S( GSL->GetEntry( (P*3+i)*6 + (Q*3+0) , 0)),
                    CD2S( GSL->GetEntry( (P*3+i)*6 + (Q*3+1) , 0)),
                    CD2S( GSL->GetEntry( (P*3+i)*6 + (Q*3+2) , 0))
                   );
          };
      };

   };
  fclose(f);
  printf("Thank you for your support.\n");

}
