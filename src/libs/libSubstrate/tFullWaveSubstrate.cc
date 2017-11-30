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
     int nrDest   = S->GetLayerIndex(XD[2]);
     int nrSource = S->GetLayerIndex(XS[2]);
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
      fprintf(ff,"rG%c%c%c%cFull(x)=($%i)\n",EM[p],EM[q],xyz[Mu],xyz[Nu],Offset+nc);
      fprintf(ff,"iG%c%c%c%cFull(x)=($%i)\n",EM[p],EM[q],xyz[Mu],xyz[Nu],Offset+nc+1);
      fprintf(ff,"mG%c%c%c%cFull(x)=(D2($%i,$%i))\n",EM[p],EM[q],xyz[Mu],xyz[Nu],Offset+nc,Offset+nc+1);

      Offset+=72;
      fprintf(ff,"rG%c%c%c%cFast(x)=($%i)\n",EM[p],EM[q],xyz[Mu],xyz[Nu],Offset+nc);
      fprintf(ff,"iG%c%c%c%cFast(x)=($%i)\n",EM[p],EM[q],xyz[Mu],xyz[Nu],Offset+nc+1);
      fprintf(ff,"mG%c%c%c%cFast(x)=(D2($%i,$%i))\n",EM[p],EM[q],xyz[Mu],xyz[Nu],Offset+nc,Offset+nc+1);

    };
  fclose(ff);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SetConsoleLogging();

  HMatrix *GFast = S->GetSubstrateDGF(Omega, XMatrix, FAST_SURFACE_CURRENT);
  HMatrix *GFull = S->GetSubstrateDGF(Omega, XMatrix, FULL_SURFACE_CURRENT);

  HMatrix *GFS = new HMatrix(36, XMatrix->NR, LHM_COMPLEX);
  AddScriptG0(S, Omega, XMatrix, GFS);

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

     GFast->GetEntries(":",nx,GVector);
     fprintVec(f,GVector,36);

     GFull->GetEntries(":",nx,GVector);
     fprintVec(f,GVector,36);

     GFS->GetEntries(":",nx,GVector);
     fprintVec(f,GVector,36);

     fprintf(f,"\n");

     //GSL->GetEntries(":",nx,GVector);
     //fprintVec(f,GVector,36);

     //GSC->GetEntries(":",nx,GVector);
     //fprintVec(f,GVector,36);

     //GFS->GetEntries(":",nx,GVector);
     //fprintVec(f,GVector,36);
     //fprintf(f,"\n");

     //if (!FreeSpace)
      //GFS->Zero();

     if (XMatrix->NR == 1)
      { 
        SetDefaultCD2SFormat("{%+.2e,%+.2e}");
        cdouble *GijFast = (cdouble *)GFast->GetColumnPointer(0);
        cdouble *GijFull = (cdouble *)GFull->GetColumnPointer(0);

        for(int P=0; P<2; P++)
         for(int Q=0; Q<2; Q++)
          { double Norm=0.0, Diff=0.0;
            for(int i=0; i<3; i++)
             for(int j=0; j<3; j++)
              { cdouble X = GijFast[ (3*P+i)*6 + (3*Q+j) ];
                cdouble Y = GijFull[ (3*P+i)*6 + (3*Q+j) ];
		Norm += norm(X);
		Diff += norm(X-Y);
              };

            printf("\n ** Quadrant %c%c: RD=%e \n\n",EM[P],EM[Q],Diff/Norm);
            for(int i=0; i<3; i++)
             printf("%s %s %s | %s %s %s | %.1e %.1e %.1e\n",
                    CD2S( GFast->GetEntry( (P*3+i)*6 + (Q*3+0) , 0)),
                    CD2S( GFast->GetEntry( (P*3+i)*6 + (Q*3+1) , 0)),
                    CD2S( GFast->GetEntry( (P*3+i)*6 + (Q*3+2) , 0)),
                    CD2S( GFull->GetEntry( (P*3+i)*6 + (Q*3+0) , 0)),
                    CD2S( GFull->GetEntry( (P*3+i)*6 + (Q*3+1) , 0)),
                    CD2S( GFull->GetEntry( (P*3+i)*6 + (Q*3+2) , 0)),
                    RD( GFast->GetEntry( (3*P+i)*6 + (3*Q+0), 0),
                        GFull->GetEntry( (3*P+i)*6 + (3*Q+0), 0)
                      ),
                    RD( GFast->GetEntry( (3*P+i)*6 + (3*Q+1), 0),
                        GFull->GetEntry( (3*P+i)*6 + (3*Q+1), 0)
                      ),
                    RD( GFast->GetEntry( (3*P+i)*6 + (3*Q+2), 0),
                        GFull->GetEntry( (3*P+i)*6 + (3*Q+2), 0)
                      )
                   );
          };
      };

   };
  fclose(f);
  printf("Thank you for your support.\n");

}
