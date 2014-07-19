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
 * scuff-test-Ewald.cc -- a test program for libscuff's routines for
 *                     -- evaluating the periodic Green's function
 * 
 * homer reid          -- 11/2005 -- 12/2012
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fenv.h>

#include <config.h>

#ifdef HAVE_LIBREADLINE
 #include <readline/readline.h>
 #include <readline/history.h>
#else
 #include "readlineReplacement.h"
#endif

#include <complex>

#include <libhrutil.h>

#include <libscuff.h>
#include <libscuffInternals.h>

#define NSUM 8
#define NFIRSTROUND 10
#define NMAX 10000

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
namespace scuff{

void GetGBarDistant(double *R, double Rho2, cdouble k, double *kBloch,
                    double Gamma[2][2], int LDim,
                    double E, int *pnCells, cdouble *Sum);

void GetGBarNearby(double *R, cdouble k, double *kBloch,
                   double *LBV[2], int LDim,
                   double E, bool ExcludeInnerCells,
                   int *pnCells, cdouble *Sum);

void AddGFull(double R[3], cdouble k, double kBloch[2],
              double Lx, double Ly, cdouble *Sum, 
              bool ValueOnly=false);

void GetRLBasis(double **L, int LDim, double Gamma[2][2],
                cdouble k, double *EOpt, double R[3], double *Rho2);

void AddGLongRealSpace(double *R, cdouble k, double *kBloch,
                       int n1, int n2, double *LBV[2], int LDim, 
                       double E, cdouble *Sum);

}


/***************************************************************/
/***************************************************************/
/***************************************************************/
void GBarVDBF(double R[3], cdouble k, double *kBloch,
              double **LBV, int LDim,
              bool ExcludeInnerCells, bool Derivatives,
              int nMax, cdouble *Sum)
{ 
  double L1[2], L2[2];
  
  L1[0]=LBV[0][0];
  L1[1]=LBV[0][1];
  if (LDim==2)
   { L2[0]=LBV[1][0];
     L2[1]=LBV[1][1];
   }
  else
   L2[0]=L2[1]=0.0;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  memset(Sum,0,NSUM*sizeof(cdouble));
  int n2Mult = (LDim==1) ? 0 : 1;
  int NumThreads = GetNumThreads();
  double rSum0=0.0, iSum0=0.0, rSum1=0.0, iSum1=0.0;
  double rSum2=0.0, iSum2=0.0, rSum3=0.0, iSum3=0.0;
  double rSum4=0.0, iSum4=0.0, rSum5=0.0, iSum5=0.0;
  double rSum6=0.0, iSum6=0.0, rSum7=0.0, iSum7=0.0;
#pragma omp parallel for schedule(dynamic,1),      \
                         num_threads(NumThreads),  \
                         reduction(+:rSum0, rSum1, rSum2, rSum3, rSum4, rSum5, rSum6, rSum7, \
                                     iSum0, iSum1, iSum2, iSum3, iSum4, iSum5, iSum6, iSum7)
  for (int n1=-nMax; n1<=nMax; n1++)
   for (int n2=-n2Mult*nMax; n2<=n2Mult*nMax; n2++)
    { 
      if ( ExcludeInnerCells && ( abs(n1)<=1 && abs(n2)<=1 ) )
       continue; // skip the innermost 9 grid cells 

      cdouble PSum[8]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      AddGFull(R, k, kBloch,
               n1*L1[0] + n2*L2[0], n1*L1[1] + n2*L2[1],
               PSum, !Derivatives);

      rSum0 += real(PSum[0]); iSum0 += imag(PSum[0]);
      rSum1 += real(PSum[1]); iSum1 += imag(PSum[1]);
      rSum2 += real(PSum[2]); iSum2 += imag(PSum[2]);
      rSum3 += real(PSum[3]); iSum3 += imag(PSum[3]);
      rSum4 += real(PSum[4]); iSum4 += imag(PSum[4]);
      rSum5 += real(PSum[5]); iSum5 += imag(PSum[5]);
      rSum6 += real(PSum[6]); iSum6 += imag(PSum[6]);
      rSum7 += real(PSum[7]); iSum7 += imag(PSum[7]);
    };

  Sum[0] = cdouble(rSum0, iSum0);
  Sum[1] = cdouble(rSum1, iSum1);
  Sum[2] = cdouble(rSum2, iSum2);
  Sum[3] = cdouble(rSum3, iSum3);
  Sum[4] = cdouble(rSum4, iSum4);
  Sum[5] = cdouble(rSum5, iSum5);
  Sum[6] = cdouble(rSum6, iSum6);
  Sum[7] = cdouble(rSum7, iSum7);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[]) 
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if ( getenv("SCUFF_ABORT_ON_FPE") )
   { feenableexcept(FE_INVALID | FE_OVERFLOW);
     Log("Enabling abort-on-floating-point-exception.");
   };

  /***************************************************************/
  /* process command-line arguments ******************************/
  /***************************************************************/
  char *GeoFileName=0;
  double LBV1[2];		int nLBV1=0;
  double LBV2[2]; 		int nLBV2=0;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"geometry",          PA_STRING,  1, 1, (void *)&GeoFileName,   0,  ".scuffgeo file"},
     {"LBV1",              PA_DOUBLE,  2, 1, (void *)LBV1,      &nLBV1,  "lattice basis vector 1"},
     {"LBV2",              PA_DOUBLE,  2, 1, (void *)LBV2,      &nLBV2,  "lattice basis vector 2"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  /***************************************************************/ 
  /* if a geometry was specified, set the lattice vectors from it*/
  /***************************************************************/
  double LBV[2][2];
  double *LBVP[2] = { LBV[0], LBV[1] };
  int LDim;
  if ( nLBV1!=0 && nLBV2==0 )
   { LBV[0][0] = LBV1[0];
     LBV[0][1] = LBV1[1];
     LDim=1;
   }
  else if (nLBV1!=0 && nLBV2!=0)
   { LBV[0][0] = LBV1[0];
     LBV[0][1] = LBV1[1];
     LBV[1][0] = LBV2[0];
     LBV[1][1] = LBV2[1];
     LDim=2;
   }
  else if (GeoFileName!=0)
   { 
     RWGGeometry *G = new RWGGeometry(GeoFileName);
     LDim = G->NumLatticeBasisVectors;
     if (LDim>=1)
      { LBV[0][0] = G->LatticeBasisVectors[0][0];
        LBV[0][1] = G->LatticeBasisVectors[0][1];
      };
     if (LDim>=2)
      { LBV[1][0] = G->LatticeBasisVectors[1][0];
        LBV[1][1] = G->LatticeBasisVectors[1][1];
      }
   }
  else // default to a square lattice with side length 1
   { LDim=2;
     LBV[0][0] = 1.0;    LBV[0][1] = 0.0;
     LBV[1][0] = 0.0;    LBV[1][1] = 1.0;
   };

  srand48(time(0));
  SetDefaultCD2SFormat("(%15.8e,%15.8e)");
  printf("Considering a %i-d lattice with ", LDim);
  if (LDim==1)
   printf(" L=(%g,%g)\n",LBV[0][0],LBV[0][1]);
  else 
   printf(" L1=(%g,%g), L2=(%g,%g)\n",
            LBV[0][0],LBV[0][1],LBV[1][0],LBV[1][1]);

  /***************************************************************/
  /* enter command loop ******************************************/
  /***************************************************************/
  using_history();
  read_history(0);
  for(;;)
   { 
     /*--------------------------------------------------------------*/
     /*- print prompt and get input string --------------------------*/
     /*--------------------------------------------------------------*/
     printf("\n");
     printf(" options: --x      xx \n");
     printf("          --y      xx \n");
     printf("          --z      xx \n");
     printf("          --k      xx \n");
     printf("          --kBloch xx xx \n");
     printf("          --E      xx\n");
     printf("          --RetainInnerCells\n");
     printf("          --SkipBF\n");
     printf("          --Derivatives\n");
     printf("          --nMax\n");
     printf("          --quit\n");
     char *p=readline("enter options: ");
     if (!p) break;
     add_history(p);
     write_history(0);

     /*--------------------------------------------------------------*/
     /*- set default variable values --------------------------------*/
     /*--------------------------------------------------------------*/
     double R[3];
     R[0] = -5.0 + 10.0*drand48();
     R[1] = -5.0 + 10.0*drand48();
     R[2] = -5.0 + 10.0*drand48();
     cdouble k = cdouble( 5.0*drand48(), 5.0*drand48() );
     double kBloch[2] = { M_PI*drand48(), 0.0 };
     if (LDim==2)
      kBloch[1] = M_PI*drand48();
     double E = -1.0;
     bool ExcludeInnerCells=true;
     bool SkipBF=false;
     bool Derivatives=false;
     int nMax=100;

     /*--------------------------------------------------------------*/
     /* parse input string                                          -*/
     /*--------------------------------------------------------------*/
     char *Tokens[50];
     int NumTokens=Tokenize(p,Tokens,50);
     for(int nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--x") )
       sscanf(Tokens[nt+1],"%le",R+0);
     for(int nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--y") )
       sscanf(Tokens[nt+1],"%le",R+1);
     for(int nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--z") )
       sscanf(Tokens[nt+1],"%le",R+2);
     for(int nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--k") )
       S2CD(Tokens[nt+1],&k);
     for(int nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--kBloch") )
       { sscanf(Tokens[nt+1],"%le",kBloch+0);
         sscanf(Tokens[nt+2],"%le",kBloch+1);
       };
     for(int nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--E") )
       sscanf(Tokens[nt+1],"%le",&E);
     for(int nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--nMax") )
       sscanf(Tokens[nt+1],"%i",&nMax);
     for(int nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--RetainInnerCells") )
       ExcludeInnerCells=false;
     for(int nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--SkipBF") )
       SkipBF=true;
     for(int nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--Derivatives") )
       Derivatives=true;
     for(int nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--quit") )
       exit(1);
       
     free(p);

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     double Gamma[2][2], EOpt, Rho2;
     GetRLBasis(LBVP, LDim, Gamma, k, &EOpt, R, &Rho2);
     if (E==-1.0) E=EOpt;

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     printf("\n** Your options:\n\n");
     printf("--x %g --y %g --z %g --k %s --kBloch %g %g --E %g --nMax %i\n",
              R[0],R[1],R[2],z2s(k),kBloch[0],kBloch[1],E,nMax);
     printf("\n");

     /*--------------------------------------------------------------*/
     /*- do the computation using the ewald summation method         */
     /*--------------------------------------------------------------*/
     cdouble GBarNearby[NSUM], GBarDistant[NSUM], GBarEwald[NSUM];
     int NearbyCells, DistantCells;
     Tic();
     GetGBarNearby(R, k, kBloch, LBVP, LDim,
                   E, ExcludeInnerCells, &NearbyCells, GBarNearby);
     GetGBarDistant(R, Rho2, k, kBloch, Gamma, LDim, 
                    E, &DistantCells, GBarDistant);
     for(int ns=0; ns<NSUM; ns++)
      GBarEwald[ns] = GBarNearby[ns] + GBarDistant[ns];
     double EwaldTime=Toc();

     printf(" Nearby:  %s (%5i cells)\n",CD2S(GBarNearby[0]), NearbyCells);
     printf(" Distant: %s (%5i cells)\n",CD2S(GBarDistant[0]), DistantCells);

     Tic();
     if (ExcludeInnerCells)
      { 
        cdouble GLongInner[NSUM];
        memset(GLongInner,0,NSUM*sizeof(cdouble));
        int n2Mult = (LDim==2) ? 1 : 0;
        for(int n1=-1; n1<=1; n1++)
         for(int n2=-1*n2Mult; n2<=1*n2Mult; n2++)
          AddGLongRealSpace(R, k, kBloch, n1, n2, LBVP, LDim, E, GLongInner);

        printf(" Inner:   %s \n",CD2S(GLongInner[0]));

        for(int ns=0; ns<NSUM; ns++)
         GBarEwald[ns] -= GLongInner[ns];
      };
     EwaldTime+=Toc();
     printf(" Ewald:   %s (%.0g us)\n",CD2S(GBarEwald[0]),EwaldTime*1e6);
   
     /*--------------------------------------------------------------*/
     /*- do the computation using the BF method if that was requested*/
     /*--------------------------------------------------------------*/
     cdouble GBarBF[8];
     if (!SkipBF)
      {
        GBarVDBF(R, k, kBloch, LBVP, LDim, 
                 ExcludeInnerCells, Derivatives, nMax, GBarBF);
        printf(" BF:      %s \n",CD2S(GBarBF[0]));

      };

    if (Derivatives)
     for(int n=1; n<8; n++)
      printf(" %i: %s | %s | %e \n",n, CD2S(GBarEwald[n]),CD2S(GBarBF[n]), RD( GBarEwald[n], GBarBF[n]));

   }; // for(;;)

}
