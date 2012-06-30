/*
 * GBarVD[0] = GBar
 * GBarVD[1] = dGBar/dX
 * GBarVD[2] = dGBar/dY
 * GBarVD[3] = dGBar/dZ
 * GBarVD[4] = d^2GBar/dXdY
 * GBarVD[5] = d^2GBar/dXdZ
 * GBarVD[6] = d^2GBar/dYdZ
 * GBarVD[7] = d^3GBar/dXdYdZ
 */
#include <stdlib.h>
#include <math.h>

#include "libhrutil.h"
#include "libMDInterp.h"

#define ABSTOL 0.0
#define RELTOL 1.0e-8

#define II cdouble(0,1)

#define NSUM 8
#define NFIRSTROUND 1
#define NMAX 10000

/***************************************************************/
/***************************************************************/
/***************************************************************/
cdouble cerfc(cdouble z);

void GBarVDBF(cdouble Beta, double *K, double *L1, double *L2, double *R,
              double AbsTol, double RelTol, int *pnCells, cdouble *GBarVD);

/***************************************************************/
/* 'convergence factor 2' **************************************/
/***************************************************************/
void GetCF2(double z, double E, cdouble kz, cdouble *CF2, cdouble *CF2Prime)
{ 
  cdouble Arg, ExpFac, dExpFac, ErfcFac, dErfcFac;
  cdouble PlusTerm, dPlusTerm, MinusTerm, dMinusTerm;

  // PlusTerm  = exp(  kz*R[2] ) * cerfc( 0.5*kz/E  + R[2]*E );
  ExpFac    = exp( kz*z );
  dExpFac   = kz*ExpFac;
  Arg       = 0.5*kz/E + z*E;
  ErfcFac   = cerfc( Arg );
  dErfcFac  = -E*exp( -Arg*Arg );
  PlusTerm  = ExpFac*ErfcFac;
  dPlusTerm = dExpFac*ErfcFac + ExpFac*dErfcFac;

  // MinusTerm  = exp( -kz*R[2] ) * cerfc( 0.5*kz/E  - R[2]*E );
  ExpFac     = exp( -kz*z );
  dExpFac    = -kz*ExpFac;
  Arg        = 0.5*kz/E - z*E;
  ErfcFac    = cerfc( Arg );
  dErfcFac   = +E*exp( -Arg*Arg );
  MinusTerm  = ExpFac*ErfcFac;
  dMinusTerm = dExpFac*ErfcFac + ExpFac*dErfcFac;
  
  *CF2      = 0.5*(PlusTerm + MinusTerm);
  *CF2Prime = 0.5*(dPlusTerm + dMinusTerm);

  if ( !isfinite(*CF2) )
   *CF2=0.0;
  if ( !isfinite(*CF2Prime) )
   *CF2Prime=0.0;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddG1Contribution(cdouble Beta, double *K, double *R, 
                       double Qx, double Qy,
                       double E, cdouble *GBarVD)
{ 
  double QpK[2];
  cdouble PhaseFactor, kz, CF2, CF2Prime;

  QpK[0] = Qx + K[0];
  QpK[1] = Qy + K[1];

  kz = sqrt ( QpK[0]*QpK[0] + QpK[1]*QpK[1] - Beta*Beta );

  PhaseFactor = exp( II * (QpK[0]*R[0] + QpK[1]*R[1]) ) / kz;

  GetCF2(R[2], E, kz, &CF2, &CF2Prime);

  GBarVD[0] += PhaseFactor * CF2;
  GBarVD[1] += II*QpK[0]*PhaseFactor*CF2;
  GBarVD[2] += II*QpK[1]*PhaseFactor*CF2;
  GBarVD[3] += PhaseFactor*CF2Prime;
  GBarVD[4] += -QpK[0]*QpK[1]*PhaseFactor*CF2;
  GBarVD[5] += II*QpK[0]*PhaseFactor*CF2Prime;
  GBarVD[6] += II*QpK[1]*PhaseFactor*CF2Prime;
  GBarVD[7] += -QpK[0]*QpK[1]*PhaseFactor*CF2Prime;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ComputeG1(cdouble Beta, double *K, double *L1, double *L2, 
               double *R, double E, int *pnCells, cdouble *Sum)
{ 
  int nx, ny;
  double RelDelta, AbsDelta;
  cdouble LastSum[NSUM];
  double MaxRelDelta, MaxAbsDelta;
  double Delta, AbsSum;
  int i, NN, ConvergedIters;
  int nCells=0;

  double Q1[2], Q2[2];

  if ( !(L1[1]==0.0 && L2[0]==0.0) )
   ErrExit("non-square lattices not supported");

  Q1[0] = 2.0*M_PI / L1[0];
  Q1[1] = 0.0;
  Q2[0] = 0.0;
  Q2[1] = 2.0*M_PI / L2[1];

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  memset(Sum,0,NSUM*sizeof(cdouble));

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for (nx=-NFIRSTROUND; nx<=NFIRSTROUND; nx++)
   for (ny=-NFIRSTROUND; ny<=NFIRSTROUND; ny++, nCells++)
    AddG1Contribution(Beta, K, R,
                      nx*Q1[0] + ny*Q2[0],
                      nx*Q1[1] + ny*Q2[1],
                      E, Sum);
         
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  memcpy(LastSum,Sum,NSUM*sizeof(cdouble));
  ConvergedIters=0;
  for(NN=NFIRSTROUND+1; ConvergedIters<3 && NN<=NMAX; NN++)
   {  
     for(nx=-NN; nx<=NN; nx++)
      for(ny=-NN; ny<=NN; ny++)
        { 
          if ( (abs(nx)<NN) && (abs(ny)<NN) )
           continue;

          AddG1Contribution(Beta, K, R, 
                            nx*Q1[0] + ny*Q2[0],
                            nx*Q1[1] + ny*Q2[1], 
                            E, Sum);

          nCells++;

        };

     /*--------------------------------------------------------------*/
     /* convergence analysis ----------------------------------------*/
     /*--------------------------------------------------------------*/
     MaxAbsDelta=MaxRelDelta=0.0;
     for(i=0; i<NSUM; i++)
      { Delta=abs(Sum[i]-LastSum[i]);
        if ( Delta>MaxAbsDelta )
         MaxAbsDelta=Delta;
        AbsSum=abs(Sum[i]);
        if ( AbsSum>0.0 && (Delta > MaxRelDelta*AbsSum) )
         MaxRelDelta=Delta/AbsSum;
      };
     if ( MaxAbsDelta<ABSTOL || MaxRelDelta<RELTOL )
      ConvergedIters++;
     else
      ConvergedIters=0;

     memcpy(LastSum,Sum,NSUM*sizeof(cdouble));

   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int ns; 
  for(ns=0; ns<NSUM; ns++)
   Sum[ns] *= 1.0 / (2.0*L1[0]*L2[1]);

  if (pnCells) 
   *pnCells=nCells;

}

/***************************************************************/
/* add the contribution of a single direct lattice vector L    */
/* to the direct-lattice sum.                                  */
/*                                                             */
/* note: the summand is:                                       */
/*                                                             */
/*  exp(i K\dot L) * g_1 *( g_{2+} g_{3+} + g_{2-}g_{3-})      */  
/*                                                             */
/* where:                                                      */
/*  g_1       = 1/(8*pi*|R|)                                   */
/*  g_{2+}    = exp(i*Beta*|R|)                                */
/*  g_{3+}    = erfc( E*R + i*Beta/(2*E) )                     */
/*  g_{2-}    = exp(-i*Beta*|R|)                               */
/*  g_{3-}    = erfc( E*R - i*Beta/(2*E) )                     */
/* (here |R| = |r+L|.)                                         */
/*                                                             */
/* if we define the following combinations:                    */
/*                                                             */
/*  ggPgg = g_{2+} g_{3+}  +  g_{2-} g_{3-}                    */
/*  ggMgg = g_{2+} g_{3+}  -  g_{2-} g_{3-}                    */
/*                                                             */
/* then we have the following derivative relationships:        */
/*                                                             */
/* d/dR (ggPgg) = i*Beta*ggMgg + g4                            */
/* d/dR (ggMgg) = i*Beta*ggPgg                                 */
/*                                                             */
/* where g4 = (-4E/sqrt(pi)) * exp( -E^2R^2 + Beta^2/(4E^2).   */
/*                                                             */
/***************************************************************/
void AddG2Contribution(cdouble Beta, double *K, double *R, 
                       double Lx, double Ly, double E, cdouble *Sum)
{ 
  cdouble PhaseFactor; 
  double RmL[3], rml2, rml, rml3, rml4, rml5, rml6, rml7;
  cdouble g2p, g3p, g2m, g3m, g4, ggPgg, ggMgg, Term;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  PhaseFactor=exp( II * (K[0]*Lx + K[1]*Ly) ) / (8.0*M_PI);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RmL[0] = (R[0]-Lx);
  RmL[1] = (R[1]-Ly);
  RmL[2] =  R[2];

  rml2=RmL[0]*RmL[0] + RmL[1]*RmL[1] + RmL[2]*RmL[2];
  rml=sqrt(rml2);
  if ( rml < 1.0e-7 ) 
   return;
  rml3=rml2*rml;
  rml4=rml3*rml;
  rml5=rml4*rml;
  rml6=rml5*rml;
  rml7=rml6*rml;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  g2p = exp( II*Beta*rml );
  g3p = cerfc( E*rml + II*Beta/(2.0*E) );

  g2m = exp( -II*Beta*rml );
  g3m = cerfc( E*rml - II*Beta/(2.0*E) );

  g4 = -2.0*M_2_SQRTPI*E*exp(-E*E*rml2 + Beta*Beta/(4.0*E*E));

  ggPgg = g2p*g3p + g2m*g3m;
  ggMgg = g2p*g3p - g2m*g3m;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Sum[0] += PhaseFactor * ggPgg / rml;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Term = -ggPgg/rml3 + (g4 + II*Beta*ggMgg)/rml2;

  Sum[1] += PhaseFactor * RmL[0] * Term;
  Sum[2] += PhaseFactor * RmL[1] * Term;
  Sum[3] += PhaseFactor * RmL[2] * Term;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Term = 3.0*ggPgg/rml5 - 3.0*(g4+II*Beta*ggMgg)/rml4 - Beta*Beta*ggPgg/rml3 - 2.0*E*E*g4/rml2;

  Sum[4] += PhaseFactor * RmL[0] * RmL[1] * Term;
  Sum[5] += PhaseFactor * RmL[0] * RmL[2] * Term;
  Sum[6] += PhaseFactor * RmL[1] * RmL[2] * Term;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Term = -15.0*ggPgg/rml7 + 15.0*(g4+II*Beta*ggMgg)/rml6 
         + 6.0*Beta*Beta*ggPgg/rml5 + 10.0*E*E*g4/rml4
         -Beta*Beta*(II*Beta*ggMgg + g4)/rml4 + 4.0*E*E*E*E*g4/rml2;

  Sum[7] += PhaseFactor * RmL[0] * RmL[1] * RmL[2] * Term;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ComputeG2(cdouble Beta, double *K, double *L1, double *L2,
               double *R, double E, int *pnCells, cdouble *Sum)
{ 
  int nx, ny;
  double RelDelta, AbsDelta;
  cdouble LastSum[NSUM];
  double MaxRelDelta, MaxAbsDelta;
  double Delta, AbsSum;
  int i, NN, ConvergedIters;

  int nCells=0;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  memset(Sum,0,NSUM*sizeof(cdouble));

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for (nx=-NFIRSTROUND; nx<=NFIRSTROUND; nx++)
   for (ny=-NFIRSTROUND; ny<=NFIRSTROUND; ny++, nCells++)
    { 
#if 0
      if ( (abs(nx)<=1) && (abs(ny)<=1) )
       continue; // skip the innermost 9 grid cells 
#endif

      AddG2Contribution(Beta, K, R, 
                        nx*L1[0] + ny*L2[0],
                        nx*L1[1] + ny*L2[1], 
                        E, Sum);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
//printf("%i %i %s \n",nx,ny,CD2S(Sum[0]));
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
    };
         
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  memcpy(LastSum,Sum,NSUM*sizeof(cdouble));
  ConvergedIters=0;
  for(NN=NFIRSTROUND+1; ConvergedIters<3 && NN<=NMAX; NN++)
   {  
     for(nx=-NN; nx<=NN; nx++)
      for(ny=-NN; ny<=NN; ny++)
        { 
          if ( (abs(nx)<NN) && (abs(ny)<NN) )
           continue;


          AddG2Contribution(Beta, K, R, 
                            nx*L1[0] + ny*L2[0],
                            nx*L1[1] + ny*L2[1], 
                            E, Sum);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
//printf("%i %i %s \n",nx,ny,CD2S(Sum[0]));
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

          nCells++;
        };

     /*--------------------------------------------------------------*/
     /* convergence analysis ----------------------------------------*/
     /*--------------------------------------------------------------*/
     MaxAbsDelta=MaxRelDelta=0.0;
     for(i=0; i<NSUM; i++)
      { Delta=abs(Sum[i]-LastSum[i]);
        if ( Delta>MaxAbsDelta )
         MaxAbsDelta=Delta;
        AbsSum=abs(Sum[i]);
        if ( AbsSum>0.0 && (Delta > MaxRelDelta*AbsSum) )
         MaxRelDelta=Delta/AbsSum;
      };
     if ( MaxAbsDelta<ABSTOL || MaxRelDelta<RELTOL )
      ConvergedIters++;
     else
      ConvergedIters=0;

     memcpy(LastSum,Sum,NSUM*sizeof(cdouble));

   };

  if (pnCells) *pnCells=nCells;

}

/***************************************************************/
/***************************************************************/
/* compute the contribution of the innermost 9 real-space      */
/* cells to the fourier-space sum                              */
/***************************************************************/
/***************************************************************/
#if 0 // this seems not to have worked
void ComputeG1First9(cdouble Beta, double *K, double *L1, double *L2, 
                     double *R, double E, cdouble *Sum)
{
  int n1, n2;
  double L[2], RmL[3], rml2, rml;
  cdouble PhaseFactor, IBR, Phi, Psi, Zeta, Upsilon;
  cdouble Part1[8], Part2[8];

  memset(Part1,0,NSUM*sizeof(cdouble));
  memset(Part2,0,NSUM*sizeof(cdouble));
  for(n1=-1; n1<=1; n1++)
   for(n2=-1; n2<=1; n2++)
    { 
      L[0] = n1*L1[0] + n2*L2[0];
      L[1] = n1*L1[1] + n2*L2[1];

      PhaseFactor = exp( II*(K[0]*L[0] + K[1]*L[1]) ) / (4.0*M_PI);

      RmL[0] = R[0]-L[0];
      RmL[1] = R[1]-L[1];
      RmL[2] = R[2];

      rml2=RmL[0]*RmL[0] + RmL[1]*RmL[1] + RmL[2]*RmL[2];
      rml=sqrt(rml2);
      if (fabs(rml)<1.0e-7) 
       continue;

      IBR=II*Beta*rml;
      Phi=exp(IBR) / (4.0*M_PI*rml);
      Psi=(-1.0 + IBR) * Phi / rml2;
      Zeta=(3.0 + IBR*(-3.0 + IBR))*Phi/(rml2*rml2);
      Upsilon=(-15.0 + IBR*(15.0 + IBR*(-6.0 + IBR)))*Phi/(rml2*rml2*rml2);
 
      Part1[0] += PhaseFactor * Phi;
      Part1[1] += PhaseFactor * RmL[0] * Psi;
      Part1[2] += PhaseFactor * RmL[1] * Psi;
      Part1[3] += PhaseFactor * RmL[2] * Psi;
      Part1[4] += PhaseFactor * RmL[0] * RmL[1] * Zeta;
      Part1[5] += PhaseFactor * RmL[0] * RmL[2] * Zeta;
      Part1[6] += PhaseFactor * RmL[1] * RmL[2] * Zeta;
      Part1[7] += PhaseFactor * RmL[1] * RmL[2] * RmL[3] * Zeta;

      AddG2Contribution(Beta, K, R, L[0], L[1], E, Part2);

    };

  int ns;
  for(ns=0; ns<NSUM; ns++)
   Sum[ns] = Part1[ns] - 0.5*Part2[ns];

};
#endif

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddGBFContribution(cdouble Beta, double K[2], double R[3],
                        double Lx, double Ly, cdouble *Sum)
{ 

  double RmL[3];
  double r2, r, kr; 
  cdouble PhaseFactor, IBR, Phi, Psi, Zeta, Upsilon;
   
  PhaseFactor=exp( II*(Lx*K[0] + Ly*K[1]) );

  RmL[0]=R[0] - Lx;
  RmL[1]=R[1] - Ly;
  RmL[2]=R[2];

  r2=RmL[0]*RmL[0] + RmL[1]*RmL[1] + RmL[2]*RmL[2];
  r=sqrt(r2);
  if (r<1.0e-7) 
   return;
  IBR=II*Beta*r;
  Phi=exp(IBR)/(4.0*M_PI*r);
  Psi=(IBR-1.0)*Phi/r2;
  Zeta=(3.0 + IBR*(-3.0 + IBR))*Phi/(r2*r2);
  Upsilon=(-15.0 + IBR*(15.0 + IBR*(-6.0 + IBR)))*Phi/(r2*r2*r2);

  Sum[0] += PhaseFactor * Phi;
  Sum[1] += PhaseFactor * RmL[0] * Psi;
  Sum[2] += PhaseFactor * RmL[1] * Psi;
  Sum[3] += PhaseFactor * RmL[2] * Psi;
  Sum[4] += PhaseFactor * RmL[0] * RmL[1] * Zeta;
  Sum[5] += PhaseFactor * RmL[0] * RmL[2] * Zeta;
  Sum[6] += PhaseFactor * RmL[1] * RmL[2] * Zeta;
  Sum[7] += PhaseFactor * RmL[0] * RmL[1] * RmL[2] * Upsilon;
 
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ComputeGBFFirst9(cdouble Beta, double *K, double *L1, double *L2,
                      double *R, cdouble *Sum)
{ 
  int nx, ny;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  memset(Sum,0,NSUM*sizeof(cdouble));
  for (nx=-1; nx<=1; nx++)
   for (ny=-1; ny<=1; ny++)
    AddGBFContribution(Beta, K, R, 
                       nx*L1[0] + ny*L2[0],
                       nx*L1[1] + ny*L2[1], 
                       Sum);
}

/***************************************************************/
/* 'GBar values and derivatives,' computed via Ewald's method **/
/***************************************************************/
double ExternalE=-1.0;
void GBarVDEwald(cdouble Beta, double *K, double *L1, double *L2, 
                 double *R, cdouble *GBarVD)
{ 
  
  if ( (L1[1]!=0.0) || (L2[0]!=0.0) ) 
   ErrExit("%s:%i: internal error",__FILE__,__LINE__);

  double E=sqrt( M_PI / (L1[0]*L2[1]) );

  /***************************************************************/
  /* FIXME *******************************************************/
  /***************************************************************/
  int OnLatticePoint=0.0;
  if ( fabs(R[2]) < 1.0e-7
      && ( (fabs(R[0]) < 1.0e-7) || (fabs(fabs(R[0])-1.0) < 1.0e-7) )
      && ( (fabs(R[1]) < 1.0e-7) || (fabs(fabs(R[1])-1.0) < 1.0e-7) )
     ) E=0.1;

if (ExternalE!=-1.0)
 E=ExternalE;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble G1[NSUM], G2[NSUM], GBFFirst9[NSUM];

  ComputeG1(Beta, K, L1, L2, R, E, 0, G1);
  ComputeG2(Beta, K, L1, L2, R, E, 0, G2);
  ComputeGBFFirst9(Beta, K, L1, L2, R, GBFFirst9);

  int ns;
  for(ns=0; ns<NSUM; ns++)
   GBarVD[ns] = G1[ns] + G2[ns] - GBFFirst9[ns];

} 
