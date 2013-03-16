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
 * FrequencyIntegral.cc -- scuff-neq module for numerical quadrature 
 *                      -- over frequencies
 *
 * homer reid           -- 5/2012
 *
 */

#include <libSGJC.h>
#include "scuff-neq.h"

// how this works: 
//  a. temperature in eV = kT = 8.6173e-5 * (T in kelvin)
//  b. temperature in our internal energy units
//     = (kT in eV) / (0.1973 eV)
//     = (8.6173e-5 / 0.1973 ) * (T in Kelvin)
#define BOLTZMANNK 4.36763e-4

// for example: suppose in the real world we have 
// omega=3e14 rad/sec, T = 300 kelvin. then
// \hbar \omega / (kT) 
//   = (6.6e-16 ev s )(3e14 s^{-1}) / (0.026 ev) 
//   = 7.6
// whereas in this code we would have Omega==1, T=300, and hence
// Omega/(BOLTZMANNK*T) = (1/(4.36763e-4*300)) = 7.6. 

/***************************************************************/
/***************************************************************/
/***************************************************************/
double Theta(double Omega, double T)
{ 
  if (T==0.0)
   return 0.0;
  return Omega / ( exp( Omega/(BOLTZMANNK*T) ) - 1.0 );
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void PutInThetaFactors(SNEQData *SNEQD, double Omega,
                       double *TSurfaces, double TEnvironment,
                       double *FluxVector)
{
  /*--------------------------------------------------------------*/
  /*- quantities arising from sources inside object nss are       */
  /*- weighted by a factor [Theta(T) - Theta(TEnv)]               */
  /*- note: nss = 'num surface, source'                           */
  /*-       nsd = 'num surface, destination'                      */
  /*--------------------------------------------------------------*/
  int NS = SNEQD->G->NumSurfaces;
  int NT = SNEQD->NumTransformations;
  int NQ = SNEQD->NQ;
  for(int nss=0; nss<NS; nss++)
   { double DeltaTheta = Theta(Omega, TSurfaces[nss]) - Theta(Omega, TEnvironment);
     for(int nt=0; nt<NT; nt++)
      for(int nsd=0; nsd<NS; nsd++)
       for(int nq=0; nq<NQ; nq++)
        FluxVector[ GetIndex(SNEQD, nt, nss, nsd, nq) ]*= DeltaTheta/M_PI;
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  FILE *f=vfopen("%s.integrand","a",SNEQD->FileBase);
  for(int nt=0; nt<NT; nt++)
   for(int nss=0; nss<NS; nss++)
    for(int nsd=0; nsd<NS; nsd++)
     { 
       fprintf(f,"%s %e %i%i ",SNEQD->GTCList[nt]->Tag,Omega,nss+1,nsd+1);
       for(int nq=0; nq<NQ; nq++)
        fprintf(f,"%.8e ",FluxVector[ GetIndex(SNEQD, nt, nss, nsd, nq) ] );
       fprintf(f,"\n");
     };
  fclose(f);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct FIData 
 {
   SNEQData *SNEQD;
   double OmegaMin;
   int Infinite;
   double *TSurfaces;
   double TEnvironment;

 } FIData;

void SGJCIntegrand(unsigned ndim, const double *x, void *params,
                   unsigned fdim, double *fval)
{
  (void) ndim; // unused
  (void) fdim;

  FIData *FID         = (FIData *)params;

  SNEQData *SNEQD     = FID->SNEQD; 
  double OmegaMin     = FID->OmegaMin;
  int Infinite        = FID->Infinite;
  double *TSurfaces    = FID->TSurfaces;
  double TEnvironment = FID->TEnvironment;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double Omega, Jacobian; 
  if(Infinite)
   { Omega    = OmegaMin + x[0] / (1.0-x[0]);
     Jacobian = 1.0 / ( (1.0-x[0]) * (1.0-x[0]) );
   }
  else
   { Omega    = x[0];
     Jacobian = 1.0;
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  GetFlux(SNEQD, Omega, fval);
  for(unsigned int nf=0; nf<fdim; nf++)
   fval[nf]*=Jacobian;
  PutInThetaFactors(SNEQD, Omega, TSurfaces, TEnvironment, fval);
}

  /*--------------------------------------------------------------*/
  /*- quantities arising from sources inside object nsp are       */
  /*- weighted by a factor [Theta(T) - Theta(TEnv)]                */
  /*--------------------------------------------------------------*/
#if 0
  int nt, nq, ns, nsp;
  int NS = SNEQD->G->NumSurfaces;
  int NT = SNEQD->NumTransformations;
  int NQ = SNEQD->NQ;
  int NSNQ = NS*NQ;
  int NS2NQ = NS*NS*NQ;
  double DeltaTheta;
  for(nsp=0; nsp<NS; nsp++)
   { DeltaTheta = Theta(Omega, TSurfaces[nsp]) - Theta(Omega, TEnvironment);
     for(nt=0; nt<NT; nt++)
      for(ns=0; ns<NS; ns++)
       for(nq=0; nq<NQ; nq++)
        fval[ nt*NS2NQ + ns*NSNQ + nsp*NQ + nq ] *= Jacobian*DeltaTheta/M_PI;
   };
#endif

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WriteDataToOutputFile(SNEQData *SNEQD, double *I, double *E)
{
  time_t MyTime;
  struct tm *MyTm;
  char TimeString[30];
  
  MyTime=time(0);
  MyTm=localtime(&MyTime);
  strftime(TimeString,30,"%D::%T",MyTm);
  FILE *f=vfopen("%s.out","a",SNEQD->FileBase);
  fprintf(f,"\n");
  fprintf(f,"# scuff-neq run on %s (%s)\n",GetHostName(),TimeString);
  fprintf(f,"# data file columns: \n");
  fprintf(f,"# 1 transform tag\n");
  fprintf(f,"# 2 (sourceObject, destObject) \n");
  int nq=3;
  if (SNEQD->QuantityFlags & QFLAG_POWER) 
   fprintf(f,"# (%i,%i) power (value,error)\n",nq++,nq++);
  if (SNEQD->QuantityFlags & QFLAG_XFORCE) 
   fprintf(f,"# (%i,%i) x-force (value,error)\n",nq++,nq++);
  if (SNEQD->QuantityFlags & QFLAG_YFORCE) 
   fprintf(f,"# (%i,%i) y-force (value,error)\n",nq++,nq++);
  if (SNEQD->QuantityFlags & QFLAG_ZFORCE) 
   fprintf(f,"# (%i,%i) z-force (value,error)\n",nq++,nq++);

  int NS = SNEQD->G->NumSurfaces;
  int NT = SNEQD->NumTransformations;
  int NQ = SNEQD->NQ;
  double TotalQuantity[MAXQUANTITIES], TotalError[MAXQUANTITIES];
  for(int nt=0; nt<NT; nt++)
   for(int nsd=0; nsd<NS; nsd++)
    { 
      memset(TotalQuantity,0,NQ*sizeof(double));
      memset(TotalError,   0,NQ*sizeof(double));
      for(int nss=0; nss<NS; nss++)
       { fprintf(f,"%s %i%i ",SNEQD->GTCList[nt]->Tag,nss+1,nsd+1);
         for(nq=0; nq<NQ; nq++)
          { int i = GetIndex(SNEQD, nt, nss, nsd, nq);
            fprintf(f,"%+16.8e %+16.8e ", I[i], E[i] );
            TotalQuantity[nq] += I[i];
            TotalError[nq] += E[i];
          };
         fprintf(f,"\n");
       };

      fprintf(f,"%s 0%i ",SNEQD->GTCList[nt]->Tag,nsd+1);
      for(nq=0; nq<NQ; nq++)
       fprintf(f,"%e %e ",TotalQuantity[nq],TotalError[nq]);
      fprintf(f,"\n");

    };
  fclose(f);   
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void EvaluateFrequencyIntegral(SNEQData *SNEQD, 
                               double OmegaMin, double OmegaMax,
                               double *TSurfaces, double TEnvironment, 
                               double AbsTol, double RelTol,
                               double *I, double *E)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  FIData MyFIData, *FID=&MyFIData;
  FID->SNEQD = SNEQD;
  FID->OmegaMin=OmegaMin;

  if (OmegaMax==-1.0)
   { FID->Infinite=1; 
     OmegaMin=0.0;
     OmegaMax=1.0;
   }
  else
   FID->Infinite=0;

  FID->TSurfaces=TSurfaces;
  FID->TEnvironment=TEnvironment;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  RWGGeometry *G = SNEQD -> G;
  int NS = G->NumSurfaces;
  int NT = SNEQD->NumTransformations;
  int NQ = SNEQD->NQ;
  int fdim = NT*NS*NS*NQ;
  adapt_integrate_log(fdim, SGJCIntegrand, (void *)FID, 1, &OmegaMin, &OmegaMax,
                      1000, AbsTol, RelTol, I, E, "scuff-neq.SGJClog",15);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  WriteDataToOutputFile(SNEQD, I, E);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
#if 0
void EvaluateFrequencyIntegral2(SNEQData *SNEQD, double OmegaMin, double OmegaMax,
                                double *TSurfaces, double TEnvironment, int NumIntervals,
                                double *I, double *E)
{ 
  bool UseVariableTransformation=true;

  double uMin, uMax;
  if ( OmegaMax == -1.0 )
   { UseVariableTransformation=true;
     uMin=0.0;
     uMax=1.0;
   }
  else
   { UseVariableTransformation=false;
     uMin=OmegaMin;
     uMax=OmegaMax;
   }

  int NS = SNEQD->G->NumSurfaces;
  int NT = SNEQD->NumTransformations;
  int NQ = SNEQD->NQ;
  int fdim = NT*NS*NS*NQ;
  double *fLeft  = new double[fdim];
  double *fMid   = new double[fdim];
  double *fRight = new double[fdim];
  double Omega, Jacobian;

  /*--------------------------------------------------------------*/
  /*- evaluate integrand at leftmost point -----------------------*/
  /*--------------------------------------------------------------*/
  if (UseVariableTransformation) 
   { Omega = OmegaMin + uMin / (1.0-uMin);
     Jacobian = 1.0/( (1.0-uMin)*(1.0-uMin) );
     GetFlux(SNEQD, Omega, fLeft);
     for(int n=0; n<fdim; n++) fLeft[n]*=Jacobian;
   }
  else
   { Omega=uMin;
     GetFlux(SNEQD, Omega, fLeft);
   };
  PutInThetaFactors(SNEQD, Omega, TSurfaces, TEnvironment, fLeft);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double u=uMin, Delta = (uMax - uMin) / NumIntervals;
  double *ISimp  = new double[fdim];
  double *ITrap  = new double[fdim];
  memset(I, 0, fdim*sizeof(double));
  memset(E, 0, fdim*sizeof(double));
  for(int nIntervals=0; nIntervals<NumIntervals; nIntervals++)
   { 
     // evaluate integrand at midpoint end of interval 
     u += 0.5*Delta;
     if (UseVariableTransformation) 
      { Omega = OmegaMin + u/ (1.0-u);
        Jacobian = 1.0/( (1.0-u)*(1.0-u) );
        GetFlux(SNEQD, Omega, fMid);
        for(int n=0; n<fdim; n++) fMid[n]*=Jacobian;
      }
     else
      { Omega=u;
        GetFlux(SNEQD, Omega, fMid);
      };
     PutInThetaFactors(SNEQD, Omega, TSurfaces, TEnvironment, fMid);

     // evaluate integrand at right end of interval 
     u += 0.5*Delta;
     if (UseVariableTransformation && nIntervals==NumIntervals-1) 
      { Omega=0.0;
        memset(fRight,0,fdim*sizeof(double)) ;
      }
     if (UseVariableTransformation) 
      { Omega = OmegaMin + u/ (1.0-u);
        Jacobian = 1.0/( (1.0-u)*(1.0-u) );
        GetFlux(SNEQD, Omega, fRight);
        for(int n=0; n<fdim; n++) fRight[n]*=Jacobian;
      }
     else
      { Omega=u;
        GetFlux(SNEQD, Omega, fRight);
      };
     PutInThetaFactors(SNEQD, Omega, TSurfaces, TEnvironment, fRight);

     // compute the simpson's rule and trapezoidal rule
     // estimates of the integral over this interval  
     // and take their difference as the error
     for(int nf=0; nf<fdim; nf++)
      { ISimp[nf] = (fLeft[nf] + 4.0*fMid[nf] + fRight[nf])*Delta/6.0;
        ITrap[nf] = (fLeft[nf] + 2.0*fMid[nf] + fRight[nf])*Delta/4.0;
        I[nf] += ISimp[nf];
        E[nf] += fabs(ISimp[nf] - ITrap[nf]);
      };

     // prepare for next iteration
     memcpy(fLeft, fRight, fdim*sizeof(double));

   };
  delete[] fLeft;
  delete[] fMid;
  delete[] fRight;
  delete[] ISimp;
  delete[] ITrap;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  WriteDataToOutputFile(SNEQD, I, E);

}
#endif

/***************************************************************/
/***************************************************************/
/***************************************************************/
void EvaluateFrequencyIntegral2(SNEQData *SNEQD, double OmegaMin, double OmegaMax,
                                double *TSurfaces, double TEnvironment, int NumIntervals,
                                double *I, double *E)
{ 
  int NS = SNEQD->G->NumSurfaces;
  int NT = SNEQD->NumTransformations;
  int NQ = SNEQD->NQ;
  int fdim = NT*NS*NS*NQ;
  double *fLeft  = new double[fdim];
  double *fMid   = new double[fdim];
  double *fRight = new double[fdim];
  double Omega;

  if ( OmegaMax == -1.0 )
   { 
     // if the user didn't specify an upper frequency bound, we choose 
     // OmegaMax to be the frequency at which Theta(OmegaMax, TMax) has 
     // decayed to 10^{-10} of the value of Theta(0,TMax), where TMax 
     // is the largest temperature of any object. 
     // note: the function x/(exp(x)-1) falls below 10^{-10} at x=26.3.

     double TMax=TEnvironment;
     for (int ns=0; ns<NS; ns++)
      if ( TSurfaces[ns]>=0.0 && TSurfaces[ns]>TMax ) 
       TMax=TSurfaces[ns];

     if (TMax==0.0) // all temperatures are zero; no point in calculating
      { memset(I,0,fdim*sizeof(double));
        memset(E,0,fdim*sizeof(double));
        return;
      };  

     OmegaMax = 26.3*BOLTZMANNK*TMax;
     Log("Integrating to a maximum frequency of k=%e um^{-1} (w=%e rad/sec)\n",OmegaMax,OmegaMax*3.0e14);

   };

  /*--------------------------------------------------------------*/
  /*- evaluate integrand at leftmost point.  ---------------------*/
  /*- If the leftmost point is less than OMEGAMIN (the smallest  -*/
  /*- frequency at which we can do reliable calculations) then   -*/
  /*- we estimate the integral from OmegaMin to OMEGAMIN by      -*/
  /*- a rectangular rule with integrand values computed at       -*/
  /*- OMEGAMIN.                                                  -*/
  /*--------------------------------------------------------------*/
  #define OMEGAMIN 0.01
  if (OmegaMin < OMEGAMIN)
   { 
     Omega=OMEGAMIN;
     GetFlux(SNEQD, Omega, fLeft);
     PutInThetaFactors(SNEQD, Omega, TSurfaces, TEnvironment, fLeft);
     for(int nf=0; nf<fdim; nf++)
      { I[nf] = fLeft[nf] * (OMEGAMIN-OmegaMin);
        E[nf] = 0.0;
      };
     OmegaMin=OMEGAMIN;
   }
  else
   { 
     Omega=OmegaMin;
     GetFlux(SNEQD, Omega, fLeft);
     PutInThetaFactors(SNEQD, Omega, TSurfaces, TEnvironment, fLeft);
     memset(I, 0, fdim*sizeof(double));
     memset(E, 0, fdim*sizeof(double));
   };


  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double Delta = (OmegaMax - OmegaMin) / NumIntervals;
  double *ISimp  = new double[fdim];
  double *ITrap  = new double[fdim];
  for(int nIntervals=0; nIntervals<NumIntervals; nIntervals++)
   { 
     // evaluate integrand at midpoint of interval 
     Omega += 0.5*Delta;
     GetFlux(SNEQD, Omega, fMid);
     PutInThetaFactors(SNEQD, Omega, TSurfaces, TEnvironment, fMid);

     // evaluate integrand at right end of interval 
     Omega += 0.5*Delta;
     GetFlux(SNEQD, Omega, fRight);
     PutInThetaFactors(SNEQD, Omega, TSurfaces, TEnvironment, fRight);

     // compute the simpson's rule and trapezoidal rule
     // estimates of the integral over this interval  
     // and take their difference as the error
     for(int nf=0; nf<fdim; nf++)
      { ISimp[nf] = (fLeft[nf] + 4.0*fMid[nf] + fRight[nf])*Delta/6.0;
        ITrap[nf] = (fLeft[nf] + 2.0*fMid[nf] + fRight[nf])*Delta/4.0;
        I[nf] += ISimp[nf];
        E[nf] += fabs(ISimp[nf] - ITrap[nf]);
      };

     // prepare for next iteration
     memcpy(fLeft, fRight, fdim*sizeof(double));

   };
  delete[] fLeft;
  delete[] fMid;
  delete[] fRight;
  delete[] ISimp;
  delete[] ITrap;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  WriteDataToOutputFile(SNEQD, I, E);

}
