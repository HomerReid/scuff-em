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
 * Quadrature.cc -- scuff-neq module for numerical quadrature
 *               -- over frequencies
 *
 * homer reid    -- 5/2012
 *
 */

static const char *SIQNames[]=
 {"PAbs", "PRad", "Fx", "Fy", "Fz", "Tx", "Ty", "Tz"};

#include <stdlib.h>
#include <libSGJC.h>
#include <libTriInt.h>
#include "scuff-neq.h"

// how this works: 
//  a. temperature in eV = kT = 8.6173e-5 * (T in kelvin)
//  b. temperature in our internal energy units
//     = (kT in eV) / (0.1973 eV)
//     = (8.6173e-5 / 0.1973 ) * (T in Kelvin)
#define BOLTZMANNK 4.36763e-4

// \hbar * \omega_0^2 in units of watts
#define HBAROMEGA02 9.491145534e-06

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
                       double *FluxVector)
{
  double TEnvironment  = SNEQD->TEnvironment;
  double *TSurfaces    = SNEQD->TSurfaces;

  /*--------------------------------------------------------------*/
  /*- quantities arising from sources inside object nss are       */
  /*- weighted by a factor of                                     */
  /*- \hbar \omega_0^2 [Theta(T) - Theta(TEnv)]                   */
  /*- note: nss = 'num surface, source'                           */
  /*-       nsd = 'num surface, destination'                      */
  /*--------------------------------------------------------------*/
  int NS = SNEQD->G->NumSurfaces;
  int NT = SNEQD->NumTransformations;
  int NQ = SNEQD->NPFT;
  int NX = SNEQD->NX;
  for(int nss=0; nss<NS; nss++)
   { 
     double DeltaTheta
      = HBAROMEGA02 * Theta(Omega, TSurfaces[nss]) - Theta(Omega, TEnvironment);

     for(int nt=0; nt<NT; nt++)
      { 
        for(int nq=0; nq<NQ; nq++)
         for(int nsd=0; nsd<NS; nsd++)
          FluxVector[ GetSIQIndex(SNEQD, nt, nss, nsd, nq) ]
           *= DeltaTheta;

        for(int nx=0; nx<NX; nx++)
         for(int nfc=0; nfc<NUMSRFLUX; nfc++)
          FluxVector[ GetSRQIndex(SNEQD, nt, nss, nx, nfc) ]
           *= DeltaTheta;
      };
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (SNEQD->NumSIQs>0)
   { 
     FILE *f=vfopen("%s.SIIntegrand","a",SNEQD->FileBase);
     for(int nt=0; nt<NT; nt++)
      for(int nss=0; nss<NS; nss++)
       for(int nsd=0; nsd<NS; nsd++)
        { 
          fprintf(f,"%s %e %i%i ",SNEQD->GTCList[nt]->Tag,Omega,nss+1,nsd+1);
          for(int nq=0; nq<NQ; nq++)
           fprintf(f,"%.8e ",FluxVector[ GetSIQIndex(SNEQD, nt, nss, nsd, nq) ] );
          fprintf(f,"\n");
        };
     fclose(f);
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (SNEQD->NumSRQs>0)
   { 
     FILE *f=vfopen("%s.SRIntegrand","a",SNEQD->FileBase);
     for(int nt=0; nt<NT; nt++)
      for(int nss=0; nss<NS; nss++)
       for(int nx=0; nx<NX; nx++)
        { 
          double X[3];
          SNEQD->SRXMatrix->GetEntriesD(nx,"0:2",X);

          fprintf(f,"%s %e %i ",SNEQD->GTCList[nt]->Tag,Omega,nss+1);
          fprintf(f,"%e %e %e ",X[0],X[1],X[2]);
          for(int nfc=0; nfc<NUMSRFLUX; nfc++)
           fprintf(f,"%.8e ",FluxVector[ GetSRQIndex(SNEQD, nt, nss, nx, nfc) ] );
          fprintf(f,"\n");
        };
     fclose(f);
   };
   
}
   
/***************************************************************/
/* data structure for GetOmegaIntegrand ************************/
/***************************************************************/
typedef struct GOIData 
 {
   SNEQData *SNEQD;
   double OmegaMin;
   bool Infinite;

 } GOIData;

int GetOmegaIntegrand(unsigned ndim, const double *x, void *params,
                      unsigned fdim, const bool *Skip, double *fval)
{
  (void) ndim; // unused

  GOIData *Data       = (GOIData *)params;

  SNEQData *SNEQD     = Data->SNEQD;
  double OmegaMin     = Data->OmegaMin;
  bool Infinite       = Data->Infinite;

  if (Skip)
   memcpy(SNEQD->OmegaConverged, Skip, fdim*sizeof(bool));
  else
   memset(SNEQD->OmegaConverged, 0, fdim*sizeof(bool));

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
#define MINOMEGA 1.0e-4
  if (Omega<MINOMEGA)
   Omega=MINOMEGA;
  GetFlux(SNEQD, Omega, fval);
  for(unsigned int nf=0; nf<fdim; nf++)
   fval[nf]*=Jacobian;
  PutInThetaFactors(SNEQD, Omega, fval);

  return 0;
}

int GetOmegaIntegrand2(unsigned ndim, const double *x, void *params,
                       unsigned fdim, double *fval)
{
  return GetOmegaIntegrand(ndim, x, params, fdim, 0, fval);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetOmegaIntegral_Adaptive(SNEQData *SNEQD,
                               double OmegaMin,
                               double OmegaMax,
                               double *I, double *E)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  GOIData MyGOIData, *Data=&MyGOIData;
  Data->SNEQD = SNEQD;
  Data->OmegaMin=OmegaMin;

  if (OmegaMax==-1.0)
   { Data->Infinite=true;
     OmegaMin=0.0;
     OmegaMax=1.0;
   }
  else
   Data->Infinite=false;


  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int fdim       = SNEQD->NumSIQs + SNEQD->NumSRQs;
  double AbsTol  = SNEQD->AbsTol;
  double RelTol  = SNEQD->RelTol;
  pcubature_log(fdim, GetOmegaIntegrand2, (void *)Data, 1,
                &OmegaMin, &OmegaMax, 1000,
                AbsTol, RelTol, ERROR_INDIVIDUAL,
                I, E, "scuff-neq.SGJClog");

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetOmegaIntegral_TrapSimp(SNEQData *SNEQD,
                               double OmegaMin,
                               double OmegaMax,
                               int NumIntervals,
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
   };

  int fdim = SNEQD->NumSIQs + SNEQD->NumSRQs;
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
  PutInThetaFactors(SNEQD, Omega, fLeft);

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
     PutInThetaFactors(SNEQD, Omega, fMid);

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
     PutInThetaFactors(SNEQD, Omega, fRight);

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

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetOmegaIntegral_Cliff(SNEQData *SNEQD,
                            double OmegaMin, double OmegaMax,
                            double *I, double *E)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  GOIData MyGOIData, *Data=&MyGOIData;
  Data->SNEQD    = SNEQD;
  Data->OmegaMin = OmegaMin;
  Data->Infinite = false;

  /***************************************************************/
  /* this really needs to be estimated automatically FIXME       */
  /***************************************************************/
  double OmegaCliff = 1.0;

  if (OmegaMax==-1.0)
   OmegaMax=OmegaMin;
  else if ( OmegaMax<=OmegaMin )
   ErrExit("degenerate frequency window in GetOmegaIntegral_Cliff");
  else
   { if ( ! (OmegaMin<OmegaCliff && OmegaCliff<OmegaMax) )
      OmegaCliff = 0.5*(OmegaMin + OmegaMax );
   };

  int fdim = SNEQD->NumSIQs + SNEQD->NumSRQs;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  char *ICFLogFile = vstrdup("%s.ICFLog",SNEQD->FileBase);
  IntegrateCliffFunction(GetOmegaIntegrand, (void *)Data, fdim,
                         OmegaMin, OmegaMax, OmegaCliff, 
                         SNEQD->AbsTol, SNEQD->RelTol, I, E,
                         ICFLogFile);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WriteSIOutputFile(SNEQData *SNEQD, double *I, double *E)
{
  /*--------------------------------------------------------------*/
  /*- open file and write preamble -------------------------------*/
  /*--------------------------------------------------------------*/
  char *TimeString=GetTimeString();
  FILE *f=vfopen("%s.NEQPFT","a",SNEQD->FileBase);
  fprintf(f,"# scuff-neq completed at %s (%s)\n",GetHostName(),TimeString);
  fprintf(f,"# data file columns: \n");
  fprintf(f,"# 1 transform tag\n");
  fprintf(f,"# 2 (sourceObject, destObject) \n");
  int nc=3;
  for (int nq=0; nq<NUMPFT; nq++)
   if (SNEQD->NeedQuantity[nq])
    { fprintf(f,"# (%i,%i) %s (value,error in frequency quadrature)\n",nc,nc+1,SIQNames[nq]);
      nc+=2; 
    };

  /*--------------------------------------------------------------*/
  /*- as we report the power/momentum exchange between all pairs  */
  /*- of bodies, we maintain running tallies of the total PFT for */
  /*- each destination body.                                      */
  /*--------------------------------------------------------------*/
  int NS = SNEQD->G->NumSurfaces;
  int NT = SNEQD->NumTransformations;
  int NQ = SNEQD->NPFT;
  double TotalQuantity[NUMPFT], TotalError[NUMPFT];
  for(int nt=0; nt<NT; nt++)
   for(int nsd=0; nsd<NS; nsd++)
    { 
      memset(TotalQuantity,0,NQ*sizeof(double));
      memset(TotalError,   0,NQ*sizeof(double));
      for(int nss=0; nss<NS; nss++)
       { fprintf(f,"%s %i%i ",SNEQD->GTCList[nt]->Tag,nss+1,nsd+1);
         for(int nq=0; nq<NQ; nq++)
          { int i = GetSIQIndex(SNEQD, nt, nss, nsd, nq);
            fprintf(f,"%+16.8e %+16.8e ", I[i], E[i] );
            TotalQuantity[nq] += I[i];
            TotalError[nq] += E[i];
          };
         fprintf(f,"\n");
       };

      fprintf(f,"%s 0%i ",SNEQD->GTCList[nt]->Tag,nsd+1);
      for(int nq=0; nq<NQ; nq++)
       fprintf(f,"%e %e ",TotalQuantity[nq],TotalError[nq]);
      fprintf(f,"\n");
    };
  fclose(f);   
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WriteSROutputFile(SNEQData *SNEQD, double *I, double *E)
{
  /*--------------------------------------------------------------*/
  /*- open file and write preamble -------------------------------*/
  /*--------------------------------------------------------------*/
  char *TimeString=GetTimeString();
  FILE *f=vfopen("%s.NEQFlux","a",SNEQD->FileBase);
  fprintf(f,"# scuff-neq complete at %s (%s)\n",GetHostName(),TimeString);
  fprintf(f,"# data file columns: \n");
  fprintf(f,"# 1 transform tag \n");
  fprintf(f,"# 2 source surface \n");
  fprintf(f,"#  3, 4, 5 (x,y,z) coordinates of evaluation point\n");
  fprintf(f,"#  6, 7, 8 (x,y,z) <P_x>, <P_y>, <P_z> (Poynting vector)\n");
  fprintf(f,"#  9,10,11 (x,y,z) <T_xx>, <T_xy>, <T_xz> (Maxwell tensor )\n");
  fprintf(f,"# 12,13,14 (x,y,z) <T_yx>, <T_yy>, <T_yz> (Maxwell tensor )\n");
  fprintf(f,"# 15,16,17 (x,y,z) <T_zx>, <T_zy>, <T_zz> (Maxwell tensor )\n");

  /*--------------------------------------------------------------*/
  /*- as we report the contributions of each source body to the   */
  /*- energy/momentum flux at each spatial point, we maintain a   */
  /*- running tally of the total fluxes at each point.            */
  /*--------------------------------------------------------------*/
  int NS = SNEQD->G->NumSurfaces;
  int NT = SNEQD->NumTransformations;
  int NX = SNEQD->SRXMatrix->NR;
  double TotalQuantity[NUMSRFLUX];
  for(int nt=0; nt<NT; nt++)
   for(int nx=0; nx<NX; nx++)
    { 
      double X[3];
      SNEQD->SRXMatrix->GetEntriesD(nx,"0:2",X);

      memset(TotalQuantity,0,NUMSRFLUX*sizeof(double));
      for(int nss=0; nss<NS; nss++)
       { fprintf(f,"%s %i ",SNEQD->GTCList[nt]->Tag,nss+1);
         fprintf(f,"%e %e %e ",X[0],X[1],X[2]);
         for(int nfc=0; nfc<NUMSRFLUX; nfc++)
          { int i = GetSRQIndex(SNEQD, nt, nss, nx, nfc);
            fprintf(f,"%+16.8e ", I[i]);
            TotalQuantity[nfc] += I[i];
          };
         fprintf(f,"\n");
       };

      fprintf(f,"%s 0 ",SNEQD->GTCList[nt]->Tag);
      fprintf(f,"%e %e %e ",X[0],X[1],X[2]);
      for(int nfc=0; nfc<NUMSRFLUX; nfc++)
       fprintf(f,"%e ",TotalQuantity[nfc]);
      fprintf(f,"\n");

    };
  fclose(f);   
}
