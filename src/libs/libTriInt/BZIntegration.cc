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
 * BZIntegration.cc -- routines for computing numerical
 *                  -- Brillouin-zone integrals in scuff-EM
 *
 * homer reid       -- 3/2015
 *
 */

/*
 *
 */

#include "libhrutil.h"
#include "libTriInt.h"
#include "BZIntegration.h"

#define MAXBZDIM 3
#define MAXSTR 1000

#ifdef HAVE_DCUTRI
 void *CreateDCUTRIWorkspace(int numfun, int maxpts);
 int DCUTRI(void *opW, double **Vertices, TriIntFun func,
            void *parms, double epsabs, double epsrel,
            double *I, double *abserr);
#else

void *CreateDCUTRIWorkspace(int numfun, int maxpts)
{ (void) numfun;
  (void) maxpts;
  ErrExit("SCUFF-EM was not compiled with DCUTRI; configure --with-dcutri");
  return 0;
}

int DCUTRI(void *opW, double **Vertices, TriIntFun func,
           void *parms, double epsabs, double epsrel,
           double *I, double *abserr)
{ 
  (void)opW; 
  (void)Vertices;
  (void)func;
  (void)parms;
  (void)epsabs;
  (void)epsrel;
  (void)I;
  (void)abserr;

  CreateDCUTRIWorkspace(0,0);
  return 0;
}
#endif

/***************************************************************/
/***************************************************************/
/***************************************************************/
int BZIntegrand_PCubature(unsigned ndim, const double *u,
                          void *pArgs, unsigned fdim,
                          double *BZIntegrand)
{
  (void) ndim; // unused
  (void) fdim; // unused

  /*--------------------------------------------------------------*/
  /*- unpack fields from user data structure ---------------------*/
  /*--------------------------------------------------------------*/
  GetBZIArgStruct *Args  = (GetBZIArgStruct *)pArgs;
  cdouble Omega          = Args->Omega;
  BZIFunction *BZIFunc   = Args->BZIFunc;
  void *UserData         = Args->UserData;
  HMatrix *RLBasis       = Args->RLBasis;
  int LDim               = RLBasis->NC;

  /*--------------------------------------------------------------*/
  /*- convert (ux, uy) variable to kBloch ------------------------*/
  /*--------------------------------------------------------------*/
  double kBloch[3]={0.0, 0.0, 0.0};
  for(int nd=0; nd<LDim; nd++)
   for(int nc=0; nc<3; nc++)
    kBloch[nc] += u[nd]*RLBasis->GetEntryD(nc,nd);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  BZIFunc(UserData, Omega, kBloch, BZIntegrand);

  Args->NumCalls++;
  return 0;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetBZIntegral_PCubature(GetBZIArgStruct *Args,
                             cdouble Omega, double *BZIntegral)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  HMatrix *RLBasis   = Args->RLBasis;
  int LDim           = RLBasis->NC;
  int FDim           = Args->FDim;
  int MaxEvals       = Args->MaxEvals;
  double *BZIError   = Args->BZIError;
  double RelTol      = Args->RelTol;
  double AbsTol      = Args->AbsTol;
  bool Reduced       = Args->Reduced;
  
  double Lower[2]={0.0, 0.0};
  double Upper[2]={1.0, 1.0};
  if (Reduced)
   Upper[0]=Upper[1]=0.5;
  Args->Omega = Omega;
  pcubature(FDim, BZIntegrand_PCubature, (void *)Args, LDim,
	    Lower, Upper, MaxEvals, AbsTol, RelTol,
	    ERROR_INDIVIDUAL, BZIntegral, BZIError);
 
}

/***************************************************************/
/* BZ integrand function passed to triangle cubature routines  */
/***************************************************************/
void BZIntegrand_TriCub(double *u, void *pArgs, double *BZIntegrand)
{
  /*--------------------------------------------------------------*/
  /*- unpack fields from user data structure ---------------------*/
  /*--------------------------------------------------------------*/
  GetBZIArgStruct *Args  = (GetBZIArgStruct *)pArgs;
  cdouble Omega          = Args->Omega;
  BZIFunction *BZIFunc   = Args->BZIFunc;
  void *UserData         = Args->UserData;
  HMatrix *RLBasis       = Args->RLBasis;
  int LDim               = RLBasis->NC;

  if (LDim!=2)
   ErrExit("Triangle-cubature BZ integrators require 2D lattices");

  /*--------------------------------------------------------------*/
  /*- convert (ux, uy) variable to kBloch ------------------------*/
  /*--------------------------------------------------------------*/
  double kBloch[3]={0.0, 0.0, 0.0};
  for(int nd=0; nd<LDim; nd++)
   for(int nc=0; nc<3; nc++)
    kBloch[nc] += u[nd]*RLBasis->GetEntryD(nc,nd);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  BZIFunc(UserData, Omega, kBloch, BZIntegrand);

  Args->NumCalls++;
}

/***************************************************************/
/* triangle cubature                                           */
/* Order = {0, 1, 2, 4, 5, 7, 9, 13, 14, 16, 20, 25}           */
/***************************************************************/
void GetBZIntegral_TC(GetBZIArgStruct *Args, cdouble Omega,
                      double *BZIntegral)
{
  int Order              = Args->Order;
  bool BZSymmetric       = Args->BZSymmetric;
  bool Reduced           = Args->Reduced;
  int FDim               = Args->FDim;
  double *BZIError       = Args->BZIError;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  static int FDimSave=0, MaxEvalsSave=0;
  static void *Workspace=0;
  if (Order==0)
   {
     if (FDimSave!=FDim || MaxEvalsSave!=Args->MaxEvals )
      { FDimSave=FDim;
        MaxEvalsSave=Args->MaxEvals;
        if (Workspace) free(Workspace);
        Workspace=CreateDCUTRIWorkspace(FDim, Args->MaxEvals);
      };
   };
  
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double V1[3]={0.0, 0.0, 0.0};
  double V2[3]={1.0, 0.0, 0.0};
  double V3[3]={1.0, 1.0, 0.0};
  double *Vertices[3]={V1,V2,V3};
  if (Reduced)
   V2[0]=V3[0]=V3[1]=0.5;
  Args->Omega=Omega; 
  memset(BZIError, 0, FDim*sizeof(double));
  if (Order==0)
   Args->NumCalls = DCUTRI(Workspace, Vertices, 
                           BZIntegrand_TriCub, (void *) Args,
                           Args->AbsTol, Args->RelTol, 
                           BZIntegral, BZIError);
  else
   Args->NumCalls = TriIntFixed(BZIntegrand_TriCub, FDim, (void *)Args,
                                V1, V2, V3, Order, BZIntegral);
  
  if (BZSymmetric)
   { for(int nf=0; nf<FDim; nf++)
      BZIntegral[nf]*=2.0;
   }
  else
   { V2[1]=V2[0]; V2[0]=0.0;
     double *F2=new double[FDim];
     if (Order==0)
      Args->NumCalls += DCUTRI(Workspace, Vertices, 
                               BZIntegrand_TriCub, (void *) Args,
                               Args->AbsTol, Args->RelTol, 
                               F2, BZIError);
     else
      Args->NumCalls += TriIntFixed(BZIntegrand_TriCub, FDim, (void *)Args,
                                    V1, V2, V3, Order, F2);

     for(int nf=0; nf<FDim; nf++)
      BZIntegral[nf] += F2[nf];
     delete[] F2;
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetBZIntegral_CC(GetBZIArgStruct *Args, cdouble Omega,
                      double *BZIntegral)
{
  BZIFunction *BZIFunc  = Args->BZIFunc;
  void *UserData        = Args->UserData;
  bool BZSymmetric      = Args->BZSymmetric;
  int FDim              = Args->FDim;
  int LDim              = Args->RLBasis->NC;
  double AbsTol         = Args->AbsTol;
  double RelTol         = Args->RelTol;
  double *BZIError      = Args->BZIError;
  bool Reduced          = Args->Reduced; 
  int Order             = Args->Order;
  HMatrix *RLBasis      = Args->RLBasis;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (Order>0)
   {
     double *CCQR=GetCCRule(Order);
     double uMin = 0.0, uMax = Reduced ? 0.5 : 1.0;
     double uAvg = 0.5*(uMax+uMin), uDelta=0.5*(uMax-uMin);

     int ncp[MAXBZDIM], NCP=Order;
     memset(ncp, 0, MAXBZDIM*sizeof(int));
     double *dBZI=BZIError;
     memset(BZIntegral, 0, FDim*sizeof(double));
     bool Done=false;
     while(!Done)
      { 
        // advance to next d-dimensional cubature point
        double u[MAXBZDIM], w=1.0;
        for(int nd=0; nd<LDim; nd++)
         { 
           int nn = ncp[nd];
           u[nd]  = uAvg - uDelta*CCQR[2*nn + 0];
           w     *=        uDelta*CCQR[2*nn + 1];

           ncp[nd] = (nn+1)%NCP;
           if(ncp[nd]) break;
           if(nd==(LDim-1)) Done=true;
         };
   
        if (LDim==2 && BZSymmetric && ncp[1]>ncp[0])
         continue;
        if (LDim==2 && BZSymmetric && ncp[1]<ncp[0])
         w*=2.0;

        double kBloch[3]={0.0, 0.0, 0.0};
        for(int nd=0; nd<LDim; nd++)
         for(int nc=0; nc<3; nc++)
          kBloch[nc] += u[nd]*RLBasis->GetEntryD(nc,nd);

         BZIFunc(UserData, Omega, kBloch, dBZI);
         VecPlusEquals(BZIntegral, w, dBZI, FDim);
         Args->NumCalls++;
      };
     return; 
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  static int FDimSave=0;
  static int LDimSave=0;
  static double *AllValues=0, *OuterValues=0;
  if (FDimSave!=FDim || LDimSave!=LDim)
   { FDimSave=FDim;
     LDimSave=LDim;
     int NP1 = (LDim==1) ? 129 : 129*129;
     int NP2 = (LDim==1) ? 65  : 65*65;
     AllValues=(double *)realloc(AllValues, NP1*FDim*sizeof(double));
     OuterValues=(double *)realloc(OuterValues, NP2*FDim*sizeof(double));
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int pMin = 2;
  int pMax = 6;
  bool Converged = false;
  double xMin[2]={0.0, 0.0};
  double xMax[2]={1.0, 1.0};
  if (Reduced)
   xMax[0]=xMax[1]=0.5;
  for(int p=pMin; p<=pMax && !Converged; p++)
   { 
     if (LDim==1)
      ECC( p, xMin[0], xMax[0], BZIntegrand_PCubature,
           (void *)Args, FDim, AllValues,
           (p==pMin ? 0 : OuterValues), BZIntegral, BZIError);
     else
      ECC2D( p, xMin, xMax, BZIntegrand_PCubature,
             (void *)Args, FDim,
             BZSymmetric, AllValues,
             (p==pMin ? 0 : OuterValues), BZIntegral, BZIError);

     int N=(1<<(p+1)+1);
     Args->NumCalls += (LDim==1 ? N : N*N);
     int OVSize = (1<<p)+1;
     if (LDim==2) OVSize*=OVSize;
     if (p<pMax) 
      memcpy(OuterValues, AllValues, OVSize*FDim*sizeof(double));

     /*--------------------------------------------------------------*/
     /*- convergence analysis ---------------------------------------*/
     /*--------------------------------------------------------------*/
     double MaxAbsError=0.0, MaxRelError=0.0; 
     for(int nf=0; nf<FDim; nf++)
      { 
        MaxAbsError = fmax(MaxAbsError, BZIError[nf]);
        if ( BZIError[nf] > AbsTol )
         MaxRelError = fmax(MaxRelError, BZIError[nf] / fabs(BZIntegral[nf]) );
      };
     if ( (MaxAbsError < AbsTol) || (MaxRelError < RelTol) )
      Converged=true;
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetBZIntegral(GetBZIArgStruct *Args, cdouble Omega, 
                   double *BZIntegral)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int FDim = Args->FDim;
  if (Args->BZIErrorSize!=FDim)
   { Args->BZIErrorSize=FDim;
     Args->BZIError
      =(double *)reallocEC(Args->BZIError, FDim*sizeof(double));
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Args->NumCalls=0;
  switch(Args->BZIMethod)
   { case BZI_ADAPTIVE:
      GetBZIntegral_PCubature(Args, Omega, BZIntegral);
      break;
     case BZI_CC:
      GetBZIntegral_CC(Args, Omega, BZIntegral);
      break;
     case BZI_TC:
      GetBZIntegral_TC(Args, Omega, BZIntegral);
      break;
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (Args->Reduced)
   { int LDim = Args->RLBasis->NC;
     double Factor = pow(2.0, LDim);
     for(int nf=0; nf<FDim; nf++)
     { BZIntegral[nf]*=Factor;
       Args->BZIError[nf]*=Factor;
     };
   };
} 

/***************************************************************/
/* get a basis for the reciprocal lattice.                     */
/***************************************************************/
double SetRLBasis(HMatrix *LBasis, HMatrix *RLBasis,
                  double *pLVolume, double *pRLVolume)
{
  if (pLVolume)   *pLVolume=0.0;
  if (pRLVolume) *pRLVolume=0.0;
 
  double LVolume=0.0, RLVolume=0.0;
  int LDim=LBasis->NC;

  double LBV[3][3], RLBV[3][3], LProduct=1.0;
  for(int nd=0; nd<LDim; nd++)
   { for(int nc=0; nc<3; nc++)
      LBV[nd][nc]=LBasis->GetEntryD(nc,nd);
     LProduct*=VecNorm(LBV[nd]);
   };

  double MinVolume = 1.0e-6*LProduct;

  switch(LDim)
   { case 1:
      LVolume = VecNorm(LBV[0]);
      if (LVolume==0.0) return 0;
      RLVolume = 2.0*M_PI/LVolume;
      VecScale(LBV[0], 2.0*M_PI/(LVolume*LVolume), RLBV[0]);
      break;
      
     case 2:
      VecCross(LBV[0],LBV[1],LBV[2]);
      LVolume = VecNormalize(LBV[2]);
      if (LVolume < MinVolume) return 0;
      RLVolume = 4.0*M_PI*M_PI/LVolume;
      VecCross(LBV[1],LBV[2],RLBV[0]);
      VecCross(LBV[2],LBV[0],RLBV[1]);
      VecScale(RLBV[0],2.0*M_PI/LVolume);
      VecScale(RLBV[1],2.0*M_PI/LVolume);
      break;

     case 3:
      double TV[3];
      LVolume = VecDot(VecCross(LBV[0],LBV[1],TV),LBV[2]);
      if (LVolume < MinVolume) return 0;
      RLVolume = 8.0*M_PI*M_PI*M_PI/LVolume;
      VecCross(LBV[1],LBV[2],RLBV[0]);
      VecCross(LBV[2],LBV[0],RLBV[1]);
      VecCross(LBV[0],LBV[1],RLBV[2]);
      VecScale(RLBV[0],2.0*M_PI/LVolume);
      VecScale(RLBV[1],2.0*M_PI/LVolume);
      VecScale(RLBV[2],2.0*M_PI/LVolume);
      break;
   };

  for(int nd=0; nd<LDim; nd++)
   for(int nc=0; nc<3; nc++)
    RLBasis->SetEntry(nc, nd, RLBV[nd][nc]);

  if (pLVolume)  *pLVolume=LVolume;
  if (pRLVolume) *pRLVolume=RLVolume;

  return RLVolume;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *GetRLBasis(HMatrix *LBasis,
                    double *pLVolume, double *pRLVolume)
{
  HMatrix *RLBasis = new HMatrix(LBasis);
  SetRLBasis(LBasis, RLBasis, pLVolume, pRLVolume);
  return RLBasis;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void BZIUsage(const char *format, ...)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (format)
   { va_list ap; 
     char buffer[MAXSTR];

     va_start(ap,format);
     vsnprintfEC(buffer,MAXSTR,format,ap);
     va_end(ap);

     fprintf(stderr,"error: %s (aborting)\n",buffer);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  printf("\n");
  printf("Options controlling Brillouin-zone integration: \n");
  printf("\n");
  printf("--BZIMethod   [TC | CC | Adaptive]\n");
  printf("--BZIOrder    xx \n");
  printf("--BZIRelTol   xx \n");
  printf("--BZIMaxEvals xx \n");
  printf("--BZSymmetric    \n");
  printf("--FullBZ         \n");
  printf("\n");

  exit(1);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
GetBZIArgStruct *InitBZIArgs(int argc, char **argv)
{
  /***************************************************************/
  /* allocate structure and fill in default values ***************/
  /***************************************************************/
  GetBZIArgStruct *BZIArgs = (GetBZIArgStruct *)mallocEC(sizeof(*BZIArgs));
  BZIArgs->BZSymmetric = false;
  BZIArgs->BZIMethod   = BZI_DEFAULT;
  BZIArgs->Order       = -1;
  BZIArgs->MaxEvals    = 1000;
  BZIArgs->RelTol      = 1.0e-2;
  BZIArgs->AbsTol      = 0.0;
  BZIArgs->Reduced     = true;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int narg=1; narg<argc; narg++)
   { 
     if (argv[narg]==0) 
      continue;

     char *Arg    = argv[narg];
     char *Option = ( narg == (argc-1) ? 0 : argv[narg+1]);

     // handle boolean options
     if ( !strcasecmp(Arg,"--BZSymmetric") )
      { BZIArgs->BZSymmetric = true; 
        argv[narg]=0;
        continue;
      };

     if ( !strcasecmp(Arg,"--FullBZ") )
      { BZIArgs->Reduced = false;
        argv[narg]=0;
        continue;
      };

     if ( !strcasecmp(Arg,"--BZIMethod") )
      { if (Option==0) 
         ErrExit("--BZIMethod requires an argument");
        if (!strcasecmp(Option,"TC")) 
         BZIArgs->BZIMethod = BZI_TC;
        else if (!strcasecmp(Option,"CC")) 
         BZIArgs->BZIMethod = BZI_CC;
        else if (!strcasecmp(Option,"Adaptive")) 
         BZIArgs->BZIMethod = BZI_ADAPTIVE;
        else
         ErrExit("unknown BZIMethod %s",Option);
        argv[narg]=argv[narg+1]=0;
        continue;
      };

     if ( !strcasecmp(Arg,"--BZIOrder") )
      { if (Option==0) 
         ErrExit("--BZIOrder requires an argument");
        if (1!=sscanf(Option,"%i",&(BZIArgs->Order)))
         ErrExit("invalid BZIOrder %s",Option);
        argv[narg]=argv[narg+1]=0;
        continue;
      };

     if ( !strcasecmp(Arg,"--BZIRelTol") )
      { if (Option==0) 
         ErrExit("--BZIRelTol requires an argument");
        if (1!=sscanf(Option,"%le",&(BZIArgs->RelTol)))
         ErrExit("invalid BZIRelTol %s",Option);
        argv[narg]=argv[narg+1]=0;
        continue;
      };

     if ( !strcasecmp(Arg,"--BZIMaxEvals") )
      { if (Option==0)
         ErrExit("--BZIMaxEvals requires an argument");
        if (1!=sscanf(Option,"%i",&(BZIArgs->MaxEvals)))
         ErrExit("invalid BZIRelTol %s",Option);
        argv[narg]=argv[narg+1]=0;
        continue;
      };

   }; // for(int narg=1; narg<argc; narg++)

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  return BZIArgs;
}

void UpdateBZIArgs(GetBZIArgStruct *Args,
                   HMatrix *RLBasis, double RLVolume)
{
  Args->RLBasis         = RLBasis;
  Args->BZVolume        = RLVolume;
  
  int LDim = RLBasis->NC;

  if (Args->BZIMethod==BZI_DEFAULT)
   Args->BZIMethod = (LDim==1 ? BZI_CC : BZI_TC);
  if (Args->Order==-1)
   Args->Order  = (LDim==1 ? 21     : 9);

  Log("Evaluating BZ integral by ");
   if (Args->BZIMethod==BZI_TC)
  LogC("triangle cubature, order %i: ",Args->Order);
   else if (Args->BZIMethod==BZI_CC)
  LogC("Clenshaw-Curtis cubature, order %i: ",Args->Order);
   else if (Args->BZIMethod==BZI_ADAPTIVE)
  LogC("adaptive cubature, {relTol, maxEvals}={%e,%i}",
        Args->RelTol, Args->MaxEvals);

}
