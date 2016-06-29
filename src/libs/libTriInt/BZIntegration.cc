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
int CCCubature(int Order, unsigned fdim, integrand f, void *fdata,
	       unsigned dim, const double *xmin, const double *xmax, 
	       size_t maxEval, double reqAbsError, double reqRelError,
               error_norm norm, double *Integral, double *Error)
{
  if (Order==0)
   return pcubature(fdim, f, fdata, dim, xmin, xmax, maxEval,
                    reqAbsError, reqRelError, norm, Integral, Error);

  double *CCQR = GetCCRule(Order);
  if (!CCQR) 
   ErrExit("invalid CCRule order (%i) in CCCubature",Order);

  if (dim>MAXBZDIM) 
   ErrExit("dimension too high in CCCubature");

  double uAvg[MAXBZDIM], uDelta[MAXBZDIM];
  for(unsigned d=0; d<dim; d++)
   { uAvg[d]   = 0.5*(xmax[d] + xmin[d]);
     uDelta[d] = 0.5*(xmax[d] - xmin[d]);
   };

  int ncp[MAXBZDIM];
  memset(ncp, 0, dim*sizeof(int));

  double *Integrand=Error;
  memset(Integral, 0, fdim*sizeof(double));
  bool Done=false;
  int nCalls=0;
  while(!Done)
   { 
     // get d-dimensional cubature point and weight
     double u[MAXBZDIM], w=1.0;
     for(unsigned nd=0; nd<dim; nd++)
      { u[nd]  = uAvg[nd] - uDelta[nd]*CCQR[2*ncp[nd] + 0];
           w  *=            uDelta[nd]*CCQR[2*ncp[nd] + 1];
      }; 

     // advance to next d-dimensional cubature point
     for(unsigned nd=0; nd<dim; nd++)
      { ncp[nd] = (ncp[nd]+1)%Order;
        if(ncp[nd]) break;
        if(nd==(dim-1)) Done=true;
      };

     f(dim, u, fdata, fdim, Integrand);
     VecPlusEquals(Integral, w, Integrand, fdim);
     nCalls++;
   };

  return nCalls;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetBZLimits(int SymmetryFactor, double Lower[2], double Upper[2])
{
  switch(SymmetryFactor)
   { case 1: Lower[0]=-0.5; Lower[1]=-0.5;
             Upper[0]=+0.5; Upper[1]=+0.5;
             return;

     case 2: Lower[0]= 0.0; Lower[1]=-0.5;
             Upper[0]=+0.5; Upper[1]=+0.5;
             return;

     case 4: Lower[0]= 0.0; Lower[1]= 0.0;
             Upper[0]=+0.5; Upper[1]=+0.5;
             return;

     case 8: Lower[0]= 0.0; Lower[1]= 0.0;
             Upper[0]=+0.5; Upper[1]=+1.0;
             return;
   };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int BZIntegrand_PCubature(unsigned ndim, const double *u,
                          void *pArgs, unsigned fdim,
                          double *BZIntegrand)
{
  (void) ndim; // unused

  /*--------------------------------------------------------------*/
  /*- unpack fields from user data structure ---------------------*/
  /*--------------------------------------------------------------*/
  GetBZIArgStruct *Args  = (GetBZIArgStruct *)pArgs;
  cdouble Omega          = Args->Omega;
  BZIFunction BZIFunc    = Args->BZIFunc;
  void *UserData         = Args->UserData;
  HMatrix *RLBasis       = Args->RLBasis;
  int SymmetryFactor     = Args->SymmetryFactor;
  int LDim               = RLBasis->NC;

  /*--------------------------------------------------------------*/
  /* if SymmetryFactor==8 we are integrating over a triangle,     */
  /* not a square, so we have to make a variable transformation   */
  /*--------------------------------------------------------------*/
  double Jacobian=1.0, uVector[2];
  uVector[0] = u[0];
  uVector[1] = u[1];
  if (SymmetryFactor==8)
   { uVector[1]*=uVector[0];
     Jacobian*=uVector[0];
   };

  /*--------------------------------------------------------------*/
  /*- convert (ux, uy) variable to kBloch ------------------------*/
  /*--------------------------------------------------------------*/
  double kBloch[3]={0.0, 0.0, 0.0};
  for(int nd=0; nd<LDim; nd++)
   for(int nc=0; nc<3; nc++)
    kBloch[nc] += uVector[nd]*RLBasis->GetEntryD(nc,nd);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  BZIFunc(UserData, Omega, kBloch, BZIntegrand);

  for(unsigned nf=0; nf<fdim; nf++)
   BZIntegrand[nf] *= Jacobian;

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
  int SymmetryFactor = Args->SymmetryFactor;
  
  double Lower[2], Upper[2];
  GetBZLimits(SymmetryFactor, Lower, Upper);
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
  BZIFunction BZIFunc    = Args->BZIFunc;
  void *UserData         = Args->UserData;
  int FDim               = Args->FDim;
  int SymmetryFactor     = Args->SymmetryFactor;
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

if (fabs(kBloch[2]>1.0e-6))
 ErrExit("%s:%i: internal error",__FILE__,__LINE__);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  static int FDimSave=0;
  static double *DeltaBZI=0;
  if (FDimSave<FDim)
   { FDimSave = FDim;
     DeltaBZI = (double *)reallocEC(DeltaBZI, FDim*sizeof(double));
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int NTheta = 8/SymmetryFactor;
  memset(BZIntegrand, 0, FDim*sizeof(double));
  for(int nTheta=0; nTheta<NTheta; nTheta++)
   { 
     double Theta = M_PI * nTheta / 4.0;
     double CosTheta = cos(Theta), SinTheta = sin(Theta);

     double RkB[3];
     RkB[0] = CosTheta*kBloch[0] - SinTheta*kBloch[1];
     RkB[1] = SinTheta*kBloch[0] + CosTheta*kBloch[1];
     RkB[2] = 0.0;

     BZIFunc(UserData, Omega, RkB, DeltaBZI);
     VecPlusEquals(BZIntegrand, 1.0, DeltaBZI, FDim);
   };

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

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetBZIntegral_CC(GetBZIArgStruct *Args, cdouble Omega,
                      double *BZIntegral)
{
  BZIFunction BZIFunc   = Args->BZIFunc;
  void *UserData        = Args->UserData;
  int FDim              = Args->FDim;
  int LDim              = Args->RLBasis->NC;
  //double AbsTol         = Args->AbsTol;
  //double RelTol         = Args->RelTol;
  double *BZIError      = Args->BZIError;
  HMatrix *RLBasis      = Args->RLBasis;
  int SymmetryFactor    = Args->SymmetryFactor;

  int FullOrder = Args->Order;
  int HalfOrder = Args->Order/2; if ( (HalfOrder%2) == 0) HalfOrder++;
  int Orders[2];
  double uMin[2], uMax[2];
  GetBZLimits(SymmetryFactor, uMin, uMax);
  Orders[0]=Orders[1]=FullOrder;
  if (SymmetryFactor>=2)
   Orders[0] = HalfOrder;
  if (SymmetryFactor>=4)
   Orders[1] = HalfOrder;
  double uAvg[2], uDelta[2];
  for(int nd=0; nd<LDim; nd++)
   { uAvg[nd]   = 0.5*(uMax[nd] + uMin[nd]);
     uDelta[nd] = 0.5*(uMax[nd] - uMin[nd]);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (FullOrder>0)
   {
     double *CCQR[2];
     CCQR[0]=GetCCRule(Orders[0]);
     CCQR[1]=GetCCRule(Orders[1]);

     int ncp[MAXBZDIM];
     memset(ncp, 0, MAXBZDIM*sizeof(int));
     double *dBZI=BZIError;
     memset(BZIntegral, 0, FDim*sizeof(double));
     bool Done=false;
     while(!Done)
      { 
        // advance to next d-dimensional cubature point
        double u[MAXBZDIM], w=1.0;
        for(int nd=0; nd<LDim; nd++)
         { u[nd]  = uAvg[nd] - uDelta[nd]*CCQR[nd][2*ncp[nd] + 0];
           w     *=            uDelta[nd]*CCQR[nd][2*ncp[nd] + 1];
         }; 

        // advance to next d-dimensional cubature point
        for(int nd=0; nd<LDim; nd++)
         { ncp[nd] = (ncp[nd]+1)%(Orders[nd]);
           if(ncp[nd]) break;
           if(nd==(LDim-1)) Done=true;
         };
   
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

#if 0
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
  double xMin[2]={-0.5, -0.5};
  double xMax[2]={+0.5, +0.5};
  if (Reduced)
   xMin[0]=xMin[1]=0.0;
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
#endif

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int BZIntegrand_Angular(unsigned ndim, const double *u, void *pArgs, 
                        unsigned fdim, double *BZIntegrand)
{ 
  (void) ndim;
  (void) fdim;

  GetBZIArgStruct *Args=(GetBZIArgStruct *)pArgs;

  BZIFunction BZIFunc = Args->BZIFunc;
  void *UserData      = Args->UserData;
  int FDim            = Args->FDim;
  double kRho         = Args->kRho;
  HMatrix *RLBasis    = Args->RLBasis; 
  int SymmetryFactor  = Args->SymmetryFactor;
  cdouble Omega       = Args->Omega;

  int NTheta=4;
  switch(SymmetryFactor)
   { case 1: NTheta=4; break; 
     case 2: NTheta=2; break; 
     case 4:
     case 8: NTheta=1; break;
   };

// FIXME
static double *Buffer=0;
if (Buffer==0)
 Buffer=(double *)mallocEC(FDim*sizeof(double));

  memset(BZIntegrand, 0, FDim*sizeof(double));
  double Gamma = RLBasis->GetEntryD(0,0);
  for(int nTheta=0; nTheta<NTheta; nTheta++)
   { double kTheta = u[0] + nTheta*0.5*M_PI;
     double kBloch[2];
     kBloch[0] = kRho * Gamma * cos(kTheta);
     kBloch[1] = kRho * Gamma * sin(kTheta);
     BZIFunc(UserData, Omega, kBloch, Buffer);
     VecPlusEquals(BZIntegrand, 1.0, Buffer, FDim);
   };
  
  return 0;

}

int BZIntegrand_Radial(unsigned ndim, const double *u, void *pArgs, 
                       unsigned fdim, double *BZIntegrand)
{
  (void) ndim;
  (void) fdim;

  GetBZIArgStruct *Args = (GetBZIArgStruct *)pArgs;

  int FDim           = Args->FDim;
  int SymmetryFactor = Args->SymmetryFactor;
  HMatrix *RLBasis   = Args->RLBasis;
  int AngularOrder   = Args->Order % 100;
  int MaxEvals       = Args->MaxEvals;
  double RelTol      = Args->RelTol;
  double AbsTol      = Args->AbsTol;
  double *BZIError   = Args->BZIError;

  double kRho = u[0];
  Args->kRho  = kRho;
  
  double Lower, Upper, kMax=0.5;
  if (kRho<=kMax)
   { Lower = 0.0;
     Upper = SymmetryFactor==8 ? 0.25*M_PI : 0.5*M_PI;
   }
  else
   { Lower = acos(kMax / kRho);
     Upper = SymmetryFactor==8 ? 0.25*M_PI : 0.5*M_PI - Lower;
   };


// FIXME
static double *Buffer=0;
if (Buffer==0)
 Buffer=(double *)mallocEC(FDim*sizeof(double));

  CCCubature(AngularOrder, FDim, BZIntegrand_Angular, pArgs, 1,
	     &Lower, &Upper, MaxEvals, AbsTol, RelTol,
	     ERROR_INDIVIDUAL, BZIntegrand, Buffer);

  for(int nf=0; nf<FDim; nf++)
   BZIntegrand[nf]*=kRho;

  return 0;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetBZIntegral_Radial(GetBZIArgStruct *Args, cdouble Omega,
                          double *BZIntegral)
{
  int FDim         = Args->FDim;
  int MaxEvals     = Args->MaxEvals;
  double AbsTol    = Args->AbsTol;
  double RelTol    = Args->RelTol;
  HMatrix *RLBasis = Args->RLBasis;
  double *BZIError = Args->BZIError;

  Args->Omega      = Omega;

  int RadialOrder  = (Args->Order) / 100;

  double Lower = 0.0;
  double Upper = 0.5*M_SQRT2;
  CCCubature(RadialOrder, FDim, BZIntegrand_Radial, (void *)Args, 1,
	     &Lower, &Upper, MaxEvals, AbsTol, RelTol,
	     ERROR_INDIVIDUAL, BZIntegral, BZIError);

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
    
  int SymmetryFactor = Args->SymmetryFactor;
  int LDim           = Args->RLBasis->NC;
  if (SymmetryFactor>=4 && LDim==1)
   { Warn("Can't have symmetry factor > 2 in 1-dimensional BZ integration! Resetting to 2.");
     SymmetryFactor=Args->SymmetryFactor=2;
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
     case BZI_RADIAL:
      GetBZIntegral_Radial(Args, Omega, BZIntegral);
      break;
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (SymmetryFactor>1)
   for(int nf=0; nf<FDim; nf++)
    { BZIntegral[nf] *= SymmetryFactor;
      Args->BZIError[nf]*=SymmetryFactor;
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
  printf("--BZIMethod   [CC | Adaptive | TC | Radial]\n");
  printf("--BZIOrder    xx \n");
  printf("--BZIRelTol   xx \n");
  printf("--BZIMaxEvals xx \n");
  printf("--BZSymmetryFactor [1|2|4|8]\n");
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
  BZIArgs->BZIMethod      = BZI_DEFAULT;
  BZIArgs->Order          = -1;
  BZIArgs->MaxEvals       = 1000;
  BZIArgs->RelTol         = 1.0e-2;
  BZIArgs->AbsTol         = 0.0;
  BZIArgs->SymmetryFactor = 1;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int narg=1; narg<argc; narg++)
   { 
     if (argv[narg]==0) 
      continue;

     char *Arg    = argv[narg];
     char *Option = ( narg == (argc-1) ? 0 : argv[narg+1]);

     if ( !strcasecmp(Arg,"--BZSymmetryFactor") )
      { 
        if (Option==0) 
         ErrExit("--BZSymmetryFactor requires an argument");
        if ( 1!=sscanf(Option,"%i",&(BZIArgs->SymmetryFactor))
            ||  (    (BZIArgs->SymmetryFactor != 1)
                  && (BZIArgs->SymmetryFactor != 2)
                  && (BZIArgs->SymmetryFactor != 4)
                  && (BZIArgs->SymmetryFactor != 8)
                )
           ) 
         ErrExit("invalid BZSymmetryFactor %s",Option);
        argv[narg]=argv[narg+1]=0;
        continue;
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
        else if (!strcasecmp(Option,"Radial"))
         BZIArgs->BZIMethod = BZI_RADIAL;
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

/***************************************************************/
/***************************************************************/
/***************************************************************/
bool LatticeIsSquare(HMatrix *LBasis)
{
  if (LBasis->NC!=2) return false;
  double LXX = LBasis->GetEntryD(0,0);
  double LXY = LBasis->GetEntryD(0,1);
  double LYX = LBasis->GetEntryD(1,0);
  double LYY = LBasis->GetEntryD(1,1);
  return ( EqualFloat(LXX,LYY) && LXY==0.0 && LYX==0.0 );
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
const char *BZIMethodNames[]=
 { "DEFAULT", "CC", "TC", "ADAPTIVE", "RADIAL" };

void UpdateBZIArgs(GetBZIArgStruct *Args,
                   HMatrix *RLBasis, double RLVolume)
{
  Args->RLBasis         = RLBasis;
  Args->BZVolume        = RLVolume;
  
  int LDim = RLBasis->NC;

  if (    (Args->BZIMethod==BZI_TC || Args->BZIMethod==BZI_RADIAL)
       && !LatticeIsSquare(RLBasis)
     ) 
   { Warn("BZ integration scheme %s is only for 2D square lattices (switching to default)",BZIMethodNames[Args->BZIMethod]);
     Args->BZIMethod = BZI_DEFAULT;
   };

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
  LogC("radial cubature, {order, relTol, maxEvals}={%i,%e,%i}",
        Args->Order, Args->RelTol, Args->MaxEvals);

}
