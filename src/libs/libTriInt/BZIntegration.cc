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
void GetOctantImage(double kBloch[2], int n, double Image[2])
{
  bool Swap    = (n%2)==1;
  double xSign = (2<= n && n<=5) ? -1.0 : 1.0;
  double ySign = (n>=4) ? -1.0:1.0;
  Image[0] = xSign * (Swap ? kBloch[1] : kBloch[0]);
  Image[1] = ySign * (Swap ? kBloch[0] : kBloch[1]);
}

/***************************************************************/
/* BZ integrand function passed to clenshaw-curtis cubature    */
/* routines                                                    */
/***************************************************************/
int BZIntegrand_CCCubature(unsigned ndim, const double *u,
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
  /*- special handling for SymmetryFactor = 8:                    */
  /*-  (a) if we are doing fixed-order CC cubature, omit          */
  /*-      points with ky > kx and halve the contributions of     */
  /*-      points with kx==ky                                     */
  /*-  (b) if we are doing adaptive cubature, use a Duffy         */
  /*-      transform to map the square integration domain into a  */
  /*-      triangle                                               */
  /*--------------------------------------------------------------*/
  double Weight=1.0;
  double uVector[3];
  memcpy(uVector, u, LDim*sizeof(double));
  if (SymmetryFactor==8)
   { if (Args->Order==0)
      { uVector[1]*=uVector[0];
        Weight=uVector[0];
      }
     else 
      { if ( EqualFloat(u[0],u[1]) ) 
         Weight=0.5;
        else if (u[1]<u[0])
         Weight=1.0;
        else // ky>kx
         { memset(BZIntegrand, 0, fdim*sizeof(double));
           return 0;
         };
      };
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
  VecScale(BZIntegrand, Weight, fdim);
  
  Args->NumCalls++;

  return 0;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetBZIntegral_CC(GetBZIArgStruct *Args,
                      cdouble Omega, double *BZIntegral)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  HMatrix *RLBasis    = Args->RLBasis;
  int LDim            = RLBasis->NC;
  int FDim            = Args->FDim;
  int Order           = Args->Order;
  int MaxEvals        = Args->MaxEvals;
  double RelTol       = Args->RelTol;
  double AbsTol       = Args->AbsTol;
  double **DataBuffer = Args->DataBuffer;
  int SymmetryFactor  = Args->SymmetryFactor;
  
  Args->Omega = Omega;

  double Lower[2], Upper[2];
  switch(SymmetryFactor)
   { case 1: Lower[0]=Lower[1]=-0.5;
             Upper[0]=Upper[1]=+0.5;
             break;
     case 2: Lower[0]=-0.5; Lower[1]=0.0;
             Upper[0]=Upper[1]=0.5;
             break;
     case 4: Lower[0]=Lower[1]=0.0;
             Upper[0]=Upper[1]=0.5;
             break;
     case 8: Lower[0]=0.0; Lower[1]=0.0;
             Upper[0]=0.5; Upper[1]=(Order==0) ? 1.0 : 0.5;
             break;
   };
  CCCubature(Order, FDim, BZIntegrand_CCCubature, (void *)Args, LDim,
	     Lower, Upper, MaxEvals, AbsTol, RelTol,
	     ERROR_INDIVIDUAL, BZIntegral, DataBuffer[0]);
  VecScale(BZIntegral, SymmetryFactor);
 
}

/***************************************************************/
/* BZ integrand function passed to triangle cubature routines  */
/***************************************************************/
void BZIntegrand_TC(double *u, void *pArgs, double *BZIntegrand)
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
  double **DataBuffer    = Args->DataBuffer;

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
  int NumOctants = 8/SymmetryFactor;
  memset(BZIntegrand, 0, FDim*sizeof(double));
  for(int n=0; n<NumOctants; n++)
   { 
     double RkB[3]={0.0, 0.0, 0.0};
     GetOctantImage(kBloch, n, RkB);
     double *DeltaBZI=DataBuffer[0];
     BZIFunc(UserData, Omega, RkB, DeltaBZI);
     VecPlusEquals(BZIntegrand, 1.0, DeltaBZI, FDim);
     Args->NumCalls++;
   };

}

/***************************************************************/
/* triangle cubature                                           */
/* Order = {0, 1, 2, 4, 5, 7, 9, 13, 14, 16, 20, 25}           */
/***************************************************************/
void GetBZIntegral_TC(GetBZIArgStruct *Args, cdouble Omega,
                      double *BZIntegral)
{
  int Order           = Args->Order;
  int FDim            = Args->FDim;
  int SymmetryFactor  = Args->SymmetryFactor;
  double **DataBuffer = Args->DataBuffer;

  Args->Omega=Omega;

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
  double V2[3]={0.5, 0.0, 0.0};
  double V3[3]={0.5, 0.5, 0.0};
  double *Vertices[3]={V1,V2,V3};
  if (Order==0)
   Args->NumCalls = DCUTRI(Workspace, Vertices,
                           BZIntegrand_TC, (void *) Args,
                           Args->AbsTol, Args->RelTol,
                           BZIntegral, DataBuffer[0]);
  else
   Args->NumCalls = TriIntFixed(BZIntegrand_TC, FDim, (void *)Args,
                                V1, V2, V3, Order, BZIntegral);

  VecScale(BZIntegral, SymmetryFactor, FDim);

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
  double *DeltaBZI    = Args->DataBuffer[2];

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double kTheta=u[0];
  double Gamma = RLBasis->GetEntryD(0,0);
  double kBloch[2];
  kBloch[0]=kRho*Gamma*cos(kTheta);
  kBloch[1]=kRho*Gamma*sin(kTheta);
  int NumOctants = 8/SymmetryFactor;
  memset(BZIntegrand, 0, FDim*sizeof(double));
  for(int n=0; n<NumOctants; n++)
   { 
     double RkB[3]={0.0, 0.0, 0.0};
     GetOctantImage(kBloch, n, RkB);
     BZIFunc(UserData, Omega, RkB, DeltaBZI);
     VecPlusEquals(BZIntegrand, 1.0, DeltaBZI, FDim);
     Args->NumCalls++;
   };  
  return 0;

}

int BZIntegrand_Radial(unsigned ndim, const double *u, void *pArgs, 
                       unsigned fdim, double *BZIntegrand)
{
  (void) ndim;
  (void) fdim;

  GetBZIArgStruct *Args = (GetBZIArgStruct *)pArgs;

  int FDim            = Args->FDim;
  int SymmetryFactor  = Args->SymmetryFactor;
  int AngularOrder    = Args->Order % 100;
  int MaxEvals        = Args->MaxEvals;
  double RelTol       = Args->RelTol;
  double AbsTol       = Args->AbsTol;
  double **DataBuffer = Args->DataBuffer;

  double kRho = u[0];
  Args->kRho  = kRho;
  
  double kMax  = 0.5;
  double Lower = (kRho<=kMax) ? 0.0 : acos(kMax / kRho);
  double Upper = 0.25*M_PI;

  if (AngularOrder==0)
   { 
      CCCubature(AngularOrder, FDim, BZIntegrand_Angular, pArgs, 1,
	         &Lower, &Upper, MaxEvals, AbsTol, RelTol,
	         ERROR_INDIVIDUAL, BZIntegrand, DataBuffer[1]);
      VecScale(BZIntegrand, SymmetryFactor, FDim);
   }
  else if ( (AngularOrder%2)==0 )
   { 
     BZIFunction BZIFunc = Args->BZIFunc;
     void *UserData      = Args->UserData;
     cdouble Omega       = Args->Omega;
     double Gamma        = Args->RLBasis->GetEntryD(0,0);
     double kBloch[3]={0.0, 0.0, 0.0};
     switch(AngularOrder)
      { case 2:  kBloch[0] = kRho*Gamma; break;
        case 4:  kBloch[1] = kRho*Gamma; break;
        case 6: 
        default: kBloch[0] = kBloch[1] = kRho*Gamma/(M_SQRT2); break;
      };
     BZIFunc(UserData, Omega, kBloch, BZIntegrand);
     VecScale(BZIntegrand, 2.0*M_PI, FDim);
   }
  else
   { memset(BZIntegrand, 0, FDim*sizeof(double));

     double DeltaTheta0 = 0.25*M_PI / AngularOrder;
     int NTheta = ceil( (Upper-Lower)/DeltaTheta0 );
     if (NTheta==0) NTheta=1;
     double DeltaTheta = (Upper-Lower) / NTheta;
     double *DeltaBZI=DataBuffer[1];
     for(int nTheta=0; nTheta<=NTheta; nTheta++)
      { double Theta = 0.25*M_PI - nTheta*DeltaTheta;
        double Weight = (nTheta==0 || nTheta==NTheta) ? 0.5 : 1.0;
        BZIntegrand_Angular(1, &Theta, pArgs, FDim, DeltaBZI);
        VecPlusEquals(BZIntegrand, Weight*SymmetryFactor*DeltaTheta, DeltaBZI, FDim);
      };
   };

  VecScale(BZIntegrand, kRho, FDim);
  return 0;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetBZIntegral_Radial(GetBZIArgStruct *Args, cdouble Omega,
                          double *BZIntegral)
{
  int FDim            = Args->FDim;
  int MaxEvals        = Args->MaxEvals;
  double AbsTol       = Args->AbsTol;
  double RelTol       = Args->RelTol;
  double **DataBuffer = Args->DataBuffer;

  Args->Omega      = Omega;

  int RadialOrder  = Args->Order / 100;

  double Lower = 0.0;
  double Upper = 0.5*M_SQRT2;
  CCCubature(RadialOrder, FDim, BZIntegrand_Radial, (void *)Args, 1,
	     &Lower, &Upper, MaxEvals, AbsTol, RelTol,
	     ERROR_INDIVIDUAL, BZIntegral, DataBuffer[0]);
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
  if (Args->BufSize<FDim)
   { Args->BufSize=FDim;
     Args->DataBuffer[0]
      =(double *)reallocEC(Args->DataBuffer[0], 3*FDim*sizeof(double));
     Args->DataBuffer[1]=Args->DataBuffer[0] + FDim;
     Args->DataBuffer[2]=Args->DataBuffer[0] + 2*FDim;
     Args->BZIError=Args->DataBuffer[0];
   };
    
  int SymmetryFactor = Args->SymmetryFactor;
  int LDim           = Args->RLBasis->NC;
  if (    SymmetryFactor!=1 && SymmetryFactor!=2 
       && SymmetryFactor!=4 && SymmetryFactor!=8
     )
   { Warn("invalid symmetry factor %i (resetting to 1)",SymmetryFactor);
     SymmetryFactor=1;
   };

  if (SymmetryFactor>=4 && LDim==1)
   { Warn("Can't have symmetry factor > 2 in 1-dimensional BZ integration! Resetting to 2.");
     SymmetryFactor=Args->SymmetryFactor=2;
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Args->NumCalls=0;
  switch(Args->BZIMethod)
   { 
     case BZI_CC:
      GetBZIntegral_CC(Args, Omega, BZIntegral);
      break;

     case BZI_TC:
      GetBZIntegral_TC(Args, Omega, BZIntegral);
      break;

     case BZI_RADIAL:
      GetBZIntegral_Radial(Args, Omega, BZIntegral);
      break;

     default:
      ErrExit("unknown BZIMethod in GetBZIntegral");
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
  printf("--BZIMethod   [CC | TC | Radial]\n");
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

  BZIArgs->BZIFunc=0;
  BZIArgs->UserData=0;
  BZIArgs->FDim=0;
  BZIArgs->RLBasis=0;
  BZIArgs->BZVolume=0.0;

  BZIArgs->BZIMethod      = BZI_DEFAULT;
  BZIArgs->Order          = -1;
  BZIArgs->MaxEvals       = DEF_BZIMAXEVALS;
  BZIArgs->RelTol         = DEF_BZIRELTOL;
  BZIArgs->AbsTol         = DEF_BZIABSTOL;
  BZIArgs->SymmetryFactor = 1;

  BZIArgs->BufSize = 0;
  BZIArgs->DataBuffer[0] = BZIArgs->DataBuffer[1] = 0;

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
 { "DEFAULT", "CC", "TC", "RADIAL" };

void UpdateBZIArgs(GetBZIArgStruct *Args,
                   HMatrix *RLBasis, double RLVolume)
{
  Args->RLBasis         = RLBasis;
  Args->BZVolume        = RLVolume;
  
  if (    (Args->BZIMethod==BZI_TC || Args->BZIMethod==BZI_RADIAL)
       && !LatticeIsSquare(RLBasis)
     ) 
   { Warn("BZ integration scheme %s is only for 2D square lattices (switching to default)",BZIMethodNames[Args->BZIMethod]);
     Args->BZIMethod = BZI_DEFAULT;
   };

  if (Args->BZIMethod==BZI_DEFAULT)
   Args->BZIMethod = BZI_CC;
  if (Args->Order==-1)
   Args->Order  = 21;

  Log("Evaluating BZ integral by ");
   if (Args->BZIMethod==BZI_TC)
  LogC("triangle cubature, order %i: ",Args->Order);
   else if (Args->BZIMethod==BZI_CC)
  LogC("Clenshaw-Curtis cubature, order %i: ",Args->Order);
   else if (Args->BZIMethod==BZI_RADIAL)
  LogC("radial cubature, radial order %i, angular order %i}",
        Args->Order/100, Args->Order%100);

}
