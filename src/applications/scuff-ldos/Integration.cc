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
 * CreateSLDData -- utility routine to initialize an Data
 *                 -- structure containing all information passed
 *                 -- around between the various scuff-ldos routines
 *
 * homer reid      -- 3/2015
 *
 */

#include "libhrutil.h"
#include "libTriInt.h"
#include "scuff-ldos.h"

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WriteLDOS(SLDData *Data, cdouble Omega,
               double *Result, double *Error)
{
  FILE *f=vfopen(Data->OutFileName,"a");
  HMatrix *XMatrix=Data->XMatrix;
  for(int nx=0; nx<XMatrix->NR; nx++)
   { 
     double X[3];
     XMatrix->GetEntriesD(nx,":",X);
     int NFun = (Data->LDOSOnly ? 2 : 20);
     fprintf(f,"%e %e %e %s ", X[0],X[1],X[2], z2s(Omega));
     for(int nf=0; nf<NFun; nf++) 
      fprintf(f,"%e %e ",Result[NFun*nx+nf], Error ? Error[NFun*nx+nf] : 0.0);
     fprintf(f,"\n");
   };
  fclose(f);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int BZIntegrand_PCubature(unsigned ndim, const double *u, void *UserData,
                          unsigned fdim, double *fval)
{
  (void) ndim; // unused

  /*--------------------------------------------------------------*/
  /*- unpack fields from user data structure ---------------------*/
  /*--------------------------------------------------------------*/
  SLDData *Data = (SLDData *)UserData;
  RWGGeometry *G    = Data->G;
  cdouble Omega     = Data->Omega;

  /*--------------------------------------------------------------*/
  /*- convert (ux, uy) variable to kBloch ------------------------*/
  /*--------------------------------------------------------------*/
  int LDim=G->LDim;
  double kBloch[2]={0.0, 0.0};
  double Jacobian;
  if (LDim==1)
   { Jacobian=2.0;
     kBloch[0] = u[0]*(Data->RLBasis[0][0]);
     kBloch[1] = 0.0; 
   }
  else // (LDim==2)
   { Jacobian=4.0;
     kBloch[0] = u[0]*(Data->RLBasis[0][0]) + u[1]*(Data->RLBasis[1][0]);
     kBloch[1] = u[0]*(Data->RLBasis[0][1]) + u[1]*(Data->RLBasis[1][1]);
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  GetLDOS(Data, Omega, kBloch, fval);

  /*--------------------------------------------------------------*/
  /*- put in the Jacobian factor; this is not really a jacobian, -*/
  /*- but just a compensation for the fact that we integrate over-*/
  /*- only half or one-fourth the full Brillouin zone.           -*/
  /*--------------------------------------------------------------*/
  for(int nf=0; nf<fdim; nf++)
   fval[nf]*=Jacobian;
 
  return 0;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetBZIntegral_PCubature(SLDData *Data, cdouble Omega)
{
  int LDim = Data->G->LDim;
  int NFun = (Data->LDOSOnly ? 2 : 20)*Data->XMatrix->NR;

  double Lower[2]={0.0, 0.0};
  double Upper[2]={0.5, 0.5};
  double *Result = new double[NFun];
  double *Error  = new double[NFun];
  Data->Omega=Omega;
  pcubature(NFun, BZIntegrand_PCubature, (void *)Data, LDim,
	    Lower, Upper, Data->MaxEvals, 0.0, Data->RelTol,
	    ERROR_INDIVIDUAL, Result, Error);

  WriteLDOS(Data, Omega, Result, Error);

  delete[] Result;
  delete[] Error;

}

/***************************************************************/
/* BZ integrand function passed to triangle cubature routines  */
/***************************************************************/
void BZIntegrand_TriCub(double *u, void *UserData, double *fval)
{
  /*--------------------------------------------------------------*/
  /*- unpack fields from user data structure ---------------------*/
  /*--------------------------------------------------------------*/
  SLDData *Data = (SLDData *)UserData;
  RWGGeometry *G    = Data->G;
  cdouble Omega     = Data->Omega;

  /*--------------------------------------------------------------*/
  /*- convert (ux, uy) variable to kBloch ------------------------*/
  /*--------------------------------------------------------------*/
  double kBloch[2]={0.0, 0.0};
  kBloch[0] = u[0]*(Data->RLBasis[0][0]) + u[1]*(Data->RLBasis[1][0]);
  kBloch[1] = u[0]*(Data->RLBasis[0][1]) + u[1]*(Data->RLBasis[1][1]);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  GetLDOS(Data, Omega, kBloch, fval);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
#ifdef HAVE_DCUTRI

void *CreateDCUTRIWorkspace(int numfun, int maxpts);
typedef void (*dcutri_func)(double *x, void *parms, double *f);
void DCUTRI(void *opW, double **Vertices, dcutri_func func, 
            void *parms, double epsabs, double epsrel, 
            double *I, double *abserr);

void GetBZIntegral_DCUTRI(SLDData *Data, cdouble Omega)
{
  int LDim = Data->G->LDim;
  int NFun = (Data->LDOSOnly ? 2 : 20)*(Data->XMatrix->NR);
  
  static void *Workspace=0;
  if (Workspace==0)
    Workspace=CreateDCUTRIWorkspace(NFun, Data->MaxEvals);

  Data->Omega=Omega;

  double V1[3]={0.0, 0.0, 0.0};
  double V2[3]={0.5, 0.0, 0.0};
  double V3[3]={0.5, 0.5, 0.0};
  double *Vertices[3]={V1, V2, V3};
  double *Result=new double[NFun];
  double *Error=new double[NFun];
  DCUTRI(Workspace, Vertices, BZIntegrand_TriCub, (void *)Data, 
         0.0, Data->RelTol, Result, Error);

  V2[0]=0.0; V2[1]=0.5;
  double *Result2=new double[NFun];
  double *Error2=new double[NFun];
  DCUTRI(Workspace, Vertices, BZIntegrand_TriCub, (void *)Data, 
         0.0, Data->RelTol, Result2, Error2);
  for(int nf=0; nf<NFun; nf++)
   { Result[nf] += Result2[nf];
     Error[nf]  += Error2[nf];
   };

  WriteLDOS(Data, Omega, Result, Error);

  delete[] Result;
  delete[] Result2;
  delete[] Error;
  delete[] Error2;
}
#else
void GetBZIntegral_DCUTRI(SLDData *Data, cdouble Omega)
 { ErrExit("scuff-ldos was not compiled with DCUTRI support; configure --with-dcutri"); }
#endif

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetBZIntegral_FOTC(SLDData *Data, cdouble Omega, char *BZIString)
{
  int LDim        = Data->G->LDim;
  int NFun        = (Data->LDOSOnly ? 2 : 20)*(Data->XMatrix->NR);
  bool BZSymmetry = Data->BZSymmetry;

  int Order;
  sscanf(BZIString+4,"%i",&Order);
  int NumPts;
  double *TCR=GetTCR(Order, &NumPts);
  if (TCR==0)
   ErrExit("unsupported FOTC order %i",Order);

  Data->Omega=Omega;

  double V1[3]={0.0, 0.0, 0.0};
  double V2[3]={0.5, 0.0, 0.0};
  double V3[3]={0.5, 0.5, 0.0};
  double *F=new double[NFun];
  TriIntFixed(BZIntegrand_TriCub, NFun, (void *)Data,
              V1, V2, V3, Order, F);
  
  if (BZSymmetry)
   { for(int nf=0; nf<NFun; nf++)
      F[nf]*=2.0;
   }
  else
   { V2[0]=0.0; V2[1]=0.5;
     double *F2=new double[NFun];
     TriIntFixed(BZIntegrand_TriCub, NFun, (void *)Data,
                 V1, V2, V3, Order, F2);
     for(int nf=0; nf<NFun; nf++)
      F[nf] += F2[nf];
     delete[] F2;
   };

  WriteLDOS(Data, Omega, F, 0);

  delete[] F;
}
