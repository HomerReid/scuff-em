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
 * SSSolver.cc -- implementation of some methods in the SSSolver 
 *               -- class for electrostatics in scuff-EM
 *
 * homer reid    -- 5/2013 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <ctype.h>

#include <libhrutil.h>
#include <libSpherical.h>
#include <libTriInt.h>

#include "libscuff.h"
#include "SSSolver.h"
#include "StaticSubstrate.h"

namespace scuff {

#define MAXSTR 1000
#define MAXTOK 50

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
SSSolver::SSSolver(const char *GeoFileName, const char *SubstrateFile, int LogLevel)
{
  Substrate=0;
  TransformLabel=0;
  FileBase=0;

  G=new RWGGeometry(GeoFileName, LogLevel);
  if (G->LDim>0)
   ErrExit("periodic geometries not yet supported for electrostatics in SCUFF-EM");

  if (SubstrateFile)
   { 
     if (Substrate)
      DestroySubstrateData(Substrate);

     char *ErrMsg=0;
     Substrate=CreateSubstrateData(SubstrateFile, &ErrMsg);
     if (ErrMsg) 
      ErrExit(ErrMsg);
   };
}

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
SSSolver::~SSSolver()
{
  delete G;
}

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
HMatrix *SSSolver::AllocateBEMMatrix()
{
  int Dim = G->TotalPanels;
  return new HMatrix(Dim,Dim);
}

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
HMatrix *SSSolver::AssembleBEMMatrix(HMatrix *M)
{
  int Dim = G->TotalPanels;

  /***************************************************************/
  /* (re)allocate the matrix as necessary                        */
  /***************************************************************/
  if ( (M->NR!=Dim) || (M->NC!=Dim) )
   { Warn("wrong-size M matrix passed to AssembleBEMMatrix (resizing...)");
     delete M;
     M=0;
   };
  if (!M)
   M = new HMatrix(Dim, Dim);

  /***************************************************************/
  /* assemble the matrix one subblock at a time                  */
  /***************************************************************/
  for (int ns=0; ns<G->NumSurfaces; ns++)
   for (int nsp=0; nsp<G->NumSurfaces; nsp++)
    AssembleBEMMatrixBlock(ns, nsp, M, 
                           G->PanelIndexOffset[ns],
                           G->PanelIndexOffset[nsp]);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  return M;

}

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
void SSSolver::AssembleBEMMatrixBlock(int nsa, int nsb,
                                      HMatrix *M,
                                      int RowOffset, int ColOffset)
{
  Log("Assembling BEM matrix block (%i,%i)...",nsa,nsb);

  RWGSurface *Sa = G->Surfaces[nsa];
  RWGSurface *Sb = G->Surfaces[nsb];

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SurfType SurfaceType;
  double Delta=0.0, Lambda=0.0;
  if (Sa->IsPEC)
   SurfaceType = PEC;
  else
   { double EpsR  = real( G->RegionMPs[ Sa->RegionIndices[0] ] -> GetEps(0.0) );
     cdouble EpsRP = G->RegionMPs[ Sa->RegionIndices[1] ] -> GetEps(0.0);

     if ( real(EpsRP)==0.0 && imag(EpsRP)<=0.0 )
      { SurfaceType = LAMBDASURFACE;
        Lambda = -imag(EpsRP);
      }
     else
      { SurfaceType = DIELECTRIC;
        Delta = 2.0*(EpsR - real(EpsRP)) / (EpsR + real(EpsRP));
      };
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
#ifdef USE_OPENMP
  int NumThreads = GetNumThreads();
  Log("OpenMP multithreading (%i threads)",NumThreads);
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int npa=0; npa<Sa->NumPanels; npa++)
   for(int npb=0; npb<Sb->NumPanels; npb++)
    { 
      if (npb==0) LogPercent(npa, Sa->NumPanels);
      double MatrixEntry=0.0;
      switch(SurfaceType)
       {
         case PEC:
           MatrixEntry = GetPPI(Sa,npa,Sb,npb,0);
           break;

         case LAMBDASURFACE:
           MatrixEntry = -1.0*GetPPI(Sa,npa,Sb,npb,0);
           if (Sa==Sb && npa==npb) MatrixEntry -= Lambda*Sa->Panels[npa]->Area;
           break;

         case DIELECTRIC:
           if (Sa==Sb && npa==npb)
            MatrixEntry = Sa->Panels[npa]->Area;
           else 
            MatrixEntry = Delta * GetPPI(Sa,npa,Sb,npb,1);
           break;
       };
      M->SetEntry(RowOffset + npa, ColOffset + npb, MatrixEntry); 
    };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (Substrate)
   AddSubstrateContributionsToBEMMatrixBlock(nsa, nsb, M,
                                             RowOffset, ColOffset);

}

/***********************************************************************/
/* Computes the integral of Phi (IntType==PHIINTEGRAL) or of nHat.E    */
/* (IntType==ENORMALINTEGRAL) over the given panel, where Phi, E are   */
/* the potential and field computed by the user's function.            */
/***********************************************************************/
double GetRHSIntegral(RWGSurface *S, int np, 
                      StaticField* SF, void *UserData,
                      IntegralType IntType)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGPanel *P = S->Panels[np];
  double *V0 = S->Vertices + 3*P->VI[0];
  double *V1 = S->Vertices + 3*P->VI[1];
  double *V2 = S->Vertices + 3*P->VI[2];
  double A[3], B[3];
  VecSub(V1, V0, A);
  VecSub(V2, V0, B);

  double *nHat = P->ZHat;
  double Area  = P->Area;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int NumPts;
  double *TCR=GetTCR(20, &NumPts);
 
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double Result=0.0;
  for(int n=0, ncp=0; n<NumPts; n++)
   { 
     double u=TCR[ncp++]; 
     double v=TCR[ncp++]; 
     double w=TCR[ncp++];

     double X[3];
     X[0] = V0[0] + u*A[0] + v*B[0];
     X[1] = V0[1] + u*A[1] + v*B[1];
     X[2] = V0[2] + u*A[2] + v*B[2];

     double PhiE[4];
     SF(X, UserData, PhiE);
 
     if (IntType==PHIINTEGRAL)
      Result += w * PhiE[0];
     else
      Result += w * (nHat[0]*PhiE[1] + nHat[1]*PhiE[2] + nHat[2]*PhiE[3]);
   };

  return Result * 2.0*Area;
  
}

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
HVector *SSSolver::AllocateRHSVector()
{ 
  int Dim = G->TotalPanels;
  return new HVector(Dim);
}

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
HVector *SSSolver::AssembleRHSVector(double *Potentials, 
                                     StaticField *SF, void *UD, 
                                     HVector *RHS)
{
  int Dim = G->TotalPanels;

  /***************************************************************/
  /* (re)allocate the vector as necessary                        */
  /***************************************************************/
  if ( (RHS->N!=Dim) )
   { Warn("wrong-size vector passed to AssembleRHSVector (resizing...)");
     delete RHS;
     RHS=0;
   };
  if (!RHS)
   RHS = new HVector(Dim);

  RHS->Zero();
  Log("Computing RHS vector...");

  /***************************************************************/
  /* add contributions of conductor potentials if present        */
  /***************************************************************/
#if 0
  if(Potentials) 
   { for(int ns=0; ns<G->NumSurfaces; ns++)
      { RWGSurface *S=G->Surfaces[ns];
        if (!S->IsPEC) continue;
        int Offset = G->PanelIndexOffset[ns];
        double V = Potentials[ns];
        for(int np=0; np<S->NumPanels; np++)
         RHS->SetEntry(Offset + np, V * S->Panels[np]->Area);
      };
   };
#endif

  /***************************************************************/
  /* add contributions of external field if present **************/
  /***************************************************************/
  Log(" ...adding external field contributions to RHS vector...");
#ifdef USE_OPENMP
  int NumThreads = GetNumThreads();
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int ns=0; ns<G->NumSurfaces; ns++)
   { 
     RWGSurface *S=G->Surfaces[ns];
     int Offset = G->PanelIndexOffset[ns];

     /*--------------------------------------------------------------*/
     /*- get prefactor for this surface -----------------------------*/
     /*--------------------------------------------------------------*/
     IntegralType IntType;
     double Delta=0.0, PotentialPreFactor=0.0, IntegralPreFactor=0.0;
     if (S->IsPEC)
      { // PEC surface
        PotentialPreFactor =  1.0;
        IntegralPreFactor  = -1.0;
        IntType = PHIINTEGRAL;
      }
     else
      { double EpsR  = real( G->RegionMPs[ S->RegionIndices[0] ] -> GetEps(0.0) );
        cdouble EpsRP = G->RegionMPs[ S->RegionIndices[1] ] -> GetEps(0.0);
   
        if ( real(EpsRP)==0.0 && imag(EpsRP)<=0.0 )
         { // \lambda-surface
           PotentialPreFactor = 0.0;
           IntegralPreFactor  = 1.0;
           IntType = PHIINTEGRAL;
         }
        else
         { // dielectric surface
           Delta = 2.0*(EpsR - real(EpsRP)) / (EpsR + real(EpsRP));
           IntegralPreFactor = -Delta;
           IntType = ENORMALINTEGRAL;
         };
      };

     /*--------------------------------------------------------------*/
     /*- contributions of fixed potentials --------------------------*/
     /*--------------------------------------------------------------*/
     if ( Potentials && IntType!=ENORMALINTEGRAL )
      for(int np=0; np<S->NumPanels; np++)
       RHS->AddEntry(Offset+np, PotentialPreFactor*(S->Panels[np]->Area)*Potentials[ns]);

     /*--------------------------------------------------------------*/
     /*- contributions of external field ----------------------------*/
     /*--------------------------------------------------------------*/
     if (SF)
      for(int np=0; np<S->NumPanels; np++)
       RHS->AddEntry(Offset+np, IntegralPreFactor*GetRHSIntegral(S,np,SF,UD,IntType) );

   };

  return RHS;

}

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
HMatrix *SSSolver::GetCartesianMoments(HVector *Sigma, HMatrix *Moments)
{
  int NS = G->NumSurfaces;
  int NM = 4;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (Moments && (Moments->NR != NS || Moments->NC < NM ) )
   { Warn("invalid matrix passed to GetCartesianMoments (reallocating)");
     delete Moments;
     Moments=0;
   };
  if (Moments==0) 
   Moments = new HMatrix(NS, NM);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  for(int ns=0; ns<NS; ns++)
   { 
     RWGSurface *S=G->Surfaces[ns];
     int Offset = G->PanelIndexOffset[ns];
     double QP[4]={0.0, 0.0, 0.0, 0.0};
     for(int np=0; np<S->NumPanels; np++)
      { 
        double Charge = Sigma->GetEntryD(Offset + np) * S->Panels[np]->Area;
        double *X0    = S->Panels[np]->Centroid;

        QP[0] += Charge;
        QP[1] += Charge * X0[0];
        QP[2] += Charge * X0[1];
        QP[3] += Charge * X0[2];
      };
     Moments->SetEntry(ns, 0, QP[0]);
     Moments->SetEntry(ns, 1, QP[1]);
     Moments->SetEntry(ns, 2, QP[2]);
     Moments->SetEntry(ns, 3, QP[3]);
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  return Moments;
  
}

/***************************************************************/
/* get the spherical moments of an individual surface **********/
/***************************************************************/
HVector *SSSolver::GetSphericalMoments(HVector *Sigma, int WhichSurface, 
                                       int lMax, HVector *Moments)
{
  int NAlpha = (lMax+1)*(lMax+1);

  /*--------------------------------------------------------------*/
  /*- (re)allocate output vector as necessary --------------------*/
  /*--------------------------------------------------------------*/
  if (Moments && Moments->N != NAlpha)
   { Warn("incorrect Moments matrix passed to GetSphericalMoments (reallocating...)");
     delete Moments;
     Moments=0;
   };
  if (!Moments)
   Moments=new HVector(NAlpha);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGSurface *S=G->Surfaces[WhichSurface];
  int Offset = G->PanelIndexOffset[WhichSurface];
  double *Ylm = new double[NAlpha];
  double *rlArray= new double[NAlpha];
  Moments->Zero();
  for(int np=0; np<S->NumPanels; np++)
   { 
     double Charge = Sigma->GetEntryD(Offset + np) * S->Panels[np]->Area;

     double *X0    = S->Panels[np]->Centroid;
     double r, Theta, Phi;
     CoordinateC2S(X0, &r, &Theta, &Phi);
     GetRealYlmArray(lMax, Theta, Phi, Ylm);

     double rl=1.0;
     for(int l=0, Alpha=0; l<=lMax; l++, rl*=r)
      for(int m=-l; m<=l; m++, Alpha++)
       rlArray[Alpha] = rl / (2.0*l+1.0);

     for(int Alpha=0; Alpha<NAlpha; Alpha++)
      Moments->AddEntry(Alpha, Charge*rlArray[Alpha]*Ylm[Alpha]);
   };
  delete[] Ylm;
  delete[] rlArray;

  return Moments;
  
}

/***************************************************************/
/* get spherical moments of entire geometry                    */
/***************************************************************/
HVector *SSSolver::GetSphericalMoments(HVector *Sigma, 
                                       int lMax, 
                                       HVector *Moments)
{
  int NAlpha = (lMax+1)*(lMax+1);

  /*--------------------------------------------------------------*/
  /*- (re)allocate output vector as necessary --------------------*/
  /*--------------------------------------------------------------*/
  if (Moments && Moments->N != NAlpha)
   { Warn("incorrect Moments matrix passed to GetSphericalMoments (reallocating...)");
     delete Moments;
     Moments=0;
   };
  if (!Moments)
   Moments=new HVector(NAlpha);
  Moments->Zero();

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  HVector *PartialMoments = new HVector(NAlpha);
  for(int ns=0; ns<G->NumSurfaces; ns++)
   { GetSphericalMoments(Sigma, ns, lMax, PartialMoments);
     for(int Alpha=0; Alpha<NAlpha; Alpha++)
      Moments->AddEntry(Alpha,PartialMoments->GetEntry(Alpha));
   };

  delete PartialMoments;
  return Moments;
}

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
void SSSolver::PlotChargeDensity(HVector *Sigma, const char *format, ...)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  va_list ap;
  char FileName[100];
  va_start(ap,format);
  vsnprintfEC(FileName,100,format,ap);
  va_end(ap);
  FILE *f=fopen(FileName,"a");
  if (!f) return;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGSurface *S;
  int np;
  fprintf(f,"View \"%s\" {\n","Surface Charge Density");
  for(int ns=0, nbf=0; ns<G->NumSurfaces; ns++)
   for(S=G->Surfaces[ns], np=0; np<S->NumPanels; np++, nbf++)
    { 
      RWGPanel *P=S->Panels[np];
      double *PV[3];
      PV[0]=S->Vertices + 3*P->VI[0];
      PV[1]=S->Vertices + 3*P->VI[1];
      PV[2]=S->Vertices + 3*P->VI[2];
      double Val = Sigma->GetEntryD( nbf );
      fprintf(f,"ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n",
                   PV[0][0], PV[0][1], PV[0][2],
                   PV[1][0], PV[1][1], PV[1][2],
                   PV[2][0], PV[2][1], PV[2][2],
                   Val, Val, Val);
   }; 
  fprintf(f,"};\n");
  fclose(f);
 
}

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
HMatrix *SSSolver::GetFields(StaticField *SF, void *UserData, HVector *Sigma, HMatrix *X, HMatrix *PhiE)
{

  int NX = X->NR;

  /*--------------------------------------------------------------*/
  /*- (re)allocate PhiE matrix as necessary ----------------------*/
  /*--------------------------------------------------------------*/
  if ( PhiE && (PhiE->NR!=NX || PhiE->NC!=4) )
   { Warn("GetFields() called with incorrect PhiE matrix (reallocating)");
     delete PhiE;
     PhiE=0;
   };
  if (PhiE==0)
   PhiE = new HMatrix(NX, 4);
  
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  PhiE->Zero();
#ifdef USE_OPENMP
  int NumThreads = GetNumThreads();
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int nx=0; nx<NX; nx++)
   { 
     double R[3];
     R[0]=X->GetEntryD(nx,0);
     R[1]=X->GetEntryD(nx,1);
     R[2]=X->GetEntryD(nx,2);

     /*--------------------------------------------------------------*/
     /*- contribution of external field if present ------------------*/
     /*--------------------------------------------------------------*/
     if (SF)
      { 
        double DeltaPhiE[4];
        SF(R, UserData, DeltaPhiE);
        PhiE->SetEntry(nx, 0, DeltaPhiE[0] );
        PhiE->SetEntry(nx, 1, DeltaPhiE[1] );
        PhiE->SetEntry(nx, 2, DeltaPhiE[2] );
        PhiE->SetEntry(nx, 3, DeltaPhiE[3] );
      };

     /*--------------------------------------------------------------*/
     /*- contribution of charges on object surfaces -----------------*/
     /*--------------------------------------------------------------*/
     for(int nbf=0, ns=0; ns<G->NumSurfaces; ns++)
      for(int np=0; np<G->Surfaces[ns]->NumPanels; np++, nbf++)
       {  
         double DeltaPhiE[4];
         GetPhiE(ns,np,R,DeltaPhiE);
         double ChargeDensity=Sigma->GetEntryD(nbf) / (4.0*M_PI);
         PhiE->AddEntry(nx, 0, ChargeDensity * DeltaPhiE[0] );
         PhiE->AddEntry(nx, 1, ChargeDensity * DeltaPhiE[1] );
         PhiE->AddEntry(nx, 2, ChargeDensity * DeltaPhiE[2] );
         PhiE->AddEntry(nx, 3, ChargeDensity * DeltaPhiE[3] );
       };

     };

  return PhiE;
  
}

} // namespace scuff
