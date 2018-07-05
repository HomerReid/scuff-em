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
 * libBeyn.cc      -- implementation of Beyn's method for
 *                 -- nonlinear eigenvalue problems
 *
 * Homer Reid      -- 6/2016
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libhrutil.h>
#include <libhmat.h>

#include <libBeyn.h>

#define II cdouble(0.0,1.0)

cdouble zrandN(double Sigma=1.0, double Mu=0.0)
{ return cdouble(randN(Sigma,Mu), randN(Sigma,Mu)); }

/***************************************************************/
/***************************************************************/
/***************************************************************/
BeynSolver *CreateBeynSolver(int M, int L)
{
  BeynSolver *Solver= (BeynSolver *)mallocEC(sizeof(*Solver));

  Solver->M = M;
  Solver->L = L;

  int MLMax = (M>L) ? M : L;
  int MLMin = (M<L) ? M : L;

  // storage for eigenvalues and eigenvectors
  Solver->Eigenvalues  = new HVector(L, LHM_COMPLEX);
  Solver->EVErrors     = new HVector(L, LHM_COMPLEX);
  Solver->Residuals    = new HVector(L, LHM_COMPLEX);
  Solver->Eigenvectors = new HMatrix(M, L, LHM_COMPLEX);

  // storage for singular values, random VHat matrix, etc. used in algorithm
  Solver->A0           = new HMatrix(M,L,LHM_COMPLEX);
  Solver->A1           = new HMatrix(M,L,LHM_COMPLEX);
  Solver->A0Coarse     = new HMatrix(M,L,LHM_COMPLEX);
  Solver->A1Coarse     = new HMatrix(M,L,LHM_COMPLEX);
  Solver->MInvVHat     = new HMatrix(M,L,LHM_COMPLEX);
  Solver->VHat         = new HMatrix(M,L,LHM_COMPLEX);
  Solver->Sigma        = new HVector(MLMin);
  ReRandomize(Solver);

  // internal workspace: need storage for 2 MxL matrices
  // plus 3 LxL matrices
  #define MLBUFFERS 2
  #define LLBUFFERS 3
  size_t ML = MLMax*L, LL = L*L;
  size_t WorkspaceSize = (MLBUFFERS*ML + LLBUFFERS*LL)*sizeof(cdouble);
 
  Solver->Workspace = (cdouble *)mallocEC( WorkspaceSize );

  return Solver;
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void DestroyBeynSolver(BeynSolver *Solver)
{
  delete  Solver->Eigenvalues;
  delete  Solver->EVErrors;  
  delete  Solver->Eigenvectors;

  delete  Solver->A0;
  delete  Solver->A1;
  delete  Solver->A0Coarse;
  delete  Solver->A1Coarse;
  delete  Solver->MInvVHat;
  delete  Solver->Sigma;
  delete  Solver->VHat;

  free(Solver->Workspace);

  delete Solver;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ReRandomize(BeynSolver *Solver, unsigned int RandSeed)
{
  if (RandSeed==0) 
   RandSeed=time(0);
  srandom(RandSeed);
  HMatrix *VHat=Solver->VHat;
  for(int nr=0; nr<VHat->NR; nr++)
   for(int nc=0; nc<VHat->NC; nc++)
    VHat->SetEntry(nr,nc,zrandN());

}

/***************************************************************/
/* perform linear-algebra manipulations on the A0 and A1       */
/* matrices (obtained via numerical quadrature) to extract     */
/* eigenvalues and eigenvectors                                */
/***************************************************************/
int ProcessAMatrices(BeynSolver *Solver, BeynFunction UserFunc, void *UserData,
                     HMatrix *A0, HMatrix *A1, cdouble z0,
                     HVector *Eigenvalues, HMatrix *Eigenvectors=0)
{
  int M          = Solver->M;
  int L          = Solver->L;
  HVector *Sigma = Solver->Sigma;

  cdouble *MLBuffers[2], *LLBuffers[3];
  int ML = M*L, LL=L*L;
  MLBuffers[0] = Solver->Workspace;
  MLBuffers[1] = MLBuffers[0] + ML;
  LLBuffers[0] = MLBuffers[1] + ML;
  LLBuffers[1] = LLBuffers[0] + LL;
  LLBuffers[2] = LLBuffers[1] + LL;

  bool Verbose=          CheckEnv("SCUFF_BEYN_VERBOSE");
  double RankTol=1.0e-4; CheckEnv("SCUFF_BEYN_RANK_TOL",&RankTol);
  double ResTol=0.0;     CheckEnv("SCUFF_BEYN_RES_TOL",&ResTol);
 
  // A0 -> V0Full * Sigma * W0TFull' 
  Log(" Beyn: computing SVD...");
  HMatrix V0Full(M,L,LHM_COMPLEX,(void *)MLBuffers[0]);
  HMatrix W0TFull(L,L,LHM_COMPLEX,(void *)LLBuffers[0]);
  A0->SVD(Sigma, &V0Full, &W0TFull);

  // compute effective rank K (number of eigenvalue candidates)
  int K=0;
  for(int k=0; k<Sigma->N; k++)
   { if (Verbose) Log("Beyn: SV(%i)=%e",k,Sigma->GetEntryD(k));
     if (Sigma->GetEntryD(k) > RankTol )
      K++;
   }
  Log(" Beyn: %i/%i relevant singular values",K,L);
  if (K==0)
   { Warn("no singular values found in Beyn eigensolver");
     return 0;
   }

  // set V0, W0T = matrices of first K right, left singular vectors
  HMatrix V0(M,K,LHM_COMPLEX,(void *)MLBuffers[0]);
  HMatrix W0T(K,L,LHM_COMPLEX,(void *)LLBuffers[1]);
  for(int k=0; k<K; k++)
   { for(int m=0; m<M; m++) V0.SetEntry(m,k,V0Full.GetEntry(m,k));
     for(int l=0; l<L; l++) W0T.SetEntry(k,l,W0TFull.GetEntry(k,l));
   }

  // B <- V0' * A1 * W0 * Sigma^-1
  HMatrix TM2(K,L,LHM_COMPLEX,(void *)LLBuffers[0]);
  HMatrix B(K,K,LHM_COMPLEX,(void *)LLBuffers[2]);
  Log(" Multiplying V0*A1->TM...");
  V0.Multiply(A1, &TM2, "--transA C");   // TM2 <- V0' * A1
  Log(" Multiplying TM*W0T->B...");
  TM2.Multiply(&W0T, &B, "--transB C");   //  B <- TM2 * W0
  for(int m=0; m<K; m++)                  //  B <- B * Sigma^{-1}
   for(int n=0; n<K; n++)
    B.ScaleEntry(m,n,1.0/Sigma->GetEntry(n));

  // B -> S*Lambda*S'
  Log(" Eigensolving (%i,%i)",K,K);
  HVector Lambda(K,LLBuffers[0]);
  HMatrix S(K,K,LLBuffers[1]);
  B.NSEig(&Lambda, &S);

  // V0S <- V0*S
  Log(" Multiplying V0*S...");
  HMatrix V0S(M,K,MLBuffers[1]);
  V0.Multiply(&S, &V0S);
 
  Eigenvalues->Zero();
  if (Eigenvectors) Eigenvectors->Zero();
  int KRetained=0;
  for(int k=0; k<K; k++)
   { 
     cdouble  z = z0 + Lambda.GetEntry(k);
     cdouble *V = (cdouble *)V0S.GetColumnPointer(k);

     double Residual=0.0;
     if (ResTol>0.0)
      { HMatrix Vk(M,1,V);
        HMatrix MVk(M,1,MLBuffers[0]);
        UserFunc(z, UserData, &Vk, &MVk);
        Residual=VecNorm(MVk.ZM, M);
        if (Verbose) Log("Beyn: Residual(%i)=%e",k,Residual);
      }
     if (ResTol>0.0 && Residual>ResTol) continue;

    Eigenvalues->SetEntry(KRetained, z);
    if (Eigenvectors) 
     { Eigenvectors->SetEntries(":", KRetained, V);
       Solver->Residuals->SetEntry(KRetained,Residual);
     }
    KRetained++;
   }
  return KRetained;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int BeynSolve(BeynSolver *Solver,
              BeynFunction UserFunc, void *UserData,
              cdouble z0, double Rx, double Ry, int N)
{  
  /***************************************************************/
  /* force N to be even so we can simultaneously evaluate        */
  /* the integral with N/2 quadrature points                     */
  /***************************************************************/
  if ( (N%2)==1 ) N++;

  if (Rx==Ry)
   Log("Applying Beyn method with z0=%s,R=%e,N=%i...",z2s(z0),Rx,N);
  else
   Log("Applying Beyn method with z0=%s,Rx=%e,Ry=%e,N=%i...",z2s(z0),Rx,Ry,N);

  int M                 = Solver->M;
  int L                 = Solver->L;
  HMatrix *A0           = Solver->A0;
  HMatrix *A1           = Solver->A1;
  HMatrix *A0Coarse     = Solver->A0Coarse;
  HMatrix *A1Coarse     = Solver->A1Coarse;
  HMatrix *MInvVHat     = Solver->MInvVHat;
  HMatrix *VHat         = Solver->VHat;

  /***************************************************************/
  /* evaluate contour integrals by numerical quadrature to get   */
  /* A0 and A1 matrices                                          */
  /***************************************************************/
  A0->Zero();
  A1->Zero();
  A0Coarse->Zero();
  A1Coarse->Zero();
  double DeltaTheta = 2.0*M_PI / ((double)N);
  Log(" Evaluating contour integral (%i points)...",N);
  for(int n=0; n<N; n++)
   { 
     double Theta = ((double)n)*DeltaTheta;
     double CT    = cos(Theta), ST=sin(Theta);
     cdouble z1   = Rx*CT + II*Ry*ST;
     cdouble dz   = (II*Rx*ST + Ry*CT)/((double)N);

     MInvVHat->Copy(VHat);
     UserFunc(z0+z1, UserData, MInvVHat, 0);

     VecPlusEquals(A0->ZM, dz,    MInvVHat->ZM, M*L);
     VecPlusEquals(A1->ZM, z1*dz, MInvVHat->ZM, M*L);

     if ( (n%2)==0 )
      { VecPlusEquals(A0Coarse->ZM, 2.0*dz,    MInvVHat->ZM, M*L);
        VecPlusEquals(A1Coarse->ZM, 2.0*z1*dz, MInvVHat->ZM, M*L);
      }
   }

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  HVector *Eigenvalues  = Solver->Eigenvalues;
  HVector *EVErrors     = Solver->EVErrors;
  HMatrix *Eigenvectors = Solver->Eigenvectors;
  
  int K       = ProcessAMatrices(Solver, UserFunc, UserData, A0,       A1,       z0, Eigenvalues, Eigenvectors);
  int KCoarse = ProcessAMatrices(Solver, UserFunc, UserData, A0Coarse, A1Coarse, z0, EVErrors);
  Log("{K,KCoarse}={%i,%i}",K,KCoarse);
  for(int k=0; k<EVErrors->N && k<Eigenvalues->N; k++)
   { EVErrors->ZV[k] -= Eigenvalues->ZV[k];
     EVErrors->ZV[k] = cdouble( fabs(real(EVErrors->ZV[k])),
                                fabs(imag(EVErrors->ZV[k]))
                              );
   }

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  return K;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int BeynSolve(BeynSolver *Solver,
              BeynFunction UserFunction, void *UserData,
              cdouble z0, double R, int N)
{ return BeynSolve(Solver, UserFunction, UserData, z0, R, R, N); }
