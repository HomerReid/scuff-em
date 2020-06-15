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
 * LBWrappers.cc  -- HMatrix class methods that amount to calls to 
 *                -- LAPACK / BLAS routines, possibly with some pre- or 
 *                -- post-processing 
 *
 * homer reid     -- 12/2009 -- 9/2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

#include <libhrutil.h>

extern "C" {
 #include "lapack.h" 
}

#include "libhmat.h"

/***************************************************************/
/* multiply matrix by B on the right to yield C                */
/* in other words, if this matrix is A, then this operation    */
/* sets C=A*B.                                                 */
/*                                                             */
/* Options is an optional string argument that allows you to   */
/* specify values for the options "--transA" and '--transB" to */
/* request that the matrix operands be transposed / conjugated,*/
/* like this:                                                  */
/*                                                             */
/*  Options="--transA T"    // use the transpose of A          */
/*  Options="--transB C"    // use the adjoint of B            */
/*  Options="--transA T --transB C"                            */
/***************************************************************/
void HMatrix::Multiply(HMatrix *B, HMatrix *C, const char *Options)
{
  double dOne=1.0, dZero=0.0;
  cdouble zOne=1.0, zZero=0.0;

  if ( StorageType!=B->StorageType || StorageType!=C->StorageType )
   ErrExit("%s:%i: storage type mismatch",__FILE__,__LINE__);
  if ( RealComplex!=B->RealComplex || RealComplex!=C->RealComplex)
   ErrExit("%s:%i: data type mismatch",__FILE__,__LINE__);

  /***************************************************************/ 
  /* parse Options string ****************************************/ 
  /***************************************************************/ 
  char TransA[2]="N", TransB[2]="N";
  if (Options)
   { int MaxTokens=4, NumTokens;
     char *Tokens[4];
     char Line[100];
     strncpy(Line, Options, 100);
     NumTokens=Tokenize(Line, Tokens, MaxTokens);
     for (int nt=0; nt<NumTokens; nt++)
      { if (!strcasecmp(Tokens[nt],"--transa"))    
         { if ( (nt+1) >= NumTokens ) 
            ErrExit("invalid options string in HMatrix::Multiply");
           strncpy(TransA,Tokens[nt+1],1);
           nt++;
         }
        else if (!strcasecmp(Tokens[nt],"--transb"))    
         { if ( (nt+1) >= NumTokens ) 
            ErrExit("invalid options string in HMatrix::Multiply");
           strncpy(TransB,Tokens[nt+1],1);
           nt++;
         }
        else 
         ErrExit("invalid options string in HMatrix::Multiply");
      };
   };

  /***************************************************************/ 
  /* sanity check on matrix sizes: if operand dimensions are     */ 
  /* (NRAxNCA), (NRBxNCB), (NRCxNCC), then we must have          */
  /* NCA=NRB, NRA=NRC, NCB=NCC.                                  */
  /***************************************************************/ 
  int NRA = (toupper(TransA[0])=='N') ? NR : NC;
  int NCA = (toupper(TransA[0])=='N') ? NC : NR;
  int NRB = (toupper(TransB[0])=='N') ? B->NR : B->NC;
  int NCB = (toupper(TransB[0])=='N') ? B->NC : B->NR;
  int NRC = C->NR;
  int NCC = C->NC;
  if ( (NCA!=NRB) || (NRA!=NRC) || (NCB!=NCC) )
   ErrExit("%s:%i: dimension mismatch",__FILE__,__LINE__);

  /***************************************************************/ 
  /***************************************************************/ 
  /***************************************************************/ 
  if ( RealComplex==LHM_REAL && StorageType==LHM_NORMAL )
   dgemm_(TransA,TransB,&NRC,&NCC,&NCA,&dOne,DM,&NR,B->DM,&(B->NR),&dZero,C->DM,&(C->NR));
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_NORMAL )
   zgemm_(TransA,TransB,&NRC,&NCC,&NCA,&zOne,ZM,&NR,B->ZM,&(B->NR),&zZero,C->ZM,&(C->NR));
  else
   ErrExit("%s:%i: multiplication not implemented for storage type",__FILE__,__LINE__);

}

/***************************************************************/
/* this routine is somewhat similar to Multiply(), but it only */
/* computes the diagonal of the product matrix, and it does    */
/* this in N^2 time (instead of N^3 for a usual matrix-matrix  */
/* multiplication) by dotting the rows of A into the columns   */
/* of B.                                                       */
/***************************************************************/
void HMatrix::GetMatrixProductDiagonal(HMatrix *B, HVector *DAB)
{
  if ( NC != B->NR )
   ErrExit("%s:%i: dimension mismatch in GetMatrixProductDiagonal",__FILE__,__LINE__);

  if ( RealComplex != B->RealComplex )
   ErrExit("%s:%i: real/complex mismatch in GetMatrixProductDiagonal",__FILE__,__LINE__);

  int n;
  int Length=NC;
  int AStride=NR;
  int BStride=1;

  if ( RealComplex==LHM_REAL )
   for(n=0; n<DAB->N; n++) 
    DAB->SetEntry( n, ddot_( &Length, DM + n, &AStride, B->DM + n*Length, &BStride) );
  else
   for(n=0; n<DAB->N; n++)
    DAB->SetEntry( n, zdotu_( &Length, ZM + n, &AStride, B->ZM + n*Length, &BStride) );

}

/***************************************************************/
/* replace the matrix with its LU factorization ****************/
/***************************************************************/
int HMatrix::LUFactorize()
{ 
  int info;

  if (ipiv==0)
   ipiv=(int *)mallocEC(NR*sizeof(int));

  if ( RealComplex==LHM_REAL && StorageType==LHM_NORMAL )
   dgetrf_(&NR, &NC, DM, &NR, ipiv, &info); 
  else if ( RealComplex==LHM_REAL && StorageType==LHM_SYMMETRIC )
   dsptrf_("U", &NR, DM, ipiv, &info);
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_NORMAL )
   zgetrf_(&NR, &NC, ZM, &NR, ipiv, &info);
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_HERMITIAN )
   zhptrf_("U", &NR, ZM, ipiv, &info);
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_SYMMETRIC ) 
   zsptrf_("U", &NR, ZM, ipiv, &info);

  return info;
}

/***************************************************************/
/* solve linear system using LU factorization ******************/
/***************************************************************/
int HMatrix::LUSolve(HVector *X)
{ 
  int info;
  int iOne=1;

  if ( RealComplex != X->RealComplex )
   ErrExit("type mismatch in LUSolve");
  if ( NR!=NC || NR!=X->N )
   ErrExit("dimension mismatch in LUSolve");

  if (ipiv==0)  
   ErrExit("LUFactorize() must be called before LUSolve()");

  if ( RealComplex==LHM_REAL && StorageType==LHM_NORMAL )
   dgetrs_("N", &NR, &iOne, DM, &NR, ipiv, X->DV, &NR, &info);
  else if ( RealComplex==LHM_REAL && StorageType==LHM_SYMMETRIC )
   dsptrs_("U", &NR, &iOne, DM, ipiv, X->DV, &NR, &info);
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_NORMAL )
   zgetrs_("N", &NR, &iOne, ZM, &NR, ipiv, X->ZV, &NR, &info);
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_HERMITIAN )
   zhptrs_("U", &NR, &iOne, ZM, ipiv, X->ZV, &NR, &info);
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_SYMMETRIC )
   zsptrs_("U", &NR, &iOne, ZM, ipiv, X->ZV, &NR, &info);

  return info;
}

/***************************************************************/
/* solve linear systems using LU factorization                 */
/*                                                             */
/* Trans should be 'N', 'T', or 'C' for no transpose,          */
/* transpose, or conjugate transpose.                         */
/***************************************************************/
int HMatrix::LUSolve(HMatrix *X, char Trans, int nrhs)
{ 
  int info;

  if ( RealComplex != X->RealComplex )
   ErrExit("type mismatch in LUSolve");  // problematic type mismatch
  if ( NR!=NC || NR!=X->NR )
   ErrExit("dimension mismatch in LUSolve");
  if ( nrhs > X->NC )
   ErrExit("too many RHSs requested in LUSolve");
  if (ipiv==0)  
   ErrExit("LUFactorize() must be called before LUSolve()");
  if ( Trans!='N' && StorageType!=LHM_NORMAL )
   ErrExit("transposed LU-solves not available for packed matrices");
if ( RealComplex==LHM_REAL && StorageType==LHM_NORMAL )
   dgetrs_(&Trans, &NR, &nrhs, DM, &NR, ipiv, X->DM, &NR, &info);
  else if ( RealComplex==LHM_REAL && StorageType==LHM_SYMMETRIC )
   dsptrs_("U", &NR, &nrhs, DM, ipiv, X->DM, &NR, &info);
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_NORMAL )
   zgetrs_(&Trans, &NR, &nrhs, ZM, &NR, ipiv, X->ZM, &NR, &info);
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_HERMITIAN )
   zhptrs_("U", &NR, &nrhs, ZM, ipiv, X->ZM, &NR, &info);
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_SYMMETRIC )
   zsptrs_("U", &NR, &nrhs, ZM, ipiv, X->ZM, &NR, &info);

  return info;
}

/* entry points with one or more parameters set to default */

int HMatrix::LUSolve(HMatrix *X, char Trans) 
 { return LUSolve(X,Trans,X->NC); }

int HMatrix::LUSolve(HMatrix *X, int nrhs) 
 { return LUSolve(X,'N',nrhs); }

int HMatrix::LUSolve(HMatrix *X) // problematic call, can hardly traceback where it is called 
 { return LUSolve(X,'N',X->NC); }

/***************************************************************/
/* replace the matrix with its inverse, assuming LUFactorize() */
/* has already been called                                     */
/***************************************************************/
int HMatrix::LUInvert()
{ 
  int info;
  double *dwork;
  cdouble *zwork;
  int lworkOptimal;
  double dlworkOptimal;
  cdouble zlworkOptimal;

  if ( NR!=NC )
   ErrExit("dimension mismatch in LUSolve");

  if (ipiv==0)  
   ErrExit("LUFactorize() must be called before LUInvert()");

  int MinusOne=-1;
  if ( RealComplex==LHM_REAL && StorageType==LHM_NORMAL )
   {
     // workspace size query
     dgetri_(&NR, DM, &NR, ipiv, &dlworkOptimal, &MinusOne, &info);
     lworkOptimal=(int)dlworkOptimal;

     if (lworkOptimal > lwork)
      { if (work) free(work);
        work=(double *)mallocEC(lworkOptimal*sizeof(double));
        lwork=lworkOptimal;
      };
     dwork=(double *)work;

     dgetri_(&NR, DM, &NR, ipiv, dwork, &lwork, &info);
   }
  else if ( RealComplex==LHM_REAL && StorageType==LHM_SYMMETRIC )
   { 
     lworkOptimal = NR;
     if (lworkOptimal > lwork)
      { if (work) free(work);
        work=(double *)mallocEC(lworkOptimal*sizeof(double));
        lwork=lworkOptimal;
      };
     dwork=(double *)work;

     dsptri_("U", &NR, DM, ipiv, dwork, &info);
 
   }
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_NORMAL )
   { 
     zgetri_(&NR, ZM, &NR, ipiv, &zlworkOptimal, &MinusOne, &info);
     lworkOptimal=(int) real(zlworkOptimal);

     if (lworkOptimal > lwork)
      { if (work) free(work);
        work=(cdouble *)mallocEC(lworkOptimal*sizeof(cdouble));
        lwork=lworkOptimal;
      };
     zwork=(cdouble *)work;

     zgetri_(&NR, ZM, &NR, ipiv, zwork, &lwork, &info);
   }
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_HERMITIAN )
   { 
     lworkOptimal=NR;
     if (lworkOptimal > lwork)
      { if (work) free(work);
        work=(cdouble *)mallocEC(lworkOptimal*sizeof(cdouble));
        lwork=lworkOptimal;
      };
     zwork=(cdouble *)work;
 
     zhptri_("U", &NR, ZM, ipiv, zwork, &info);
   }
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_SYMMETRIC )
   { 
     lworkOptimal=NR;
     if (lworkOptimal > lwork)
      { if (work) free(work);
        work=(cdouble *)mallocEC(lworkOptimal*sizeof(cdouble));
        lwork=lworkOptimal;
      };
     zwork=(cdouble *)work;
 
     zsptri_("U", &NR, ZM, ipiv, zwork, &info);
   };

  return info;
}

/***************************************************************/
/* replace the matrix with its cholesky factorization **********/
/***************************************************************/
int HMatrix::CholFactorize()
{ 
  int info;

  if ( RealComplex==LHM_REAL && StorageType==LHM_NORMAL )
   dpotrf_("U", &NR, DM, &NR, &info);
  else if ( RealComplex==LHM_REAL && StorageType==LHM_SYMMETRIC )
   dpptrf_("U", &NR, DM, &info);
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_NORMAL )
   zpotrf_("U", &NR, ZM, &NR, &info);
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_HERMITIAN )
   zpptrf_("U", &NR, ZM, &info);
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_SYMMETRIC )
   { fprintf(stderr,"\n*\n* ERROR: cannot cholesky-factorize a complex symmetric matrix \n*\n");
     return -1;
   };

  return info;
}

/***************************************************************/
/* solve linear system using cholesky factorization ************/
/***************************************************************/
int HMatrix::CholSolve(HVector *X)
{ 
  int info;
  int iOne=1;

  if ( NR!=NC || NR!=X->N )
   ErrExit("dimension mismatch in CholSolve");

  if ( RealComplex==LHM_REAL && StorageType==LHM_NORMAL )
   dpotrs_("U", &NR, &iOne, DM, &NR, X->DV, &NR, &info);
  else if ( RealComplex==LHM_REAL && StorageType==LHM_SYMMETRIC )
   dpptrs_("U", &NR, &iOne, DM, X->DV, &NR, &info);
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_NORMAL )
   zpotrs_("U", &NR, &iOne, ZM, &NR, X->ZV, &NR, &info);
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_HERMITIAN )
   zpptrs_("U", &NR, &iOne, ZM, X->ZV, &NR, &info);
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_SYMMETRIC )
   { fprintf(stderr,"\n*\n* ERROR: cannot cholesky-factorize a complex symmetric matrix \n*\n");
     return -1;
   };

  return info;
}

/***************************************************************/
/* solve linear systems using cholesky factorization ***********/
/***************************************************************/
int HMatrix::CholSolve(HMatrix *X, int nrhs)
{ 
  int info;

  if ( NR!=NC || NR!=X->NR )
   ErrExit("dimension mismatch in CholSolve");
  if ( nrhs > X->NC )
   ErrExit("too many RHSs requested in CholSolve");

  if ( RealComplex==LHM_REAL && StorageType==LHM_NORMAL )
   dpotrs_("U", &NR, &nrhs, DM, &NR, X->DM, &NR, &info);
  else if ( RealComplex==LHM_REAL && StorageType==LHM_SYMMETRIC )
   dpptrs_("U", &NR, &nrhs, DM, X->DM, &NR, &info);
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_NORMAL )
   zpotrs_("U", &NR, &nrhs, ZM, &NR, X->ZM, &NR, &info);
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_HERMITIAN )
   zpptrs_("U", &NR, &nrhs, ZM, X->ZM, &NR, &info);
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_SYMMETRIC )
   { fprintf(stderr,"\n*\n* ERROR: cannot cholesky-factorize a complex symmetric matrix \n*\n");
     return -1;
   };

  return info;
}

int HMatrix::CholSolve(HMatrix *X) { return CholSolve(X,X->NC); }
 
/***************************************************************/
/* QR-factorization in which 'this' is replaced by its own     */
/* QR-factorization. If *pQ and *pR are NULL on entry, or if   */
/* they point to matrices of the wrong sizes, then new         */
/* HMatrices of the appropriate size are allocated and returned*/
/* in those slots.                                             */
/***************************************************************/
int HMatrix::QR(HMatrix **pQ, HMatrix **pR)
{
  if (pQ==0 || pR==0)
   ErrExit("HMatrix::QR called with null pointers");

  HMatrix *Q = *pQ; 
  HMatrix *R = *pR; 

  // make sure Q and R are the right sizes, namely Q is NRxNR and R is NRxNC
  if ( Q==0 || Q->NR!=NR || Q->NC!=NR )
   { if (Q!=0) Warn("invalid Q matrix on entry to HMatrix::QR; reallocating...");
     Q=new HMatrix(NR, NR, RealComplex);
   };

  if ( R==0 || R->NR!=NR || R->NC!=NC )
   { if (R!=0) Warn("invalid Q matrix on entry to HMatrix::QR; reallocating...");
     R=new HMatrix(NR, NC, RealComplex);
   };

  /*--------------------------------------------------------------*/
  /*- set Q equal to the unit matrix -----------------------------*/
  /*--------------------------------------------------------------*/
  Q->Zero();
  for(int nr=0; nr<NR; nr++)
   Q->SetEntry(nr, nr, 1.0);

  /*--------------------------------------------------------------*/
  /*- call lapack's QR routine to get the R matrix, plus the     -*/
  /*- sequence of reflectors that define Q, in the body of the   -*/
  /*- R HMatrix                                                  -*/
  /*--------------------------------------------------------------*/
  R->Copy(this);
  int info, MinusOne=-1;
  int TauLength = NR<NC ? NR:NC;
  if ( RealComplex==LHM_REAL )
   { 
     double *Tau=R->DM, *dwork;
     double dlworkOptimal;
     int lworkOptimal;

     // workspace size query
     dgeqrf_(&NR, &NC, R->DM, &NR, Tau, &dlworkOptimal, &MinusOne, &info);
     lworkOptimal = (int) dlworkOptimal;

     // make sure our internally-stored workspace is large enough
     // to store TauLength + lworkOptimal doubles
     if ( lwork < (TauLength+lworkOptimal) )
      { if (work) free(work);
        work = (double *)mallocEC( (TauLength+lworkOptimal)*sizeof(double) );
        lwork=(TauLength+lworkOptimal);
      };
     Tau   = (double *)work;
     dwork = Tau + TauLength;

     // do the QR factorization 
     dgeqrf_(&NR, &NC, R->DM, &NR, Tau, dwork, &lworkOptimal, &info);

     // reconstruct the Q matrix from the Householder reflectors
     dormqr_("L", "N", &NR, &NR, &TauLength, R->DM, &NR, Tau, Q->DM, &NR,
             dwork, &lworkOptimal, &info);

   }
  else //( RealComplex==LHM_COMPLEX
   { 
     cdouble *Tau=R->ZM, *zwork;
     cdouble zlworkOptimal;
     int lworkOptimal;

     // workspace size query
     zgeqrf_(&NR, &NC, R->ZM, &NR, Tau, &zlworkOptimal, &MinusOne, &info);
     lworkOptimal = (int) real(zlworkOptimal);

     // make sure our internally-stored workspace is large enough
     // to store TauLength + lworkOptimal cdoubles
     if ( lwork < (TauLength+lworkOptimal) )
      { if (work) free(work);
        work = (cdouble *)mallocEC( (TauLength+lworkOptimal)*sizeof(cdouble) );
        lwork=(TauLength+lworkOptimal);
      };

     Tau   = (cdouble *)work;
     zwork = Tau + TauLength;

     // do the QR factorization 
     zgeqrf_(&NR, &NC, R->ZM, &NR, Tau, zwork, &lworkOptimal, &info);

     // reconstruct the Q matrix from the Householder reflectors
     zunmqr_("L", "N", &NR, &NR, &TauLength, R->ZM, &NR, Tau, Q->ZM, &NR,
             zwork, &lworkOptimal, &info);

   };

  // zero out the lower diagonal of R 
  for(int nc=0; nc<NC; nc++)
   for(int nr=nc+1; nr<NR; nr++)
    R->SetEntry(nr,nc,0.0);

  *pQ = Q;
  *pR = R;
  return info;
}

/***************************************************************/
/* routine for eigenvalues and optionally vectors of           */
/* symmetric/hermitian matrices (dsyevr / zheevr)              */
/***************************************************************/
HVector *HMatrix::Eig(HVector *Lambda, HMatrix *U)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (NR!=NC)
   ErrExit("non-square matrix passed to Eig");
  if ( Lambda && ( (Lambda->N!=NR) || (Lambda->RealComplex!=LHM_REAL)) )
   { Warn("Incorrect Lambda vector passed to Eig (reallocating)");
     delete Lambda;
     Lambda=0;
   };
  if (Lambda==0)
   Lambda=new HVector(NR,LHM_REAL);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if ( U && ( (U->NR!=NR) || (U->NC!=NC)) )
   { Warn("Incorrect U matrix passed to Eig (reallocating)");
     delete U;
     U = new HMatrix(NR, NC, RealComplex);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  const char *jobz  = (U==0) ? "N" : "V";
  int idum1, idum2, idum3, lworkOptimal, liworkOptimal, info, MinusOne=-1;
  double ddum1, ddum2, AbsTol=0.0, dlworkOptimal;

  int *isuppz = new int[2*NR]; 
 
  if ( RealComplex==LHM_REAL && StorageType==LHM_NORMAL )
   { 
     /*--------------------------------------------------------------*/
     /*- query optimal workspace sizes and (re)allocate internally-  */
     /*- stored workspaces as necessary                              */
     /*--------------------------------------------------------------*/
     dsyevr_(jobz, "A", "U", &NR, DM, &NR, &ddum1, &ddum2, &idum1, &idum2, &AbsTol, &idum3,
             Lambda->DV, U ? U->DM : 0, &NR, isuppz, &dlworkOptimal, &MinusOne,
             &liworkOptimal, &MinusOne, &info);
     lworkOptimal = (int)dlworkOptimal;
     if (lworkOptimal > lwork)
      { work=realloc(work,lworkOptimal*sizeof(double));
        lwork=lworkOptimal;
      };
     if (liworkOptimal > liwork)
      { iwork=(int *)realloc(iwork,liworkOptimal*sizeof(int));
        liwork=liworkOptimal;
      };

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     dsyevr_(jobz, "A", "U", &NR, DM, &NR, &ddum1, &ddum2, &idum1, &idum2, &AbsTol, &idum3,
             Lambda->DV, U ? U->DM : 0, &NR, isuppz, (double *)work, &lwork,
             iwork, &liwork, &info);
   }
  else if ( RealComplex==LHM_REAL && StorageType==LHM_SYMMETRIC )
   {
      ErrExit("Eig() not yet implemented for packed matrices");
   }
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_NORMAL )
   {
     /*--------------------------------------------------------------*/
     /*- query optimal workspace sizes and (re)allocate internally-  */
     /*- stored workspaces as necessary                              */
     /*--------------------------------------------------------------*/
     cdouble zlworkOptimal;
     double dlrworkOptimal;
     zheevr_(jobz, "A", "U", &NR, ZM, &NR, &ddum1, &ddum2, &idum1, &idum2,
             &AbsTol, &idum3, Lambda->DV, U ? U->ZM : 0, &NR, isuppz,
             &zlworkOptimal, &MinusOne, &dlrworkOptimal, &MinusOne,
             &liworkOptimal, &MinusOne, &info);

     int workSize  = (int)real(zlworkOptimal);
     int rworkSize = (int)dlrworkOptimal;
     lworkOptimal = 2*workSize + rworkSize;
     if (lworkOptimal > lwork)
      { work=realloc(work,lworkOptimal*sizeof(double));
        lwork=lworkOptimal;
      };

     cdouble *zwork = (cdouble *)work;
     double *rwork = (double *)(zwork + workSize);

     if (liworkOptimal > liwork)
      { iwork=(int *)realloc(iwork,liworkOptimal*sizeof(int));
        liwork=liworkOptimal;
      };

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     zheevr_(jobz, "A", "U", &NR, ZM, &NR, &ddum1, &ddum2, &idum1, &idum2, &AbsTol, &idum3,
             Lambda->DV, U ? U->ZM : 0, &NR, isuppz, (cdouble *)work, &workSize,
             (double *)rwork, &rworkSize, iwork, &liwork, &info);

   }
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_HERMITIAN )
   {
      ErrExit("Eig() not yet implemented for packed matrices");
   }
  else if ( RealComplex==LHM_COMPLEX && StorageType==LHM_SYMMETRIC )
   {
      ErrExit("use NSEig() for non-hermitian eigenproblems");
   }

  delete[] isuppz;
  return Lambda;
}

/***************************************************************/
/* routine for eigenvalues and optionally vectors of           */
/* non-symmetric matrices (dsyevr / zheevr)                    */
/***************************************************************/
HVector *HMatrix::NSEig(HVector *Lambda, HMatrix *U)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (NR!=NC)
   ErrExit("non-square matrix passed to NSEig");
  if (StorageType!=LHM_NORMAL)
   ErrExit("packed-storage matrix passed to NSEig()");
  if ( Lambda && ( (Lambda->N!=NR) || (Lambda->RealComplex!=RealComplex)) )
   { Warn("Incorrect Lambda vector passed to NSEig (reallocating)");
     delete Lambda;
     Lambda=0;
   }
  if (Lambda==0)
   Lambda=new HVector(NR,LHM_COMPLEX);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if ( U && ( (U->NR!=NR) || (U->NC!=NC)) )
   { Warn("Incorrect U matrix passed to NSEig (reallocating)");
     delete U;
     U = new HMatrix(NR, NC, LHM_COMPLEX);
   }

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  const char *jobvl = "N";
  const char *jobvr = ((U==0) ? "N" : "V");
  int lworkOptimal, info, MinusOne=-1;

  if (RealComplex==LHM_REAL)
   { 
     /*--------------------------------------------------------------*/
     /*- query optimal workspace sizes and (re)allocate internally-  */
     /*- stored workspaces as necessary                              */
     /*--------------------------------------------------------------*/
     double *wr = new double[NR];
     double *wi = new double[NR];
     double dlworkOptimal;
     dgeev_(jobvl, jobvr, &NR, DM, &NR, wr, wi, 0, &NR, 
            U ? U->DM : 0, &NR, &dlworkOptimal, &MinusOne, &info);

     lworkOptimal = (int)dlworkOptimal;
     if (lworkOptimal > lwork)
      { work=realloc(work,lworkOptimal*sizeof(double));
        lwork=lworkOptimal;
      };

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     dgeev_(jobvl, jobvr, &NR, DM, &NR, wr, wi, 0, &NR, 
            U ? U->DM : 0, &NR, (double *)work, &lwork, &info);

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     for(int n=0; n<NR; n++)
      Lambda->SetEntry( n, cdouble(wr[n], wi[n]));

     delete[] wr;
     delete[] wi;
   }
  else // (RealComplex==LHM_COMPLEX)
   {
     /*--------------------------------------------------------------*/
     /*- query optimal workspace sizes and allocate workspaces       */
     /*--------------------------------------------------------------*/
     double *rwork = (double *)mallocEC(2*NR*sizeof(double));

     cdouble zlworkOptimal;
     zgeev_(jobvl, jobvr, &NR, ZM, &NR, Lambda->ZV, 0, &NR, 
            (U ? U->ZM : 0), &NR, 
            &zlworkOptimal, &MinusOne, rwork, &info);
     lworkOptimal = 2*(int)(ceil(real(zlworkOptimal)));
     cdouble *zwork = (cdouble *)mallocEC(lworkOptimal * sizeof(cdouble));

     zgeev_(jobvl, jobvr, &NR, ZM, &NR, Lambda->ZV, 0, &NR,
            U ? U->ZM : 0, &NR, zwork, &lworkOptimal, rwork, &info);

     Log("zgeev(%i): info=%i",NR,NR,info);
     free(zwork);
     free(rwork);
   }

  return Lambda;

}

/***************************************************************/
/* routine for singular-value decomposition (dgesvd / zgesvd)  */
/*                                                             */
/* If Sigma is NULL on entry, a new vector is allocated and    */
/* returned.                                                   */
/*                                                             */
/* If U is:                                                    */
/*  -- NULL:     no right singular vectors computed            */
/*  -- NR x NR:  U <-- matrix of full RSVs (JOBU='A')          */
/*  -- NR x NS:  U <-- matrix of reduced RSVs (JOBU='S')       */
/*                                                             */
/* If VT is:                                                   */
/*  -- NULL:     no left singular vectors computed             */
/*  -- NC x NC:  VT <-- (adjoint of) matrix of full LSVs (JOBV='A')  */
/*  -- NS x NC:  VT <-- (adjoint of) matrix of reduced LSVs (JOBV='S') */
/***************************************************************/
HVector *HMatrix::SVD(HVector *Sigma, HMatrix *U, HMatrix *VT)
{
  if (StorageType!=LHM_NORMAL)
   ErrExit("SVD() not supported for packed-storage matrices");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NRCMin = (NR<NC ? NR : NC);
  if ( Sigma && ( (Sigma->N!=NRCMin) ) )
   { Warn("Incorrect Sigma vector passed to SVD (reallocating)");
     delete Sigma;
     Sigma=0;
   };
  if (Sigma==0)
   Sigma=new HVector(NRCMin,LHM_REAL);
  
  char jobu[2]="";
  if ( U==0 ) 
   jobu[0]='N';
  else if ( U->NR==NR && U->NC==NR     && U->RealComplex==RealComplex) 
   jobu[0]='A';
  else if ( U->NR==NR && U->NC==NRCMin && U->RealComplex==RealComplex)
   jobu[0]='S';
  else 
   { Warn("Incorrect U matrix passed to SVD (omitting RSV calculation)");
     jobu[0]='N';
     U=0;
   };
  int LDU = U ? U->NR : 0;

  char jobvt[2]="";
  if ( VT==0 ) 
   jobvt[0]='N';
  else if ( VT->NR==NC && VT->NC==NC     && VT->RealComplex==RealComplex)
   jobvt[0]='A';
  else if ( VT->NR==NRCMin && VT->NC==NC && VT->RealComplex==RealComplex)
   jobvt[0]='S';
  else 
   { Warn("Incorrect VT matrix passed to SVD (omitting LSV calculation)");
     jobvt[0]='N';
     VT=0;
   };
  int LDVT = VT ? VT->NR : 0;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int lworkOptimal, info, MinusOne=-1;
  if (RealComplex==LHM_REAL)
   { 
     /*--------------------------------------------------------------*/
     /*- query optimal workspace sizes and (re)allocate internally-  */
     /*- stored workspaces as necessary                              */
     /*--------------------------------------------------------------*/
     double dlworkOptimal;
     dgesvd_(jobu, jobvt, &NR, &NC, DM, &NR, Sigma->DV,
             U ? U->DM : 0, &LDU, VT ? VT->DM : 0, &LDVT,
             &dlworkOptimal, &MinusOne, &info);

     lworkOptimal = (int)dlworkOptimal;
     if (lworkOptimal > lwork)
      { work=realloc(work,lworkOptimal*sizeof(double));
        lwork=lworkOptimal;
      };

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     dgesvd_(jobu, jobvt, &NR, &NC, DM, &NR, Sigma->DV,
             U ? U->DM : 0, &LDU, VT ? VT->DM : 0, &LDVT,
             (double *)work, &lwork, &info);
   }
  else // (RealComplex==LHM_COMPLEX)
   {
     /*--------------------------------------------------------------*/
     /*- query optimal workspace sizes and (re)allocate internally-  */
     /*- stored workspaces as necessary                              */
     /*--------------------------------------------------------------*/
     double *rwork = new double[5*NRCMin];
     cdouble zlworkOptimal;
     zgesvd_(jobu, jobvt, &NR, &NC, ZM, &NR, Sigma->DV, 
             U ? U->ZM : 0, &LDU, VT ? VT->ZM : 0, &LDVT, 
             &zlworkOptimal, &MinusOne, rwork, &info);

     lworkOptimal = 2*(int)(real(zlworkOptimal));
     if (lworkOptimal > lwork)
      { work=realloc(work,lworkOptimal*sizeof(double));
        lwork=lworkOptimal;
      };

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     zgesvd_(jobu, jobvt, &NR, &NC, ZM, &NR, Sigma->DV, 
             U ? U->ZM : 0, &LDU, VT ? VT->ZM : 0, &LDVT,
             (cdouble *)work, &lwork, rwork, &info);
     delete[] rwork;

   };

  return Sigma;

}

/***************************************************************/
/***************************************************************/ 
/***************************************************************/
double HMatrix::GetNorm(bool UseInfinityNorm)
{
  char *Norm = const_cast <char *>(UseInfinityNorm ? "I" : "1");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (RealComplex==LHM_REAL)
   { 
     /*--------------------------------------------------------------*/
     /*- (re)allocate internally--stored workspaces as necessary     */
     /*--------------------------------------------------------------*/
     int lworkOptimal = UseInfinityNorm ? NR : 0;
     if (lworkOptimal > lwork)
      { work=realloc(work,lworkOptimal*sizeof(double));
        lwork=lworkOptimal;
      };

     double *rwork=(double *)work;
     return dlange_(Norm, &NR, &NC, DM, &NR, rwork);
   }
  else // (RealComplex==LHM_COMPLEX)
   {
     /*--------------------------------------------------------------*/
     /*- (re)allocate internally--stored workspaces as necessary     */
     /*--------------------------------------------------------------*/
     int lworkOptimal = UseInfinityNorm ? NR : 0;
     if (lworkOptimal > lwork)
      { work=realloc(work,lworkOptimal*sizeof(double));
        lwork=lworkOptimal;
      };
     double *rwork  = (double *)work;
     return zlange_(Norm, &NR, &NC, ZM, &NR, rwork);
   };

  return 0.0;

}

/***************************************************************/
/***************************************************************/ 
/***************************************************************/
double HMatrix::GetRCond(double ANorm, bool UseInfinityNorm)
{
  char *Norm = const_cast<char *> (UseInfinityNorm ? "I" : "1");
  double RCond;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (RealComplex==LHM_REAL)
   { 
     /*--------------------------------------------------------------*/
     /*- (re)allocate internally--stored workspaces as necessary     */
     /*--------------------------------------------------------------*/
     int lworkOptimal = 4*NR;
     if (lworkOptimal > lwork)
      { work=realloc(work,lworkOptimal*sizeof(double));
        lwork=lworkOptimal;
      };
     int iworkOptimal = NR;
     if (iworkOptimal > liwork)
      { iwork=(int *)realloc(iwork,iworkOptimal*sizeof(int));
        liwork=iworkOptimal;
      };

     double *rwork=(double *)work;
     int info;
     dgecon_(Norm, &NR, DM, &NR, &ANorm, &RCond, rwork, (int *)iwork, &info);
   }
  else // (RealComplex==LHM_COMPLEX)
   {
     /*--------------------------------------------------------------*/
     /*- (re)allocate internally--stored workspaces as necessary     */
     /*--------------------------------------------------------------*/
     int lworkOptimal = 6*NR; // 2N doubles plus 2N cdoubles
     if (lworkOptimal > lwork)
      { work=realloc(work,lworkOptimal*sizeof(double));
        lwork=lworkOptimal;
      };
  
     double *rwork  = (double *)work;
     double *work2  = rwork + 2*NR;
     cdouble *zwork = (cdouble *)work2;
     int info;
     zgecon_(Norm, &NR, ZM, &NR, &ANorm, &RCond, zwork, rwork, &info);

   };

  return RCond;

}

/***************************************************************/
/* dot product of vectors **************************************/
/***************************************************************/
cdouble DoDot(HVector *A, HVector *B, bool Hermitian=true)
{
  if ( (A->N != A->N) || (A->RealComplex != B->RealComplex) )
   { Warn("attempt to dot-product incompatible vectors (returning 0)");
     return 0.0;
   };

  int N=A->N, IncX=1, IncY=1;
  if (A->RealComplex==LHM_REAL)
   return ddot_(&N, A->DV, &IncX, B->DV, &IncY);
  else if (Hermitian)
   return zdotc_(&N, A->ZV, &IncX, B->ZV, &IncY);
  else 
   return zdotu_(&N, A->ZV, &IncX, B->ZV, &IncY);
}

cdouble HVector::Dot(HVector *B)
 { return DoDot(this, B, true); }

cdouble HVector::DotU(HVector *B)
 { return DoDot(this, B, false); }

double HVector::DotD(HVector *B)
 { return real(DoDot(this, B, true)); }

/***************************************************************/
/* matrix-vector multiplication, A*X = Y (where this=A)        */
/*                                                             */
/*  Transpose=='C'   // use the hermitian transpose            */
/*  Transpose=='T'   // use the non-hermitian transpose        */
/*  Transpose==anything else for no transpose                  */
/***************************************************************/
void HMatrix::Apply(HVector *X, HVector *Y, char TransChar)
{
  if ( StorageType!=LHM_NORMAL )
   ErrExit("%s:%i: HMatrix::Apply not yet implemented for packed matrices");
  if ( X==0 || Y==0  )
   ErrExit("%s:%i: invalid operand in HMatrix::Apply",__FILE__,__LINE__);
  if ( (RealComplex!=X->RealComplex) || (RealComplex!=Y->RealComplex) )
   ErrExit("%s:%i: data type mismatch in HMatrix::Apply",__FILE__,__LINE__);

  char Trans[2]="N";
  int NREffective=NR;
  int NCEffective=NC;
  if (TransChar=='T' || TransChar=='C')
   { Trans[0]=TransChar;
     NREffective=NC;
     NCEffective=NR;
   };

  if ( (X->N != NCEffective) || (Y->N != NREffective) )
   ErrExit("%s:%i: dimension mismatch in HMatrix::Apply",__FILE__,__LINE__);

  int IncX=1, IncY=1;
  double dOne=1.0, dZero=0.0;
  cdouble zOne=1.0, zZero=0.0;

  if (RealComplex==LHM_REAL)
   dgemv_(Trans, &NR, &NC, &dOne, DM, &NR, X->DV, &IncX, &dZero, Y->DV, &IncY);
  else 
   zgemv_(Trans, &NR, &NC, &zOne, ZM, &NR, X->ZV, &IncX, &zZero, Y->ZV, &IncY);

}
