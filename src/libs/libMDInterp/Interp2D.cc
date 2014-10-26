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
 *  Interp2D.cc -- implementation of Interp2D class for 2D interpolation
 *
 *  homer reid  -- 3/2011
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>
#include <libhmat.h>

#include "libMDInterp.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#ifdef USE_PTHREAD
#  include <pthread.h>
#endif

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define NEQUATIONS 16  // number of matrix rows
#define NCOEFF     16  // number of matrix columns
#define NDATA      4   // number of values stored at each grid point

/****************************************************************/
/****************************************************************/
/****************************************************************/
static int GetPhiVDTableOffset(int nf, int nFun, int n1, int N1, int n2, int N2)
{
  (void) N1; // unused;
  return NDATA*(nf + nFun*(n2 + N2*(n1)));
}

/****************************************************************/
/****************************************************************/
/****************************************************************/
static int GetCTableOffset(int nf, int nFun, int n1, int N1, int n2, int N2)
{
  (void) N1; // unused;
  return NCOEFF*(nf + nFun*(n2 + (N2-1)*(n1)));
}

/****************************************************************/
/* thread routine used to compute values and derivatives of     */
/* user's function at grid points                               */
/****************************************************************/
typedef struct ThreadData
 {
   double *X1Points, *X2Points;
   int N1, N2;
   double X1Min, X2Min;
   double DX1, DX2;

   int nFun;

   Phi2D PhiFunc;
   void *UserData;
   double *PhiVDTable;

   int nt, nThread;

 } ThreadData;

static void *GetPhiVD_Thread(void *data)
{
  
  ThreadData *TD=(ThreadData *)data;

  /***************************************************************/
  /* variables unpacked from thread data structure ***************/
  /***************************************************************/
  double *X1Points=TD->X1Points;
  double *X2Points=TD->X2Points;
  int N1   = TD->N1;
  int N2   = TD->N2;
  double X1Min=TD->X1Min;
  double X2Min=TD->X2Min;
  double DX1=TD->DX1;
  double DX2=TD->DX2;
  int nFun = TD->nFun;
  double *PhiVDTable = TD->PhiVDTable;

  /***************************************************************/
  /* other local variables ***************************************/
  /***************************************************************/
  int nt, nf;
  double *PhiVD = new double[nFun*NDATA];
  double *P;

  /*--------------------------------------------------------------*/
  /*- hack to force scheduler on gnu/linux systems to assign each-*/
  /*- thread to a different core                                 -*/
  /*--------------------------------------------------------------*/
#if defined(_GNU_SOURCE) && defined(USE_PTHREAD)
  cpu_set_t cpuset;
  CPU_ZERO(&cpuset);
  CPU_SET(TD->nt,&cpuset);
  pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
#endif
  /*--------------------------------------------------------------*/
  /*- end gnu/linux hack -----------------------------------------*/
  /*--------------------------------------------------------------*/

  /***************************************************************/
  /* loop over all grid points, calling the user's function to   */ 
  /* get function values and derivatives at each point, and      */ 
  /* inserting those values into the table                       */
  /***************************************************************/
  nt=0;
  int nPoints=0;
  int n1, n2;
  for(n1=0; n1<N1; n1++)
   for(n2=0; n2<N2; n2++)
    {
      /*--------------------------------------------------------------*/
      /*- status logging ---------------------------------------------*/
      /*--------------------------------------------------------------*/
      if (TD->nt==0)
       { nPoints++;
        // for(pc=10; pc<=90; pc+=10)
        //  if (nPoints == (pc*TotalPoints/100) )
        //   Log("...%i %%",pc);
       };
 
      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      nt++;
      if (nt==TD->nThread) nt=0;
      if (nt!=TD->nt) continue;

      if (X1Points)
       TD->PhiFunc(X1Points[n1], X2Points[n2], TD->UserData, PhiVD);
      else
       TD->PhiFunc(X1Min + n1*DX1, X2Min + n2*DX2, TD->UserData, PhiVD);

      for(nf=0; nf<nFun; nf++)
       { 
         P = PhiVDTable + GetPhiVDTableOffset(nf, nFun, n1, N1, n2, N2);
         memcpy(P, PhiVD + NDATA*nf, NDATA*sizeof(double));
       };
   
     };

  delete[] PhiVD;

  return 0;
  
}

/****************************************************************/
/* class constructor 1: construct the class from a user-supplied*/
/* function and grid                                            */
/****************************************************************/
Interp2D::Interp2D(double *pX1Points, int pN1,
                   double *pX2Points, int pN2,
                   int pnFun, Phi2D PhiFunc, void *UserData) 
{
   /*--------------------------------------------------------------*/ 
   /*- simple setup: initialize class data                ---------*/ 
   /*--------------------------------------------------------------*/ 
   if (pN1<2 || pN2<2)
    ErrExit("Interp2D: grid must have 2 or more points in every dimension");
   N1=pN1; 
   N2=pN2;
   X1Points=(double *)memdup(pX1Points, N1*sizeof(double));
   X2Points=(double *)memdup(pX2Points, N2*sizeof(double));
   nFun=pnFun;

   /*--------------------------------------------------------------*/ 
   /*- allocate space for the CTable, which stores all             */ 
   /*- coefficients for all functions and all grid cells.          */
   /*-                                                             */
   /*- how it works:                                               */
   /*-                                                             */
   /*- CTable + GetCTableOffset(nf,nFun,n1,N1,n2,N2)               */
   /*-                                                             */
   /*- is a pointer to an array of NCOEFF doubles which are the    */
   /*- coefficients for function #nf within grid cell (n1,n2).     */
   /*-                                                             */
   /*- grid cell (n1,n2) is defined to be the region               */
   /*-  X1Points[n1] <= X1 < X1Points[n1+1]                        */
   /*-  X2Points[n2] <= X2 < X2Points[n2+1]                        */
   /*--------------------------------------------------------------*/ 
   CTable=(double *)mallocEC((N1-1)*(N2-1)*nFun*NCOEFF*sizeof(double));
   if (!CTable)
    ErrExit("%s:%i:out of memory",__FILE__,__LINE__);

   /*--------------------------------------------------------------*/ 
   /*--------------------------------------------------------------*/ 
   /*--------------------------------------------------------------*/ 
   if (PhiFunc==0)
    return;

   /*--------------------------------------------------------------*/ 
   /*--------------------------------------------------------------*/ 
   /*--------------------------------------------------------------*/ 
   InitInterp2D(PhiFunc, UserData);

}

/****************************************************************/
/* class constructor 2: construct the class from user-supplied  */
/* uniform grids                                                */
/****************************************************************/
Interp2D::Interp2D(double pX1Min, double X1Max, int pN1,
                   double pX2Min, double X2Max, int pN2,
                   int pnFun, Phi2D PhiFunc, void *UserData)
{
   /*--------------------------------------------------------------*/ 
   /*- simple setup: initialize class data                ---------*/ 
   /*--------------------------------------------------------------*/ 
   if (pN1<2 || pN2<2)
    ErrExit("Interp2D: grid must have 2 or more points in every dimension");

   N1=pN1; 
   X1Points=0;
   X1Min=pX1Min;
   DX1 = (X1Max - X1Min) / ((double)(N1-1));

   N2=pN2; 
   X2Points=0;
   X2Min=pX2Min;
   DX2 = (X2Max - X2Min) / ((double)(N2-1));

   if ( DX1==0.0 || DX2==0.0 )
    ErrExit("Interp2D: grid has zero size in one or more dimensions");

   nFun=pnFun;

   /*--------------------------------------------------------------*/ 
   /*- allocate space for the CTable (see comment above on how it  */ 
   /*- works)                                                      */
   /*--------------------------------------------------------------*/ 
   CTable=(double *)mallocEC((N1-1)*(N2-1)*nFun*NCOEFF*sizeof(double));
   if (!CTable)
    ErrExit("%s:%i:out of memory",__FILE__,__LINE__);

   /*--------------------------------------------------------------*/ 
   /*--------------------------------------------------------------*/ 
   /*--------------------------------------------------------------*/ 
   if (PhiFunc==0)
    return;

   /*--------------------------------------------------------------*/ 
   /*--------------------------------------------------------------*/ 
   /*--------------------------------------------------------------*/ 
   InitInterp2D(PhiFunc, UserData);

}

/****************************************************************/
/* main body of class constructor for the above two entry points*/
/****************************************************************/
void Interp2D::InitInterp2D(Phi2D PhiFunc, void *UserData)
{

 //  Log("Creating interpolation grid with %i grid points...",N1*N2);

   /*--------------------------------------------------------------*/ 
   /*- temporarily allocate a big array in which we will store     */ 
   /*- values and derivatives of the user's function. allocating   */ 
   /*- this entire table all at once is totally storage-inefficient*/ 
   /*- and could be much optimized since we only ever need 4       */ 
   /*- grid points worth of data at a time, not the full grid.     */ 
   /*-                                                             */
   /*- ('PVD' = 'phi values and derivatives')                      */
   /*-                                                             */
   /*- how it works:                                               */
   /*-                                                             */
   /*-  PVDTable + GetPVDTableOffset(nf,nFun,n1,N1,n2,N2)          */
   /*-                                                             */
   /*- is a pointer to an array of NDATA doubles which are the     */
   /*- value and derivatives of function #nf at grid point (n1,n2).*/
   /*--------------------------------------------------------------*/
   double *PhiVDTable=(double *)mallocEC(N1*N2*nFun*NDATA*sizeof(double));
   if (!PhiVDTable)
    ErrExit("%s:%i:out of memory",__FILE__,__LINE__);

   /*--------------------------------------------------------------*/
   /*- fire off threads that will call the user's function to     -*/
   /*- populate the PVDTable table.                               -*/
   /*--------------------------------------------------------------*/
 //  Log("Computing user function at grid points...");
   int nThread=GetNumThreads();

#ifdef USE_PTHREAD
   ThreadData *TDs = new ThreadData[nThread], *TD;
   pthread_t *Threads = new pthread_t[nThread];
#endif
   ThreadData TD1;
   int nt;

   TD1.X1Points=X1Points;
   TD1.X2Points=X2Points;
   TD1.N1=N1;
   TD1.N2=N2;
   TD1.X1Min=X1Min;
   TD1.X2Min=X2Min;
   TD1.DX1=DX1;
   TD1.DX2=DX2;
   TD1.nFun=nFun;
   TD1.PhiFunc=PhiFunc;
   TD1.UserData=UserData;
   TD1.PhiVDTable=PhiVDTable;
   TD1.nThread=nThread;
#ifdef USE_OPENMP
#pragma omp parallel for firstprivate(TD1), schedule(static,1), num_threads(nThread)
#endif
   for(nt=0; nt<nThread; nt++)
    {
#ifdef USE_PTHREAD
      TD=&(TDs[nt]);
      *TD = TD1;
#else
      ThreadData *TD=&TD1;
#endif
      TD->nt=nt;
#ifdef USE_PTHREAD
      if (nt+1 == nThread)
	GetPhiVD_Thread((void *)TD);
      else
	pthread_create( &(Threads[nt]), 0, GetPhiVD_Thread, (void *)TD);
#else
      GetPhiVD_Thread((void *)TD);
#endif
    };

#ifdef USE_PTHREAD
   /*--------------------------------------------------------------*/
   /*- wait for threads to terminate ------------------------------*/
   /*--------------------------------------------------------------*/ 
   for(nt=0; nt<nThread-1; nt++)
    pthread_join(Threads[nt],0);

   delete[] Threads;
   delete[] TDs;
#endif

   /*--------------------------------------------------------------*/
   /*- construct the matrix whose inverse maps a vector of        -*/
   /*- phi values and derivatives into polynomial coefficients.   -*/
   /*- (since this matrix is a fixed constant entity i could just -*/
   /*-  assemble and invert it once and hard-code the inverse     -*/
   /*-  here, but i find it just as easy to do it this way.)      -*/
   /*-                                                            -*/
   /*- the outer loop (actually two loops) loops over the 4       -*/
   /*- corners of the grid cell.                                  -*/
   /*-                                                            -*/
   /*- the inner loop (actually two loops) loops over the NCOEFF  -*/
   /*- monomials in the interpolating polynomial.                 -*/
   /*--------------------------------------------------------------*/
   double X1Bar, X2Bar;
   double X1Powers[5], X2Powers[5];
   int nm, ncp, p, q;
   HMatrix *M=new HMatrix(NEQUATIONS, NCOEFF);

   M->Zero();

   /* note: X1Powers[p] = 0        ,  p==0 */
   /*                   = X1^{p-1} ,  p>0  */
   X1Powers[0]=0.0; X1Powers[1]=1.0; 
   X2Powers[0]=0.0; X2Powers[1]=1.0; 

   for(ncp=0, X1Bar=0.0; X1Bar<=1.0; X1Bar+=1.0)
    for(X2Bar=0.0; X2Bar<=1.0; X2Bar+=1.0, ncp++)
     { 
       X1Powers[2]=X1Powers[3]=X1Powers[4]=X1Bar;
       X2Powers[2]=X2Powers[3]=X2Powers[4]=X2Bar;

       for(nm=0, p=0; p<4; p++)
        for(q=0; q<4; q++, nm++)
         { 
           /* the nmth monomial is x1^p x2^q                 */
           /* its x1 derivative is p x1^{p-1} x2^q ,     etc */
           /* since xi={0,1}, we have xi^n = xi for all n    */
           M->SetEntry(NDATA*ncp + 0, nm, X1Powers[p+1]*X2Powers[q+1]);
           M->SetEntry(NDATA*ncp + 1, nm, p*X1Powers[p]*X2Powers[q+1]);
           M->SetEntry(NDATA*ncp + 2, nm, q*X1Powers[p+1]*X2Powers[q]);
           M->SetEntry(NDATA*ncp + 3, nm, p*q*X1Powers[p]*X2Powers[q]);
         };
     }; 
   M->LUFactorize();

   /*--------------------------------------------------------------*/
   /*- now go through and solve a linear system for each grid      */
   /*- cell to compute the coefficients of the interpolating       */
   /*- polynomial in that cell                                     */
   /*--------------------------------------------------------------*/
   // int nPoints=0, TotalPoints=(N1-1)*(N2-1), pc;
   int n1, dn1, n2, dn2;
   int nf;
   double *P;
   HVector *C = new HVector(NCOEFF);

   //Log("Computing coefficients of interpolating polynomials...");
   for(n1=0; n1<(N1-1); n1++)
    for(n2=0; n2<(N2-1); n2++)
     { 
       if (X1Points)
        { DX1=X1Points[n1+1]-X1Points[n1];
          DX2=X2Points[n2+1]-X2Points[n2];
        };

       /* separately compute interpolation coefficients for each function  */
       for(nf=0; nf<nFun; nf++)
        { 
          /* construct the RHS vector by extracting from PhiVDTable the 4 data */
          /* values for each of the 4 corners of this grid cell.               */
          for(ncp=0, dn1=0; dn1<=1; dn1++)
           for(dn2=0; dn2<=1; dn2++, ncp++)
            { P=PhiVDTable + GetPhiVDTableOffset(nf, nFun, n1+dn1, N1, n2+dn2, N2);
              C->SetEntry( NDATA*ncp + 0, P[0]);
              C->SetEntry( NDATA*ncp + 1, DX1*P[1]);
              C->SetEntry( NDATA*ncp + 2, DX2*P[2]);
              C->SetEntry( NDATA*ncp + 3, DX1*DX2*P[3]);
            };

          /* solve the 16x16 system */
          M->LUSolve(C);

          /* store the coefficients in the appropriate place in CTable */
          P=CTable + GetCTableOffset(nf, nFun, n1, N1, n2, N2);
          memcpy(P, C->DV, NCOEFF*sizeof(double));
        };

     };

   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   delete C;
   delete M;
   free(PhiVDTable);
   //Log("...done!");

}  

/****************************************************************/
/****************************************************************/
/****************************************************************/
void Interp2D::ReInitialize(Phi2D PhiFunc, void *UserData)
{ 
  InitInterp2D(PhiFunc, UserData);
}

/****************************************************************/
/* class constructor 3: construct the class from a data file    */
/* generated by a previous call to Interp2D::WriteToFile()      */
/****************************************************************/
Interp2D::Interp2D(const char *FileName)
{
   Log("Attempting to read interpolation table from file %s...",FileName);

   FILE *f=fopen(FileName,"r");
   if (!f)
    ErrExit("could not open file %s");
   
   freadEC(&N1  , sizeof(int), 1, f, FileName);
   freadEC(&N2  , sizeof(int), 1, f, FileName);
   freadEC(&nFun, sizeof(int), 1, f, FileName);

   int NonUniform;
   freadEC(&NonUniform, sizeof(int), 1, f, FileName);

   if (NonUniform)
    { X1Points=(double *)mallocEC(N1*sizeof(double));
      freadEC(X1Points, sizeof(double), N1, f, FileName);
      X2Points=(double *)mallocEC(N2*sizeof(double));
      freadEC(X2Points, sizeof(double), N2, f, FileName);
    }
   else
    { X1Points=X2Points=0;
      freadEC(&X1Min, sizeof(double), 1, f, FileName);
      freadEC(&DX1,   sizeof(double), 1, f, FileName);
      freadEC(&X2Min, sizeof(double), 1, f, FileName);
      freadEC(&DX2,   sizeof(double), 1, f, FileName);
    };

   int CTableSize=(N1-1)*(N2-1)*nFun*NCOEFF;
   int fCTableSize;
   freadEC(&fCTableSize, sizeof(int), 1, f, FileName);
   if (fCTableSize!=CTableSize)
    ErrExit("%s: data file is invalid",FileName);

   CTable=(double *)mallocEC(CTableSize*sizeof(double));
   freadEC(CTable, sizeof(double), CTableSize, f, FileName);

   fclose(f);

   Log("...success!");

}

/****************************************************************/
/* write all class data to a binary file for subsequent recovery*/
/* by the above constructor routine                             */
/****************************************************************/
void Interp2D::WriteToFile(const char *FileName)
{
   Log("Writing interpolation table to file %s...",FileName);

   FILE *f=fopen(FileName,"w");
   if (!f)
    ErrExit("could not open file %s");

   fwrite(&N1  , sizeof(int), 1, f);
   fwrite(&N2  , sizeof(int), 1, f);
   fwrite(&nFun, sizeof(int), 1, f);

   int NonUniform=(X1Points ? 1 : 0);
   fwrite(&NonUniform, sizeof(int), 1, f);

   if (NonUniform)
    { fwrite(X1Points, sizeof(double), N1, f);
      fwrite(X2Points, sizeof(double), N2, f);
    }
   else 
    { fwrite(&X1Min, sizeof(double), 1, f);
      fwrite(&DX1,   sizeof(double), 1, f);
      fwrite(&X2Min, sizeof(double), 1, f);
      fwrite(&DX2,   sizeof(double), 1, f);
    };

   int CTableSize=(N1-1)*(N2-1)*nFun*NCOEFF;
   fwrite(&CTableSize, sizeof(int), 1, f);
   fwrite(CTable, sizeof(double), CTableSize, f);

   fclose(f);
   Log("...success!");
}

/****************************************************************/
/* class destructor *********************************************/
/****************************************************************/
Interp2D::~Interp2D()
{ 
  if (X1Points) free(X1Points);
  if (X2Points) free(X2Points);
  free(CTable);
}

/****************************************************************/
/* the routine that actually does the interpolation.            */
/* TODO: Consolidate Evaluate, EvaluatePlus, EvaluatePlusPlus   */
/* into a single routine.                                       */
/****************************************************************/
void Interp2D::Evaluate(double X1, double X2, double *Phi)
{
  int n1, n2;
  double X1Bar, X2Bar;

  FindInterval(X1, X1Points, N1, X1Min, DX1, &n1, &X1Bar);
  FindInterval(X2, X2Points, N2, X2Min, DX2, &n2, &X2Bar);

  /****************************************************************/
  /* tabulate powers of the scaled/shifted coordinates            */
  /****************************************************************/
  int np;
  double X1Powers[4], X2Powers[4];
  X1Powers[0]=X2Powers[0]=1.0;
  for(np=1; np<4; np++)
   { X1Powers[np]=X1Powers[np-1]*X1Bar;
     X2Powers[np]=X2Powers[np-1]*X2Bar;
   };

  /****************************************************************/
  /* construct the NCOEFF-element vector of monomial values that  */
  /* I dot with the NCOEFF-element vector of coefficients to get  */
  /* the value of the interpolating polynomial                    */
  /****************************************************************/
  double MV[NCOEFF];
  int nm, p, q;
  for(nm=p=0; p<4; p++)
   for(q=0; q<4; q++, nm++)
    MV[nm]=X1Powers[p]*X2Powers[q];

  /****************************************************************/
  /* for each component of the Phi function vector, look up the   */
  /* coefficients of the interpolating polynomial for that        */
  /* component and use them to get the value of the interpolant   */
  /****************************************************************/
  double P;
  double *C;
  int nf;
  for(nf=0; nf<nFun; nf++)
   { 
     C=CTable + GetCTableOffset(nf, nFun, n1, N1, n2, N2);

     for(P=0.0, nm=0; nm<NCOEFF; nm++)
      P+=C[nm]*MV[nm];

     Phi[nf]=P;

   };

}

/****************************************************************/
/* Like Evaluate but also computes derivatives.                 */
/* Phi[ 4*nf + 0 ] = Phi_{nf}                                   */
/* Phi[ 4*nf + 1 ] = d/dx Phi_{nf}                              */
/* Phi[ 4*nf + 2 ] = d/dy Phi_{nf}                              */
/* Phi[ 4*nf + 3 ] = d^2/dxdy Phi_{nf}                          */
/****************************************************************/
void Interp2D::EvaluatePlus(double X1, double X2, double *Phi)
{
  int n1, n2;
  double X1Bar, X2Bar;

  /****************************************************************/
  /* construct the scaled and shifted versions of the input       */
  /* coordinates.                                                 */
  /* the scaling and shifting operation is constructed so that    */
  /* the corners of the grid cell are at {(0,0),(0,1),(1,0),(1,1)}*/
  /****************************************************************/
  FindInterval(X1, X1Points, N1, X1Min, DX1, &n1, &X1Bar);
  FindInterval(X2, X2Points, N2, X2Min, DX2, &n2, &X2Bar);

  /****************************************************************/
  /* tabulate powers of the scaled/shifted coordinates            */
  /****************************************************************/
  double X1Powers[4], X2Powers[4];
  X1Powers[0]=X2Powers[0]=1.0;
  for(int np=1; np<4; np++)
   { X1Powers[np]=X1Powers[np-1]*X1Bar;
     X2Powers[np]=X2Powers[np-1]*X2Bar;
   };

  /****************************************************************/
  /* get the widths of this grid cell *****************************/
  /****************************************************************/
  double L1, L2;
  if (X1Points)
   { L1=X1Points[n1+1]-X1Points[n1];
     L2=X2Points[n2+1]-X2Points[n2];
   }
  else 
   { L1 = DX1;
     L2 = DX2;
   };

  /****************************************************************/
  /* construct the NCOEFF-element vector of monomial values that  */
  /* I dot with the NCOEFF-element vector of coefficients to get  */
  /* the value of the interpolating polynomial                    */
  /****************************************************************/
  double MVP[NCOEFF], MVPX[NCOEFF], MVPY[NCOEFF], MVPXY[NCOEFF];
  int nm, p, q;
  for(nm=p=0; p<4; p++)
   for(q=0; q<4; q++, nm++)
    { MVP[nm]   = X1Powers[p]*X2Powers[q];
      MVPX[nm]  = p==0 ? 0.0 : p*X1Powers[p-1]*X2Powers[q]/L1;
      MVPY[nm]  = q==0 ? 0.0 : q*X1Powers[p]*X2Powers[q-1]/L2;
      MVPXY[nm] = (p==0 || q==0) ? 0.0 : p*q*X1Powers[p-1]*X2Powers[q-1]/(L1*L2);
    };

  /****************************************************************/
  /* for each component of the Phi function vector, look up the   */
  /* coefficients of the interpolating polynomial for that        */
  /* component and use them to get the value of the interpolant   */
  /****************************************************************/
  double P, PX, PY, PXY;
  double *C;
  int nf;
  for(nf=0; nf<nFun; nf++)
   { 
     C=CTable + GetCTableOffset(nf, nFun, n1, N1, n2, N2);

     for(P=PX=PY=PXY=0.0, nm=0; nm<NCOEFF; nm++)
      { P+=C[nm]*MVP[nm];
        PX+=C[nm]*MVPX[nm];
        PY+=C[nm]*MVPY[nm];
        PXY+=C[nm]*MVPXY[nm];
      };

     Phi[4*nf+0]=P;
     Phi[4*nf+1]=PX;
     Phi[4*nf+2]=PY;
     Phi[4*nf+3]=PXY;

   };

}

/****************************************************************/
/* Like EvaluatePlus but also computes unmixed second derivatives.*/
/* Phi[ 6*nf + 0 ] = Phi_{nf}                                   */
/* Phi[ 6*nf + 1 ] = d/dx Phi_{nf}                              */
/* Phi[ 6*nf + 2 ] = d/dy Phi_{nf}                              */
/* Phi[ 6*nf + 3 ] = d^2/dx^2 Phi_{nf}                          */
/* Phi[ 6*nf + 4 ] = d^2/dxdy Phi_{nf}                          */
/* Phi[ 6*nf + 5 ] = d^2/dy^2 Phi_{nf}                          */
/****************************************************************/
void Interp2D::EvaluatePlusPlus(double X1, double X2, double *Phi)
{
  /****************************************************************/
  /* construct the scaled and shifted versions of the input       */
  /* coordinates.                                                 */
  /* the scaling and shifting operation is constructed so that    */
  /* the corners of the grid cell are at {(0,0),(0,1),(1,0),(1,1)}*/
  /****************************************************************/
  int n1, n2;
  double X1Bar, X2Bar;
  FindInterval(X1, X1Points, N1, X1Min, DX1, &n1, &X1Bar);
  FindInterval(X2, X2Points, N2, X2Min, DX2, &n2, &X2Bar);

  /****************************************************************/
  /* tabulate powers of the scaled/shifted coordinates            */
  /****************************************************************/
  double X1Powers[4], X2Powers[4];
  int np;
  X1Powers[0]=X2Powers[0]=1.0;
  for(np=1; np<4; np++)
   { X1Powers[np]=X1Powers[np-1]*X1Bar;
     X2Powers[np]=X2Powers[np-1]*X2Bar;
   };

  /****************************************************************/
  /* get the widths of this grid cell *****************************/
  /****************************************************************/
  double L1, L2;
  if (X1Points)
   { L1=X1Points[n1+1]-X1Points[n1];
     L2=X2Points[n2+1]-X2Points[n2];
   }
  else 
   { L1 = DX1;
     L2 = DX2;
   };

  /****************************************************************/
  /* construct the NCOEFF-element vector of monomial values that  */
  /* I dot with the NCOEFF-element vector of coefficients to get  */
  /* the value of the interpolating polynomial                    */
  /****************************************************************/
  double MVP[NCOEFF], MVPX[NCOEFF], MVPY[NCOEFF]; 
  double MVPXX[NCOEFF], MVPXY[NCOEFF], MVPYY[NCOEFF];
  int nm, p, q;
  for(nm=p=0; p<4; p++)
   for(q=0; q<4; q++, nm++)
    { MVP[nm]   = X1Powers[p]*X2Powers[q];
      MVPX[nm]  = p==0 ? 0.0 : p*X1Powers[p-1]*X2Powers[q]/L1;
      MVPY[nm]  = q==0 ? 0.0 : q*X1Powers[p]*X2Powers[q-1]/L2;
      MVPXX[nm] = (p<=1)         ? 0.0 : p*(p-1)*X1Powers[p-2]*X2Powers[q]/(L1*L1);
      MVPXY[nm] = (p==0 || q==0) ? 0.0 : p*q*X1Powers[p-1]*X2Powers[q-1]/(L1*L2);
      MVPYY[nm] = (q<=1)         ? 0.0 : q*(q-1)*X1Powers[p]*X2Powers[q-2]/(L2*L2);
    };

  /****************************************************************/
  /* for each component of the Phi function vector, look up the   */
  /* coefficients of the interpolating polynomial for that        */
  /* component and use them to get the value of the interpolant   */
  /****************************************************************/
  for(int nf=0; nf<nFun; nf++)
   { 
     double *C=CTable + GetCTableOffset(nf, nFun, n1, N1, n2, N2);

     double P, PX, PY, PXX, PXY, PYY;
     for(P=PX=PY=PXX=PXY=PYY=0.0, nm=0; nm<NCOEFF; nm++)
      { P   += C[nm]*MVP[nm];
        PX  += C[nm]*MVPX[nm];
        PY  += C[nm]*MVPY[nm];
        PXX += C[nm]*MVPXX[nm];
        PXY += C[nm]*MVPXY[nm];
        PYY += C[nm]*MVPYY[nm];
      };

     Phi[6*nf+0]=P;
     Phi[6*nf+1]=PX;
     Phi[6*nf+2]=PY;
     Phi[6*nf+3]=PXX;
     Phi[6*nf+4]=PXY;
     Phi[6*nf+5]=PYY;

   };

}
