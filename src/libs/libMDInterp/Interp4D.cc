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
 *  Interp4D.cc -- implementation of Interp4D class
 *
 *  homer reid     -- 3/2011
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
#define NEQUATIONS 256 // number of matrix rows
#define NCOEFF     256 // number of matrix columns
#define NDATA      16  // number of values stored at each grid point

/****************************************************************/
/****************************************************************/
/****************************************************************/
static int GetPhiVDTableOffset(int nf, int nFun, int n1, int N1, 
                               int n2, int N2, int n3, int N3, 
                               int n4, int N4)
{
  (void) N1; // unused;
  return NDATA*(nf + nFun*(n4 + N4*(n3 + N3*(n2 + N2*n1))));
}

/****************************************************************/
/****************************************************************/
/****************************************************************/
static int GetCTableOffset(int nf, int nFun, int n1, int N1,
                           int n2, int N2, int n3, int N3, 
                           int n4, int N4)
{
  (void) N1; // unused;
  return NCOEFF*(nf + nFun*(n4 + (N4-1)*(n3 + (N3-1)*(n2 + (N2-1)*n1))));
}

/****************************************************************/
/* thread routine used to compute values and derivatives of     */
/* user's function at grid points                               */
/****************************************************************/
typedef struct ThreadData
 {
   double *X1Points, *X2Points, *X3Points, *X4Points;
   int N1, N2, N3, N4;
   double X1Min, X2Min, X3Min, X4Min;
   double DX1, DX2, DX3, DX4;

   int nFun;

   Phi4D PhiFunc; 
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
  double *X3Points=TD->X3Points;
  double *X4Points=TD->X4Points;
  int N1   = TD->N1;
  int N2   = TD->N2;
  int N3   = TD->N3;
  int N4   = TD->N4;
  double X1Min=TD->X1Min;
  double X2Min=TD->X2Min;
  double X3Min=TD->X3Min;
  double X4Min=TD->X4Min;
  double DX1=TD->DX1;
  double DX2=TD->DX2;
  double DX3=TD->DX3;
  double DX4=TD->DX4;
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
  int nPoints=0, TotalPoints=N1*N2*N3*N4, pc;
  int n1, n2, n3, n4;
  for(n1=0; n1<N1; n1++)
   for(n2=0; n2<N2; n2++)
    for(n3=0; n3<N3; n3++)
     for(n4=0; n4<N4; n4++)
      {

        /*--------------------------------------------------------------*/
        /*- status logging ---------------------------------------------*/
        /*--------------------------------------------------------------*/
        if (TD->nt==0)
         { nPoints++;
           for(pc=10; pc<=90; pc+=10)
            if (nPoints == (pc*TotalPoints/100) )
             Log("...%i %%",pc);
         };
        /*--------------------------------------------------------------*/
        /*--------------------------------------------------------------*/
        /*--------------------------------------------------------------*/

        nt++;
        if (nt==TD->nThread) nt=0;
        if (nt!=TD->nt) continue;

       if (X1Points)
        TD->PhiFunc(X1Points[n1], X2Points[n2], X3Points[n3], X4Points[n4], 
                    TD->UserData, PhiVD);
       else
        TD->PhiFunc(X1Min + n1*DX1, X2Min + n2*DX2, X3Min + n3*DX3, X4Min + n4*DX4,
                    TD->UserData, PhiVD);

        for(nf=0; nf<nFun; nf++)
         { 
           P = PhiVDTable + GetPhiVDTableOffset(nf, nFun, n1, N1, n2, N2, n3, N3, n4, N4);
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
Interp4D::Interp4D(double *pX1Points, int pN1,
                   double *pX2Points, int pN2,
                   double *pX3Points, int pN3,
                   double *pX4Points, int pN4,
                   int pnFun, Phi4D PhiFunc, void *UserData)
{
   /*--------------------------------------------------------------*/ 
   /*- simple setup: initialize class data                ---------*/ 
   /*--------------------------------------------------------------*/ 
   if (pN1<2 || pN2<2 || pN3<2 || pN4<2)
    ErrExit("Interp4D: grid must have 2 or more points in every dimension");

   N1=pN1; 
   N2=pN2;
   N3=pN3; 
   N4=pN4; 
   X1Points=(double *)memdup(pX1Points, N1*sizeof(double));
   X2Points=(double *)memdup(pX2Points, N2*sizeof(double));
   X3Points=(double *)memdup(pX3Points, N3*sizeof(double));
   X4Points=(double *)memdup(pX4Points, N4*sizeof(double));

   nFun=pnFun;

   /*--------------------------------------------------------------*/ 
   /*- allocate space for the CTable, which stores all             */ 
   /*- coefficients for all functions and all grid cells.          */
   /*-                                                             */
   /*- how it works:                                               */
   /*-                                                             */
   /*- CTable + GetCTableOffset(nf,nFun,n1,N1,n2,N2,n3,N3,n4,N4)   */
   /*-                                                             */
   /*- is a pointer to an array of NCOEFF doubles which are the    */
   /*- coefficients for function #nf within grid cell (n1,n2,n3,n4)*/
   /*-                                                             */
   /*- grid cell (n1,n2,n3,n4) is defined to be the region         */
   /*-  X1Points[n1] <= X1 < X1Points[n1+1]                        */
   /*-  X2Points[n2] <= X2 < X2Points[n2+1]                        */
   /*-  X3Points[n3] <= X3 < X3Points[n3+1].                       */
   /*-  X4Points[n4] <= X4 < X4Points[n4+1].                       */
   /*--------------------------------------------------------------*/ 
   CTable=(double *)mallocEC((N1-1)*(N2-1)*(N3-1)*(N4-1)*nFun*NCOEFF*sizeof(double));
   if (!CTable)
    ErrExit("%s:%i:out of memory",__FILE__,__LINE__);

   /*--------------------------------------------------------------*/ 
   /*--------------------------------------------------------------*/ 
   /*--------------------------------------------------------------*/ 
   InitInterp4D(PhiFunc, UserData);
}

/****************************************************************/
/* class constructor 1: construct the class from a user-supplied*/
/* function and grid                                            */
/****************************************************************/
Interp4D::Interp4D(double pX1Min, double X1Max, int pN1,
                   double pX2Min, double X2Max, int pN2,
                   double pX3Min, double X3Max, int pN3,
                   double pX4Min, double X4Max, int pN4,
                   int pnFun, Phi4D PhiFunc, void *UserData)
{
   /*--------------------------------------------------------------*/ 
   /*- simple setup: initialize class data                ---------*/ 
   /*--------------------------------------------------------------*/ 
   if (pN1<=2 || pN2<=2 || pN3<=2 || pN4<=2)
    ErrExit("Interp4D: grid must have 2 or more points in every dimension");

   N1=pN1; 
   X1Points=0;
   X1Min=pX1Min;
   DX1 = (X1Max - X1Min) / ((double)(N1-1));

   N2=pN2; 
   X2Points=0;
   X2Min=pX2Min;
   DX2 = (X2Max - X2Min) / ((double)(N2-1));

   N3=pN3; 
   X3Points=0;
   X3Min=pX3Min;
   DX3 = (X3Max - X3Min) / ((double)(N3-1));

   N4=pN4; 
   X4Points=0;
   X4Min=pX4Min;
   DX4 = (X4Max - X4Min) / ((double)(N4-1));

   if ( DX1==0.0 || DX2==0.0 || DX3==0.0 || DX4==0.0 )
    ErrExit("Interp4D: grid has zero size in one or more dimensions");

   nFun=pnFun;

   /*--------------------------------------------------------------*/ 
   /*- allocate space for the CTable (see previous routine)        */ 
   /*--------------------------------------------------------------*/ 
   CTable=(double *)mallocEC((N1-1)*(N2-1)*(N3-1)*(N4-1)*nFun*NCOEFF*sizeof(double));
   if (!CTable)
    ErrExit("%s:%i:out of memory",__FILE__,__LINE__);

   /*--------------------------------------------------------------*/ 
   /*--------------------------------------------------------------*/ 
   /*--------------------------------------------------------------*/ 
   InitInterp4D(PhiFunc, UserData);
}

/****************************************************************/
/****************************************************************/
/****************************************************************/
void Interp4D::InitInterp4D(Phi4D PhiFunc, void *UserData)
{
   Log("Creating interpolation grid with %i grid points...",N1*N2*N3*N4);

   /*--------------------------------------------------------------*/ 
   /*- temporarily allocate a big array in which we will store     */ 
   /*- values and derivatives of the user's function. allocating   */ 
   /*- this entire table all at once is totally storage-inefficient*/ 
   /*- and could be much optimized since we only ever need 8 grid  */ 
   /*- points worth of data at a time.                             */ 
   /*-                                                             */
   /*- ('PVD' = 'phi values and derivatives')                      */
   /*-                                                             */
   /*- how it works:                                               */
   /*-                                                             */
   /*-  PVDTable + GetPVDTableOffset(nf,nFun,n1,N1,n2,N2,n3,N3,n4,N4) */
   /*-                                                             */
   /*- is a pointer to an array of NDATA doubles which are the     */
   /*- value and derivatives of function #nf at grid point         */
   /*- (n1,n2,n3,n4).                                              */
   /*--------------------------------------------------------------*/
   double *PhiVDTable=(double *)mallocEC(N1*N2*N3*N4*nFun*NDATA*sizeof(double));
   if (!PhiVDTable)
    ErrExit("%s:%i:out of memory",__FILE__,__LINE__);

   /*--------------------------------------------------------------*/
   /*- fire off threads that will call the user's function to     -*/
   /*- populate the PVDTable table.                               -*/
   /*--------------------------------------------------------------*/
   Log("Computing user function at grid points...");
   int nThread=GetNumThreads();
   
#ifdef USE_PTHREAD
   ThreadData *TDs = new ThreadData[nThread], *TD;
   pthread_t *Threads = new pthread_t[nThread];
#endif
   ThreadData TD1;
   int nt;

   TD1.X1Points=X1Points;
   TD1.X2Points=X2Points;
   TD1.X3Points=X3Points;
   TD1.X4Points=X4Points;
   TD1.N1=N1;
   TD1.N2=N2;
   TD1.N3=N3;
   TD1.N4=N4;
   TD1.X1Min=X1Min;
   TD1.X2Min=X2Min;
   TD1.X3Min=X3Min;
   TD1.X4Min=X4Min;
   TD1.DX1=DX1;
   TD1.DX2=DX2;
   TD1.DX3=DX3;
   TD1.DX4=DX4;
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
   /*- the outer loop (actually four loops) loops over the 16     -*/
   /*- corners of the grid cell.                                  -*/
   /*-                                                            -*/
   /*- the inner loop (actually four loops) loops over the NCOEFF -*/
   /*- monomials in the interpolating polynomial.                 -*/
   /*--------------------------------------------------------------*/
   double X1Bar, X2Bar, X3Bar, X4Bar;
   double X1Powers[5], X2Powers[5], X3Powers[5], X4Powers[5];
   int nm, ncp, p, q, r, s;
   HMatrix *M=new HMatrix(NEQUATIONS, NCOEFF);

   M->Zero();

   /* note: X1Powers[p] = 0        ,  p==0 */
   /*                   = X1^{p-1} ,  p>0  */
   X1Powers[0]=0.0; X1Powers[1]=1.0; 
   X2Powers[0]=0.0; X2Powers[1]=1.0; 
   X3Powers[0]=0.0; X3Powers[1]=1.0; 
   X4Powers[0]=0.0; X4Powers[1]=1.0; 

   for(ncp=0, X1Bar=0.0; X1Bar<=1.0; X1Bar+=1.0)
    for(X2Bar=0.0; X2Bar<=1.0; X2Bar+=1.0)
     for(X3Bar=0.0; X3Bar<=1.0; X3Bar+=1.0)
      for(X4Bar=0.0; X4Bar<=1.0; X4Bar+=1.0, ncp++)
       { 
         X1Powers[2]=X1Powers[3]=X1Powers[4]=X1Bar;
         X2Powers[2]=X2Powers[3]=X2Powers[4]=X2Bar;
         X3Powers[2]=X3Powers[3]=X3Powers[4]=X3Bar;
         X4Powers[2]=X4Powers[3]=X4Powers[4]=X4Bar;

         for(nm=0,p=0; p<4; p++)
          for(q=0; q<4; q++)
           for(r=0; r<4; r++)
            for(s=0; s<4; s++, nm++)
             { 
               /* the nmth monomial is x1^p x2^q x3^r x4^s            */
               /* its x1 derivative is p x1^{p-1} x2^q x3^r x4^s, etc */
               M->SetEntry(NDATA*ncp + 0, nm, X1Powers[p+1]*X2Powers[q+1]*X3Powers[r+1]*X4Powers[s+1]);
               M->SetEntry(NDATA*ncp + 1, nm, p*X1Powers[p]*X2Powers[q+1]*X3Powers[r+1]*X4Powers[s+1]);
               M->SetEntry(NDATA*ncp + 2, nm, q*X1Powers[p+1]*X2Powers[q]*X3Powers[r+1]*X4Powers[s+1]);
               M->SetEntry(NDATA*ncp + 3, nm, r*X1Powers[p+1]*X2Powers[q+1]*X3Powers[r]*X4Powers[s+1]);
               M->SetEntry(NDATA*ncp + 4, nm, s*X1Powers[p+1]*X2Powers[q+1]*X3Powers[r+1]*X4Powers[s]);
               M->SetEntry(NDATA*ncp + 5, nm, p*q*X1Powers[p]*X2Powers[q]*X3Powers[r+1]*X4Powers[s+1]);
               M->SetEntry(NDATA*ncp + 6, nm, p*r*X1Powers[p]*X2Powers[q+1]*X3Powers[r]*X4Powers[s+1]);
               M->SetEntry(NDATA*ncp + 7, nm, p*s*X1Powers[p]*X2Powers[q+1]*X3Powers[r+1]*X4Powers[s]);
               M->SetEntry(NDATA*ncp + 8, nm, q*r*X1Powers[p+1]*X2Powers[q]*X3Powers[r]*X4Powers[s+1]);
               M->SetEntry(NDATA*ncp + 9, nm, q*s*X1Powers[p+1]*X2Powers[q]*X3Powers[r+1]*X4Powers[s]);
               M->SetEntry(NDATA*ncp + 10, nm, r*s*X1Powers[p+1]*X2Powers[q+1]*X3Powers[r]*X4Powers[s]);
               M->SetEntry(NDATA*ncp + 11, nm, p*q*r*X1Powers[p]*X2Powers[q]*X3Powers[r]*X4Powers[s+1]);
               M->SetEntry(NDATA*ncp + 12, nm, p*q*s*X1Powers[p]*X2Powers[q]*X3Powers[r+1]*X4Powers[s]);
               M->SetEntry(NDATA*ncp + 13, nm, p*r*s*X1Powers[p]*X2Powers[q+1]*X3Powers[r]*X4Powers[s]);
               M->SetEntry(NDATA*ncp + 14, nm, q*r*s*X1Powers[p+1]*X2Powers[q]*X3Powers[r]*X4Powers[s]);
               M->SetEntry(NDATA*ncp + 15, nm, p*q*r*s*X1Powers[p]*X2Powers[q]*X3Powers[r]*X4Powers[s]);
            };
       }; 
   M->LUFactorize();

   /*--------------------------------------------------------------*/
   /*- now go through and solve a linear system for each grid      */
   /*- cell to compute the coefficients of the interpolating       */
   /*- polynomial in that cell                                     */
   /*--------------------------------------------------------------*/
   int nPoints=0, TotalPoints=(N1-1)*(N2-1)*(N3-1)*(N4-1), pc;
   int n1, dn1, n2, dn2, n3, dn3, n4, dn4;
   int nf;
   double *P;
   HVector *C = new HVector(NCOEFF);

   Log("Computing coefficients of interpolating polynomials...");
   for(n1=0; n1<(N1-1); n1++)
    for(n2=0; n2<(N2-1); n2++)
     for(n3=0; n3<(N3-1); n3++)
      for(n4=0; n4<(N4-1); n4++)
       { 
          for(pc=10; pc<=90; pc++)
           if ( nPoints++ == pc*TotalPoints/100 )
            Log("...%i %%...\n",pc);

          if (X1Points)
           { DX1=X1Points[n1+1]-X1Points[n1];
             DX2=X2Points[n2+1]-X2Points[n2];
             DX3=X3Points[n3+1]-X3Points[n3];
             DX4=X4Points[n4+1]-X4Points[n4];
           };

          /* separately compute interpolation coefficients for each function  */
          for(nf=0; nf<nFun; nf++)
           { 
             /* construct the RHS vector by extracting from PhiVDTable the 16 data */
             /* values for each of the 16 corners of this grid cell                */
             for(ncp=0, dn1=0; dn1<=1; dn1++)
              for(dn2=0; dn2<=1; dn2++)
               for(dn3=0; dn3<=1; dn3++)
                for(dn4=0; dn4<=1; dn4++, ncp++)
                 { 
                   P=PhiVDTable + GetPhiVDTableOffset(nf, nFun, n1+dn1, N1, n2+dn2, N2, n3+dn3, N3, n4+dn4, N4);

                   C->SetEntry( NDATA*ncp +  0, P[0]);
                   C->SetEntry( NDATA*ncp +  1, DX1*P[1]);
                   C->SetEntry( NDATA*ncp +  2, DX2*P[2]);
                   C->SetEntry( NDATA*ncp +  3, DX3*P[3]);
                   C->SetEntry( NDATA*ncp +  4, DX4*P[4]);
                   C->SetEntry( NDATA*ncp +  5, DX1*DX2*P[5]);
                   C->SetEntry( NDATA*ncp +  6, DX1*DX3*P[6]);
                   C->SetEntry( NDATA*ncp +  7, DX1*DX4*P[7]);
                   C->SetEntry( NDATA*ncp +  8, DX2*DX3*P[8]);
                   C->SetEntry( NDATA*ncp +  9, DX2*DX4*P[9]);
                   C->SetEntry( NDATA*ncp + 10, DX3*DX4*P[10]);
                   C->SetEntry( NDATA*ncp + 11, DX1*DX2*DX3*P[11]);
                   C->SetEntry( NDATA*ncp + 12, DX1*DX2*DX4*P[12]);
                   C->SetEntry( NDATA*ncp + 13, DX1*DX3*DX4*P[13]);
                   C->SetEntry( NDATA*ncp + 14, DX2*DX3*DX4*P[14]);
                   C->SetEntry( NDATA*ncp + 15, DX1*DX2*DX3*DX4*P[15]);

                 };

             /* solve the 256x256 system */
             M->LUSolve(C);

             /* store the coefficients in the appropriate place in CTable */
             P=CTable + GetCTableOffset(nf, nFun, n1, N1, n2, N2, n3, N3, n4, N4);
             memcpy(P,C->DV, NCOEFF*sizeof(double));
           };

        };

   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   delete C;
   delete M;
   free(PhiVDTable);
   Log("...done!");

}  

/****************************************************************/
/****************************************************************/
/****************************************************************/
void Interp4D::ReInitialize(Phi4D PhiFunc, void *UserData)
{ 
  InitInterp4D(PhiFunc, UserData);
}

/****************************************************************/
/* class constructor 2: construct the class from a data file    */
/* generated by a previous call to Interp4D::WriteToFile()      */
/****************************************************************/
Interp4D::Interp4D(const char *FileName)
{
   Log("Attempting to read interpolation table from file %s...",FileName);

   FILE *f=fopen(FileName,"r");
   if (!f)
    ErrExit("could not open file %s");
   
   freadEC(&N1  , sizeof(int), 1, f, FileName);
   freadEC(&N2  , sizeof(int), 1, f, FileName);
   freadEC(&N3  , sizeof(int), 1, f, FileName);
   freadEC(&N4  , sizeof(int), 1, f, FileName);
   freadEC(&nFun, sizeof(int), 1, f, FileName);
  
   int NonUniform;
   freadEC(&NonUniform, sizeof(int), 1, f, FileName);

   if (NonUniform)
    { X1Points=(double *)mallocEC(N1*sizeof(double));
      freadEC(X1Points, sizeof(double), N1, f, FileName);
      X2Points=(double *)mallocEC(N2*sizeof(double));
      freadEC(X2Points, sizeof(double), N2, f, FileName);
      X3Points=(double *)mallocEC(N3*sizeof(double));
      freadEC(X3Points, sizeof(double), N3, f, FileName);
      X4Points=(double *)mallocEC(N4*sizeof(double));
      freadEC(X4Points, sizeof(double), N4, f, FileName);
    }
   else
    { X1Points=X2Points=X3Points=X4Points=0;
      freadEC(&X1Min, sizeof(double), 1, f, FileName);
      freadEC(&DX1,   sizeof(double), 1, f, FileName);
      freadEC(&X2Min, sizeof(double), 1, f, FileName);
      freadEC(&DX2,   sizeof(double), 1, f, FileName);
      freadEC(&X3Min, sizeof(double), 1, f, FileName);
      freadEC(&DX3,   sizeof(double), 1, f, FileName);
      freadEC(&X4Min, sizeof(double), 1, f, FileName);
      freadEC(&DX4,   sizeof(double), 1, f, FileName);
    };

   int CTableSize=(N1-1)*(N2-1)*(N3-1)*(N4-1)*nFun*NCOEFF;
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
void Interp4D::WriteToFile(const char *FileName)
{
   Log("Writing interpolation table to file %s...",FileName);

   FILE *f=fopen(FileName,"w");
   if (!f)
    ErrExit("could not open file %s");

   fwrite(&N1  , sizeof(int), 1, f);
   fwrite(&N2  , sizeof(int), 1, f);
   fwrite(&N3  , sizeof(int), 1, f);
   fwrite(&N4  , sizeof(int), 1, f);
   fwrite(&nFun, sizeof(int), 1, f);

   int NonUniform=(X1Points ? 1 : 0);
   fwrite(&NonUniform, sizeof(int), 1, f);
    
   if (NonUniform)
    { fwrite(X1Points, sizeof(double), N1, f);
      fwrite(X2Points, sizeof(double), N2, f);
      fwrite(X3Points, sizeof(double), N3, f);
      fwrite(X4Points, sizeof(double), N4, f);
    }
   else
    { fwrite(&X1Min, sizeof(double), 1, f);
      fwrite(&DX1,   sizeof(double), 1, f);
      fwrite(&X2Min, sizeof(double), 1, f);
      fwrite(&DX2,   sizeof(double), 1, f);
      fwrite(&X3Min, sizeof(double), 1, f);
      fwrite(&DX3,   sizeof(double), 1, f);
      fwrite(&X4Min, sizeof(double), 1, f);
      fwrite(&DX4,   sizeof(double), 1, f);
    };

   int CTableSize=(N1-1)*(N2-1)*(N3-1)*(N4-1)*nFun*NCOEFF;
   fwrite(&CTableSize, sizeof(int), 1, f);
   fwrite(CTable, sizeof(double), CTableSize, f);

   fclose(f);
}

/****************************************************************/
/* class destructor *********************************************/
/****************************************************************/
Interp4D::~Interp4D()
{ 
  if (X1Points) free(X1Points);
  if (X2Points) free(X2Points);
  if (X3Points) free(X3Points);
  if (X4Points) free(X4Points);
  free(CTable);
}

/****************************************************************/
/* the routine that actually does the interpolation             */
/****************************************************************/
void Interp4D::Evaluate(double X1, double X2, double X3, double X4, double *Phi)
{
  int n1, n2, n3, n4;
  double X1Bar, X2Bar, X3Bar, X4Bar;

  FindInterval(X1, X1Points, N1, X1Min, DX1, &n1, &X1Bar);
  FindInterval(X2, X2Points, N2, X2Min, DX2, &n2, &X2Bar);
  FindInterval(X3, X3Points, N3, X3Min, DX3, &n3, &X3Bar);
  FindInterval(X4, X4Points, N4, X4Min, DX4, &n4, &X4Bar);

  /****************************************************************/
  /* tabulate powers of the scaled/shifted coordinates            */
  /****************************************************************/
  int np;
  double X1Powers[4], X2Powers[4], X3Powers[4], X4Powers[4];
  X1Powers[0]=X2Powers[0]=X3Powers[0]=X4Powers[0]=1.0;
  for(np=1; np<4; np++)
   { X1Powers[np]=X1Powers[np-1]*X1Bar;
     X2Powers[np]=X2Powers[np-1]*X2Bar;
     X3Powers[np]=X3Powers[np-1]*X3Bar;
     X4Powers[np]=X4Powers[np-1]*X4Bar;
   };

  /****************************************************************/
  /* construct the NCOEFF-element vector of monomial values that  */
  /* we dot with the NCOEFF-element vector of coefficients to get */
  /* the value of the interpolating polynomial                    */
  /****************************************************************/
  double MV[NCOEFF];
  int nm, p, q, r, s;
  for(nm=p=0; p<4; p++)
   for(q=0; q<4; q++)
    for(r=0; r<4; r++)
     for(s=0; s<4; s++, nm++)
      MV[nm]=X1Powers[p]*X2Powers[q]*X3Powers[r]*X4Powers[s];

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
     C=CTable + GetCTableOffset(nf, nFun, n1, N1, n2, N2, n3, N3, n4, N4);

     for(P=0.0, nm=0; nm<NCOEFF; nm++)
      P+=C[nm]*MV[nm];

     Phi[nf]=P;

   };

}
