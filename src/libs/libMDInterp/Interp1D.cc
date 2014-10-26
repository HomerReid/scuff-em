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
 *  Interp1D.cc -- implementation of Interp1D class for 1D interpolation
 *
 *  homer reid  -- 5/2011
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
#define NEQUATIONS 4  // number of matrix rows
#define NCOEFF     4  // number of matrix columns
#define NDATA      2  // number of values stored at each grid point

/****************************************************************/
/****************************************************************/
/****************************************************************/
static int GetPhiVDTableOffset(int nf, int nFun, int n)
{
  return NDATA*(nf + nFun*n);
}

/****************************************************************/
/****************************************************************/
/****************************************************************/
static int GetCTableOffset(int nf, int nFun, int n)
{
  return NCOEFF*(nf + nFun*n);
}

/****************************************************************/
/* thread routine used to compute values and derivatives of     */
/* user's function at grid points                               */
/****************************************************************/
typedef struct ThreadData
 {
   double *XPoints;
   double XMin, DX;
   int N;

   int nFun;

   Phi1D PhiFunc;
   void *UserData;
   double *PhiVDTable;

   int nt, nThread;

 } ThreadData;

void *GetPhiVD_Thread(void *data)
{
  
  ThreadData *TD=(ThreadData *)data;

  /***************************************************************/
  /* variables unpacked from thread data structure ***************/
  /***************************************************************/
  double *XPoints    = TD->XPoints;
  double XMin        = TD->XMin;
  double DX          = TD->DX;
  int N              = TD->N;
  int nFun           = TD->nFun;
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
  int n, nPoints=0, pc;
  double X;
  for(n=0; n<N; n++)
   {
     /*--------------------------------------------------------------*/
     /*- status logging ---------------------------------------------*/
     /*--------------------------------------------------------------*/
     if (TD->nt==0)
      { nPoints++;
        for(pc=10; pc<=90; pc+=10)
         if (nPoints == (pc*N/100) )
          Log("...%i %%",pc);
      };

      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      nt++;
      if (nt==TD->nThread) nt=0;
      if (nt!=TD->nt) continue;
 
      if (XPoints)
       X=XPoints[n];
      else
       X=XMin + ((double)n)*DX;
      TD->PhiFunc(X, TD->UserData, PhiVD);

      for(nf=0; nf<nFun; nf++)
       { 
         P = PhiVDTable + GetPhiVDTableOffset(nf, nFun, n);
         memcpy(P, PhiVD + NDATA*nf, NDATA*sizeof(double));
       };
   
     };

  delete[] PhiVD;

  return 0;
  
}

/****************************************************************/
/* class constructor 1: construct the class from a user-supplied*/
/* function and nonuniform grid                                 */
/****************************************************************/
Interp1D::Interp1D(double *pXPoints, int pN, int pnFun,
                   Phi1D PhiFunc, void *UserData)
{ 
  if ( pN < 2 )
   ErrExit("Interp1D: grid must have 2 or more data points");
  if (pnFun<1)
   ErrExit("Interp1D: grid must have 1 or more data values");

  N=pN;
  XPoints=(double *)memdup(pXPoints, N*sizeof(double));
  nFun=pnFun;
  InitInterp1D(PhiFunc, UserData);
}

/****************************************************************/
/* class constructor 2: construct the class from a user-supplied*/
/* uniform and grid                                             */
/****************************************************************/
Interp1D::Interp1D(double pXMin, double XMax, int pN,
                   int pnFun, Phi1D PhiFunc, void *UserData)
{ 
  if ( pN < 2 )
   ErrExit("Interp1D: grid must have 2 or more data points");
  if (pnFun<1)
   ErrExit("Interp1D: grid must have 1 or more data values");


  XPoints=0;
  XMin=pXMin;
  N=pN;
  DX = (XMax-XMin) / ( (double)(N-1) );
  nFun=pnFun;
  InitInterp1D(PhiFunc, UserData);
}

/****************************************************************/
/* body of the class constructor for the above 2 cases **********/
/****************************************************************/
void Interp1D::InitInterp1D(Phi1D PhiFunc, void *UserData)
{
   Log("Creating interpolation grid with %i grid points...",N);

   /*--------------------------------------------------------------*/ 
   /*- allocate space for the CTable, which stores all             */ 
   /*- coefficients for all functions and all grid cells.          */
   /*-                                                             */
   /*- how it works:                                               */
   /*-                                                             */
   /*- CTable + GetCTableOffset(nf,nFun,n)                         */
   /*-                                                             */
   /*- is a pointer to an array of NCOEFF doubles which are the    */
   /*- coefficients for function #nf within grid cell #n.          */
   /*-                                                             */
   /*- grid cell #n is defined to be the region                    */
   /*-  XPoints[n] <= X < XPoints[n+1]                             */
   /*--------------------------------------------------------------*/ 
   CTable=(double *)mallocEC((N-1)*nFun*NCOEFF*sizeof(double));
   if (!CTable)
    ErrExit("%s:%i:out of memory",__FILE__,__LINE__);

   /*--------------------------------------------------------------*/ 
   /*- temporarily allocate a big array in which we will store     */ 
   /*- values and derivatives of the user's function. allocating   */ 
   /*- this entire table all at once is totally storage-inefficient*/ 
   /*- and could be much optimized since we only ever need 2       */ 
   /*- grid points worth of data at a time, not the full grid.     */
   /*-                                                             */
   /*- ('PVD' = 'phi values and derivatives')                      */
   /*-                                                             */
   /*- how it works:                                               */
   /*-                                                             */
   /*-  PVDTable + GetPVDTableOffset(nf,nFun,n)                    */
   /*-                                                             */
   /*- is a pointer to an array of NDATA doubles which are the     */
   /*- value and derivatives of function #nf at grid point #n.     */
   /*--------------------------------------------------------------*/
   double *PhiVDTable=(double *)mallocEC(N*nFun*NDATA*sizeof(double));
   if (!PhiVDTable)
    ErrExit("%s:%i:out of memory",__FILE__,__LINE__);

   /*--------------------------------------------------------------*/
   /*- fire off threads that will call the user's function to     -*/
   /*- populate the PVDTable table.                               -*/
   /*--------------------------------------------------------------*/
   Log("Computing user function at grid points...");
   int nThread = GetNumThreads();

#ifdef USE_PTHREAD
   ThreadData *TDs = new ThreadData[nThread], *TD;
   pthread_t *Threads = new pthread_t[nThread];
#endif
   ThreadData TD1;
   int nt;

   TD1.XPoints=XPoints;
   TD1.XMin=XMin;
   TD1.DX=DX;
   TD1.N=N;
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
   /*- the outer loop runs over the 2 corners of the grid cell.   -*/
   /*-                                                            -*/
   /*- the inner loop runs over the NCOEFF monomials in the       -*/
   /*- interpolating polynomial.                                  -*/
   /*--------------------------------------------------------------*/
   double XBar;
   double XPowers[5];
   int ncp, p;
   HMatrix *M=new HMatrix(NEQUATIONS, NCOEFF);

   M->Zero();

   /* note: XPowers[p]  = 0        ,  p==0 */
   /*                   = X^{p-1}  ,  p>0  */
   XPowers[0]=0.0; XPowers[1]=1.0; 

   for(ncp=0, XBar=0.0; XBar<=1.0; XBar+=1.0, ncp++)
    { 
      XPowers[2]=XPowers[3]=XPowers[4]=XBar;

      for(p=0; p<4; p++)
       { 
         /* the pth monomial is x^p                   */
         /* its x derivative is p x^{p-1}             */
         /* since x={0,1}, we have x^p = x for all p  */
         M->SetEntry(NDATA*ncp + 0, p, XPowers[p+1]);
         M->SetEntry(NDATA*ncp + 1, p, p*XPowers[p]);
       };
    }; 
   M->LUFactorize();

   /*--------------------------------------------------------------*/
   /*- now go through and solve a linear system for each grid      */
   /*- cell to compute the coefficients of the interpolating       */
   /*- polynomial in that cell                                     */
   /*--------------------------------------------------------------*/
   int nPoints=0, TotalPoints=(N-1), pc;
   int n, dn;
   int nf;
   double *P;
   HVector *C = new HVector(NCOEFF);

   Log("Computing coefficients of interpolating polynomials...");
   for(n=0; n<(N-1); n++)
    { 
      for(pc=10; pc<=90; pc+=10)
       if ( nPoints++ == pc*TotalPoints/100 )
        Log("...%i %%...\n",pc);

      if (XPoints)
       DX=XPoints[n+1]-XPoints[n];

      /* separately compute interpolation coefficients for each function  */
      for(nf=0; nf<nFun; nf++)
       { 
         /* construct the RHS vector by extracting from PhiVDTable the 2 data */
         /* values for each of the 2 corners of this grid cell.               */
         for(ncp=0, dn=0; dn<=1; dn++, ncp++)
          { P=PhiVDTable + GetPhiVDTableOffset(nf, nFun, n+dn);
            C->SetEntry( NDATA*ncp + 0, P[0]);
            C->SetEntry( NDATA*ncp + 1, DX*P[1]);
          };

         /* solve the 4x4 system */
         M->LUSolve(C);

         /* store the coefficients in the appropriate place in CTable */
         P=CTable + GetCTableOffset(nf, nFun, n);
         memcpy(P, C->DV, NCOEFF*sizeof(double));
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

/*--------------------------------------------------------------*/
/*- class constructor 2: construct from a user-specified set of */
/*- X and Y data values.                                        */
/*-                                                             */
/*- The YPoints array should be laid out as follows:            */
/*-                                                             */
/*-  YPoints[0]      = Phi_0    at XPoints[0]                   */
/*-  YPoints[1]      = Phi_1    at XPoints[0]                   */
/*-  ...                                                        */
/*-  YPoints[nfun-1] = Phi_nfun at XPoints[0]                   */
/*-                                                             */
/*-  YPoints[nfun]   = Phi_0    at XPoints[1]                   */
/*-  YPoints[nfun+1] = Phi_1    at XPoints[1]                   */
/*-  ...                                                        */
/*-                                                             */
/*- and in general                                              */
/*-  YPoints[ n*nfun + nf ] = Phi_nf at XPoints[n]              */
/*--------------------------------------------------------------*/
Interp1D::Interp1D(double *pXPoints, double *YPoints, int pN, int pnFun)
{
   if (pN<2)
    ErrExit("Interp1D: grid must have 2 or more data points");
   if (pnFun<1)
    ErrExit("Interp1D: grid must have 1 or more data values");

   /*--------------------------------------------------------------*/ 
   /*- simple setup: initialize class data                ---------*/ 
   /*--------------------------------------------------------------*/ 
   N=pN; 
   XPoints=(double *)memdup(pXPoints, N*sizeof(double));
   nFun=pnFun;

   /*--------------------------------------------------------------*/ 
   /*- allocate space for the CTable, which stores all             */ 
   /*- coefficients for all functions and all grid cells.          */
   /*-                                                             */
   /*- how it works:                                               */
   /*-                                                             */
   /*- CTable + GetCTableOffset(nf,nFun,n)                         */
   /*-                                                             */
   /*- is a pointer to an array of NCOEFF doubles which are the    */
   /*- coefficients for function #nf within grid cell #n.          */
   /*-                                                             */
   /*- grid cell #n is defined to be the region                    */
   /*-  XPoints[n] <= X < XPoints[n+1]                             */
   /*--------------------------------------------------------------*/ 
   CTable=(double *)mallocEC((N-1)*nFun*NCOEFF*sizeof(double));
   if (!CTable)
    ErrExit("%s:%i:out of memory",__FILE__,__LINE__);

   /*--------------------------------------------------------------*/ 
   /*- now just go through and initialize the coefficients for each*/ 
   /*- grid cell by hand.                                          */ 
   /*-                                                             */ 
   /*- how it works:                                               */ 
   /*-                                                             */ 
   /*- 1. for the first and last grid cells, i do a simple linear  */ 
   /*-    interpolation between the function values at the cell    */ 
   /*-    boundary points                                          */ 
   /*- 2. for interior grid cells (X_n < X < X_{n+1}), i compute   */ 
   /*-    a cubic polynomial in X that matches the function values */ 
   /*-    at points X_{n-1}, X_n, X_{n+1}, and X_{n+2}.            */ 
   /*--------------------------------------------------------------*/ 
   HMatrix *M = new HMatrix(NEQUATIONS,NCOEFF);
   HVector *C = new HVector(NCOEFF);
   int n, nf;
   double Xm1, Xp2, *P;
   for(n=0; n<N-1; n++)
    for(nf=0; nf<nFun; nf++)
     { 
       DX = XPoints[n+1] - XPoints[n];
       P = CTable + GetCTableOffset(nf, nFun, n);

       // handle first and last grid cells separately
       if (n==0 || n==N-2)
        { P[0]=YPoints[ n*nFun + nf ];
          P[1]=YPoints[ (n+1)*nFun + nf ] - P[0];
          P[2]=P[3]=0.0;
        }
       else
        { 
          Xm1=(XPoints[n-1]-XPoints[n])/DX;
          Xp2=(XPoints[n+2]-XPoints[n])/DX;
 
          /* form the matrix and RHS vector of the 4x4 linear system    */
          /* that we solve for the coefficients of the cubic polynomial */
          /* for this grid cell                                         */
          /* the system is:                                             */
          /*  C_0 + C_1 x_{n-1} + C_2 x_{n-1}^2 + C_3 x_{n-1}^3 = y_{n-1} */
          /*  C_0 + C_1 x_{n+0} + C_2 x_{n+0}^2 + C_3 x_{n+0}^3 = y_{n+0} */
          /*  C_0 + C_1 x_{n+1} + C_2 x_{n+1}^2 + C_3 x_{n+1}^3 = y_{n+1} */
          /*  C_0 + C_1 x_{n+2} + C_2 x_{n+2}^2 + C_3 x_{n+2}^3 = y_{n+2} */
          /* (note that in our units we have x_{n+0} = 0, x_{n+1}=1)      */
          M->SetEntry(0,0,1.0);
          M->SetEntry(0,1,Xm1);
          M->SetEntry(0,2,Xm1*Xm1);
          M->SetEntry(0,3,Xm1*Xm1*Xm1);
          M->SetEntry(1,0,1.0);
          M->SetEntry(1,1,0.0);
          M->SetEntry(1,2,0.0);
          M->SetEntry(1,3,0.0);
          M->SetEntry(2,0,1.0);
          M->SetEntry(2,1,1.0);
          M->SetEntry(2,2,1.0);
          M->SetEntry(2,3,1.0);
          M->SetEntry(3,0,1.0);
          M->SetEntry(3,1,Xp2);
          M->SetEntry(3,2,Xp2*Xp2);
          M->SetEntry(3,3,Xp2*Xp2*Xp2);

          /* */
          C->SetEntry(0,YPoints[ (n-1)*nFun + nf ]);
          C->SetEntry(1,YPoints[ (n+0)*nFun + nf ]);
          C->SetEntry(2,YPoints[ (n+1)*nFun + nf ]);
          C->SetEntry(3,YPoints[ (n+2)*nFun + nf ]);

          /* */
          M->LUFactorize();
          M->LUSolve(C);

          /* store the coefficients in the appropriate place in CTable */
          memcpy(P, C->DV, NCOEFF*sizeof(double));
        };
     };
   delete M;
   delete C;

}

/****************************************************************/
/* class constructor 3: construct the class from a data file    */
/* generated by a previous call to Interp1::WriteToFile()      */
/****************************************************************/
Interp1D::Interp1D(const char *FileName)
{
  Log("Attempting to read interpolation table from file %s...",FileName);

  FILE *f=fopen(FileName,"r");
  if (!f)
   ErrExit("could not open file %s");

  freadEC(&Type, sizeof(int), 1, f, FileName);
  freadEC(&N   , sizeof(int), 1, f, FileName);
  freadEC(&nFun, sizeof(int), 1, f, FileName);

  int XPointsIsNULL;
  freadEC(&XPointsIsNULL, sizeof(int), 1, f, FileName);
  if (XPointsIsNULL)
   { XPoints=0;
     freadEC(&XMin, sizeof(double), 1, f, FileName);
     freadEC(&DX  , sizeof(double), 1, f, FileName);
   }
  else
   { XPoints=(double *)mallocEC(N*sizeof(double));
     freadEC(XPoints, sizeof(double), N, f, FileName);
   };

  int CTableSize = (N-1)*nFun*NCOEFF;
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
/* copy constructor                                             */
/****************************************************************/
Interp1D::Interp1D(Interp1D *Original)
{
  XMin = Original->XMin;
  DX   = Original->DX;
  N    = Original->N;
  nFun = Original->nFun;
  Type = Original->Type;

  if (Original->XPoints)
   { XPoints=(double *)mallocEC(N*sizeof(double));
     memcpy(XPoints, Original->XPoints, N*sizeof(double));
   }
  else
   XPoints=0;
     
  int CTableSize = (N-1)*nFun*NCOEFF*sizeof(double);
  CTable=(double *)mallocEC(CTableSize);
  memcpy(CTable, Original->CTable, CTableSize);

}

/****************************************************************/
/* write all class data to a binary file for subsequent recovery*/
/* by the above constructor routine                             */
/****************************************************************/
void Interp1D::WriteToFile(const char *FileName)
{
  Log("Writing interpolation table to file %s...",FileName);

  FILE *f=fopen(FileName,"w");
  if (!f)
   ErrExit("could not open file %s");

  fwrite(&Type, sizeof(int), 1, f); 
  fwrite(&N   , sizeof(int), 1, f); 
  fwrite(&nFun, sizeof(int), 1, f);

  int XPointsIsNULL = (XPoints==0 ? 1 : 0);
  fwrite(&XPointsIsNULL, 1, sizeof(int), f);
  if (XPointsIsNULL)
   { fwrite(&XMin, sizeof(double), 1, f);
     fwrite(&DX  , sizeof(double), 1, f);
   }
  else
   fwrite(XPoints, sizeof(double), N, f);

  int CTableSize = (N-1)*nFun*NCOEFF;
  fwrite(&CTableSize, sizeof(int), 1, f);
  fwrite(CTable, sizeof(double), CTableSize, f);

  fclose(f);
  Log("...success!");

}

/****************************************************************/
/* class destructor *********************************************/
/****************************************************************/
Interp1D::~Interp1D()
{ 
  if (XPoints)
   free(XPoints);
  free(CTable);
}

/****************************************************************/
/* the routine that actually does the interpolation             */
/****************************************************************/
void Interp1D::Evaluate(double X, double *Phi)
{
  /****************************************************************/
  /* find the grid cell in which the given point lies, and set    */
  /* XBar = scaled and shifted value of the coordinate as         */
  /* appropriate for that grid cell.                              */
  /****************************************************************/
  int n;
  double XBar;
  FindInterval(X, XPoints, N, XMin, DX, &n, &XBar);

  /****************************************************************/
  /* tabulate powers of the scaled/shifted coordinates            */
  /****************************************************************/
  int p;
  double XPowers[4];
  XPowers[0]=1.0;
  for(p=1; p<4; p++)
   XPowers[p]=XPowers[p-1]*XBar;

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
     C=CTable + GetCTableOffset(nf, nFun, n);

     for(P=0.0, p=0; p<NCOEFF; p++)
      P+=C[p]*XPowers[p];

     Phi[nf]=P;

   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
double Interp1D::Evaluate(double X)
{ double *Phi = new double[nFun];
  Evaluate(X,Phi);
  double Phi0 = Phi[0];
  delete[] Phi;
  return Phi0;
}
