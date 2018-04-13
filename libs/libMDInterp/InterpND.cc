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
 *  InterpND.cc -- N-dimensional polynomial interpolation
 *
 *  homer reid  -- 3/2011 - 1/2018
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <vector>

#include <libhrutil.h>
#include <libhmat.h>
#include "libMDInterp.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

/***************************************************************/
/***************************************************************/
/***************************************************************/
bool Increment(iVec &n, iVec &N)
{ 
  for(unsigned d=0; d<n.size(); d++)
   { if ( ++n[d] < N[d] )
      return false;
     n[d]=0;
   };
  return true;
}

double Monomial(dVec xVec, iVec pVec)
{ double Value=1.0;
  for(unsigned d=0; d<pVec.size(); d++)
   for(int p=0; p<pVec[d]; p++)
    Value*=xVec[d];
  return Value;
}

size_t InterpND::GetPointIndex(iVec nVec)
{ size_t Index=0;
  for(int d=0; d<D; d++) 
   Index += nVec[d]*PointStride[d];
  return Index;
}

iVec InterpND::GetPoint(size_t Index)
{ iVec nVec(D);
  for(int d=0; d<D; d++)
   { nVec[d] = Index % (NVec[d]+1);
     Index /= (NVec[d]+1);
   };
  return nVec;
}

size_t InterpND::GetCellIndex(iVec nVec)
{ size_t Index=0;
  for(int d=0; d<D; d++) 
   Index += nVec[d]*CellStride[d];
  return Index;
}

iVec InterpND::GetCell(size_t Index)
{ iVec nVec(D);
  for(int d=0; d<D; d++)
   { nVec[d] = Index % NVec[d];
     Index /= NVec[d];
   };
  return nVec;
}

dVec InterpND::n2x(iVec nVec)
{ dVec xVec(D);
  for(int d=0, n=nVec[0]; d<D; n=nVec[++d])
   xVec[d] = xPoints[d] ? xPoints[d][n] : XMin[d] + n*DX[d];
  return xVec;
}

/****************************************************************/
/* layout of PhiVDTable:                                        */
/*  Phi values, derivatives for function #1  at point #1        */
/*  Phi values, derivatives for function #2  at point #1        */
/*  ...                                                         */
/*  Phi values, derivatives for function #NF at point #1        */
/*  Phi values, derivatives for function #1  at point #2        */
/*  Phi values, derivatives for function #2  at point #2        */
/*  ...                                                         */
/*  Phi values, derivatives for function #NF at point #NPoint   */
/* layout of PhiVDTable:                                        */
/*  Phi values, derivatives for function #1  at point #1        */
/****************************************************************/
size_t InterpND::GetPhiVDTableOffset(int nPoint, int nFun)
{ return nPoint*NF*NVD + nFun*NVD; }

size_t InterpND::GetPhiVDTableOffset(iVec nVec, iVec tauVec, int nFun)
{
  iVec npVec(D);
  for(int d=0; d<D; d++)
   npVec[d]=nVec[d] + tauVec[d];
  return GetPhiVDTableOffset(GetPointIndex(npVec),nFun);
}

/****************************************************************/
/* layout of CTable:                                            */
/*  Coefficients for function #1 in grid cell #1                */
/*  Coefficients for function #2 in grid cell #1                */
/*  ...                                                         */
/*  Coefficients for function #NF in grid cell #1               */
/*  Coefficients for function #1  in grid cell #2               */
/*  ...                                                         */
/*  Coefficients for function #NF in grid cell #NCell           */
/****************************************************************/
size_t InterpND::GetCTableOffset(int nCell, int nFun)
{ return nCell*NF*NCoeff+ nFun*NCoeff; }

size_t InterpND::GetCTableOffset(iVec nVec, int nFun)
{ return GetCTableOffset(GetCellIndex(nVec),nFun); }

/****************************************************************/
/* class constructor 1: construct the class from a user-supplied*/
/* function and uniform grids                                   */
/****************************************************************/
InterpND::InterpND(dVec _XMin, dVec XMax, iVec _NVec, int _NF,
                   PhiVDFunc UserFunc, void *UserData, bool Verbose)
 : NVec(_NVec), NF(_NF), XMin(_XMin), DX(XMax)
{   
  memset(xPoints, 0, MAXDIM*sizeof(double *));
  D = XMin.size();
  for(int d=0; d<D; d++)
   { if (NVec[d]<2)
      ErrExit("%s:%i: grid must have 2 or more points in all dimensions (N[%i]=%i)",__FILE__,__LINE__,d,NVec[d]);
     DX[d] = (XMax[d] - XMin[d]) / ((double)(NVec[d]));
   };

  Initialize(UserFunc, UserData, Verbose);
}

/****************************************************************/
/****************************************************************/
/****************************************************************/
InterpND::~InterpND()
{ 
  if (CTable)
   free(CTable);
}

/****************************************************************/
/****************************************************************/
/****************************************************************/
void InterpND::Initialize(PhiVDFunc UserFunc, void *UserData, bool Verbose)
{
   /*--------------------------------------------------------------*/
   /*- compute some statistics ------------------------------------*/
   /*--------------------------------------------------------------*/
   D = NVec.size();
   CellStride.reserve(D);
   PointStride.reserve(D);
   CellStride[0]=PointStride[0]=1;
   size_t NCell  = NVec[0];        // # grid cells
   size_t NPoint = NVec[0]+1;      // # grid points
   for(int d=1; d<D; d++)
    { CellStride[d] = NCell;
      PointStride[d] = NPoint;
      NCell  *= NVec[d];
      NPoint *= NVec[d]+1;
    };

   NVD    = (1<<D);     // # function vals, derivs per grid point
   NCoeff = NVD*NVD;    // # polynomial coefficients per grid cell

   size_t PhiVDTableSize = (NPoint * NF * NVD) * sizeof(double);
   double *PhiVDTable=(double *)mallocEC(PhiVDTableSize);
   size_t CTableSize = (NCell * NF * NCoeff) * sizeof(double);
   CTable=(double *)mallocEC(CTableSize);
   if (!PhiVDTable || !CTable)
    ErrExit("%s:%i:out of memory",__FILE__,__LINE__);

   /*--------------------------------------------------------------*/
   /*- call user's function to populate table of function values   */
   /*- and derivatives at grid points                              */
   /*--------------------------------------------------------------*/
   int NumThreads = GetNumThreads();
   Log("Evaluating interpolation function at %i points (%i threads)",NPoint,NumThreads);
#ifdef USE_OPENMP
#pragma omp parallel for num_threads(NumThreads)
#endif
   for(size_t nPoint=0; nPoint<NPoint; nPoint++)
    {
      if (Verbose) LogPercent(nPoint,NPoint);
      size_t Offset = GetPhiVDTableOffset(nPoint,0);
      dVec xVec = n2x(GetPoint(nPoint));
      UserFunc(&(xVec[0]), UserData, PhiVDTable + Offset);
    };

   /*--------------------------------------------------------------*/
   /*- construct the NCoeff x NCoeff matrix that operates on a     */
   /*- vector of polynomial coefficients to yield a vector of      */
   /*- function values and derivatives at grid points bounding     */
   /*- a single grid cell                                          */
   /*--------------------------------------------------------------*/
   int NEq    = NCoeff;     // # matching equations per grid cell
   HMatrix *M = new HMatrix(NEq, NCoeff);
   M->Zero();
   iVec Twos(D,2), Fours(D,4);
   // yPowers[0,1][p] = (\pm 1)^p
   double yPowers[2][4]={ {1.0,-1.0,1.0,-1.0}, {1.0, 1.0, 1.0, 1.0} };
   LOOP_OVER_IVECS(nTau, tauVec, Twos)
    { LOOP_OVER_IVECS(nVD, sigmaVec, Twos)
       { int nEq = nTau*NVD + nVD; // index of equation [0,4^D-1]
         LOOP_OVER_IVECS(nCoeff, pVec, Fours)
          { double Entry=1.0;
            for(int d=0; d<D; d++)
             { int tau=tauVec[d], sigma=sigmaVec[d], p=pVec[d];
               if (sigma)
                Entry *= (p==0 ? 0.0 : p*yPowers[tau][p-1]);
               else
                Entry *= yPowers[tau][p];
             };
            M->SetEntry(nEq, nCoeff, Entry);
          }
       }
    }
   M->LUFactorize();

   /*--------------------------------------------------------------*/
   /*- loop over grid cells; for each grid cell and each component */
   /*- of the function vector, populate a vector of function values*/
   /*- and derivatives at the grid-cell corners, then solve linear */
   /*- system to compute polynomial coefficients for this grid cell*/
   /*--------------------------------------------------------------*/
   LOOP_OVER_IVECS(nCell, nVec, NVec)
    for(int nf=0; nf<NF; nf++)
     { 
       // vector of scaling factors to accomodate grid-cell dimensions
       dVec D2Vec(D); // DeltaOver2 vector
       for(int d=0, n=nVec[d]; d<D; d++, n=(d==D ? 0 : nVec[d]))
        D2Vec[d] = 0.5*( xPoints[d] ? (xPoints[d][n+1]-xPoints[d][n])
                                    : DX[d]
                       );

       // populate RHS vector with function values and derivatives
       // at all corners of grid cell
       HVector RHS(NCoeff, LHM_REAL, CTable + GetCTableOffset(nCell,nf));
       LOOP_OVER_IVECS(nTau, tauVec, Twos)
        { double *PhiVD = PhiVDTable + GetPhiVDTableOffset(nVec, tauVec, nf);
          LOOP_OVER_IVECS(nSigma, sigmaVec, Twos)
           RHS.SetEntry(nTau*NVD + nSigma, Monomial(D2Vec, sigmaVec)*PhiVD[nSigma]);
        };

       // operate with inverse M matrix to yield C coefficients
       M->LUSolve(&RHS);
     };

   free(PhiVDTable);
   delete M;
}

/****************************************************************/
/****************************************************************/
/****************************************************************/
void InterpND::Evaluate(double *xVec, double *Phi)
{
  /****************************************************************/
  /****************************************************************/
  /****************************************************************/
  iVec nVec(D);
  dVec xBarVec(D);
  for(int d=0; d<D; d++)
   { int n;
     double xBar;
     FindInterval(xVec[d], xPoints[d], NVec[d], XMin[d], DX[d], &n, &xBar);
     nVec[d]=n;
     xBarVec[d] = 2.0*(xBar-0.5);
   };
  
  /****************************************************************/
  /* tabulate powers of the scaled/shifted coordinates            */
  /****************************************************************/
  double xPowers[MAXDIM][4];
  for(int d=0; d<D; d++)
   { xPowers[d][0]=1.0;
     xPowers[d][1]=xBarVec[d];
     xPowers[d][2]=xPowers[d][1]*xBarVec[d];
     xPowers[d][3]=xPowers[d][2]*xBarVec[d];
   };

  /****************************************************************/
  /* construct the NCOEFF-element vector of monomial values that  */
  /* I dot with the NCOEFF-element vector of coefficients to get  */
  /* the value of the interpolating polynomial                    */
  /****************************************************************/
  memset(Phi, 0, NF*sizeof(double));
  iVec Fours(D,4);
  LOOP_OVER_IVECS(nCoeff,pVec,Fours)
   { double xFactor=1.0;
     for(int d=0; d<D; d++)
      xFactor*=xPowers[d][pVec[d]];
     for(int nf=0; nf<NF; nf++)
      { double *C = CTable + GetCTableOffset(nVec, nf);
        Phi[nf] += C[nCoeff]*xFactor;
      };
   };
}

/****************************************************************/
/****************************************************************/
/****************************************************************/
double InterpND::PlotInterpolationError(PhiVDFunc UserFunc,
                                        void *UserData,
                                        char *OutFileName,
                                        bool CentersOnly)
{
  FILE *f=fopen(OutFileName,"w");
  if (!f) return 0.0;

  int NVD = 1 << D;
  double *PhiExact  = new double[NVD*NF];
  double *PhiInterp = new double[NF];

  iVec Threes(D, 3);
  double MaxRelErr=0.0;
  LOOP_OVER_IVECS(nCell, nVec, NVec)
   { LOOP_OVER_IVECS(nTau, tauVec, Threes)
      {
        //
        bool Skip=false;
        if (CentersOnly)
         for(int d=0; !Skip && d<D; d++)
          if (tauVec[d]!=1) Skip=true;
        for(int d=0; !Skip && d<D; d++)
         if (tauVec[d]==2 && nVec[d]!=(NVec[d]-1) )
          Skip=true;
        if (Skip) continue;

        dVec X(D);
        for(int d=0; d<D; d++)
         X[d] = XMin[d] + (nVec[d] + 0.5*tauVec[d])*DX[d];

        UserFunc(&(X[0]), UserData, PhiExact);
        Evaluate(&(X[0]), PhiInterp);

        for(int nf=0; nf<NF; nf++)
         { double Exact  = PhiExact[nf*NVD+0];
           double Interp = PhiInterp[nf];
           if (fabs(Exact) > 0.0)
            MaxRelErr=fmax(MaxRelErr, fabs(Exact-Interp)/fabs(Exact));
         }

        fprintVec(f,&(X[0]),D);
        for(int d=0; d<D; d++)
         fprintf(f,"%e ",nVec[d] + 0.5*tauVec[d]);
        for(int nf=0; nf<NF; nf++)
         fprintf(f,"%e ",PhiExact[nf*NVD+0]);
        fprintVecCR(f,PhiInterp,NF);
      }
   }
  fclose(f);
  delete[] PhiExact;
  delete[] PhiInterp;
  return MaxRelErr;
}

/****************************************************************/
/* If PhiVEMatrix is non-null, it should be an NFx2 matrix      */
/* which will be filled in as follows on return:                */
/* Column 1: mean value of Phi at sample points                 */
/* Column 2: mean absolute error in Phi at sample points        */
/****************************************************************/
double GetInterpolationError(PhiVDFunc UserFunc, void *UserData, int NF,
                             dVec XVec, dVec dXVec,
                             double *MeanRelError, double *MeanAbsError)
{ 
  int D = XVec.size();
  dVec XMax = XVec;
  for(int d=0; d<D; d++) XMax[d]+=dXVec[d];
  iVec NVec(D,2);

  InterpND Interp(XVec, XMax, NVec, NF, UserFunc, UserData);

  int NVD = (1<<D);     // # function vals, derivs per grid point

  double *PhiExact  = new double[NVD*NF];
  double *PhiInterp = new double[NF];

  bool OwnBuffers=false;
  if (MeanRelError==0 || MeanAbsError==0)
   { MeanRelError = new double[D*NF];
     MeanAbsError = new double[D*NF];
     OwnBuffers=true;
   }
  memset(MeanRelError, 0, D*NF*sizeof(double));
  memset(MeanAbsError, 0, D*NF*sizeof(double));

  double NumSamples = (double)(1<<(D-1));
  iVec Twos(D,2);
  for(int d=0; d<D; d++)
   {
     LOOP_OVER_IVECS(nTau, tauVec, Twos)
      {
        if (tauVec[d]!=0) continue;

        dVec XSample(D);
        for(int dd=0; dd<D; dd++)
         XSample[dd] = XVec[dd] + (dd==d ? 0.5 : tauVec[dd])*dXVec[dd];

        UserFunc(&(XSample[0]), UserData, PhiExact);
        Interp.Evaluate(&(XSample[0]), PhiInterp);

        for(int nf=0; nf<NF; nf++)
         { double Error = fabs(PhiExact[nf*NVD+0] - PhiInterp[nf])/NumSamples;
           MeanAbsError[d*NF + nf] += Error;
           if ( fabs(PhiExact[nf*NVD+0])!=0.0 )
            MeanRelError[d*NF + nf] += Error/fabs(PhiExact[nf*NVD+0]);
         }
      }
   }

  double MaxRE = 0.0;
  for(int nfd=0; nfd<NF*D; nfd++)
   MaxRE = fmax(MaxRE, MeanRelError[nfd]);

  delete[] PhiExact;
  delete[] PhiInterp;
  if (OwnBuffers)
   { delete[] MeanAbsError;
     delete[] MeanRelError;
   }

  return MaxRE;
}
