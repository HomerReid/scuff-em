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
   }
  return true;
}

double Monomial(dVec XVec, iVec pVec)
{ double Value=1.0;
  for(unsigned d=0; d<pVec.size(); d++)
   for(int p=0; p<pVec[d]; p++)
    Value*=XVec[d];
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
   }
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
   }
  return nVec;
}

dVec InterpND::X2X0(dVec XVec)
{ dVec X0Vec = FixedCoordinates;
  for(int d0=0, d=0; d0<D0; d0++)
   if (isinf(X0Vec[d0])) 
    X0Vec[d0] = XVec[d++];
  return X0Vec;
}

dVec InterpND::X02X(dVec X0Vec)
{
  dVec XVec;
  for(int d0=0; d0<D0; d0++)
   if (isinf(FixedCoordinates[d0]))
    XVec.push_back(X0Vec[d0]);
   else if (!EqualFloat(X0Vec[d0],FixedCoordinates[d0]))
    ErrExit("%s:%i: X0Vec[%i]=%e != %e",__FILE__,__LINE__,d0,X0Vec[d0],FixedCoordinates[d0]);
  return XVec;
}

dVec InterpND::n2X(iVec nVec)
{ dVec XVec(D);
  for(int d=0; d<D; d++)
   { int n=nVec[d];
     XVec[d] = (XGrids.size()>0 ? XGrids[d][n] : XMin[d] + n*DX[d]);
   }
  return XVec;
}

dVec InterpND::n2X0(iVec nVec)
{ return X2X0(n2X(nVec)); }
 
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
InterpND::InterpND(PhiVDFunc UserFunc, void *UserData, int _NF,
                   dVec X0Min, dVec X0Max, iVec N0Vec, bool Verbose)
 : NF(_NF)
{ 
  D0=X0Min.size();
  for(int d0=0; d0<D0; d0++)
   if (N0Vec[d0] < 2 || EqualFloat(X0Min[d0],X0Max[d0]) )
    FixedCoordinates.push_back( X0Min[d0] );
   else
    { FixedCoordinates.push_back( HUGE_VAL );
      NVec.push_back( N0Vec[d0] );
      XMin.push_back( X0Min[d0] );
      DX.push_back(  (X0Max[d0] - X0Min[d0]) / N0Vec[d0] );
    }

  Initialize(UserFunc, UserData, Verbose);
}

/****************************************************************/
/* class constructor 2: construct the class from a user-supplied*/
/* function and user-supplied grids for each variable           */
/****************************************************************/
InterpND::InterpND(PhiVDFunc UserFunc, void *UserData, int _NF,
                   vector<dVec> &X0Grids, bool Verbose)
  : NF(_NF)
 
{ 
  D0=X0Grids.size();
  for(int d0=0; d0<D0; d0++)
   if ( X0Grids[d0].size() == 1 )
    FixedCoordinates.push_back( X0Grids[d0][0] );
   else
    { FixedCoordinates.push_back( HUGE_VAL );
      XGrids.push_back(X0Grids[d0]);
      NVec.push_back( X0Grids[d0].size() - 1 );
    }

  Initialize(UserFunc, UserData, Verbose);
}

/****************************************************************/
/* class constructor 3: construct the class from a user-supplied*/
/* function, ranges for the variables, and an error tolerance;  */
/* in this case, non-uniform grids are determined automatically */
/****************************************************************/
InterpND::InterpND(PhiVDFunc UserFunc, void *UserData, int _NF,
                   dVec X0Min, dVec X0Max, double MaxRelError, bool Verbose)
  : NF(_NF)
{ 
  D0=X0Min.size();

  for(int d0=0; d0<D0; d0++)
   if ( EqualFloat(X0Min[d0], X0Max[d0]) )
    FixedCoordinates.push_back( X0Min[d0] );
   else
    { FixedCoordinates.push_back( HUGE_VAL );
      dVec XGrid = GetXGrid(UserFunc, UserData, NF, X0Min, X0Max, d0, MaxRelError);
      XGrids.push_back(XGrid);
      NVec.push_back( XGrid.size() - 1 );
    }
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
      if (NVec[d]<1) 
       ErrExit("%s:%i: invalid number of points (%i) in dimension %i",__FILE__,__LINE__,NVec[d],d);
    }

   NVD    = (1<<D);     // # function vals, derivs per grid point
   NCoeff = NVD*NVD;    // # polynomial coefficients per grid cell

   D0     = FixedCoordinates.size();
   NVD0   = (1<<D0);

   size_t PhiVDTableSize = (NPoint * NF * NVD0) * sizeof(double);
   double *PhiVDTable = (double *)mallocEC(PhiVDTableSize);
   size_t CTableSize = (NCell * NF * NCoeff) * sizeof(double);
   CTable=(double *)mallocEC(CTableSize);
   if (!PhiVDTable || !CTable)
    ErrExit("%s:%i:out of memory",__FILE__,__LINE__);

   iVec dxMax(D0);
   for(int d0=0; d0<D0; d0++)
    dxMax[d0] = (isinf(FixedCoordinates[d0]) ? 2 : 1);

   /*--------------------------------------------------------------*/
   /*- call user's function to populate table of function values   */
   /*- and derivatives at grid points                              */
   /*--------------------------------------------------------------*/
   int ThreadTaskThreshold = 16; CheckEnv("SCUFF_INTERP_MIN_TASKS",&ThreadTaskThreshold);
   int NumThreads = ( ((int)NPoint) <= ThreadTaskThreshold ? 1 :  GetNumThreads() );
   if (Verbose)
    Log("Evaluating interpolation function at %i points (%i threads)",NPoint,NumThreads);
#ifdef USE_OPENMP
#pragma omp parallel for num_threads(NumThreads)
#endif
   for(size_t nPoint=0; nPoint<NPoint; nPoint++)
    {
      if (Verbose) LogPercent(nPoint,NPoint);
      size_t Offset = GetPhiVDTableOffset(nPoint,0);
      dVec X0Vec  = n2X0(GetPoint(nPoint));
      UserFunc(X0Vec, UserData, PhiVDTable + Offset, dxMax);
    }

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
  {
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
             }
            M->SetEntry(nEq, nCoeff, Entry);
          }
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
       for(int d=0; d<D; d++)
        if (XGrids.size()==0)
         D2Vec[d] = 0.5*DX[d];
        else 
         { size_t n=nVec[d], np1=n+1;
           if (np1 >= XGrids[d].size() )
            { np1 = XGrids[d].size() - 1;
              n   = np1-1;
            }
           D2Vec[d] = 0.5*(XGrids[d][np1] - XGrids[d][n]);
         }

       // populate RHS vector with function values and derivatives
       // at all corners of grid cell
       HVector RHS(NCoeff, LHM_REAL, CTable + GetCTableOffset(nCell,nf));
       LOOP_OVER_IVECS(nTau, tauVec, Twos)
        { 
          double *PhiVD = PhiVDTable + GetPhiVDTableOffset(nVec, tauVec, nf);

          LOOP_OVER_IVECS(nSigma, sigmaVec, Twos)
           RHS.SetEntry(nTau*NVD + nSigma, Monomial(D2Vec, sigmaVec)*PhiVD[nSigma]);
        }

       // operate with inverse M matrix to yield C coefficients
       M->LUSolve(&RHS);
     }

   free(PhiVDTable);
   delete M;
}

/****************************************************************/
/****************************************************************/
/****************************************************************/
bool InterpND::PointInGrid(double *X0Vec, int *nVec, double *XBarVec)
{ 
  for(int d0=0, d=0; d0<D0; d0++)
   if ( !isinf(FixedCoordinates[d0]) )
    {
      if (!EqualFloat(X0Vec[d0], FixedCoordinates[d0])) return false;
    }
   else
    { 
      double *XdGrid = (XGrids.size() > 0 ? &(XGrids[d][0]) : 0);
      double XdMin   = (XGrids.size() > 0 ? 0.0              : XMin[d] );
      double DXd     = (XGrids.size() > 0 ? 0.0              : DX[d] );
      int nd;
      double XdBar;
      if (!FindInterval(X0Vec[d0], XdGrid, NVec[d], XdMin, DXd, &nd, &XdBar))
       return false;
      if (nVec) nVec[d]=nd;
      if (XBarVec) XBarVec[d] = 2.0*(XdBar-0.5);
      d++;
    }
  return true;
}

/****************************************************************/
/****************************************************************/
/****************************************************************/
bool InterpND::Evaluate(double *X0, double *Phi)
{
  /****************************************************************/
  /****************************************************************/
  /****************************************************************/
  iVec nVec(D);
  dVec XBarVec(D);
  if ( !PointInGrid(X0, &(nVec[0]), &(XBarVec[0])) ) return false;
  
  /****************************************************************/
  /* tabulate powers of the scaled/shifted coordinates            */
  /****************************************************************/
  double XBarPowers[MAXDIM][4]; 
  for(int d=0; d<D; d++)
   { XBarPowers[d][0]=1.0;
     XBarPowers[d][1]=XBarVec[d];
     XBarPowers[d][2]=XBarPowers[d][1]*XBarVec[d];
     XBarPowers[d][3]=XBarPowers[d][2]*XBarVec[d];
   }

  /****************************************************************/
  /* construct the NCOEFF-element vector of monomial values that  */
  /* I dot with the NCOEFF-element vector of coefficients to get  */
  /* the value of the interpolating polynomial                    */
  /****************************************************************/
  memset(Phi, 0, NF*sizeof(double));
  iVec Fours(D,4);
  LOOP_OVER_IVECS(nCoeff,pVec,Fours)
   { 
     double Monomial=1.0;
     for(int d=0; d<D; d++)
      Monomial*= XBarPowers[d][pVec[d]];
     for(int nf=0; nf<NF; nf++)
      Phi[nf] += CTable[ GetCTableOffset(nVec,nf) + nCoeff ]*Monomial;
   }
  return true;
}

/****************************************************************/
/****************************************************************/
/****************************************************************/
bool InterpND::EvaluateVD(double *X0, double *PhiVD)
{
  /****************************************************************/
  /****************************************************************/
  /****************************************************************/
  iVec nVec(D);
  dVec XBarVec(D);
  if ( !PointInGrid(X0, &(nVec[0]), &(XBarVec[0])) ) return false;
  
  /****************************************************************/
  /* tabulate powers of the scaled/shifted coordinates            */
  /****************************************************************/
  double XBarPowers[MAXDIM][4]; 
  for(int d=0; d<D; d++)
   { XBarPowers[d][0]=1.0;
     XBarPowers[d][1]=XBarVec[d];
     XBarPowers[d][2]=XBarPowers[d][1]*XBarVec[d];
     XBarPowers[d][3]=XBarPowers[d][2]*XBarVec[d];
   }

  double dXBarPowers[MAXDIM][4];
  for(int d=0; d<D; d++)
   { double LO2 = 0.5*( XGrids[d].size()==0 ? DX[d] : (XGrids[d][nVec[d]+1]-XGrids[d][nVec[d]]));
     dXBarPowers[d][0]=0.0;
     dXBarPowers[d][1]=1.0/LO2;
     dXBarPowers[d][2]=2.0*XBarVec[d]/LO2;
     dXBarPowers[d][3]=3.0*XBarVec[d]*XBarVec[d]/LO2;
   }

  /****************************************************************/
  /* construct the NCOEFF-element vector of monomial values that  */
  /* I dot with the NCOEFF-element vector of coefficients to get  */
  /* the value of the interpolating polynomial                    */
  /****************************************************************/
  memset(PhiVD, 0, NVD*NF*sizeof(double));
  iVec Fours(D,4), Twos(D,2);
  LOOP_OVER_IVECS(nCoeff,pVec,Fours)
   { LOOP_OVER_IVECS(nVD, tauVec, Twos)
      { 
        double Monomial=1.0;
        for(int d=0; d<D; d++)
         Monomial *= (tauVec[d] ? dXBarPowers[d][pVec[d]] : XBarPowers[d][pVec[d]]);
        for(int nf=0; nf<NF; nf++)
         PhiVD[nf*NVD + nVD] += CTable[GetCTableOffset(nVec,nf) + nCoeff]*Monomial;
      }
   }
  return true;
}

/****************************************************************/
/****************************************************************/
/****************************************************************/
bool InterpND::EvaluateVDD(double *X0, double *PhiVD)
{
(void) X0; 
(void) PhiVD;
ErrExit("%s:%i: internal error",__FILE__,__LINE__);
#if 0
  /****************************************************************/
  /****************************************************************/
  /****************************************************************/
  iVec nVec(D);
  dVec XBarVec(D);
  if ( !PointInGrid(X0, &(nVec[0]), &(XBarVec[0])) ) return false;
  
  /****************************************************************/
  /* tabulate powers of the scaled/shifted coordinates            */
  /****************************************************************/
  double XBarPowers[MAXDIM][4]; 
  for(int d=0; d<D; d++)
   { XBarPowers[d][0]=1.0;
     XBarPowers[d][1]=XBarVec[d];
     XBarPowers[d][2]=XBarPowers[d][1]*XBarVec[d];
     XBarPowers[d][3]=XBarPowers[d][2]*XBarVec[d];
   }

  // sizes of grid cells in all dimensions
  dVec LO2(D);
  for(int d=0; d<D; d++)
   LO2[d] = 0.5*( XGrids[d].size()==0 ? DX[d] : (XGrids[d][nVec[d]+1]-XGrids[d][nVec[d]]));

  // number of values and derivatives per function:
  //  1 value, D first derivatives, D*(D+1)/2 second derivatives.
  int NumVDD = 1 + D + D*(D+1)/2;

  // MonomialV
  double MonomialFactors[NumVDD][D][4];
  for(int d=0; d<D; d++)
   for(int p=0; p<4; p++)
    { 
      int nmf=0;
      MonomialFactors[nmf++][d][p] = XBarPowers[d][p];

      for(int dd=0; dd<D; dd++)
       MonomialFactors[nmf++][d][p] = (dd==d ? (p==0 ? 0.0 : XBarPowers[d][p-1]/L02[d]) : XBarPowers[dp][p]);

      for(int dd=0; dd<D; dd++)
       for(int ddd=dd; ddd<D; ddd++)
        MonomialFactors[nmf++][d][p] = (dd=ddd ? (p<=1 ? 0.0 : p*(p-1)*
 
     dXBarPowers[d][0]=0.0;
     dXBarPowers[d][1]=1.0/LO2;
     dXBarPowers[d][2]=2.0*XBarVec[d]/LO2;
     dXBarPowers[d][3]=3.0*XBarVec[d]*XBarVec[d]/LO2;
   }

  /****************************************************************/
  /* construct the NCOEFF-element vector of monomial values that  */
  /* I dot with the NCOEFF-element vector of coefficients to get  */
  /* the value of the interpolating polynomial                    */
  /****************************************************************/
  memset(PhiVD, 0, NumVDD*NF*sizeof(double));
  iVec Fours(D,4);
  LOOP_OVER_IVECS(nCoeff,pVec,Fours)
   { for(int nvdd=0; nvdd<NumVDD; nvdd++)
      { double Monomial=1.0;
        for(int d=0; d<D; d++)
         Monomial *= MonomialFactors[nvdd][d][pVec[d]];
        for(int nf=0; nf<NF; nf++)
         PhiVD[nf*NVD + nVD] += CTable[GetCTableOffset(nVec,nf) + nCoeff]*Monomial;
      }
   }
#endif
  return true;
}

/****************************************************************/
/****************************************************************/
/****************************************************************/
double InterpND::PlotInterpolationError(PhiVDFunc UserFunc, void *UserData, char *OutFileName, bool CentersOnly)
{
  FILE *f=fopen(OutFileName,"w");
  if (!f) return 0.0;

  double *PhiExact  = new double[NVD0*NF];
  double *PhiInterp = new double[NF];

  iVec Threes(D, 3);
  double MaxRelErr=0.0;
  iVec Nm1Vec(D);
  for(size_t d=0; d<NVec.size(); d++) Nm1Vec[d]=NVec[d]-1;
  LOOP_OVER_IVECS(nCell, nVec, Nm1Vec)
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

        dVec X0 = FixedCoordinates;
        for(int d0=0, d=0; d0<D0; d0++)
         if (isinf(X0[d0]))
          { if (XGrids.size() == 0 )
             X0[d0] = XMin[d] + (nVec[d] + 0.4999*tauVec[d])*DX[d]; //EXPLAIN OR FIX ME
            else
             { size_t n = nVec[d], np1=n+1;
               if (np1 >= XGrids[d].size() )
                { n   = XGrids[d].size() - 2;
                  np1 = XGrids[d].size() - 1;
                }
               double Delta = XGrids[d][np1] - XGrids[d][n];
               X0[d0] = XGrids[d][n] + 0.4999*tauVec[d]*Delta; //EXPLAIN OR FIX ME
             }
            d++;
          }

        iVec dXMax(D0,1);
        UserFunc(X0, UserData, PhiExact, dXMax);
        Evaluate(&(X0[0]), PhiInterp);

        for(int nf=0; nf<NF; nf++)
         { double Exact  = PhiExact[nf];
           double Interp = PhiInterp[nf];
           if (fabs(Exact) > 0.0)
            MaxRelErr=fmax(MaxRelErr, fabs(Exact-Interp)/fabs(Exact));
         }

        fprintVec(f,&(X0[0]),D0);
        for(int nf=0; nf<NF; nf++)
         fprintf(f,"%e ",PhiExact[nf]);
        fprintVecCR(f,PhiInterp,NF);
      }
   }
  fclose(f);
  delete[] PhiExact;
  delete[] PhiInterp;
  return MaxRelErr;
}

/****************************************************************/
/****************************************************************/
/****************************************************************/
double GetInterpolationError(PhiVDFunc UserFunc, void *UserData, int NF,
                             dVec X0Vec, dVec DeltaVec,
                             double AbsTol, char *LogFileName)
{ 
  int D0 = X0Vec.size();
  dVec X0Max = X0Vec;
  iVec N0Vec(D0,1);
  for(int d0=0; d0<D0; d0++) 
   if (DeltaVec[d0]!=0.0)
    { N0Vec[d0]=2;
      X0Max[d0]+=2.0*DeltaVec[d0];
    }

  InterpND Interp(UserFunc, UserData, NF, X0Vec, X0Max, N0Vec);

  int NVD = (1<<D0);     // # function vals, derivs per grid point

  double *PhiExact  = new double[NVD*NF];
  double *PhiInterp = new double[NF];

  FILE *LogFile = LogFileName ? fopen(LogFileName,"a") : 0;

  iVec TauMax(D0);
  for(int d0=0; d0<D0; d0++)
   TauMax[d0] = (N0Vec[d0]==1 ? 1 : 3);
  double MaxRelError=0.0;
  LOOP_OVER_IVECS(nTau, tauVec, TauMax)
   {
     dVec XSample=X0Vec;
     for(int d0=0; d0<D0; d0++)
      XSample[d0] += 0.4999*tauVec[d0]*DeltaVec[d0]; // EXPLAIN OR FIX ME

     iVec dXMax(D0,1);
     UserFunc(XSample, UserData, PhiExact, dXMax);
     bool Status=Interp.Evaluate(&(XSample[0]), PhiInterp);
     if (!Status) ErrExit("%s:%i: internal error",__FILE__,__LINE__);
   
     double ThisMaxRelError=0.0;

     for(int nf=0; nf<NF; nf++)  
      { double AbsError = fabs(PhiExact[nf] - PhiInterp[nf]);
        if (AbsError<AbsTol) continue;
        double RelError = AbsError / ( fabs(PhiExact[nf])==0.0 ? 1.0 : fabs(PhiExact[nf]));
        ThisMaxRelError = fmax(ThisMaxRelError, RelError);
      }

     if(LogFile) 
      { fprintVec(LogFile,&(XSample[0]),D0);
        for(int nf=0; nf<NF; nf++)
         fprintf(LogFile,"%e %e ",fabs(PhiExact[nf]),fabs(PhiExact[nf]-PhiInterp[nf]));
        fprintf(LogFile,"%e\n",ThisMaxRelError);
      }

     MaxRelError=fmax(MaxRelError,ThisMaxRelError);
   }
  if(LogFile) { fprintf(LogFile,"\n\n"); fclose(LogFile); }

  delete[] PhiExact;
  delete[] PhiInterp;
  return MaxRelError;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
double GetInterpolationError(PhiVDFunc UserFunc, void *UserData, int NF,
                             int d0, double X, double Delta, dVec X0Min, dVec X0Max, 
                             double AbsTol, char *LogFileName)
{
  int D0=X0Min.size();
  dVec X0Vec(D0), DeltaVec(D0,0.0);
  X0Vec[d0]    = X;
  DeltaVec[d0] = Delta;
  
  int PointsPerDimension=5;
  iVec N0Vec(D0, PointsPerDimension);
  N0Vec[d0]=1;
  for(int d0Prime=0; d0Prime<D0; d0Prime++)
   if (EqualFloat(X0Min[d0Prime],X0Max[d0Prime]))
    N0Vec[d0Prime]=1;
  
  double MaxError=0.0;
  LOOP_OVER_IVECS(n, n0Vec, N0Vec)
   { 
     for(int d0Prime=0; d0Prime<D0; d0Prime++)
      if (d0Prime!=d0)
       X0Vec[d0Prime] 
       = X0Min[d0Prime] + n0Vec[d0Prime]*(X0Max[d0Prime]-X0Min[d0Prime])/(PointsPerDimension-1);

     MaxError=fmax(MaxError,GetInterpolationError(UserFunc, UserData, NF, X0Vec, DeltaVec, AbsTol, LogFileName));
   }
  return MaxError;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
#if 0
dVec GetXdGrid(PhiVDFunc UserFunc, void *UserData, int NF, dVec X0, int d,
               double XdMin, double XdMax, double DesiredMaxRE)
{ 
  bool ExtraVerbose = CheckEnv("SCUFF_INTERPOLATION_EXTRA_VERBOSE");
  Log("Autotuning interpolation grid for coordinate %i, range {%e,%e}",d,XdMin,XdMax);
  if (X0.size() != 1)
   { Log(" X0 = { ");
     for(size_t d0=0; d0<X0.size(); d0++)
      { if ( ((int)d0) == d ) 
         LogC("xx");
        else 
         LogC("%e",X0[d0]);
        LogC("%c",d0==X0.size()-1 ? '}' : ',');
      }
   }
  
  double DeltaMin = (XdMax - XdMin) / (1000.0);
  double DeltaMax = (XdMax - XdMin) / 2.0;
  double Delta    = (XdMax - XdMin) / (10.0); // initial guess
 
  double Xd = XdMin;
  dVec XdGrid(1,Xd);
  while( Xd<XdMax )
   { 
     dVec X0Vec = X0;  X0Vec[d] = Xd;
     bool TooBig=false, JustRight=false;
     while(!JustRight)
      { dVec DeltaVec(X0.size(), 0.0); DeltaVec[d]=Delta;
        double Err=GetInterpolationError(UserFunc, UserData, NF, X0Vec, DeltaVec);
        if (ExtraVerbose) Log("   x%i=%+f, Delta=%f: MRE %.2e",d,Xd,Delta,Err);
        if ( Err > 10.0*DesiredMaxRE )
         { TooBig=true; Delta *= 0.2; }
        else if ( Err > DesiredMaxRE )
         { TooBig=true; Delta *= 0.75; }
        else if ( Err<0.1*DesiredMaxRE )
         { if (TooBig) 
            JustRight=true; // prevent oscillatory behavior
           else
            Delta *= 2.0; 
         }
        else
         JustRight = true;

        if (Delta<DeltaMin || Delta>DeltaMax)
         { Delta = (Delta<DeltaMin ? DeltaMin : DeltaMax);
           if (ExtraVerbose) Log("   setting Delta=%f",Delta);
           JustRight = true;
         }
      }
     Xd+=Delta;
     if (Xd>XdMax) Xd=XdMax;
     XdGrid.push_back(Xd);
   }
  Log(" ...%i grid points",XdGrid.size()); 
  return XdGrid;
}
#endif

/***************************************************************/
/***************************************************************/
/***************************************************************/
dVec GetXGrid(PhiVDFunc UserFunc, void *UserData, int NF, dVec X0Min, dVec X0Max, int d0,
              double DesiredMaxRE)
{ 
  bool ExtraVerbose = CheckEnv("SCUFF_INTERPOLATION_EXTRA_VERBOSE");
  Log("Autotuning interpolation grid for coordinate %i, range {%e,%e}",d0,X0Min[d0],X0Max[d0]);

  int NMin=2;    CheckEnv("SCUFF_INTERPOLATION_NMIN",&NMin);
  int NMax=100;  CheckEnv("SCUFF_INTERPOLATION_NMAX",&NMax);
  
  double DeltaMin = (X0Max[d0] - X0Min[d0]) / NMax;
  double DeltaMax = (X0Max[d0] - X0Min[d0]) / NMin;
  double Delta    = (X0Max[d0] - X0Min[d0]) / 10;   // initial guess
 
  double X = X0Min[d0];
  dVec XGrid(1,X);
  while( X<X0Max[d0] )
   { 
     bool TooBig=false, JustRight=false;
     while(!JustRight)
      { double Err=GetInterpolationError(UserFunc, UserData, NF, d0, X, Delta, X0Min, X0Max);
        if (ExtraVerbose) Log("   x%i=%+f, Delta=%f: MRE %.2e",d0,X,Delta,Err);
        if ( Err > 10.0*DesiredMaxRE )
         { TooBig=true; Delta *= 0.2; }
        else if ( Err > DesiredMaxRE )
         { TooBig=true; Delta *= 0.75; }
        else if ( Err<0.1*DesiredMaxRE )
         { if (TooBig) 
            JustRight=true; // prevent oscillatory behavior
           else
            Delta *= 2.0; 
         }
        else
         JustRight = true;

        if (Delta<DeltaMin || Delta>DeltaMax)
         { Delta = (Delta<DeltaMin ? DeltaMin : DeltaMax);
           if (ExtraVerbose) Log("   setting Delta=%f",Delta);
           JustRight = true;
         }
      }
     X+=Delta;
     if (X>X0Max[d0]) X=X0Max[d0];
     XGrid.push_back(X);
   }
  Log(" ...%i grid points",XGrid.size()); 
  return XGrid;
}
