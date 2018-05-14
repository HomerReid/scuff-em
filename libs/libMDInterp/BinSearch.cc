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
 * BinSearch -- binary search to find the proper slot for a number in 
 *           -- a sorted but possibly nonuniformly spaced array
 */

#include <stdlib.h>
#include <math.h>
#include <libhrutil.h> // needed for EqualFloat()

/***************************************************************/
/* given an array of N points, which is ordered (i.e. for all  */
/* i we have XPoints[i]<XPoints[i+1]) but not necessarily      */
/* uniformly spaced, return n such that                        */
/*  XPoints[n] <= X < XPoints[n+1].                            */
/* note: if X < XPoints[0] then we return 0, but the caller has*/ 
/* no way to distinguish this return value from the case       */ 
/* XPoints[0] <= X < XPoints[1], and so the caller should      */ 
/* test separately for that case.                              */ 
/***************************************************************/
int BinSearch(double X, double *XPoints, int N)
{
   if ( X<XPoints[0] ) return 0;
   if ( X>XPoints[N-1] ) return N-1;

   int nMin=0; 
   int nMax=N-1;

   for(;;)
    { 
      int nHalf=(nMin + nMax)/2;

      if ( X>=XPoints[nHalf] )
       { if ( X<XPoints[nHalf+1] ) 
          return nHalf;
         nMin=nHalf;
       }
      else
       nMax=nHalf;
    }
}

/***************************************************************/
/* 1D cell-finding routine: given an array of intervals, find  */
/* the interval in which the given point lies.                 */
/*                                                             */
/* if XPoints is non-NULL, then this routine computes n such   */
/* that XPoints[n] < X < XPoints[n+1]. (if X<XPoints[0] then   */
/* we set n=0, and if X>XPoints[N-1] then we set n=N-1.)       */ 
/*                                                             */
/* if XPoints is NULL then this routine computes n such that   */
/* XMin + n*DX < X < XMin + (n+1)*DX.  (again, if X<XMin then  */
/* we put n=0, and if X>XMin + (N-1)*DX then we set n=(N-1).   */
/*                                                             */
/* in either case, XBar is a scaled version of X within the    */
/* interval defined such that XBar=0(1) if X lies at the left  */
/* (right) endpoint of the interval.                           */
/*                                                             */
/* 20121129 the test case for the NOTE below is:               */
/*          X=0.75                                             */
/*          XMin = -0.72260686734410873                        */
/*          N = 30                                             */
/* in this case what happens is:                               */
/*  XMax = XMin + ((double)(N-1)*DX) = 0.75000000000000011     */
/* and hence X < XMax, so we fall through from the second      */
/* clause of the if..else to the third clause; but then        */
/* trunc( (X-XMin)/DX ) = 29 and XBar gets set to 0, i.e.      */
/* we wind up at the far left edge of the (nonexistent) Nth    */
/* bin, whereas where we want to be is at the far right edge   */
/* of the (N-1)th bin.                                         */
/*                                                             */
/* 20140812 the return value is true if the point lies in      */
/* the interval, false otherwise.                              */
/***************************************************************/
bool FindInterval(double X, double *XPoints, int N, double XMin, double DX,
                  int *pn, double *pXBar)
{
  int n;
  double XBar;
  bool Inside=true; // true if point lies within interval

  if ( XPoints )  // nonuniform grid
   { if ( X<=XPoints[0] )
      { n=0;
        XBar=0.0;
        Inside=EqualFloat(X,XPoints[0]);
        X=XPoints[0]; 
      }
    else if ( X>=XPoints[N-1] )
      { n=N-2;
        XBar=1.0;
        Inside = EqualFloat(X,XPoints[N-1]);
        X=XPoints[N-1]; 
      }
    else
      { n=BinSearch(X,XPoints,N);
        XBar = (X - XPoints[n]) / (XPoints[n+1]-XPoints[n]);
      }
   }
  else // uniform grid
   { if ( X <= XMin )
      { n=0;
        XBar=0.0;
        Inside = EqualFloat(X,XMin);
        X=XMin;
      }
     else if ( X >= (XMin + ((double)(N-1))*DX) )
      { n=N-2;
        XBar=1.0;
        Inside = EqualFloat(X, (XMin + ((double)N-1))*DX);
        X=XMin + ((double)(N-1))*DX;
      }
     else
      { 
        XBar=(X-XMin)/DX;
        n=(int)trunc(XBar);
        XBar -= (double)n;

        // see NOTE above
        if ( n==(N-1) )
         { n=N-2;
           XBar=1.0;
         }
      }
   }

  *pn=n;
  *pXBar=XBar;
  return Inside;

}
