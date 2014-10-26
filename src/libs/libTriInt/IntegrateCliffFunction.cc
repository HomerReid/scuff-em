#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#include <libhrutil.h>
#include <libTriInt.h>
#include <libSGJC.h>

/***************************************************************/
/* Estimate the integral of fCliff over interval [x0,x0+Delta] */
/* using (a) the composite trapezoidal rule with two           */
/* subintervals of width Delta/2, and (b) the (non-composite)  */
/* Simpson's rule with a single interval of width Delta.       */
/*                                                             */
/* Estimate (b) is taken as our estimate of the integral, and  */
/* the difference between the estimates is taken as a rough    */
/* estimate of the error.                                      */
/*                                                             */
/* Each of Buffer[0..2] must point to caller-allocated buffers */
/*  with room to store nFun doubles. The user's integrand      */
/*  function will be called to fill these buffers with the     */
/*  values of f at the following three points:                 */
/*   Buffer[0] = f[x0 + 0.0*Delta]                             */
/*   Buffer[1] = f[x0 + 0.5*Delta]                             */
/*   Buffer[2] = f[x0 + 1.0*Delta]                             */
/*                                                             */
/* If BufferFull[i] is true (i=0,1,2) then Buffer[i] is assumed*/
/*  to contain already the appropriate values of f and the call*/
/*  to the user's function is skipped for that point.          */
/***************************************************************/
int TrapSimp(double x0, double Delta,
             CliffFunction fCliff, void *UserData, int nFun, bool *Skip,
             double *fBuffer[3], bool BufferFull[3],
             double *Result, double *Error, double *AbsResult,
             double *pMaxAbsError, double *pMinAbsError,
             double *pMaxRelError, double *pMinRelError)
{
  double *fLeft  = fBuffer[0];
  double *fMid   = fBuffer[1];
  double *fRight = fBuffer[2];
  int fEvals=0;

  if (!BufferFull[0])
   { double x=x0;
     fCliff(1, &x, UserData, nFun, Skip, fLeft);
     fEvals++;
   };

  if (!BufferFull[1])
   { double x=x0+0.5*Delta;
     fCliff(1, &x, UserData, nFun, Skip, fMid);
     fEvals++;
   };

  if (!BufferFull[2])
   { double x=x0+1.0*Delta;
     fCliff(1, &x, UserData, nFun, Skip, fRight);
     fEvals++;
   };

  double MaxAbsError=0.0, MinAbsError=1.23e45;
  double MaxRelError=0.0, MinRelError=1.23e45;
  for(int nf=0; nf<nFun; nf++)
   { 
     if (Skip[nf])
      continue;

     double Trap = Delta*(fLeft[nf] + 2.0*fMid[nf] + fRight[nf]) / 4.0;
     double Simp = Delta*(fLeft[nf] + 4.0*fMid[nf] + fRight[nf]) / 6.0;
     Result[nf] = Simp;
     Error[nf]  = fabs(Simp-Trap);

     MaxAbsError = fmax(MaxAbsError, Error[nf]);
     MinAbsError = fmin(MinAbsError, Error[nf]);

     double RelError = Error[nf] / fabs(Result[nf]);
     MaxRelError = fmax(MaxRelError, RelError);
     MinRelError = fmin(MinRelError, RelError);

     AbsResult[nf] = Delta*( fabs(fLeft[nf]) + 4.0*fabs(fMid[nf]) + fabs(fRight[nf])) / 6.0;

   };
  if (MinAbsError==1.23e45) MinAbsError=0.0;
  if (MinRelError==1.23e45) MinRelError=0.0;

  if (pMaxAbsError) *pMaxAbsError = MaxAbsError;
  if (pMinAbsError) *pMinAbsError = MinAbsError;
  if (pMaxRelError) *pMaxRelError = MaxRelError;
  if (pMinRelError) *pMinRelError = MinRelError;
  return fEvals;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int IntegrateCliffFunction(CliffFunction fCliff, void *UserData, int nFun,
                           double xMin, double xMax, double xCliff,
                           double AbsTol, double RelTol, 
                           double *Integral, double *Error,
                           char *LogFileName)
{
  if (LogFileName)
   { FILE *f=fopen(LogFileName,"a");
     fprintf(f,"\n## Cliff integrator running on %s\n",GetHostName());
     fprintf(f,"## Columns: \n");
     fprintf(f,"## 1 integrand index\n");
     fprintf(f,"## 2 x \n");
     fprintf(f,"## 3 integral so far \n");
     fprintf(f,"## 4 error so far \n");
     fprintf(f,"## 5 contribution of this subinterval\n");
     fprintf(f,"## 6 error in this subinterval\n");
     fprintf(f,"## 7 integral so far of abs(f) \n");
     fprintf(f,"## 8 relative slope\n");
     fclose(f);
   };

  /***************************************************************/
  /* sanity check on input parameter values **********************/
  /***************************************************************/
  // if xCliff<=xMin we assume the user wants us to guess xCliff automatically.
  if (xCliff<=xMin)
   { if (xMax > xMin)
      xCliff = 0.5*(xMax-xMin);
     else if (xMin>0.0)
      xCliff = 1.5*xMin;
     else // interval is [0,infinity] so we have no way to guess a length scale
      xCliff = 1.0;
   };

  // if xMax<=xMin we assume the user intended xMax=infinity.
  if (xMax<=xMin)
   xMax=1.234e89;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  bool *Converged = new bool[nFun];
  memset(Converged, 0, nFun*sizeof(bool));

  int *ConvergedIters = new int[nFun];
  memset(ConvergedIters, 0, nFun*sizeof(int));

  double *Workspace = new double[7*nFun];

  double *fBuffer[3];
  double *fLeft  = fBuffer[0] = Workspace + 0*nFun;
  double *fMid   = fBuffer[1] = Workspace + 1*nFun;
  double *fRight = fBuffer[2] = Workspace + 2*nFun;
  bool BufferFull[3] = {false, false, false};

  /***************************************************************/
  /* We stop integrating when the contributions of the most      */
  /* recent MinConvergedIters iterations have failed to change   */
  /* the integral value by more than the tolerances. We don't    */
  /* begin to assess this criterion until we have accumulated    */
  /* the contributions of at least the first MinIters intervals. */
  /*                                                             */
  /* The value of MinConvergedIters is initially set to 3, but   */
  /* is increased later if the integrand is observed to be       */
  /* non-monotonic.                                              */
  /***************************************************************/
  int MinConvergedIters=2;
  int MinIters=3;

  /***************************************************************/
  /*- Starting at xMin, we inch along towards xMax in steps of   */
  /*- width Delta. We estimate the integral over each interval   */
  /*- and the error in our estimate. If the error estimate is too*/
  /*- large, we shrink Delta and try again. Otherwise we add the */
  /*- contributions of this interval to the running tally of the */
  /*- full integral value and move on, possibly enlarging Delta  */
  /*- if the error estimate was ``too small.''                   */
  /***************************************************************/
  double *dIntegral    = Workspace + 3*nFun;
  double *dError       = Workspace + 4*nFun;
  double *dAbsIntegral = Workspace + 5*nFun;
  double *AbsIntegral  = Workspace + 6*nFun;
  double Delta = 0.25*xCliff;
  int Iters=0;
  int fEvals=0;
  memset(Integral,    0, nFun*sizeof(double));
  memset(Error,       0, nFun*sizeof(double));
  memset(AbsIntegral, 0, nFun*sizeof(double));
  double x=xMin;
  bool AllConverged=false;
  while( (AllConverged==false) && (x<xMax) )
   { 
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     if (LogFileName)
      { FILE *f=fopen(LogFileName,"a");
        time_t MyTime=time(0);
        struct tm *MyTm=localtime(&MyTime);
        char TimeStr[30];
        strftime(TimeStr,30,"%D::%T",MyTm);
        fprintf(f,"# %s: iter %i: x=%e, Delta=%e\n",TimeStr,Iters+1,x,Delta);
        fclose(f);
      };

     /*--------------------------------------------------------------*/
     /* estimate the integral over the subinterval [x, x+Delta]      */
     /*--------------------------------------------------------------*/
     if (x+Delta > xMax) // avoid overshoot at right endpoint 
      Delta=(xMax-x);
     double MaxAbsError, MinAbsError;
     double MaxRelError, MinRelError;
     fEvals += TrapSimp(x, Delta, fCliff, UserData, nFun, Converged,
                        fBuffer, BufferFull,
                        dIntegral, dError, dAbsIntegral,
                        &MaxAbsError, &MinAbsError,
                        &MaxRelError, &MinRelError);

     /*--------------------------------------------------------------*/
     /* if the errors on this subinterval were too large,            */
     /* shrink the width and try again.                              */
     /*--------------------------------------------------------------*/
     if (    (Delta>1.0e-8)
          && ( (MaxRelError > 2.0*RelTol) && (MaxAbsError > AbsTol) ) 
        )
      { 
        Delta*=0.5;
        memcpy(fRight,fMid,nFun*sizeof(double));
        BufferFull[0]=BufferFull[2]=true; BufferFull[1]=false;
        continue;
      }

     /*--------------------------------------------------------------*/
     /* accumulate contributions to \int f dx and  \int |f| dx       */
     /*--------------------------------------------------------------*/
     x+=Delta;
     AllConverged=true;
     for(int nf=0; nf<nFun; nf++)
      { 
        if (Converged[nf]) 
         continue;

        Integral[nf]    += dIntegral[nf];
        Error[nf]       += dError[nf];
        AbsIntegral[nf] += dAbsIntegral[nf];

        // compute the 'relative slope' of this integrand component.
        // if it is small and has been small for at least MinConvergedIters,
        // then mark this integrand component as having converged. 
        if ( Iters<MinIters )
         AllConverged=false;
        else
         { 
           double FullSlope   = AbsIntegral[nf] / (x-xMin);
           double LocalSlope = dAbsIntegral[nf] / Delta;

           if ( LocalSlope <= RelTol*FullSlope )
            ConvergedIters[nf]++;
           else
            ConvergedIters[nf]=0;

           if (ConvergedIters[nf]==MinConvergedIters)
            Converged[nf] = true;
           else  
            AllConverged = false;
        };
      };

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     if (LogFileName)
      { FILE *f=fopen(LogFileName,"a");
        for(int nf=0; nf<nFun; nf++)
         { 
            double FullSlope  = (Iters==0) ? 1.0 : AbsIntegral[nf] / (x-xMin);
            double LocalSlope = dAbsIntegral[nf] / Delta;
            if (Converged[nf])
             fprintf(f,"  %2i %.4e  %+.8e %.1e  (converged)\n",
                          nf, x, Integral[nf], Error[nf]);
            else
             fprintf(f,"  %2i %.4e  %+.4e %.1e   %+.4e  %+.1e %.4e  %.1e\n",
                         nf, x, 
                         Integral[nf], Error[nf], dIntegral[nf], dError[nf],
                         AbsIntegral[nf], LocalSlope/FullSlope);
         };
        fclose(f);
      };
    
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     Iters++;
     memcpy(fLeft,fRight,nFun*sizeof(double));
     BufferFull[0]=true; BufferFull[1]=BufferFull[2]=false;

     // if the errors on this subinterval were ``too small,''
     // enlarge the width but keep the integral estimates
     if (MaxRelError < RelTol)
      Delta*=2.0;

   };

  delete[] Workspace;
  delete[] ConvergedIters;
  delete[] Converged;
  return fEvals;
  
}
