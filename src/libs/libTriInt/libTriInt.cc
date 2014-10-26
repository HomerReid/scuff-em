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
 * libTriInt.cc -- implementation of libTriInt 
 *
 * homer reid   -- 12/2008
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>

#include "libTriInt.h"
#include "libhrutil.h"

extern FILE *LogFile;

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/* 1. data structure definitions.                               */
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
typedef struct Triangle
 { 
   double *V1, *V2, *V3;
   double *I, *E;
   double MaxError;
 } Triangle;

typedef struct TIWorkspace
 { 
   int nFun;
   int MaxIters, MaxNodes, MaxTriangles;
   int NumNodes, NumTriangles;
   double *Nodes;
   Triangle **Triangles;
   Triangle *TBuffer;
   double *IBuffer, *EBuffer;
   char *LogFileName;

} TIWorkspace;

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/* 2A. user-callable functions that manipulate TIWorkspaces.    */
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/***************************************************************/
/* allocate and return an opaque pointer to a new TIWorkspace  */
/* structure.                                                  */
/***************************************************************/
void *CreateTIWorkspace(int nFun, int MaxIters)
{ 
  int nt;
  TIWorkspace *TIW;
  Triangle *T;

  TIW=(TIWorkspace *)malloc(sizeof *TIW);

  TIW->nFun=nFun;
  TIW->MaxIters=MaxIters;

  /* each iteration adds 1 new vertex and 2 new triangles */ 
  TIW->MaxNodes=4 + MaxIters;
  TIW->MaxTriangles=2 + 2*MaxIters;

  TIW->Nodes=(double *)malloc(3*TIW->MaxNodes*sizeof(double));
  TIW->TBuffer=(Triangle *)malloc(TIW->MaxTriangles*sizeof(Triangle));
  TIW->IBuffer=(double *)malloc(TIW->MaxTriangles*nFun*sizeof(Triangle));
  TIW->EBuffer=(double *)malloc(TIW->MaxTriangles*nFun*sizeof(Triangle));
  TIW->Triangles=(Triangle **)malloc(TIW->MaxTriangles*sizeof(Triangle *));
  for(nt=0; nt<TIW->MaxTriangles; nt++)
   { T=TIW->Triangles[nt]=TIW->TBuffer + nt;
     T->I=TIW->IBuffer + nt*nFun;
     T->E=TIW->EBuffer + nt*nFun;
   };

  TIW->LogFileName=strdupEC("libTriInt.log");
  return TIW;

} 

/***************************************************************/
/* deallocate a previously allocated TIWorkspace structure. ****/
/***************************************************************/
void DeleteTIWorkspace(void *pTIW)
{ 
  if (!pTIW) return;
 
  TIWorkspace *TIW=(TIWorkspace *)pTIW;
  free(TIW->Nodes);
  free(TIW->IBuffer);
  free(TIW->EBuffer);
  free(TIW->TBuffer);
  free(TIW->Triangles);
  free(TIW->LogFileName);
  free(TIW);
}

/***************************************************************/
/* set the name of the log file. *******************************/
/***************************************************************/
void SetTriIntLogFileName(void *pTIW, const char *format, ...)
{ 
  TIWorkspace *TIW=(TIWorkspace *)pTIW;
  char buffer[200];
  va_list ap;
  
  va_start(ap,format);
  vsnprintfEC(buffer,200,format,ap);
  va_end(ap);

  free(TIW->LogFileName);

  TIW->LogFileName=strdupEC(buffer);
}

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/* 2B. non-user-callable functions that manipulate TIWorkspaces.*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/***************************************************************/
/***************************************************************/
/***************************************************************/
static void ResetTIWorkspace(TIWorkspace *TIW)
{ TIW->NumNodes=TIW->NumTriangles=0;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
static double *AddNode(TIWorkspace *TIW, double *V)
{ 
  double *Slot;

  Slot=TIW->Nodes + 3*TIW->NumNodes;

  memcpy(Slot, V, 3*sizeof(double));
  TIW->NumNodes++;
  return Slot;
  
} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
static void AddTriangle(TIWorkspace *TIW, double *V1, double *V2, double *V3, 
                        double *I, double *E)
{ 
  Triangle *T;
  int nf=0;

  T=TIW->Triangles[TIW->NumTriangles++];
  T->V1=V1;
  T->V2=V2;
  T->V3=V3;  
  memcpy(T->I,I,TIW->nFun*sizeof(double));
  memcpy(T->E,E,TIW->nFun*sizeof(double));

  T->MaxError=T->E[0];
  for(nf=1; nf<TIW->nFun; nf++)
   if( T->E[nf] > T->MaxError ) 
    T->MaxError=T->E[nf];
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
static void DeleteFirstTriangle(TIWorkspace *TIW)
{ 
  Triangle *TEnd;

  if (TIW->NumTriangles==0) 
   return;
  
  TEnd=TIW->Triangles[ TIW->NumTriangles-1 ];
  TIW->Triangles[ TIW->NumTriangles-1 ] = TIW->Triangles[0]; 
  TIW->Triangles[ 0 ] = TEnd;
  TIW->NumTriangles--;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
static int tcmp(const void *pT1, const void *pT2)
{ 
  Triangle *T1=*(Triangle **)pT1;
  Triangle *T2=*(Triangle **)pT2;

  if (T1->MaxError > T2->MaxError) 
   return -1;
  return +1;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
static void SortTriangles(TIWorkspace *TIW)
{ 
  qsort(TIW->Triangles, TIW->NumTriangles, sizeof(Triangle *), tcmp);
}
                 
/***************************************************************/
/***************************************************************/
/***************************************************************/
static void Log(TIWorkspace *TIW, const char *format, ...)
{
 char buffer[200], TimeString[30];
 struct tm *MyTm;
 time_t MyTime;
 va_list ap;
 FILE *f;

 if ( !(f=fopen(TIW->LogFileName,"a")) )
  return;

 va_start(ap,format);
 vsnprintfEC(buffer,200,format,ap);
 va_end(ap);

 MyTime=time(0);
 MyTm=localtime(&MyTime);
 strftime(TimeString,30,"%D::%T",MyTm);

 if(buffer[strlen(buffer)-1]=='\n')
  buffer[strlen(buffer)-1]=0;
 fprintf(f,"\n%s %s",TimeString,buffer);
 fclose(f);

}

/* somewhat more portable replacement for gcc-style variable-length arrays */
#ifdef __GNUC__
#  define DALLOCA(n) (double *) __builtin_alloca(sizeof(double) * (n))
#else
#  include <alloca.h>
#  define DALLOCA(n) (double *) alloca(sizeof(double) * (n))
#endif

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/* 3. Apply fixed-order quadrature rules to estimate the        */
/*    the integral and its error over the specified triangle.   */
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
void TriIntFixed(TriIntFun F, int nFun, void *UserData, 
                 double *V1, double *V2, double *V3, 
                 int Order, double *I)
{ 
  int i, nf, nqr;
  double u, v, Weight, C, J;
  double X[3], A[3], B[3];
  double *dI = DALLOCA(nFun);
  int NumPts;
  double *TCR;

  /* get the triangle cubature rule for the specified order */
  TCR=GetTCR(Order, &NumPts); 

  for(i=0; i<3; i++)
   { A[i]=V2[i]-V1[i];
     B[i]=V3[i]-V1[i];
   };

  /* jacobian = 2*area = norm(A\cross B) */
  C=A[1]*B[2]-A[2]*B[1]; J=C*C;
  C=A[2]*B[0]-A[0]*B[2]; J+=C*C;
  C=A[0]*B[1]-A[1]*B[0]; J+=C*C;
  J=sqrt(J);

  memset(I,0,nFun*sizeof(double));
  for(nqr=0; nqr<3*NumPts; )
   { 
     u=TCR[nqr++];
     v=TCR[nqr++];
     Weight=J*TCR[nqr++];

     for(i=0; i<3; i++)
      X[i] = V1[i] + u*A[i] + v*B[i];

     F(X,UserData,dI);

     for(nf=0; nf<nFun; nf++)
      I[nf]+=Weight*dI[nf];
   };

}

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
void TriIntEmbedded(TriIntFun F, int nFun, void *UserData, 
                    double *V1, double *V2, double *V3, 
                    double *I, double *E)
{ 
  int i, nf, ncp;
  double u, v, Weight, C, J;
  double X[3], A[3], B[3];
  double *TCR5, *TCR7;
  int NumPts5, NumPts7;
  double *dI = DALLOCA(nFun), *I5 = DALLOCA(nFun);

  TCR5=GetTCR(5, &NumPts5); 
  TCR7=GetTCR(7, &NumPts7); 

  for(i=0; i<3; i++)
   { A[i]=V2[i]-V1[i];
     B[i]=V3[i]-V1[i];
   };

  /* jacobian = 2*area = norm(A\cross B) */
  C=A[1]*B[2]-A[2]*B[1]; J=C*C;
  C=A[2]*B[0]-A[0]*B[2]; J+=C*C;
  C=A[0]*B[1]-A[1]*B[0]; J+=C*C;
  J=sqrt(J);

  memset(I5,0,nFun*sizeof(double));
  memset(I,0,nFun*sizeof(double));
  for(ncp=0; ncp<NumPts7; ncp++)
   { 
     u=TCR7[3*ncp+0]; v=TCR7[3*ncp+1]; Weight=J*TCR7[3*ncp+2];

     for(i=0; i<3; i++)
      X[i] = V1[i] + u*A[i] + v*B[i];

     F(X,UserData,dI);

     for(nf=0; nf<nFun; nf++)
      I[nf]+=Weight*dI[nf];

     if (ncp < NumPts5)
      for(nf=0; nf<nFun; nf++)
       I5[nf]+=J*TCR5[3*ncp+2]*dI[nf];

   };

 for(nf=0; nf<nFun; nf++)
  E[nf]=fabs(I[nf]-I5[nf]);

}

/***************************************************************/
/* get integral and error estimate over a triangle.     ********/
/***************************************************************/
void GetIandE(TriIntFun F, int nFun, void *UserData,
              double *V1, double *V2, double *V3, double *I, double *E)
{ 
  int nf;
  double *I1 = DALLOCA(nFun);

  TriIntFixed(F, nFun, UserData, V1, V2, V3, 7, I1);
  TriIntFixed(F, nFun, UserData, V1, V2, V3, 13, I);
  for(nf=0; nf<nFun; nf++)
   E[nf]=fabs(I[nf]-I1[nf]);
}

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- 4. user-callable entry point for TriIntAdaptive routine.  --*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
int TriIntAdaptive(void *pTIW, TriIntFun F, void *UserData, 
                   double *V1, double *V2, double *V3,
                   double AbsTol, double RelTol, double *I, double *E)
{
  TIWorkspace *TIW=(TIWorkspace *)pTIW; 
  int nFun=TIW->nFun;
  int i, nf, NumIters;
  double MaxAbsErr, MaxRelErr, RelErr;
  double *IA = DALLOCA(nFun), *IB = DALLOCA(nFun), *IC = DALLOCA(nFun);
  double *EA = DALLOCA(nFun), *EB = DALLOCA(nFun), *EC = DALLOCA(nFun);
  double *V1t, *V2t, *V3t;
  Triangle *T;
  double Centroid[3], *VC;
  
  /***************************************************************/
  /* begin by zeroing out the TIW structure.   *******************/
  /***************************************************************/
  ResetTIWorkspace(TIW);

  /***************************************************************/
  /* get integral and error estimate for entire triangle and add */
  /* it to list of triangles                                    */
  /***************************************************************/
  GetIandE(F, nFun, UserData, V1, V2, V3, I, E);
  V1t=AddNode(TIW, V1);
  V2t=AddNode(TIW, V2);
  V3t=AddNode(TIW, V3);
  AddTriangle(TIW, V1t, V2t, V3t, I, E);

  /***************************************************************/
  /* main loop. on each iteration, we trisect the triangle with  */
  /* the largest error estimate, evaluate the integral and the   */
  /* error estimate for each subtriangle, and add the three new  */  
  /* subtriangles to our full list of triangles. (we also        */  
  /* remove the triangle we trisected).                          */  
  /***************************************************************/
  for(NumIters=0; NumIters<TIW->MaxIters; NumIters++)
   { 
    /*--------------------------------------------------------------*/
    /*- Convergence analysis ---------------------------------------*/
    /*--------------------------------------------------------------*/
     MaxAbsErr=0.0;
     MaxRelErr=0.0;
     Log(TIW,"Iter %i: %i nodes, %i triangles.",
              NumIters,TIW->NumNodes,TIW->NumTriangles);
     for(nf=0; nf<nFun; nf++)
      { if (E[nf] > MaxAbsErr) 
         MaxAbsErr=E[nf]; 
        if (E[nf] > AbsTol) 
         { RelErr=E[nf]/fabs(I[nf]);
           if (RelErr > MaxRelErr) 
            MaxRelErr=RelErr;
         };
        Log(TIW," %2i: (%+12.9e,%7.3e,%7.3e) (%7.3e,%7.3e)",
                  nf, I[nf],E[nf],E[nf]/fabs(I[nf]),AbsTol,RelTol);
      };
    
     if ( MaxAbsErr<AbsTol && MaxRelErr < RelTol )
      return 0;

     /*--------------------------------------------------------------*/
     /*- get centroid of first triangle -----------------------------*/
     /*--------------------------------------------------------------*/
     T=TIW->Triangles[0];
     for(i=0; i<3; i++)
      Centroid[i]=(T->V1[i] + T->V2[i] + T->V3[i]) / 3.0;
     VC=AddNode(TIW, Centroid);

     /*--------------------------------------------------------------*/
     /*- for each subtriangle, get integral and error estimate and  -*/
     /*- add subtriangle to list                                    -*/
     /*--------------------------------------------------------------*/
     GetIandE(F, nFun, UserData, T->V1, T->V2, VC, IA, EA);
     AddTriangle(TIW, T->V1, T->V2, VC, IA, EA);

     GetIandE(F, nFun, UserData, T->V2, VC, T->V3, IB, EB);
     AddTriangle(TIW, T->V2, VC, T->V3, IB, EB);

     GetIandE(F, nFun, UserData, VC, T->V3, T->V1, IC, EC);
     AddTriangle(TIW, VC, T->V3, T->V1, IC, EC);

     for(nf=0; nf<nFun; nf++)
      { I[nf] += (IA[nf] + IB[nf] + IC[nf] - T->I[nf]);
        E[nf] += (EA[nf] + EB[nf] + EC[nf] - T->E[nf]);
      };

     DeleteFirstTriangle(TIW);
  
     SortTriangles(TIW);

   }; 

  return -1;

}
