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
 * tlibTriInt2.cc
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "libTriInt.h"
#include "libhrutil.h"
#include "libDCUTRI.h"

#define LENGTH 1.0
#define NFUN 5

int nCalls;

#define AD(x,y) fabs(x-y)

/***************************************************************/
/* polynomial in three variables of degree Degree that is      */
/* initialized on first invocation (at each degree) to have    */
/* random coefficients                                         */
/***************************************************************/
#define MAXDEGREE 12
double Poly(double *X, int Degree)
{ 
  static double *Coefficients[MAXDEGREE]={0,0,0,0,0,0,0,0,0,0,0,0};
  
  double f;
  int XPow, YPow, ZPow;
  double XFac, YFac, ZFac;
  int nc;

  if (Degree>=MAXDEGREE) 
   ErrExit("%s:%i: internal error",__FILE__,__LINE__);

  if ( Coefficients[Degree]==0 )
   { 
     int NumCoefficients=Degree*Degree*Degree;
     double NormFac = 1.0 /( (double) NumCoefficients);

     Coefficients[Degree]=(double *)malloc(NumCoefficients*sizeof(double));

     if ( !Coefficients[Degree] ) 
      ErrExit("out of memory");

     for(nc=0; nc<NumCoefficients; nc++)
      Coefficients[Degree][nc++]=NormFac*(drand48()-0.5);
   };

  f=0.0; 
  nc=0;
  for(XPow=0, XFac=1.0; XPow<Degree; XFac*=X[0], XPow++)
   for(YPow=0, YFac=1.0; YPow<Degree; YFac*=X[1], YPow++)
    for(ZPow=0, ZFac=1.0; ZPow<Degree; ZFac*=X[2], ZPow++)
     f+=Coefficients[Degree][nc++]*XFac*YFac*ZFac;

  return f;
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void Integrand(double *X, void *UserData, double *F) 
{
  nCalls++;

  F[0]=Poly(X,3);
  F[1]=Poly(X,5);
  F[2]=Poly(X,7);
  F[3]=Poly(X,9);
  F[4]=Poly(X,11);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{ 
  int Orders[]={1,4,5,7,9,13,14,16,20,25};
  int NumOrders = (sizeof(Orders)/sizeof(Orders[0]));
  int i, j, no;
  double I[NFUN], IRef[NFUN], ERef[NFUN], E[NFUN];
  void *DCW;
  double V[9];
  double *Vertices[3]={V,V+3,V+6};
  double Centroid[3];
  int nCallsDCUTRI;

  srand48(time(0));
  DCW=CreateDCUTRIWorkspace(NFUN, 100000);

  /*--------------------------------------------------------------*/
  /*- generate random triangle -----------------------------------*/
  /*--------------------------------------------------------------*/
  for(i=0; i<3; i++)
   Centroid[i]=LENGTH*(drand48()-0.5);
  for(i=0; i<3; i++)
   memcpy(Vertices[i],Centroid,3*sizeof(double));
  for(i=0; i<3; i++)
   for(j=0; j<3; j++)
    Vertices[i][j]+=LENGTH*(drand48()-0.5);

  /*--------------------------------------------------------------*/
  /*- use DCUTRI to get an 'exact' value for the integrals       -*/
  /*--------------------------------------------------------------*/
  nCalls=0;
  DCUTRI(DCW, Vertices, Integrand, 0, 0.0, 1.0e-20, IRef, ERef);
  nCallsDCUTRI=nCalls;

  /*--------------------------------------------------------------*/
  /*- use TriIntFixed to get estimates of the integrals ----------*/
  /*--------------------------------------------------------------*/
  printf("%3s (%4s): %8s  %8s  %8s  %8s  %8s\n","ORD","NFUN",
         " DEG 3  "," DEG 5  "," DEG 7  "," DEG 9  "," DEG 11 ");
  printf("%3s--%4s---%8s--%8s--%8s--%8s--%8s\n","---","----",
         "--------","--------","--------","--------","--------");

  for(no=0; no<NumOrders; no++)
   { 
     nCalls=0;

     TriIntFixed(Integrand, NFUN, 0, Vertices[0], Vertices[1], Vertices[2], Orders[no], I);
        
     printf("%3i (%4i): %8.2e  %8.2e  %8.2e  %8.2e  %8.2e\n",
             Orders[no], nCalls,
             AD(I[0],IRef[0]), AD(I[1],IRef[1]), AD(I[2],IRef[2]), AD(I[3],IRef[3]), AD(I[4],IRef[4]));
   };
  printf("INF (%4i): %8.2e  %8.2e  %8.2e  %8.2e  %8.2e\n",
          nCallsDCUTRI, ERef[0], ERef[1], ERef[2], ERef[3], ERef[4]);

  /*--------------------------------------------------------------*/
  /*- use TriIntEmbedded to estimate errors                       */
  /*--------------------------------------------------------------*/
  TriIntEmbedded(Integrand, NFUN, 0, Vertices[0], Vertices[1], Vertices[2], I, E);
  printf("(Absolute) Errors reported by embedded scheme: \n");
  printf(" Degree 3:  %.3e  (actual %.3e)\n",E[0],fabs(I[0]-IRef[0]));
  printf(" Degree 5:  %.3e  (actual %.3e)\n",E[1],fabs(I[1]-IRef[1]));
  printf(" Degree 7:  %.3e  (actual %.3e)\n",E[2],fabs(I[2]-IRef[2]));
  printf(" Degree 9:  %.3e  (actual %.3e)\n",E[3],fabs(I[3]-IRef[3]));
  printf(" Degree 11: %.3e  (actual %.3e)\n",E[4],fabs(I[4]-IRef[4]));

}
