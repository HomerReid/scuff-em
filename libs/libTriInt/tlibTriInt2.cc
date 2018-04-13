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
#include <libhrutil.h>
#include <libDCUTRI.h>

double Kappa=1.0;
int nCalls;
#define NFUN 3

/***************************************************************/
/***************************************************************/
/***************************************************************/
void Integrand(double *X, void *UserData, double *F) 
{
  double r=sqrt(X[0]*X[0]+X[1]*X[1]+X[2]*X[2]);
  nCalls++;

  //F[0]=(exp(-Kappa*r)-1+(Kappa*r)-0.5*(Kappa*Kappa*r*r))/r;
  //F[1]=(cos(Kappa*r)-1+0.5*(Kappa*Kappa*r*r))/r;
  F[0]=exp(-Kappa*r) /(r+1.0e-6);
  F[1]=cos(Kappa*r) /(r+1.0e-6);
  F[2]=cos(Kappa*X[0])*sin(Kappa*X[1]);
}

/***************************************************************/
/* utility vector routines *************************************/
/***************************************************************/
double *VecCross(double *v1, double *v2, double *v3)
{ v3[0]=v1[1]*v2[2] - v1[2]*v2[1];
  v3[1]=v1[2]*v2[0] - v1[0]*v2[2];
  v3[2]=v1[0]*v2[1] - v1[1]*v2[0];
  return v3;
}

/* v3 = v1 - v2 */
double *VecSub(double *v1, double *v2, double *v3)
{ v3[0]=v1[0] - v2[0];
  v3[1]=v1[1] - v2[1];
  v3[2]=v1[2] - v2[2];
  return v3;
}

double VecDot(double *v1, double *v2)
{ return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]; }

double VecNorm2(double *v)
{ return VecDot(v,v); }

double VecNorm(double *v)
{ return sqrt(VecDot(v,v)); }

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{ 
  int Orders[]={1,2,4,5,7,9,13,14,16,20,25};
  int NumOrders = (sizeof(Orders)/sizeof(Orders[0]));
  int i, j, no, nPts;
  double I[NFUN], IRef[NFUN];
  double E[NFUN], ERef[NFUN];
  void *DCW;
  double V[9];
  double *Vertices[3]={V,V+3,V+6};
  double Centroid[3];
  FILE *f;

  /*--------------------------------------------------------------*/
  /*- process command-line arguments -----------------------------*/
  /*--------------------------------------------------------------*/
  ArgStruct ASArray[]=
   { {"Kappa",     PA_DOUBLE, (void *)&Kappa,     "1.0",  "Kappa"},
     {0,0,0,0,0}
   };
  ProcessArguments(argc, argv, ASArray);

  /*--------------------------------------------------------------*/
  /*- generate a random triangle ---------------------------------*/
  /*--------------------------------------------------------------*/
  for(i=0; i<3; i++)
   Centroid[i]=(drand48()-0.5) + 1.0;
  for(i=0; i<3; i++)
   memcpy(Vertices[i],Centroid,3*sizeof(double));
  for(i=0; i<3; i++)
   for(j=0; j<3; j++)
    Vertices[i][j]+=(drand48()-0.5);
  
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  DCW=CreateDCUTRIWorkspace(NFUN, 100000);
  nCalls=0;
  DCUTRI(DCW, Vertices, Integrand, 0, 0.0, 1.0e-10, IRef, ERef);
  int nCallsDCUTRI=nCalls;

  /*--------------------------------------------------------------*/
  /*- use TriIntFixed to get estimates of the integrals ----------*/
  /*--------------------------------------------------------------*/
  printf("%3s (%4s): %14s  %8s  %14s  %8s  %14s  %8s\n","ORD","NFUN",
         "   VAL 1    "," AD  1  ", "   VAL 2    "," AD  2  ", "   VAL 3    "," AD  3  ");
  printf("%3s--%4s---%14s--%8s--%14s--%8s--%14s--%8s\n","---","----",
         "--------------","--------",
         "--------------","--------",
         "--------------","--------");

  for(no=0; no<NumOrders; no++)
   { 
     nCalls=0;

     TriIntFixed(Integrand, NFUN, 0, Vertices[0], Vertices[1], Vertices[2], Orders[no], I);
        
     printf("%3i (%4i): %14.6e  %8.2e  %14.6e  %8.2e  %14.6e  %8.2e\n",
             Orders[no], nCalls, 
             I[0], fabs(I[0]-IRef[0]), I[1], fabs(I[1]-IRef[1]), I[2], fabs(I[2]-IRef[2]));
   };
  printf("INF (%4i): %14.6e  %8.2e  %14.6e  %8.2e  %14.6e  %8.2e\n", 
          nCallsDCUTRI, IRef[0], ERef[0], IRef[1], ERef[1], IRef[2], ERef[2]);

  /*--------------------------------------------------------------*/
  /*- use TriIntEmbedded to estimate errors                       */
  /*--------------------------------------------------------------*/
#if 0
  TriIntEmbedded(Integrand, NFUN, 0, Vertices[0], Vertices[1], Vertices[2], I, E);
  printf("\n");
  printf("(Absolute) Errors reported by embedded scheme: \n");
  printf(" F1:        %.3e  (actual %.3e)\n",E[0],fabs(I[0]-IRef[0]));
  printf(" F2:        %.3e  (actual %.3e)\n",E[1],fabs(I[1]-IRef[1]));
  printf(" F3:        %.3e  (actual %.3e)\n",E[2],fabs(I[2]-IRef[2]));
#endif

}
