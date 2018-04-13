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
 * tlibTriInt.cc
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "libTriInt.h"
#include "libdcutri.h"
#include "libhrutil.h"

int nCalls=0;
FILE *LogFile;
int NFUN

void Integrand(double *X, void *UserData, double *F) 
{
  nCalls++;

  F[0]=1.0;
  F[1]=exp( -(X[0]*X[0] + X[1]*X[1]) );
  F[2]=Poly(X,5);
  F[3]=Poly(X,7);
  F[4]=Poly(X,10);

 // fprintf(LogFile,"%e %e %e %e\n",X[0],X[1],F[0],F[1]);

}

int main(int argc, char *argv[])
{ 
  int nf;
  void *pTIW;
  double AbsTol, RelTol;
  double F[2], E[2];
  DCWorkspace *DCW;
  double V[9]={  0.0,  0.0,  0.0, 
                10.0,  0.0,  0.0,
                 0.0, 10.0,  0.0};
  double *Vertices[3]={V,V+3,V+6};


  if (argc!=3)
   { fprintf(stdout,"usage: %s AbsTol RelTol\n",argv[0]);
     exit(1);
   };
  AbsTol=strtod(argv[1],0);
  RelTol=strtod(argv[2],0);

  LogFile=fopen("doomage","w");

  /***************************************************************/
  /* compute integral using libTriInt ****************************/
  /***************************************************************/
  pTIW=CreateTIWorkspace(2, 10000);
  nCalls=0;
  TriIntAdaptive(pTIW, Integrand, 0, V+0, V+3, V+6, AbsTol, RelTol, F, E);
    
  printf("libTriInt: \n");
  printf("%i calls.\n",nCalls);
  for(nf=0; nf<2; nf++)
   printf("%i: %+12.9e %+12.9e (rel %e)\n",nf,F[nf],E[nf],E[nf]/fabs(F[nf]));

  printf("actual error=%e (rel %e)\n",
          fabs(F[0]-0.25*M_PI),fabs(F[0]-0.25*M_PI)/(0.25*M_PI));

  /***************************************************************/
  /* compute integral using libdcutri ****************************/
  /***************************************************************/
  DCW=CreateDCWorkspace(2, 10000);
  nCalls=0;
  DCUTRI(DCW, Vertices, Integrand, 0, AbsTol, RelTol, F, E);
    
  printf("libDCUTRI: \n");
  printf("%i calls.\n",nCalls);
  for(nf=0; nf<2; nf++)
   printf("%i: %+12.9e %+12.9e (rel %e)\n",nf,F[nf],E[nf],E[nf]/fabs(F[nf]));

  printf("actual error=%e (rel %e)\n",
          fabs(F[0]-0.25*M_PI),fabs(F[0]-0.25*M_PI)/(0.25*M_PI));

  fclose(LogFile);

} 
