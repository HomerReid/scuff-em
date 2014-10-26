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
 * tLUSolve.cc 
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <libhrutil.h>
#include "libhmat.h"

#if defined(_WIN32)
#  define srand48 srand
#  define drand48 my_drand48
static double my_drand48(void) {
  return rand() * 1.0 / RAND_MAX;
}
#endif

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{ 
  /*--------------------------------------------------------------*/
  /*- process options  -------------------------------------------*/
  /*--------------------------------------------------------------*/
  int N=1000;
  int Complex=0;
  char *Flag=0;
  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { {"N",       PA_INT,     1, 1, (void *)&N,       0, "dimension "},
     {"Complex", PA_BOOL,    0, 1, (void *)&Complex, 0, "complex-valued matrix"},
     {"Flag",    PA_STRING,  1, 1, (void *)&Flag,    0, "either N, C, or T"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  if (Flag==0)
   Flag=strdupEC("N");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  printf("Creating random %ix%i matrices ... \n",N,N);
  HMatrix *M1=new HMatrix(N, N, Complex ? LHM_COMPLEX : LHM_REAL);
  HMatrix *M2=new HMatrix(N, N, Complex ? LHM_COMPLEX : LHM_REAL);

  cdouble II = Complex ? cdouble(0.0,1.0) : cdouble(0.0,0.0);
  srand48(time(0));

  int m,n;
  double Elapsed;
  Tic();
  for(m=0; m<N; m++)
   for(n=0; n<N; n++)
    { M1->SetEntry(m, n, drand48() + II*drand48());
      M2->SetEntry(m, n, drand48() + II*drand48());
    };
  Elapsed=Toc();
  printf("...%.3f s\n",Elapsed);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  printf("LU-factorizing M1...");
  Tic();
  M1->LUFactorize();
  Elapsed=Toc();
  printf("...%.3f s\n",Elapsed);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  printf("LU-solving M2 (Flag=%s)...",Flag);
  Tic();
  M1->LUSolve(M2,Flag[0]);
  Elapsed=Toc();
  printf("...%.3f s\n",Elapsed);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  HMatrix *M3=new HMatrix(N, N, Complex ? LHM_COMPLEX : LHM_REAL );
  if (M3)
   { 
     printf("Multiplying M1*M2 ...\n");
     Tic();
     M1->Multiply(M2,M3);
     Elapsed=Toc();
     printf("...%.3f s\n",Elapsed);
   };

}
