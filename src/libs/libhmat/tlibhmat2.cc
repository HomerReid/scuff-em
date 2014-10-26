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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "libhrutil.h"
#include "libhmat.h"

#if defined(_WIN32)
#  define srand48 srand
#  define drand48 my_drand48
static double my_drand48(void) {
  return rand() * 1.0 / RAND_MAX;
}
#endif

int main(int argc, char *argv[])
{ 
  HMatrix *M, *MSymm, *dM, *dMCopy;
  double Trace1, Trace2, Elapsed;
  int m, n, N;
  
  /*--------------------------------------------------------------*/
  /*- process command-line arguments -----------------------------*/
  /*--------------------------------------------------------------*/
  ArgStruct ASArray[]=
   { {"N",         PA_INT,    (void *)&N,            "0",  "dimension"},
     {0,0,0,0,0}
   };
  ProcessArguments(argc, argv, ASArray);
  if (N==0)
   ASUsage(argv[0],ASArray,"--N option is mandatory");

  /*--------------------------------------------------------------*/
  /*- allocate matrices ------------------------------------------*/
  /*--------------------------------------------------------------*/
  M=new HMatrix(N,N,LHM_REAL);
  MSymm=new HMatrix(N,N,LHM_REAL,LHM_SYMMETRIC);
  dM=new HMatrix(N,N/2);
  dMCopy=new HMatrix(N,N/2);

  /*--------------------------------------------------------------*/
  /*- initialize matrices to random entries ----------------------*/
  /*--------------------------------------------------------------*/
  srand48(time(0));
  for(m=0; m<N; m++)
   for(n=m; n<N; n++)
    { M->SetEntry(m,n,drand48());
      M->SetEntry(n,m,M->GetEntryD(m,n));
      MSymm->SetEntry(m,n,M->GetEntryD(m,n));
    };

  for(m=0; m<N; m++)
   for(n=0; n<N/2; n++)
    dM->SetEntry(m,n,drand48());
  dMCopy->Copy(dM);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Tic();
  M->LUFactorize();
  Elapsed=Toc();
  printf("LU(M)        : %9.3e s \n",Elapsed);

  Tic();
  MSymm->LUFactorize();
  Elapsed=Toc();
  printf("LU(MSymm)    : %9.3e s \n",Elapsed);

  Tic();
  M->LUSolve(dM);
  Elapsed=Toc();
  for(Trace1=0.0, m=0; m<N/2; m++)
   Trace1+=dM->GetEntryD(m,m);
  printf("Solve(M)     : %9.3e s (%+10.3e)\n",Elapsed,Trace1);

  Tic();
  MSymm->LUSolve(dMCopy);
  Elapsed=Toc();
  for(Trace2=0.0, m=0; m<N/2; m++)
   Trace2+=dMCopy->GetEntryD(m,m);
  printf("Solve(Msymm) : %9.3e s (%+10.3e)\n",Elapsed,Trace2);

}
